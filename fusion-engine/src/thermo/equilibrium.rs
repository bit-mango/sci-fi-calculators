// At some fixed T, P you want the mole numbers n_j of each species that minmize total Gibbs energY:
// G/RT = Σ_j n_j * µ_j/RT
// Where the chemical potential of species j is:
// µ_j/RT = g_j + ln(n_j/n) + ln(P)
// g_j = G°_j/RT: The species standard state reduced Gibbs energy(use NASA polynomials)
// n = Σ n_j: total moles of the gas
// ln(n_j/n): The mixing/entropy term which is a mole fraction(n_j/n)
// ln(P): pressure correction (P in atm, standard state 1 atm)
//
// The atoms have to add up to what you started with
// Σ_j a_ij * n_j = b_i for each element i
// a_ij: atoms of element i in species j(a fixed known matrix). b_i: total moles of element i (fixed by reactants).
// So if i is rows and j is columns, then for the reactants: H2O we would have something like.
//    j  H2O  H2  O2
// i
// H      2   2   0
// O      1   0   2
//
// Go from Langrangian -> element potentials
// Introduce 1 Lagrange multiplier π_i per element( literally called element potentials)
// The stationarity condition (∂L/∂n_j = 0) gives, for every species:
// g_j + ln(n_j/n) + ln(P) = Σ a_ij * π_i: The chemical potential of species j equals the sum of the potentials of the elements its built from
// When you rearrange it, it tells you n_j directly:
// ln(n_j) = -g_j +ln(P) + Σ a_ij*π_i + ln(n)
//
// Newton's method is actually solving for π_i(one per element) and ln(n) (total moles of products)
// You need 3 equations:
// Mass balance for H: 2*n_H2O + 2*n_H2 = b_H
// Mass balance for O: 1*n_H2O + 2*n_O2 = b_O
// Definition of total moles = n_H2O + n_H2 + n_O2 = n
// n_j is an is nonlinear(exponential) function of π_H, π_O, ln(n) these three equations are nonlinear, so we use Newton's method.
// The Jacobian entries are just derivatives of n_j = exp(...)
// ∂n_j/∂π_k = a_kj * n_j
// ∂n_j/∂ln(n) = n_j
//
use nalgebra::{DMatrix, DVector};
use std::collections::HashSet;
use thermo_species::{Constituent, Species};

/// Temporary window into solve_for_products' internal state, just enough to
/// validate D.1-D.5 against the guide's numbers/tests. Not a final API —
/// Section G will need a real Result/error type once singular-matrix
/// recovery and the outer driver exist.
struct TpSolveDebug {
    /// The assembled (pre-solve) matrix from the very first Newton
    /// iteration, at the D.6 canonical initial guess.
    assembled_first_iteration: DMatrix<f64>,
    converged: bool,
    iterations_used: usize,
    /// Final gas amounts, RAW (exp(ln_nj[j])), length ng.
    nj: Vec<f64>,
    n: f64,
    /// Final temperature — equal to state1 for TP, solved-for otherwise.
    temperature_k: f64,
    /// F.1: which condensed species (indices into the condensed slice,
    /// i.e. species[ng+c]) ended up active, and their final linear amounts.
    condensed_active: Vec<bool>,
    condensed_nj: Vec<f64>,
    /// The actual gas/condensed species lists solve_for_products settled
    /// on (order matches `nj`/`condensed_active`/`condensed_nj`) — needed
    /// by callers/tests when `only` wasn't given, so the product set was
    /// auto-generated rather than known in advance.
    gas_species: Vec<Species>,
    condensed_species: Vec<Species>,

    // --- Section I: post-processing & thermodynamic derivatives ---
    // Computed once at the final state (converged or not — matches the
    // guide's "never a hard failure" framing for I.2/I.3); never gates
    // `converged` itself.
    pressure_bar: f64,
    volume_m3_per_kg: f64,
    enthalpy_kj_per_kg: f64,
    energy_kj_per_kg: f64, // U
    entropy_kj_per_kg_k: f64,
    gibbs_kj_per_kg: f64,
    molar_mass_kg_per_kmol: f64,             // M = 1/n
    mean_molecular_weight_kg_per_kmol: f64,  // MW, gas-only (Eq. 2.4)
    cp_frozen_kj_per_kg_k: f64,
    cv_frozen_kj_per_kg_k: f64,
    // I.2: equilibrium (reacting) specific heats and isentropic exponent.
    // Fall back to I.1/I.3's degrade-gracefully constants if the partials
    // system (which shares the converged matrix's LHS) turns out singular
    // — rare, since that same LHS structure already solved successfully.
    cp_eq_kj_per_kg_k: f64,
    cv_eq_kj_per_kg_k: f64,
    gamma_s: f64,
    // Mole/mass fractions, split gas/condensed like gas_species/condensed_species.
    gas_mole_fractions: Vec<f64>,
    condensed_mole_fractions: Vec<f64>,
    gas_mass_fractions: Vec<f64>,
    condensed_mass_fractions: Vec<f64>,
}

struct MixtureThermo {
    // All length NG+NC, reduced units(num_gas + num_condensed)
    cp: Vec<f64>,       // Cp/R
    cv: Vec<f64>,       // Cp/ R -1
    enthalpy: Vec<f64>, // H/(R*T)
    entropy: Vec<f64>,  // S/R
    energy: Vec<f64>,   // U/(R*T) = enthalpy - 1
}

fn calc_thermo(
    species: &[Species],
    temperature_k: f64,
    ng: usize, // Number of gas species.
    nc: usize, // Number of condensed species.
    include_condensed: bool,
    out: &mut MixtureThermo,
) {
    for j in 0..ng {
        out.cp[j] = species[j].data().cp_over_r(temperature_k);
        out.cv[j] = species[j].data().cp_over_r(temperature_k) - 1.0;
        // h_over_r/u_over_r return H/R and U/R (Kelvin units, guide A.1);
        // the solver needs the reduced H/(RT) and U/(RT), hence the /T here.
        out.enthalpy[j] = species[j].data().h_over_r(temperature_k) / temperature_k;
        out.entropy[j] = species[j].data().s_over_r(temperature_k);
        out.energy[j] = species[j].data().u_over_r(temperature_k) / temperature_k;
    }
    if include_condensed {
        for j in ng..(ng + nc) {
            out.cp[j] = species[j].data().cp_over_r(temperature_k);
            out.cv[j] = species[j].data().cp_over_r(temperature_k) - 1.0;
            out.enthalpy[j] = species[j].data().h_over_r(temperature_k) / temperature_k;
            out.entropy[j] = species[j].data().s_over_r(temperature_k);
            out.energy[j] = species[j].data().u_over_r(temperature_k) / temperature_k;
        }
    }
}

// Solve G*x = rhs in place. G is (n, n+1) augmented. Returns Err(k)
// with the 0-based column index that could not be pivoted.
fn gauss(g: &mut DMatrix<f64>) -> Result<(), usize> {
    let n = g.nrows();
    for k in 0..(n - 1) {
        //  Custom pivot search
        let mut pivot = None;
        let mut min_ratio = 1e25;
        for i in k..n {
            let gk = g[(i, k)].abs();
            let ratio = if gk != 0.0 {
                let mut m = 0.0_f64;
                for j in (k + 1)..=n {
                    m = m.max(g[(i, j)].abs());
                }
                m / gk
            } else {
                1e25
            };
            if ratio < min_ratio {
                min_ratio = ratio;
                pivot = Some(i);
            }
        }

        let Some(p) = pivot else { return Err(k) };
        if min_ratio >= 1e25 {
            return Err(k);
        }
        g.swap_rows(p, k);
        // Eliminate
        for j in (k + 1)..=n {
            g[(k, j)] /= g[(k, k)];
        }
        for i in (k + 1)..n {
            for j in (k + 1)..=n {
                g[(i, j)] -= g[(i, k)] * g[(k, j)];
            }
        }
    }
    // Back substitution
    g[(n - 1, n)] /= g[(n - 1, n - 1)];
    for k in (0..n - 1).rev() {
        for i in (k + 1)..n {
            let t = g[(k, i)] * g[(i, n)];
            g[(k, n)] -= t;
        }
    }

    Ok(())
}

// D.2 (RP-1311 Eqs. 3.1-3.3): caps the Newton step so no gas species jumps by
// more than a couple orders of magnitude in one iteration. dln_t is 0.0 for
// TP problems (no ΔlnT unknown until Section E adds variable-T problems).
fn damping_lambda(ln_nj: &[f64], dln_nj: &[f64], dln_n: f64, dln_t: f64, n: f64, size: f64) -> f64 {
    const FACTOR: f64 = -9.2103404; // ln(1e-4)
    let log_n = n.ln();

    // l1_denom starts as the "5x the T/n step" floor, then becomes the
    // largest Δln n_j among species not already near the trace floor.
    let mut l1_denom = f64::max(5.0 * dln_t.abs(), 5.0 * dln_n.abs());
    let mut lambda1 = 1.0_f64;
    let mut lambda2 = 1.0_f64;

    for j in 0..ln_nj.len() {
        let d = dln_nj[j];
        if d <= 0.0 {
            continue; // only growth needs damping
        }
        if ln_nj[j] - log_n + size <= 0.0 {
            // "small" species: below the SIZE threshold (n_j < n * 1e-8),
            // trying to grow. Cap how far it can jump toward the 1e-4*n
            // reseed level in this step.
            let l2_denom = (d - dln_n).abs();
            if l2_denom >= size + FACTOR {
                lambda2 = lambda2.min((FACTOR - ln_nj[j] + log_n).abs() / l2_denom);
            }
        } else if d > l1_denom {
            // "large" species: track the biggest positive step.
            l1_denom = d;
        }
    }
    if l1_denom > 2.0 {
        lambda1 = 2.0 / l1_denom;
    }

    1.0_f64.min(lambda1).min(lambda2)
}

// AI.3.9: condensed add/remove/phase-swap and singular-recovery events are
// the primary debugging tool when a condensed set oscillates. Set CEA_DEBUG=1
// to mirror CEA's logging of them.
macro_rules! cea_debug {
    ($($arg:tt)*) => {
        if std::env::var_os("CEA_DEBUG").is_some() {
            eprintln!($($arg)*);
        }
    };
}

// F.4 (`check_condensed_phases`): phase-validity and phase-transition logic
// for condensed species with sibling phases (H2O(cr) vs H2O(L), Al2O3(a) vs
// Al2O3(L), ...). Sibling discovery uses Species::phases() (the build-time
// equivalent of CEA's trailing-suffix stripping), sorted by max fit
// temperature ("melting point") low→high; `phase()` is thermo.inp's raw ifaz,
// where adjacent phases of one substance differ by 1. At most ONE change per
// call, reported by the return value so the caller restarts its epoch.
//
// One deviation from the Fortran: the melting-point coexistence path pins T
// to the transition temperature, which only makes sense when T is a Newton
// unknown — for const-T problems the user's T stays untouched (the doomed
// phase then dies through F.3.1's negative-amount removal instead).
#[allow(clippy::too_many_arguments)]
fn check_condensed_phases(
    species: &[Species],
    ng: usize,
    const_t: bool,
    temperature_k: &mut f64,
    nj_c: &mut [f64],
    is_active: &mut [bool],
    active_ranked: &mut Vec<usize>,
    j_switch: &mut Option<usize>,
    j_sol: &mut Option<usize>,
    j_liq: &mut Option<usize>,
) -> bool {
    const DT_PHASE: f64 = 50.0; // K, "near the transition" half-width
    let cond_index = |sp: Species| species[ng..].iter().position(|&s| s == sp);

    let snapshot = active_ranked.clone();
    for &i in &snapshot {
        // The pinned solid+liquid pair is exempt from further checks.
        if Some(i) == *j_sol || Some(i) == *j_liq {
            continue;
        }
        let sc_i = ng + i;
        let (t_low_i, t_high_i) = species[sc_i].data().valid_temperature_range();

        // Sibling phases (self included), sorted by melting point low→high.
        let mut siblings: Vec<usize> = species[sc_i]
            .phases()
            .iter()
            .filter_map(|&sp| cond_index(sp))
            .collect();
        siblings.sort_by(|&a, &b| {
            let ta = species[ng + a].data().valid_temperature_range().1;
            let tb = species[ng + b].data().valid_temperature_range().1;
            ta.partial_cmp(&tb).unwrap()
        });
        let max_t_sibling = siblings
            .iter()
            .map(|&j| species[ng + j].data().valid_temperature_range().1)
            .fold(0.0_f64, f64::max);

        for &j in &siblings {
            let (_, t_high_j) = species[ng + j].data().valid_temperature_range();
            if *temperature_k > t_high_j {
                continue; // sibling j can't exist at T
            }
            let d_phase =
                species[sc_i].data().phase() as i32 - species[ng + j].data().phase() as i32;
            if d_phase == 0 {
                break; // reached i itself (same phase id) — stop scanning
            }
            if is_active[j] {
                continue; // both already active; nothing to transition
            }
            let out_of_range =
                *temperature_k < t_low_i - DT_PHASE || *temperature_k > t_high_i + DT_PHASE;
            if d_phase.abs() > 1 || (Some(j) != *j_switch && out_of_range) {
                cea_debug!(
                    "[F.4] replacing {} with {} at T={temperature_k}",
                    species[sc_i].data().symbol(),
                    species[ng + j].data().symbol()
                );
                // REPLACE i by j: transfer the amount, preserve i's rank.
                let pos = active_ranked.iter().position(|&x| x == i).unwrap();
                active_ranked[pos] = j;
                is_active[i] = false;
                is_active[j] = true;
                nj_c[j] = nj_c[i];
                nj_c[i] = 0.0;
                *j_switch = Some(i);
                *j_sol = None;
                *j_liq = None;
            } else {
                // Melting-point coexistence (either we just swapped i in for
                // j and the iterate is bouncing back — j == j_switch — or T
                // sits within ±50 K of i's own range edge): keep BOTH
                // phases, pin T to the transition temperature (variable-T
                // problems only), split the amount 50/50.
                cea_debug!(
                    "[F.4] coexistence {} + {} at T={temperature_k}",
                    species[sc_i].data().symbol(),
                    species[ng + j].data().symbol()
                );
                if !const_t {
                    *temperature_k = t_high_i.min(t_high_j);
                }
                if t_high_i > t_high_j {
                    *j_sol = Some(j);
                    *j_liq = Some(i);
                } else {
                    *j_sol = Some(i);
                    *j_liq = Some(j);
                }
                is_active[j] = true;
                active_ranked.insert(0, j);
                nj_c[i] /= 2.0;
                nj_c[j] = nj_c[i];
            }
            return true;
        }

        if *temperature_k > 1.2 * max_t_sibling {
            cea_debug!(
                "[F.4] removing {} — T={temperature_k} above 1.2*{max_t_sibling}",
                species[sc_i].data().symbol()
            );
            // Way above every sibling's range: remove i entirely.
            is_active[i] = false;
            nj_c[i] = 0.0;
            active_ranked.retain(|&x| x != i);
            *j_switch = None;
            *j_sol = None;
            *j_liq = None;
            return true;
        }
    }
    false
}

// unused_assignments: pi_e's zeroing/accumulation sites mirror CEA's control
// flow literally (guide H.2's ⚠ note); some are dead under H.3's overwrite,
// exactly as in the Fortran.
#[allow(unused_assignments)]
fn solve_for_products(
    reactants: &[(f64, Species)],
    // Section E.1: which state1/state2 pair this problem fixes. Only 't'
    // (TP), 'h' (HP), and 's' (SP) are implemented so far — TV/UV/SV need
    // const_v support, which isn't written yet.
    const_t: bool,
    // Only meaningful when !const_t: true = SP (state1 = S⁰/R), false = HP
    // (state1 = H⁰/R). Selects the 's·' vs 'h·' branch of E.2's tmp_j table.
    const_s: bool,
    // true = P fixed (state2 = P [bar], TP/HP/SP); false = V fixed
    // (state2 = V [m³/kg], TV — UV/SV aren't wired up yet). When !const_p,
    // pressure is recomputed every iteration from the ideal-gas closure
    // (Global conventions: P = 1e-5 * R * n * T / V) instead of being fixed.
    const_p: bool,
    // state1: T [K] if const_t; H⁰/R [kmol·K/kg] if HP; S⁰/R [kmol/kg] if SP.
    state1: f64,
    // state2: P [bar] if const_p, else V [m³/kg].
    state2: f64,
    include_condensed: bool,
    include_ions: bool,
    // Restrict products to exactly this list (mirrors CEA's `only` keyword),
    // instead of auto-generating every species compatible with the element
    // pool. Mainly for tests, where the ground truth was computed for a
    // specific, curated product list.
    only: Option<&[Species]>,
) -> TpSolveDebug {
    let mut constituent_set: HashSet<Constituent> = HashSet::new();
    reactants.iter().for_each(|(_, parent_species)| {
        let children = parent_species.data().constituents();
        children.iter().for_each(|&(_, child_species)| {
            constituent_set.insert(child_species);
        });
    });
    if !include_ions && constituent_set.contains(&Constituent::E) {
        panic!("Ions disabled but reactants have ions!");
    }
    if include_ions && !constituent_set.contains(&Constituent::E) {
        // Could be introduced by solver.
        constituent_set.insert(Constituent::E);
    }

    let possible_species: Vec<Species> = Species::all()
        .iter()
        .filter(|&&species_i| {
            // In order for a species to be possible every one if its constituents must be present.
            !species_i
                .data()
                .constituents()
                .iter()
                .any(|(_, constituent_i)| !constituent_set.contains(constituent_i))
        })
        .filter(|&&species_i| only.is_none_or(|list| list.contains(&species_i)))
        .map(|&species_i| species_i)
        .collect();
    let mut species: Vec<Species> = possible_species
        .iter()
        .filter(|species_i| species_i.data().phase() == 0)
        .map(|&species_i| species_i)
        .collect();
    // Append condensed to the END of the species
    let ng = species.len();
    let nc;
    if include_condensed {
        species.extend(
            possible_species
                .iter()
                .filter(|species_i| species_i.data().phase() > 0)
                .map(|&species_i| species_i),
        );
        nc = species.len() - ng;
    } else {
        nc = 0;
    }
    let constituent_pool = if constituent_set.contains(&Constituent::E) {
        // Remove it temporarily.
        constituent_set.remove(&Constituent::E);
        let mut pool: Vec<Constituent> = constituent_set.drain().collect();
        // Sort it for determinism.
        pool.sort();
        // Add it back so its at the end.
        pool.push(Constituent::E);
        pool
    } else {
        let mut pool: Vec<Constituent> = constituent_set.drain().collect();
        // Sort it for determinism.
        pool.sort();
        pool
    };

    // H.1: E (the electron pseudo-element) is always last when present —
    // already guaranteed by the construction above. `stoich[j][ne-1]` is
    // then species j's charge count (D.1/H.1: -1 * actual charge).
    let has_electron = constituent_pool.last() == Some(&Constituent::E);

    let ne = constituent_pool.len();
    let ns = species.len();

    let mut stoich = vec![vec![0.0; ne]; ns];
    for (j, sp) in species.iter().enumerate() {
        for (i, constituent) in constituent_pool.iter().enumerate() {
            for (moles, c) in sp.data().constituents() {
                if c == constituent {
                    // We have a match!
                    stoich[j][i] = *moles;
                }
            }
        }
    }

    let b0 = element_amounts(reactants, &constituent_pool);

    // D.4 gas truncation threshold: SIZE = ln(1e8). Below n/1e8 a species
    // contributes 0 to every sum, though its ln_nj keeps evolving. Mutable
    // because G.2's correct_singular sets it to 80 (disabling truncation)
    // as its first action on any singular matrix.
    let mut tsize: f64 = 18.420681;
    // D.4/D.5: the looser threshold tsize relaxes to on each convergence
    // (default ≈1e-11; would be -ln(trace) if a user trace level existed) —
    // this second stage is how CEA resolves trace species accurately after
    // the main composition settles.
    let mut xsize: f64 = 25.328436;
    // H.1: charged species get the tighter esize threshold instead of tsize
    // (ions live at much smaller mole fractions) = min(80, xsize + ln 1000).
    let mut esize: f64 = 32.236191;
    // D.5 tolerances (guide's Global conventions).
    const NJ_TOL: f64 = 0.5e-5;
    const B_TOL: f64 = 1e-6;
    const T_TOL: f64 = 1e-4;
    const S_TOL: f64 = 0.5e-4;
    // CEA's legacy R (Global conventions) — only needed here, for the P/V
    // closure; every reduced thermo quantity elsewhere is already R-free.
    const R_CEA: f64 = 8314.51;

    // H.1/H.3: whether ions are still part of the active solve (can
    // permanently turn off if every charged species dies out), and the
    // electron multiplier from the last H.2 polish, fed back into the next
    // iteration's species recovery (E.3-style, but for charge).
    let mut active_ions = include_ions && has_electron;
    let mut pi_e: f64 = 0.0;

    // D.6: T = state1 for TP; otherwise starts at the fixed 3800K guess and
    // becomes an unknown (ΔlnT) that the Newton iteration solves for.
    let mut temperature_k = if const_t { state1 } else { 3800.0 };

    let mut ln_nj = vec![(0.1 / ng as f64).ln(); ns];
    let mut n: f64 = 0.1; // D.6: total gas moles, tracked as its own state var.
    // D.3/D.5's "stored" nj: truncated/effective gas amounts, updated once
    // per iteration and read by check_convergence's tests 1 and 3. Distinct
    // from the RAW exp(ln_nj[j]) that test 4 (element balance) uses instead.
    let mut nj_stored = vec![0.1 / ng as f64; ng];

    // F.1: condensed state. nj_c is LINEAR (not log-space), one entry per
    // condensed species (indices into 0..nc, i.e. species[ng+c]). Invariant
    // relied on below: nj_c[c] == 0.0 whenever !is_active[c] — both
    // activation and deactivation always reset it, so every "sum over all
    // species" below can just use nj_c directly without re-checking
    // is_active. No rank ordering (F.4's phase-transition logic needs it;
    // not implemented — see the module note on Section F's scope).
    let mut nj_c = vec![0.0; nc];
    let mut is_active = vec![false; nc];
    // F.1: active condensed species ordered by activation rank, NEWEST FIRST
    // — this ordering decides matrix row order and therefore how a singular
    // ierr in the condensed block maps back to a species (G.2c / AI.1.2).
    let mut active_ranked: Vec<usize> = Vec::new();
    // F.4 phase-transition bookkeeping: the phase most recently swapped out
    // (re-encountering it means melting-point coexistence), and the pinned
    // solid+liquid pair at a melting point.
    let mut j_switch: Option<usize> = None;
    let mut j_sol: Option<usize> = None;
    let mut j_liq: Option<usize> = None;
    // F.3.3/G.2: a condensed species removed because it made the matrix
    // singular; the insertion test must never re-add it (and per the guide,
    // its delta-G still being negative afterwards is a hard failure).
    let mut singular_blocklist: Option<usize> = None;
    // F.3.4/G.1: how many times the inner Newton loop has converged this
    // solve — repeated re-convergence means the condensed set is thrashing.
    let mut times_converged: usize = 0;
    // H.2/H.3: whether the electron-balance polish has converged, and whether
    // the last convergence check's base tests passed (the pair gates the
    // pi_e feedback: it fires only on the Newton sweep right after a
    // converged-but-ions-unpolished check).
    let mut ions_converged = false;
    let mut base_converged_prev = false;
    // G.2(a): the condensed species most recently activated by F.3.3 —
    // needed because a singular matrix on the very first attempt after an
    // insertion usually means the new species conflicts with an existing
    // one over some element, not that the row indicated by `ierr` itself is
    // individually degenerate.
    let mut last_cond_idx: Option<usize> = None;
    // pi_prev doubles as both "the pi from the last successful solve" (what
    // F.3.3's insertion test calls pi_prev) and "the current pi" used to
    // recover dln_nj (D.1/E.3) — the guide sets pi_prev = pi after every
    // successful solve regardless of convergence, so within one solve these
    // are the same array. G.2(e)'s element reduction zeroes an entry here
    // when it excludes that element from the matrix.
    let mut pi_prev = vec![0.0; ne];
    // G.2(e): elements whose row/column has been dropped from the matrix
    // because they caused an unrecoverable-any-other-way singularity
    // (linearly dependent element balance). Persists for the rest of this
    // solve — restoring it is only about bookkeeping/output naming in the
    // guide, which we don't expose, so there's nothing to undo here.
    let mut element_excluded = vec![false; ne];

    // D.7 default budget, now per "condensed epoch" (F.3 resets it whenever
    // the active set changes) — MAX_TOTAL_ITERATIONS is the hard safety cap
    // across all epochs combined, standing in for Section G's thrash guard.
    const MAX_ITERATIONS: usize = 50;
    const MAX_TOTAL_ITERATIONS: usize = 400;
    const MAX_SINGULAR_RETRIES: usize = 8; // G.1: give up after this many.
    let mut assembled_first_iteration: Option<DMatrix<f64>> = None;
    let mut converged = false;
    let mut iterations_used = 0;
    let mut iteration = 0usize;
    let mut times_singular = 0usize;
    // G.3: fires at most once per solve.
    let mut g3_fallback_used = false;
    // HP is const_t==false, const_p==true, const_s==false (the only
    // problem type among TP/HP/SP/TV/UV/SV that matches that combination).
    let is_hp = !const_t && const_p && !const_s;

    for total_iteration in 0..MAX_TOTAL_ITERATIONS {
        if iteration >= MAX_ITERATIONS {
            // G.3: max-iteration high-T condensed fallback. If the iterate
            // has collapsed into a single condensed phase with essentially
            // no gas left (n <= 1e-4), the Newton state has nothing left to
            // work with — reset to the uniform initial guess and drop that
            // condensed species, then keep going instead of giving up.
            let single_active_condensed = (0..nc).filter(|&c| is_active[c]).count() == 1;
            if !g3_fallback_used && (!is_hp || temperature_k > 100.0) && single_active_condensed && n <= 1e-4
            {
                g3_fallback_used = true;
                for c in 0..nc {
                    if is_active[c] {
                        is_active[c] = false;
                        nj_c[c] = 0.0;
                        j_switch = Some(c); // guide G.3: mark the dropped phase
                    }
                }
                active_ranked.clear();
                j_sol = None;
                j_liq = None;
                times_converged = 0;
                n = 0.1;
                for j in 0..ng {
                    ln_nj[j] = (0.1 / ng as f64).ln();
                }
                iteration = 0;
                continue;
            }
            break; // this epoch's budget is exhausted, no fallback applies — give up.
        }
        iteration += 1;
        iterations_used = total_iteration + 1;

        // Divergence guard (not in the guide): the phase-rule-degenerate
        // family (a lone condensed phase + nearly-pure gas of the same
        // substance) can drift n upward unboundedly instead of failing
        // cleanly — n above ~500 kmol/kg (1 kg of free electrons) is
        // physically impossible, so stop wasting the iteration budget.
        if n > 1e4 {
            cea_debug!("[guard] n={n} diverged — giving up");
            break;
        }

        let active_condensed: Vec<usize> = active_ranked.clone();
        let na = active_condensed.len();
        // H.3/G.2(d): once ions die, the E pseudo-element's row leaves the
        // matrix (num_active_elements shrinks from the end). This is only
        // safe because the dead charged species are simultaneously parked
        // below LOG_MIN and PINNED (dln_nj forced to 0 below) — an unpinned
        // charged species with no charge-balance row and no other element
        // (e-) would grow without bound, and an unparked one re-zeroes the
        // E row into an endless singular-recovery loop.
        let ions_dead = has_electron && !active_ions;
        let active_elements: Vec<usize> = (0..ne)
            .filter(|&k| !element_excluded[k] && !(ions_dead && k == ne - 1))
            .collect();
        let ne_active = active_elements.len();

        // H.1: is_ion[j] = ions requested, still active, and species j
        // actually carries a nonzero electron count (charge). Snapshotted
        // once per iteration (rather than a closure over `active_ions`) so
        // H.3's later flip of `active_ions` doesn't retroactively change
        // truncation decisions already made earlier this same iteration.
        let is_ion: Vec<bool> =
            (0..ng).map(|j| active_ions && stoich[j][ne - 1] != 0.0).collect();

        let mut mixture_thermo: MixtureThermo = MixtureThermo {
            cp: vec![0.0; ns],
            cv: vec![0.0; ns],
            enthalpy: vec![0.0; ns],
            entropy: vec![0.0; ns],
            energy: vec![0.0; ns],
        };
        calc_thermo(
            &species,
            temperature_k,
            ng,
            nc,
            include_condensed,
            &mut mixture_thermo,
        );

        // Pressure/volume closure (Global conventions): fixed for const_p
        // problems, otherwise recomputed each iteration from n/T. n counts
        // gas only (AI.3.4) — condensed never enters this.
        let pressure_bar = if const_p {
            state2
        } else {
            1e-5 * R_CEA * n * temperature_k / state2
        };

        // mu_g[j] (D.1): uses raw ln_nj, not the truncated nj_eff.
        let ln_p_over_n = (pressure_bar / n).ln();
        let mu_g: Vec<f64> = (0..ng)
            .map(|j| {
                mixture_thermo.enthalpy[j] - mixture_thermo.entropy[j] + ln_nj[j] + ln_p_over_n
            })
            .collect();

        // nj_eff[j] (D.4/H.1): truncated/effective gas amount used in every
        // sum. Charged species use the tighter esize threshold.
        let nj_eff: Vec<f64> = (0..ng)
            .map(|j| {
                let threshold = n.ln() - if is_ion[j] { esize } else { tsize };
                if ln_nj[j] > threshold {
                    ln_nj[j].exp()
                } else {
                    0.0
                }
            })
            .collect();

        // b_delta[k] = b0[k] - Σ_j a_kj*n_j, over ALL species — gas
        // (truncated) plus active condensed (F.2's explicit list of where
        // condensed amounts enter).
        let b_delta: Vec<f64> = (0..ne)
            .map(|k| {
                let gas: f64 = (0..ng).map(|j| stoich[j][k] * nj_eff[j]).sum();
                let cond: f64 =
                    active_condensed.iter().map(|&c| stoich[ng + c][k] * nj_c[c]).sum();
                b0[k] - gas - cond
            })
            .collect();
        let n_delta = n - nj_eff.iter().sum::<f64>();

        // "h_j if const P, u_j if const V" (E.2a, E.3) — used wherever the
        // guide's formula depends on P/V rather than on the energy type.
        let energy_weight: Vec<f64> = if const_p {
            mixture_thermo.enthalpy[..ng].to_vec()
        } else {
            mixture_thermo.energy[..ng].to_vec()
        };
        // Same choice, but for one condensed species: F.2's condensed row's
        // own ΔlnT entry (pressure/volume type — not the energy type).
        let condensed_pressure_weight = |sc: usize| -> f64 {
            if const_p { mixture_thermo.enthalpy[sc] } else { mixture_thermo.energy[sc] }
        };
        // F.2's "condensed phase column entry" (h_c/s_c/u_c per the ENERGY
        // row's own type — h· for HP, s· for SP/SV, u· for UV).
        let condensed_energy_type_value = |sc: usize| -> f64 {
            if const_s {
                mixture_thermo.entropy[sc]
            } else if const_p {
                mixture_thermo.enthalpy[sc]
            } else {
                mixture_thermo.energy[sc]
            }
        };

        // E.1/F.2: unknown layout is [pi(ne) | Δn_c(na) | Δln n (iff
        // const_p) | ΔlnT (iff !const_t)]. For TV (const_t && !const_p) with
        // no condensed active, neq == ne: no Δln n and no ΔlnT unknown at
        // all — n is just Σnj after each update (E.1), not solved for here.
        // G.2(e): matrix rows/columns only exist for ACTIVE elements —
        // excluded ones (element_excluded) keep their last pi (pi_prev) but
        // get no equation of their own this solve. b_delta/test4 still
        // check every element regardless (the constraint itself doesn't go
        // away, only its dedicated row in the linearized system).
        let cond_col = |idx: usize| ne_active + idx;
        let lnn_col = ne_active + na; // only valid to read/write when const_p
        let lnt_col = ne_active + na + (const_p as usize); // only valid when !const_t
        let neq = ne_active + na + (const_p as usize) + (!const_t as usize);
        let mut g: DMatrix<f64> = DMatrix::zeros(neq, neq + 1);

        for (row_idx, &k) in active_elements.iter().enumerate() {
            // tmp[j] = a_kj * n_j, reused for every column of this row.
            let tmp: Vec<f64> = (0..ng).map(|j| nj_eff[j] * stoich[j][k]).collect();
            for (col_idx, &i) in active_elements.iter().enumerate() {
                // Σ_j a_kj * a_ij * n_j
                g[(row_idx, col_idx)] = (0..ng).map(|j| tmp[j] * stoich[j][i]).sum::<f64>();
            }
            // F.2: element row's condensed column is the bare stoichiometric
            // coefficient a_kc — Δn_c is additive (not log-space), so unlike
            // gas's Δln n_j substitution there's no n_c weighting here.
            // Symmetric with the condensed row's own pi columns below.
            for (idx, &c) in active_condensed.iter().enumerate() {
                let a_kc = stoich[ng + c][k];
                g[(row_idx, cond_col(idx))] = a_kc;
                g[(cond_col(idx), row_idx)] = a_kc;
            }
            if const_p {
                g[(row_idx, lnn_col)] = tmp.iter().sum::<f64>(); // Δln n column
                g[(lnn_col, row_idx)] = g[(row_idx, lnn_col)]; // symmetric fill of the moles row
            }
            if !const_t {
                // E.2a: element row's ΔlnT column = Σ_j a_kj * n_j * (h_j or
                // u_j). Depends only on const_p/const_v, not on the energy
                // type (const_s) — NOT generally equal to the energy row's
                // own pi-column entry below, which uses the type-specific
                // tmp_energy. They happen to coincide for HP/UV (tmp_energy
                // == n_j*energy_weight there) but not for SP/SV, so this is
                // filled on its own, not mirrored.
                g[(row_idx, lnt_col)] = (0..ng).map(|j| tmp[j] * energy_weight[j]).sum();
            }
            g[(row_idx, neq)] = b_delta[k] + (0..ng).map(|j| tmp[j] * mu_g[j]).sum::<f64>();
        }

        // F.2: condensed rows. No diagonal, no Δln n entry — n counts gas
        // only (AI.3.4), and Δn_c isn't coupled to itself.
        for (idx, &c) in active_condensed.iter().enumerate() {
            let sc = ng + c;
            let row = cond_col(idx);
            if !const_t {
                g[(row, lnt_col)] = condensed_pressure_weight(sc);
            }
            // mu_c/RT = h_c - s_c (D.1: condensed has no ln(n_c/n) or ln(P)
            // mixing term).
            g[(row, neq)] = mixture_thermo.enthalpy[sc] - mixture_thermo.entropy[sc];
        }

        // Moles row (Eq. 2.26); pi coefficients already filled symmetrically
        // above. Absent entirely for const_v — n isn't a Newton unknown
        // there, and (AI.3.4) never involves condensed even when const_p.
        if const_p {
            g[(lnn_col, lnn_col)] = -n_delta;
            g[(lnn_col, neq)] = n_delta + (0..ng).map(|j| nj_eff[j] * mu_g[j]).sum::<f64>();
        }
        if !const_t {
            // E.2b: moles row's ΔlnT column = Σ_j n_j * h_j (const_p only —
            // there's no moles row at all otherwise, so always enthalpy,
            // never energy_weight, and gas-only per AI.3.4).
            if const_p {
                g[(lnn_col, lnt_col)] =
                    (0..ng).map(|j| nj_eff[j] * mixture_thermo.enthalpy[j]).sum::<f64>();
            }

            // E.2c energy row. tmp_j/hsu_delta depend on the energy type:
            //   HP ('h'): tmp_j = n_j*h_j            (energy_weight = h_j since const_p)
            //   UV ('u'): tmp_j = n_j*u_j             (energy_weight = u_j since !const_p)
            //   SP/SV ('s'): tmp_j = n_j*(h_j - mu_j) (always enthalpy, per the guide's table)
            // hsu_delta = state1[/T] - Σ_gas tmp_j - Σ_cond n_c*(type value)
            // — the SP/SV gas residual (guide's Σ n_j(s_j - ln n_j -
            // ln(P/n))) simplifies to the same Σ tmp_j because that bracket
            // equals h_j - mu_g[j]; condensed instead uses raw s_c directly
            // (it has no ln(n_c/n) mixing term to cancel against).
            let tmp_energy: Vec<f64> = if const_s {
                (0..ng)
                    .map(|j| nj_eff[j] * (mixture_thermo.enthalpy[j] - mu_g[j]))
                    .collect()
            } else {
                (0..ng).map(|j| nj_eff[j] * energy_weight[j]).collect()
            };
            let tmp_energy_condensed: Vec<f64> = active_condensed
                .iter()
                .map(|&c| nj_c[c] * condensed_energy_type_value(ng + c))
                .collect();
            let hsu_delta = (if const_s { state1 } else { state1 / temperature_k })
                - tmp_energy.iter().sum::<f64>()
                - tmp_energy_condensed.iter().sum::<f64>();

            // Energy row's own pi columns and Δln n column — type-specific,
            // computed independently rather than mirrored (see E.2a note).
            // Per F.2's explicit list of where condensed amounts enter, this
            // stays gas-only (Δn_c gets its own dedicated column instead,
            // filled just below). SV extra term (Eq. 2.48): !const_p &&
            // const_s subtracts the plain element-row sum on top of the
            // tmp_energy-based one — also gas-only in the guide's wording.
            for (col_idx, &i) in active_elements.iter().enumerate() {
                g[(lnt_col, col_idx)] = (0..ng).map(|j| tmp_energy[j] * stoich[j][i]).sum::<f64>();
                if !const_p && const_s {
                    g[(lnt_col, col_idx)] -= (0..ng).map(|j| nj_eff[j] * stoich[j][i]).sum::<f64>();
                }
            }
            // F.2: energy row's own condensed columns — bare per-mole value
            // (like the element row's a_kc above), not weighted by nj_c.
            for (idx, &c) in active_condensed.iter().enumerate() {
                g[(lnt_col, cond_col(idx))] = condensed_energy_type_value(ng + c);
            }
            if const_p {
                g[(lnt_col, lnn_col)] = tmp_energy.iter().sum::<f64>();
            }

            g[(lnt_col, lnt_col)] = if const_p {
                // F.2: "the ΔlnT diagonal (Σ n_c cp_c)" — condensed enters
                // via the nj*cp term only, not the tmp*h cross term (gas-only
                // there, per the guide's explicit condensed-entry list).
                let cond_cp: f64 =
                    active_condensed.iter().map(|&c| nj_c[c] * mixture_thermo.cp[ng + c]).sum();
                (0..ng).map(|j| nj_eff[j] * mixture_thermo.cp[j]).sum::<f64>()
                    + cond_cp
                    + (0..ng).map(|j| tmp_energy[j] * mixture_thermo.enthalpy[j]).sum::<f64>()
            } else {
                let cond_cv: f64 =
                    active_condensed.iter().map(|&c| nj_c[c] * mixture_thermo.cv[ng + c]).sum();
                let mut v = (0..ng).map(|j| nj_eff[j] * mixture_thermo.cv[j]).sum::<f64>()
                    + cond_cv
                    + (0..ng).map(|j| tmp_energy[j] * mixture_thermo.energy[j]).sum::<f64>();
                if const_s {
                    // SV extra term (Eq. 2.48), gas-only per the guide.
                    v -= (0..ng).map(|j| nj_eff[j] * mixture_thermo.energy[j]).sum::<f64>();
                }
                v
            };
            // const_s's extra RHS term differs by const_p: SP adds n_delta,
            // SV subtracts Σ nj_eff*mu_g instead (Eq. 2.28 vs 2.48).
            let extra = if const_s {
                if const_p {
                    n_delta
                } else {
                    -(0..ng).map(|j| nj_eff[j] * mu_g[j]).sum::<f64>()
                }
            } else {
                0.0
            };
            g[(lnt_col, neq)] =
                hsu_delta + (0..ng).map(|j| tmp_energy[j] * mu_g[j]).sum::<f64>() + extra;
        }

        // Snapshot the assembled (pre-solve) matrix from the very first pass
        // only — that's the canonical-guess matrix Section D validation
        // checks against. gauss() overwrites g in place with the solution.
        if total_iteration == 0 {
            assembled_first_iteration = Some(g.clone());
        }

        // Solve for [pi_0..pi_{ne-1}, Δn_c, dln_n, (dlnT)], then recover
        // dln_nj (Eq. 2.18, gaining the ΔlnT*h_j term for variable-T
        // problems, E.3).
        if let Err(ierr) = gauss(&mut g) {
            // G.2 correct_singular(ierr): "first action always" is to
            // disable truncation, then try recovery branches in order,
            // first match wins.
            times_singular += 1;
            cea_debug!(
                "[G.2] singular ierr={ierr} times={times_singular} iter={iteration} ne_active={ne_active} na={na}"
            );
            if times_singular > MAX_SINGULAR_RETRIES {
                break; // G.1: give up: converged stays false.
            }
            // G.2's "first action always": disable truncation entirely —
            // all three thresholds, so charged species un-truncate too, and
            // xsize so a later convergence's tsize = xsize keeps it off.
            tsize = 80.0;
            xsize = 80.0;
            esize = 80.0;
            let mut changed = false;

            // (a) freshly-added condensed conflicts: a singular matrix on
            // the very first attempt after F.3.3 activated a new condensed
            // species, hitting a condensed row, usually means the new
            // species and an existing one are fighting over some shared
            // element — remove whichever OTHER active condensed species
            // shares an element with the freshly-added one and has the
            // smallest amount, rather than trusting `ierr`'s own index.
            if !changed
                && ierr >= ne_active
                && iteration <= 1
                && na > 1
                && let Some(last) = last_cond_idx
                && is_active.get(last).copied().unwrap_or(false)
            {
                let last_sc = ng + last;
                let worst = active_condensed
                    .iter()
                    .filter(|&&c| {
                        c != last && (0..ne).any(|k| stoich[ng + c][k] != 0.0 && stoich[last_sc][k] != 0.0)
                    })
                    .min_by(|&&a, &&b| nj_c[a].partial_cmp(&nj_c[b]).unwrap());
                if let Some(&c) = worst {
                    is_active[c] = false;
                    nj_c[c] = 0.0;
                    active_ranked.retain(|&x| x != c);
                    j_switch = Some(c);
                    singular_blocklist = Some(c);
                    if j_sol == Some(c) || j_liq == Some(c) {
                        j_sol = None;
                        j_liq = None;
                    }
                    changed = true;
                }
            }

            // (b) element-row singularity: remove the smallest-nj active
            // condensed species that actually touches this element (same
            // j_switch/blocklist bookkeeping as (a)).
            if !changed && ierr < ne_active && na > 0 {
                let k = active_elements[ierr];
                let worst = active_condensed
                    .iter()
                    .filter(|&&c| stoich[ng + c][k].abs() > 1e-8)
                    .min_by(|&&a, &&b| nj_c[a].partial_cmp(&nj_c[b]).unwrap());
                if let Some(&c) = worst {
                    is_active[c] = false;
                    nj_c[c] = 0.0;
                    active_ranked.retain(|&x| x != c);
                    j_switch = Some(c);
                    singular_blocklist = Some(c);
                    if j_sol == Some(c) || j_liq == Some(c) {
                        j_sol = None;
                        j_liq = None;
                    }
                    changed = true;
                }
            }

            // (c) condensed-row singularity: deactivate the condensed
            // species at that active RANK (active_condensed is rank-ordered,
            // so this index maps exactly like the guide's G.2c).
            if !changed && ierr >= ne_active && ierr < ne_active + na {
                let c = active_condensed[ierr - ne_active];
                is_active[c] = false;
                nj_c[c] = 0.0;
                active_ranked.retain(|&x| x != c);
                if j_sol == Some(c) || j_liq == Some(c) {
                    j_sol = None;
                    j_liq = None;
                }
                changed = true;
            }

            // (d) electron-row fallback: the E (electron) element row
            // itself is degenerate — the ion channel isn't solvable this
            // way for this state. Zero out every charged gas species (down
            // to SMNOL, not literally -inf, so they stay Newton-updatable)
            // and permanently disable ions for the rest of this solve.
            if !changed
                && ierr < ne_active
                && active_ions
                && has_electron
                && active_elements[ierr] == ne - 1
            {
                for j in 0..ng {
                    if stoich[j][ne - 1] != 0.0 {
                        ln_nj[j] = -87.3; // below LOG_MIN — see the H.3 parking note
                    }
                }
                active_ions = false;
                pi_e = 0.0;
                pi_prev[ne - 1] = 0.0;
                changed = true;
            }

            // (e) element reduction: this element's balance is linearly
            // dependent on the others (e.g. present in too few surviving
            // species to pin its own multiplier) — drop its row/column for
            // the rest of this solve. Gated to only the first attempt of a
            // fresh epoch, or after several other attempts have already
            // failed (guide's "iter < 1 or times_singular >= 4") — this
            // matters: firing it eagerly can permanently drop a real
            // constraint (nothing left to pin, say, total oxygen) when the
            // real fix was just to reseed a decayed species (f) instead.
            // (guide's extra gate: never element-reduce while ions are live —
            // the E row's degeneracy is (d)'s job, not a reason to drop a
            // real element constraint.)
            if !changed
                && ierr < ne_active
                && ne_active > 1
                && (iteration <= 1 || times_singular >= 4)
                && !(include_ions && active_ions)
            {
                let k = active_elements[ierr];
                cea_debug!("[G.2e] excluding element {:?}", constituent_pool[k]);
                element_excluded[k] = true;
                pi_prev[k] = 0.0;
                changed = true;
            }

            // (f) trace reseed: last resort. Re-inflate any gas species
            // that's decayed to negligible amounts, in case truncation
            // zeroed it into a singular sub-matrix.
            if !changed {
                for j in 0..ng {
                    // Dead-ion parking is exempt: re-inflating a pinned
                    // charged species would put phantom element amounts
                    // back into the balance sums (see the H.3 parking note).
                    if has_electron && !active_ions && stoich[j][ne - 1] != 0.0 {
                        continue;
                    }
                    if ln_nj[j] < -80.0 {
                        ln_nj[j] = -13.815511; // guide's SMNOL = ln(1e-6)
                    }
                }
            }

            iteration = 0; // restart this epoch's iteration budget
            continue;
        }
        // H.3 (CEA: get_solution_vars, i.e. AFTER the linear solve): if ions
        // are active but every charged species has decayed to nothing, ions
        // turn off permanently for this solve and the charged species get
        // parked at SMNOL. This must NOT run before assembly — an all-zero E
        // row must reach `gauss` while active_ions is still true, so the
        // singularity routes to G.2(d)'s ion shutdown rather than G.2(e)
        // wrongly excluding the E element (which would leave the parked
        // species with no charge balance pinning them, free to run away).
        if active_ions && (0..ng).all(|j| !(stoich[j][ne - 1] != 0.0 && nj_eff[j] > 0.0)) {
            cea_debug!("[H.3] all charged species decayed — ions off");
            active_ions = false;
            pi_prev[ne - 1] = 0.0;
            for j in 0..ng {
                if stoich[j][ne - 1] != 0.0 {
                    // Deliberate deviation from the guide's SMNOL (1e-6):
                    // park BELOW the LOG_MIN reporting cutoff. At 1e-6 a
                    // parked molecular ion still carries ~1e-6 of its real
                    // elements into the balance sums — above test 4's
                    // 1e-6·max(b0) tolerance — so convergence could never
                    // pass; and once it decays under the truncation
                    // threshold it re-zeroes the E row into an endless
                    // singular loop. Parked+pinned at 1e-38 it is inert.
                    ln_nj[j] = -87.3;
                }
            }
        }
        // pi_prev doubles as "current pi" here (see its declaration) — only
        // active positions get overwritten; excluded ones (G.2e) keep
        // whatever they were (zeroed at the moment they got excluded).
        for (row_idx, &k) in active_elements.iter().enumerate() {
            pi_prev[k] = g[(row_idx, neq)];
        }
        let dnj_c: Vec<f64> = (0..na).map(|idx| g[(cond_col(idx), neq)]).collect();
        // TV: no Δln n unknown at all — dln_n stays 0 (E.1).
        let dln_n = if const_p { g[(lnn_col, neq)] } else { 0.0 };
        let dln_t = if const_t { 0.0 } else { g[(lnt_col, neq)] };
        // H.3: the pi_e feedback fires only on the Newton sweep immediately
        // after a check that converged on the base tests but whose electron
        // polish did NOT settle — then pi_e is the E-element multiplier from
        // THIS solve, layered on top of the regular pi sum so the neutral
        // species re-equilibrate against the shifted ions. Otherwise 0.
        pi_e = if active_ions && has_electron && base_converged_prev && !ions_converged {
            pi_prev[ne - 1]
        } else {
            0.0
        };
        let dln_nj: Vec<f64> = (0..ng)
            .map(|j| {
                // Live active_ions, not the top-of-iteration ions_dead: the
                // death check just above may have parked them this same
                // iteration, and they must not take one last unpinned step.
                if has_electron && !active_ions && stoich[j][ne - 1] != 0.0 {
                    return 0.0; // pinned: see the active_elements note above
                }
                -mu_g[j]
                    + dln_n
                    + (0..ne).map(|i| stoich[j][i] * pi_prev[i]).sum::<f64>()
                    + dln_t * energy_weight[j] // h_j if const_p, u_j if const_v (E.3)
                    + stoich[j][ne - 1] * pi_e
            })
            .collect();

        // D.2: cap the step before it's applied. (Condensed Δn_c isn't part
        // of this — the guide's damping formula only ever mentions gas
        // Δln n_j, ΔlnT, and Δln n.)
        let lambda = damping_lambda(&ln_nj[..ng], &dln_nj, dln_n, dln_t, n, tsize);

        // D.3: apply the damped update. nj_stored's truncation threshold
        // uses the OLD n (not yet updated below) — matches the guide's
        // exact pseudocode order: nj[j] update happens inside the same loop
        // that bumps ln_nj[j], before n itself is reassigned.
        let applied_dln_nj: Vec<f64> = (0..ng).map(|j| lambda * dln_nj[j]).collect();
        for j in 0..ng {
            ln_nj[j] += applied_dln_nj[j];
            let threshold = n.ln() - if is_ion[j] { esize } else { tsize };
            nj_stored[j] = if ln_nj[j] > threshold { ln_nj[j].exp() } else { 0.0 };
        }
        // F.2: condensed update is LINEAR, not log-space, and can go
        // negative transiently — that's F.3.1's removal signal.
        let applied_dnj_c: Vec<f64> = (0..na).map(|idx| lambda * dnj_c[idx]).collect();
        for (idx, &c) in active_condensed.iter().enumerate() {
            nj_c[c] += applied_dnj_c[idx];
        }
        let applied_dln_n = lambda * dln_n;
        n = if const_p {
            (n.ln() + applied_dln_n).exp()
        } else {
            // TV: n isn't a Newton unknown — it's just Σ nj after the
            // species update (E.1), using the freshly truncated amounts.
            nj_stored.iter().sum()
        };
        let applied_dln_t = lambda * dln_t;
        if !const_t {
            temperature_k = (temperature_k.ln() + applied_dln_t).exp();
        }

        // D.3 (cont.): the guide's update_solution recomputes thermo at the
        // new T before returning — needed here so test 6 (entropy, SP) and
        // F.3 (condensed adjustment) read values at the state actually being
        // checked, not the pre-update one still sitting in `mixture_thermo`.
        let mut mixture_thermo_post = MixtureThermo {
            cp: vec![0.0; ns],
            cv: vec![0.0; ns],
            enthalpy: vec![0.0; ns],
            entropy: vec![0.0; ns],
            energy: vec![0.0; ns],
        };
        calc_thermo(
            &species,
            temperature_k,
            ng,
            nc,
            include_condensed,
            &mut mixture_thermo_post,
        );

        // D.5: convergence check. Tests that apply: condensed (test 2, only
        // over active species), total moles (test 3, only when const_p —
        // for TV, n is a direct sum, not a Newton unknown, so there's
        // nothing to check separately from test 1), temperature (test 5,
        // only when !const_t), entropy (test 6, only when const_s), no
        // trace/ions (test 7).
        let sum_nj: f64 =
            nj_stored.iter().sum::<f64>() + active_condensed.iter().map(|&c| nj_c[c]).sum::<f64>();

        // Test 1: gas species.
        let test1 = (0..ng).all(|j| nj_stored[j] * applied_dln_nj[j].abs() / sum_nj <= NJ_TOL);

        // Test 2: condensed species (active only).
        let test2 = (0..na).all(|idx| applied_dnj_c[idx].abs() / sum_nj <= NJ_TOL);

        // Test 3: total moles (const_p only). Denominator is Σ over GAS
        // species only, unlike tests 1/2's gas+condensed sum_nj.
        let gas_sum_nj: f64 = nj_stored.iter().sum();
        let test3 = !const_p || n * applied_dln_n.abs() / gas_sum_nj <= NJ_TOL;

        // Test 4: elements, using RAW nj (not the stored/truncated values)
        // for gas, and the linear nj_c for condensed.
        let b_max = b0.iter().cloned().fold(0.0_f64, f64::max);
        let test4 = (0..ne).all(|k| {
            if b0[k] <= B_TOL {
                return true;
            }
            let gas_b_k: f64 = (0..ng).map(|j| stoich[j][k] * ln_nj[j].exp()).sum();
            let cond_b_k: f64 =
                active_condensed.iter().map(|&c| stoich[ng + c][k] * nj_c[c]).sum();
            (b0[k] - gas_b_k - cond_b_k).abs() <= B_TOL * b_max
        });

        // Test 5: temperature (only for variable-T problems).
        let test5 = const_t || applied_dln_t.abs() <= T_TOL;

        // Test 6: entropy (SP/SV only). s_delta = state1 - Σ_gas nj*(s_j -
        // ln nj - ln(P/n)) - Σ_cond n_c*s_c, evaluated at the post-update
        // state (RAW ln_nj, like test 4 — the ln(nj) term isn't meaningful
        // for a truncated 0).
        let test6 = if const_s {
            let new_ln_p_over_n = (pressure_bar / n).ln();
            let gas_s: f64 = (0..ng)
                .map(|j| {
                    let nj_raw = ln_nj[j].exp();
                    nj_raw * (mixture_thermo_post.entropy[j] - ln_nj[j] - new_ln_p_over_n)
                })
                .sum();
            let cond_s: f64 = active_condensed
                .iter()
                .map(|&c| nj_c[c] * mixture_thermo_post.entropy[ng + c])
                .sum();
            let s_delta = state1 - gas_s - cond_s;
            s_delta.abs() <= S_TOL
        } else {
            true
        };

        let base_tests_pass = test1 && test2 && test3 && test4 && test5 && test6;
        base_converged_prev = base_tests_pass;

        if base_tests_pass {
            // D.5's on-success actions: loosen truncation to xsize (this is
            // a persistent mutation — how CEA resolves trace species after
            // the main composition settles) and recompute the effective
            // amounts at the looser threshold.
            times_converged += 1;
            tsize = xsize;
            for j in 0..ng {
                let threshold = n.ln() - if is_ion[j] { esize } else { tsize };
                nj_stored[j] = if ln_nj[j] > threshold { ln_nj[j].exp() } else { 0.0 };
            }
            cea_debug!(
                "[D.5] base converged (times={times_converged}) T={temperature_k} n={n} iter={iteration}"
            );
            // F.3.4/G.1: repeated re-convergence means the condensed set is
            // thrashing — give up (converged stays false).
            if times_converged > 3 * ne_active.max(1) {
                cea_debug!("[F.3.4] condensed-set thrash — giving up");
                break;
            }
        }

        // H.2: post-convergence electron-balance polish. Runs only once the
        // other six tests already pass — the tiny electron concentration
        // makes the main matrix's own pi_E ill-conditioned, so this does a
        // separate, well-conditioned 1D Newton iteration purely on the
        // charge-neutrality residual (Σ a_Ej*n_j == 0) instead of trusting
        // that channel of the big solve.
        let mut ion_gate_pass = true;
        if base_tests_pass && active_ions && has_electron {
            let e_idx = ne - 1;
            let mut polished = false;
            for _ in 0..80 {
                let mut sum1 = 0.0_f64;
                let mut sum2 = 0.0_f64;
                for j in 0..ng {
                    let a_e = stoich[j][e_idx];
                    if a_e == 0.0 {
                        continue;
                    }
                    let raw = if ln_nj[j] > -87.0 { ln_nj[j].exp() } else { 0.0 };
                    nj_stored[j] = if ln_nj[j] > n.ln() - esize { raw } else { 0.0 };
                    sum1 += a_e * raw; // RAW amount in the sums, not truncated
                    sum2 += a_e * a_e * raw;
                }
                let mut dpi_e = 0.0;
                if sum2 != 0.0 {
                    dpi_e = -sum1 / sum2;
                    for j in 0..ng {
                        let a_e = stoich[j][e_idx];
                        if a_e != 0.0 {
                            ln_nj[j] += a_e * dpi_e;
                        }
                    }
                }
                if dpi_e.abs() > 1e-4 {
                    pi_e += dpi_e; // keep polishing
                } else {
                    polished = true; // small enough — done, pi_e keeps its
                    break; // accumulated value (not this final increment).
                }
            }
            ions_converged = polished;
            if !polished {
                pi_e = 0.0; // 80 sub-iterations exhausted: reset and retry.
                ion_gate_pass = false;
            }
        }

        if !(base_tests_pass && ion_gate_pass) {
            continue;
        }

        // F.3: test_condensed — runs only after the inner Newton loop
        // converges, and may still flip things back to "keep going" by
        // changing the active set. Priority order per the guide: F.3.1
        // remove negatives, F.3.2/F.4 phase validity, F.3.3 insert; one
        // change per call, and any change restarts the epoch budget.
        let mut changed = false;

        // F.3.1: remove the active condensed species with the most negative
        // nj_c * cp_c (only one removal per call).
        if nc > 0 {
            let mut worst: Option<(usize, f64)> = None;
            for &c in &active_condensed {
                let val = nj_c[c] * mixture_thermo_post.cp[ng + c];
                if val <= 0.0 && worst.is_none_or(|(_, w)| val < w) {
                    worst = Some((c, val));
                }
            }
            if let Some((c, _)) = worst {
                cea_debug!("[F.3.1] removing {} (nj={})", species[ng + c].data().symbol(), nj_c[c]);
                is_active[c] = false;
                nj_c[c] = 0.0;
                active_ranked.retain(|&x| x != c);
                if j_sol == Some(c) || j_liq == Some(c) {
                    j_sol = None;
                    j_liq = None;
                }
                changed = true;
            }
        }

        // F.3.2/F.4: phase-validity and phase-transition check.
        if !changed && nc > 0 {
            changed = check_condensed_phases(
                &species,
                ng,
                const_t,
                &mut temperature_k,
                &mut nj_c,
                &mut is_active,
                &mut active_ranked,
                &mut j_switch,
                &mut j_sol,
                &mut j_liq,
            );
        }

        // F.3.3: insertion test — only if nothing was removed this pass.
        // For every inactive species whose fit range allows the current T,
        // compute the per-mass Gibbs decrease from forming it; activate the
        // most negative one (only one insertion per call).
        if !changed && nc > 0 {
            let mut best: Option<(usize, f64)> = None;
            let mut blocked_wants_in = false;
            for c in 0..nc {
                if is_active[c] {
                    continue;
                }
                let sc = ng + c;
                let (t_low, t_high) = species[sc].data().valid_temperature_range();
                let allowed = (temperature_k >= t_low || t_low == 200.0) && temperature_k <= t_high;
                if !allowed {
                    continue;
                }
                let a_pi: f64 = (0..ne).map(|i| stoich[sc][i] * pi_prev[i]).sum();
                let delta_g = (mixture_thermo_post.enthalpy[sc]
                    - mixture_thermo_post.entropy[sc]
                    - a_pi)
                    / species[sc].data().mw();
                if Some(c) == singular_blocklist {
                    // G.2 removed this species for singularizing the matrix;
                    // never re-add it. Wanting back in is CEA's hard-error
                    // "re-insertion likely to cause singular matrix".
                    if delta_g < -1e-12 {
                        blocked_wants_in = true;
                    }
                    continue;
                }
                if delta_g < 0.0 && best.is_none_or(|(_, b)| delta_g < b) {
                    best = Some((c, delta_g));
                }
            }
            match best {
                Some((c, delta_g)) if delta_g.abs() > 1e-12 => {
                    cea_debug!(
                        "[F.3.3] inserting {} (delta_g={delta_g}) T={temperature_k}",
                        species[ng + c].data().symbol()
                    );
                    // F.1: activate at rank 1 ("front"); nj stays 0 and the
                    // first Newton step sizes it.
                    is_active[c] = true;
                    nj_c[c] = 0.0;
                    active_ranked.insert(0, c);
                    last_cond_idx = Some(c); // G.2(a): remember for the next singular-matrix check.
                    changed = true;
                }
                None if blocked_wants_in => {
                    // Nothing else can be inserted, and the equilibrium
                    // genuinely wants the blocklisted species: unsolvable as
                    // configured (CEA aborts here) — return unconverged.
                    break;
                }
                _ => {}
            }
        }

        if changed {
            iteration = 0; // F.3: restart this epoch's iteration budget.
            ions_converged = false; // the changed set needs a fresh polish.
            continue;
        }

        converged = true;
        break;
    }

    // G.1: final sanity bounds — a "converged" T outside CEA's physical
    // range means the iterate ran away, not that the problem was solved.
    if converged && !(160.0..=22000.0).contains(&temperature_k) {
        converged = false;
    }

    // Section I: post-processing & thermodynamic derivatives. Runs once, at
    // the final state, independent of the main loop's per-iteration scoping.
    let active_condensed: Vec<usize> = (0..nc).filter(|&c| is_active[c]).collect();
    let na = active_condensed.len();
    let active_elements: Vec<usize> = (0..ne).filter(|&k| !element_excluded[k]).collect();
    let ne_active = active_elements.len();

    let mut final_thermo = MixtureThermo {
        cp: vec![0.0; ns],
        cv: vec![0.0; ns],
        enthalpy: vec![0.0; ns],
        entropy: vec![0.0; ns],
        energy: vec![0.0; ns],
    };
    calc_thermo(&species, temperature_k, ng, nc, true, &mut final_thermo);

    let final_pressure_bar =
        if const_p { state2 } else { 1e-5 * R_CEA * n * temperature_k / state2 };
    let final_volume_m3_per_kg = if const_p {
        R_CEA * n * temperature_k / (final_pressure_bar * 1e5)
    } else {
        state2
    };

    // I.1: reported mixture properties. Uses the LOG_MIN-mapped ("reporting")
    // gas amounts (nj_stored, D.4-truncated), except entropy, which the
    // guide explicitly calls out as using the raw log instead (below).
    let sum_all: f64 =
        nj_stored.iter().sum::<f64>() + active_condensed.iter().map(|&c| nj_c[c]).sum::<f64>();
    let mut mole_fractions = vec![0.0; ns];
    let mut mass_fractions = vec![0.0; ns];
    let mut total_mass_flow: f64 =
        (0..ng).map(|j| nj_stored[j] * species[j].data().mw()).sum();
    total_mass_flow += active_condensed.iter().map(|&c| nj_c[c] * species[ng + c].data().mw()).sum::<f64>();
    for j in 0..ng {
        mole_fractions[j] = if sum_all > 0.0 { nj_stored[j] / sum_all } else { 0.0 };
        mass_fractions[j] =
            if total_mass_flow > 0.0 { nj_stored[j] * species[j].data().mw() / total_mass_flow } else { 0.0 };
    }
    for c in 0..nc {
        mole_fractions[ng + c] = if sum_all > 0.0 { nj_c[c] / sum_all } else { 0.0 };
        mass_fractions[ng + c] = if total_mass_flow > 0.0 {
            nj_c[c] * species[ng + c].data().mw() / total_mass_flow
        } else {
            0.0
        };
    }

    let enthalpy_kj_per_kg = ((0..ng).map(|j| nj_stored[j] * final_thermo.enthalpy[j]).sum::<f64>()
        + active_condensed.iter().map(|&c| nj_c[c] * final_thermo.enthalpy[ng + c]).sum::<f64>())
        * R_CEA
        * temperature_k
        / 1000.0;
    let energy_kj_per_kg = enthalpy_kj_per_kg - n * R_CEA * temperature_k / 1000.0;

    // S uses RAW ln_nj (guide's I.1 note): any species that survived the
    // D.4 truncation map (nj_stored[j] > 0) contributes its true log.
    let ln_p_over_n_final = (final_pressure_bar / n).ln();
    let gas_s: f64 = (0..ng)
        .filter(|&j| nj_stored[j] > 0.0)
        .map(|j| nj_stored[j] * (final_thermo.entropy[j] - ln_nj[j] - ln_p_over_n_final))
        .sum();
    let cond_s: f64 =
        active_condensed.iter().map(|&c| nj_c[c] * final_thermo.entropy[ng + c]).sum();
    let entropy_kj_per_kg_k = R_CEA / 1000.0 * (gas_s + cond_s);
    let gibbs_kj_per_kg = enthalpy_kj_per_kg - temperature_k * entropy_kj_per_kg_k;

    let molar_mass_kg_per_kmol = 1.0 / n;
    let active_cond_mole_frac: f64 = active_condensed.iter().map(|&c| mole_fractions[ng + c]).sum();
    let mean_molecular_weight_kg_per_kmol = (1.0 - active_cond_mole_frac) / n;

    let cp_frozen_kj_per_kg_k = ((0..ng).map(|j| nj_stored[j] * final_thermo.cp[j]).sum::<f64>()
        + active_condensed.iter().map(|&c| nj_c[c] * final_thermo.cp[ng + c]).sum::<f64>())
        * R_CEA
        / 1000.0;
    let cv_frozen_kj_per_kg_k = cp_frozen_kj_per_kg_k - n * R_CEA / 1000.0;

    // I.2: equilibrium partials — two linear systems sharing the converged
    // matrix's LHS (element rows + condensed rows + moles row, always
    // ne_active+na+1 equations regardless of the problem's const_p/const_t
    // flags — this analysis always differentiates w.r.t. lnT at fixed P and
    // w.r.t. lnP at fixed T) built from UNTRUNCATED current nj (raw
    // exp(ln_nj), not nj_stored).
    let nj_raw: Vec<f64> = (0..ng).map(|j| ln_nj[j].exp()).collect();
    let neq2 = ne_active + na + 1;
    let cond_col2 = |idx: usize| ne_active + idx;
    let lnn_col2 = ne_active + na;

    let mut lhs = DMatrix::<f64>::zeros(neq2, neq2);
    for (row_idx, &k) in active_elements.iter().enumerate() {
        let tmp: Vec<f64> = (0..ng).map(|j| nj_raw[j] * stoich[j][k]).collect();
        for (col_idx, &i) in active_elements.iter().enumerate() {
            lhs[(row_idx, col_idx)] = (0..ng).map(|j| tmp[j] * stoich[j][i]).sum::<f64>();
        }
        for (idx, &c) in active_condensed.iter().enumerate() {
            let a_kc = stoich[ng + c][k];
            lhs[(row_idx, cond_col2(idx))] = a_kc;
            lhs[(cond_col2(idx), row_idx)] = a_kc;
        }
        lhs[(row_idx, lnn_col2)] = tmp.iter().sum::<f64>();
        lhs[(lnn_col2, row_idx)] = lhs[(row_idx, lnn_col2)];
    }
    let n_delta_raw = n - nj_raw.iter().sum::<f64>();
    lhs[(lnn_col2, lnn_col2)] = -n_delta_raw;

    let mut rhs_t = vec![0.0; neq2];
    for (row_idx, &k) in active_elements.iter().enumerate() {
        rhs_t[row_idx] =
            -(0..ng).map(|j| stoich[j][k] * nj_raw[j] * final_thermo.enthalpy[j]).sum::<f64>();
    }
    for (idx, &c) in active_condensed.iter().enumerate() {
        rhs_t[cond_col2(idx)] = -final_thermo.enthalpy[ng + c];
    }
    let sum_nj_hj_raw: f64 = (0..ng).map(|j| nj_raw[j] * final_thermo.enthalpy[j]).sum();
    rhs_t[lnn_col2] = -sum_nj_hj_raw;

    let mut rhs_p = vec![0.0; neq2];
    for (row_idx, &k) in active_elements.iter().enumerate() {
        rhs_p[row_idx] = (0..ng).map(|j| stoich[j][k] * nj_raw[j]).sum::<f64>();
    }
    rhs_p[lnn_col2] = nj_raw.iter().sum::<f64>();

    let solve_partial = |rhs: &[f64]| -> Option<DVector<f64>> {
        if neq2 == 0 {
            return Some(DVector::zeros(0));
        }
        let mut g = DMatrix::<f64>::zeros(neq2, neq2 + 1);
        g.view_mut((0, 0), (neq2, neq2)).copy_from(&lhs);
        for i in 0..neq2 {
            g[(i, neq2)] = rhs[i];
        }
        gauss(&mut g).ok()?;
        Some(DVector::from_fn(neq2, |i, _| g[(i, neq2)]))
    };

    let sol_t = solve_partial(&rhs_t);
    let sol_p = solve_partial(&rhs_p);

    let (dlnv_dlnt, dlnv_dlnp, cp_eq_over_r) = match (&sol_t, &sol_p) {
        (Some(st), Some(sp)) => {
            let dlnn_dlnt = st[lnn_col2];
            let dlnn_dlnp = sp[lnn_col2];
            let dlnv_dlnt = 1.0 + dlnn_dlnt;
            let dlnv_dlnp = -1.0 + dlnn_dlnp;

            let term1: f64 = active_elements
                .iter()
                .enumerate()
                .map(|(row_idx, &k)| {
                    let sum_akj_nj_hj: f64 = (0..ng)
                        .map(|j| stoich[j][k] * nj_raw[j] * final_thermo.enthalpy[j])
                        .sum();
                    sum_akj_nj_hj * st[row_idx]
                })
                .sum();
            let term2: f64 = active_condensed
                .iter()
                .enumerate()
                .map(|(idx, &c)| final_thermo.enthalpy[ng + c] * st[cond_col2(idx)])
                .sum();
            let term3 = sum_nj_hj_raw * dlnn_dlnt;
            let term4 = (0..ng).map(|j| nj_raw[j] * final_thermo.cp[j]).sum::<f64>()
                + active_condensed.iter().map(|&c| nj_c[c] * final_thermo.cp[ng + c]).sum::<f64>();
            let term5 =
                (0..ng).map(|j| nj_raw[j] * final_thermo.enthalpy[j].powi(2)).sum::<f64>();

            (dlnv_dlnt, dlnv_dlnp, term1 + term2 + term3 + term4 + term5)
        }
        // I.3's degrade-gracefully fallback ("never a hard failure"): a
        // singular partials system falls back to the ideal-gas-like
        // constants plus I.1's own frozen-cp fallback formula.
        _ => (1.0, -1.0, cp_frozen_kj_per_kg_k * 1000.0 / R_CEA),
    };

    let gamma_s = -1.0 / (dlnv_dlnp + dlnv_dlnt * dlnv_dlnt * n / cp_eq_over_r);
    let cv_eq_over_r = -cp_eq_over_r / (gamma_s * dlnv_dlnp);
    let cp_eq_kj_per_kg_k = cp_eq_over_r * R_CEA / 1000.0;
    let cv_eq_kj_per_kg_k = cv_eq_over_r * R_CEA / 1000.0;

    // D.7's final reporting map: species below LOG_MIN are reported as
    // exactly 0 rather than denormal-ish exp() values.
    let nj: Vec<f64> = (0..ng)
        .map(|j| if ln_nj[j] > -87.0 { ln_nj[j].exp() } else { 0.0 })
        .collect();
    TpSolveDebug {
        assembled_first_iteration: assembled_first_iteration
            .expect("MAX_TOTAL_ITERATIONS > 0, so iteration 0 always runs"),
        converged,
        iterations_used,
        nj,
        n,
        temperature_k,
        condensed_active: is_active,
        condensed_nj: nj_c,
        gas_species: species[..ng].to_vec(),
        condensed_species: species[ng..].to_vec(),
        pressure_bar: final_pressure_bar,
        volume_m3_per_kg: final_volume_m3_per_kg,
        enthalpy_kj_per_kg,
        energy_kj_per_kg,
        entropy_kj_per_kg_k,
        gibbs_kj_per_kg,
        molar_mass_kg_per_kmol,
        mean_molecular_weight_kg_per_kmol,
        cp_frozen_kj_per_kg_k,
        cv_frozen_kj_per_kg_k,
        cp_eq_kj_per_kg_k,
        cv_eq_kj_per_kg_k,
        gamma_s,
        gas_mole_fractions: mole_fractions[..ng].to_vec(),
        condensed_mole_fractions: mole_fractions[ng..].to_vec(),
        gas_mass_fractions: mass_fractions[..ng].to_vec(),
        condensed_mass_fractions: mass_fractions[ng..].to_vec(),
    }
}

// b0_i = Σ_j a_ij^reac * x_j / Σ_j x_j * M_j  (moles of element i per kg mixture)
fn element_amounts(reactants: &[(f64, Species)], constituent_pool: &[Constituent]) -> Vec<f64> {
    let total_mass: f64 = reactants.iter().map(|(x, sp)| x * sp.data().mw()).sum();
    let mut b0 = vec![0.0; constituent_pool.len()];
    for (x, sp) in reactants {
        for (i, constituent) in constituent_pool.iter().enumerate() {
            for (count, c) in sp.data().constituents() {
                if c == constituent {
                    b0[i] += count * x / total_mass;
                }
            }
        }
    }
    b0
}

// state1 for HP problems (Appendix AI.2.6): the reactants' assigned
// enthalpy at their own reference temperature, divided by R.
// state1 = Σ_j x_j * H_j(T_init) / (R * Σ_j x_j * M_j)  [kmol*K/kg]
// h_over_r already returns H_j(T)/R directly (Kelvin units, guide A.1).
fn reactant_enthalpy_over_r(reactants: &[(f64, Species)], reactant_temperature_k: f64) -> f64 {
    let total_mass: f64 = reactants.iter().map(|(x, sp)| x * sp.data().mw()).sum();
    reactants
        .iter()
        .map(|(x, sp)| x * sp.data().h_over_r(reactant_temperature_k))
        .sum::<f64>()
        / total_mass
}

#[cfg(test)]
mod gauss_validation {
    //! Section C validation (CEA_RUST_PORT_GUIDE.md): the guide's suggested
    //! checks are (1) solve the commented-out test_h2_o2 4x4 system given in
    //! Section J.1 and (2) cross-check a random SPD system against nalgebra.
    use super::*;

    /// The assembled HP matrix for test_h2_o2 at the canonical initial guess
    /// (T=3800, nj=0.1/6 each), unknown layout [pi_H, pi_O, dln_n, dlnT].
    fn test_h2_o2_initial_matrix() -> DMatrix<f64> {
        #[rustfmt::skip]
        let rows: [[f64; 5]; 4] = [
            [0.16666666666666666, 0.05,                 0.1,                 0.29041578093904313, -2.8522519562949125],
            [0.05,                 0.11666666666666667, 0.083333333333333,   0.35522308781408640, -2.5626653151659138],
            [0.1,                  0.083333333333333,   0.0,                 0.50248557985367082, -2.5915522029331699],
            [0.29041578093904313,  0.35522308781408640, 0.50248557985367082, 4.60022515802027690, -10.014870869177022],
        ];
        DMatrix::from_fn(4, 5, |r, c| rows[r][c])
    }

    /// Split an (n, n+1) augmented matrix into its (n,n) coefficient matrix
    /// and length-n RHS, for cross-checking against nalgebra's own solver.
    fn split_augmented(g: &DMatrix<f64>) -> (DMatrix<f64>, DVector<f64>) {
        let n = g.nrows();
        let a = g.view((0, 0), (n, n)).clone_owned();
        let b = DVector::from_fn(n, |r, _| g[(r, n)]);
        (a, b)
    }

    #[test]
    fn solves_test_h2_o2_initial_matrix() {
        let mut g = test_h2_o2_initial_matrix();
        let (a, b) = split_augmented(&g);
        let expected = a.lu().solve(&b).expect("nalgebra solve failed");

        gauss(&mut g).expect("gauss reported singular on a well-posed system");
        for i in 0..4 {
            let got = g[(i, 4)];
            assert!(
                (got - expected[i]).abs() < 1e-9,
                "row {i}: gauss = {got}, nalgebra = {}",
                expected[i]
            );
        }
    }

    #[test]
    fn solves_random_spd_system() {
        // A simple deterministic PRNG (xorshift) so the test has no new deps
        // and is reproducible without pulling in `rand`.
        let mut state: u64 = 0x2545F4914F6CDD1D;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64 // in [0, 1)
        };

        const N: usize = 5;
        let m = DMatrix::from_fn(N, N, |_, _| next() * 2.0 - 1.0);
        let spd = &m * m.transpose() + DMatrix::identity(N, N) * (N as f64); // diagonally dominant SPD
        let b = DVector::from_fn(N, |_, _| next() * 2.0 - 1.0);

        let expected = spd.clone().lu().solve(&b).expect("nalgebra solve failed");

        let mut g = DMatrix::zeros(N, N + 1);
        g.view_mut((0, 0), (N, N)).copy_from(&spd);
        g.view_mut((0, N), (N, 1)).copy_from(&b);

        gauss(&mut g).expect("gauss reported singular on an SPD system");
        for i in 0..N {
            let got = g[(i, N)];
            assert!(
                (got - expected[i]).abs() < 1e-6,
                "row {i}: gauss = {got}, nalgebra = {}",
                expected[i]
            );
        }
    }
}

#[cfg(test)]
mod d1_assembly_validation {
    //! Section D.1 validation. The guide's own test_h2_o2 (J.1) matrix can't
    //! be used directly here: it's an HP problem (adds a ΔlnT column/energy
    //! row we haven't built — Section E), and the guide doesn't give us
    //! test_h2_air's (J.2) full 20-species list/order either, so there's no
    //! guide-provided TP ground truth to compare against as-is.
    //!
    //! Instead: restrict products to J.1's own 6-species H2/O2 list via
    //! `only`, and solve it as a TP problem (T=3800K, P=1.01325 bar, same
    //! canonical D.6 initial guess) — then compare against a matrix computed
    //! independently in Python from the exact same D.1 formulas and the same
    //! thermo.inp fit coefficients (already validated against real CEA2 in
    //! `cea_validation`). This exercises the same indexing/assembly code
    //! path the guide's matrix would, without needing HP or the missing J.2
    //! species list.
    use super::*;

    #[test]
    fn assembles_h2_o2_tp_matrix_at_canonical_guess() {
        let of_ratio = 15.87336;
        let mw_h2 = Species::H2.data().mw();
        let mw_o2 = Species::O2.data().mw();
        let reactants = [(1.0 / mw_h2, Species::H2), (of_ratio / mw_o2, Species::O2)];
        let only = [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ];

        let result = solve_for_products(
            &reactants, true, false, true, 3800.0, 1.01325, false, false, Some(&only),
        );
        let g = &result.assembled_first_iteration;

        #[rustfmt::skip]
        let expected: [[f64; 4]; 3] = [
            [0.16666666666666666, 0.05,                 0.1,                 -2.852251937205511],
            [0.05,                 0.11666666666666667, 0.08333333333333333, -2.5626653163685194],
            [0.1,                  0.08333333333333333, -0.0,                -2.59155220293317],
        ];

        assert_eq!(g.nrows(), 3, "expected 2 elements + dln_n = 3 equations");
        for r in 0..3 {
            for c in 0..4 {
                let got = g[(r, c)];
                assert!(
                    (got - expected[r][c]).abs() < 1e-8,
                    "G[{r}][{c}] = {got}, expected {}",
                    expected[r][c]
                );
            }
        }
    }
}

#[cfg(test)]
mod d5_convergence_validation {
    //! Section D.5 validation. There's no guide-provided converged TP
    //! ground truth available to us (see d1_assembly_validation's note on
    //! why J.1/J.2 don't apply here), so instead of comparing against a
    //! borrowed number, this checks the two things that must be true at any
    //! correct equilibrium, computed independently from the public
    //! thermo_species API rather than any of solve_for_products' internal
    //! state:
    //!   1. Element conservation: Σ_j a_ij * nj[j] == b0_i.
    //!   2. Lagrangian stationarity: there's a SINGLE (pi_H, pi_O) such that
    //!      mu_j/RT == a_Hj*pi_H + a_Oj*pi_O for every one of the 6 species
    //!      (solved from H2/O2 alone, then checked against the other 4).
    use super::*;

    #[test]
    fn h2_o2_tp_converges_to_a_valid_equilibrium() {
        let of_ratio = 15.87336;
        let mw_h2 = Species::H2.data().mw();
        let mw_o2 = Species::O2.data().mw();
        let reactants = [(1.0 / mw_h2, Species::H2), (of_ratio / mw_o2, Species::O2)];
        let only = [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ];
        // species order matches `only`, confirmed against Species::all()'s
        // enumeration order (alphabetical-ish); [H, O] atoms per species:
        let stoich: [[f64; 2]; 6] = [
            [1.0, 0.0], // H
            [2.0, 0.0], // H2
            [2.0, 1.0], // H2O
            [0.0, 1.0], // O
            [0.0, 2.0], // O2
            [1.0, 1.0], // OH
        ];

        let temperature_k = 3800.0;
        let pressure_bar = 1.01325;
        let result = solve_for_products(
            &reactants,
            true,
            false,
            true,
            temperature_k,
            pressure_bar,
            false,
            false,
            Some(&only),
        );

        assert!(
            result.converged,
            "did not converge within {} iterations",
            result.iterations_used
        );
        assert!(result.iterations_used < 50);

        // --- Check 1: element conservation ---
        let b0 = element_amounts(&reactants, &[Constituent::H, Constituent::O]);
        for (k, &b0_k) in b0.iter().enumerate() {
            let b_k: f64 = (0..6).map(|j| stoich[j][k] * result.nj[j]).sum();
            assert!(
                (b0_k - b_k).abs() < 1e-8,
                "element {k}: b0 = {b0_k}, Σ a_ij*nj = {b_k}"
            );
        }

        // --- Check 2: Lagrangian stationarity, independent of any internal
        // mu_g/matrix state — recomputed here straight from the public API.
        let species = [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ];
        let mu: Vec<f64> = (0..6)
            .map(|j| {
                let d = species[j].data();
                d.h_over_r(temperature_k) / temperature_k - d.s_over_r(temperature_k)
                    + result.nj[j].ln()
                    + (pressure_bar / result.n).ln()
            })
            .collect();

        // Solve for (pi_H, pi_O) from H2 (index 1, [2,0]) and O2 (index 4, [0,2]).
        let pi_h = mu[1] / 2.0;
        let pi_o = mu[4] / 2.0;

        for j in 0..6 {
            let predicted = stoich[j][0] * pi_h + stoich[j][1] * pi_o;
            assert!(
                (mu[j] - predicted).abs() < 1e-6,
                "species {j}: mu/RT = {}, a.pi = {predicted}",
                mu[j]
            );
        }
    }
}

#[cfg(test)]
mod e_hp_validation {
    //! Section E validation, test_h2_o2 (J.1) — this is finally the real
    //! deal: an HP problem, which is what J.1 actually is (unlike D's
    //! TP-only tests, which had to fall back to independently-derived
    //! numbers because J.1 didn't apply). Both the canonical-guess matrix
    //! AND the converged T/mole fractions are the guide's actual numbers.
    use super::*;

    fn h2_o2_reactants() -> [(f64, Species); 2] {
        let of_ratio = 15.87336;
        let mw_h2 = Species::H2.data().mw();
        let mw_o2 = Species::O2.data().mw();
        [(1.0 / mw_h2, Species::H2), (of_ratio / mw_o2, Species::O2)]
    }

    fn only_h2_o2_products() -> [Species; 6] {
        [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ]
    }

    #[test]
    fn assembles_h2_o2_hp_initial_matrix() {
        let reactants = h2_o2_reactants();
        let only = only_h2_o2_products();
        let state1 = reactant_enthalpy_over_r(&reactants, 2000.0);

        let result = solve_for_products(&reactants, false, false, true, state1, 1.01325, false, false, Some(&only));
        let g = &result.assembled_first_iteration;

        #[rustfmt::skip]
        let expected: [[f64; 5]; 4] = [
            [0.16666666666666666, 0.05,                 0.1,                 0.29041578093904313, -2.8522519562949125],
            [0.05,                 0.11666666666666667, 0.083333333333333,   0.35522308781408640, -2.5626653151659138],
            [0.1,                  0.083333333333333,   0.0,                 0.50248557985367082, -2.5915522029331699],
            [0.29041578093904313,  0.35522308781408640, 0.50248557985367082, 4.60022515802027690, -10.014870869177022],
        ];

        assert_eq!(g.nrows(), 4, "expected 2 elements + dln_n + dlnT = 4 equations");
        for r in 0..4 {
            for c in 0..5 {
                let got = g[(r, c)];
                assert!(
                    (got - expected[r][c]).abs() < 1e-6,
                    "G[{r}][{c}] = {got}, expected {}",
                    expected[r][c]
                );
            }
        }
    }

    #[test]
    fn converges_h2_o2_hp_to_guide_numbers() {
        let reactants = h2_o2_reactants();
        let only = only_h2_o2_products();
        let state1 = reactant_enthalpy_over_r(&reactants, 2000.0);

        let result = solve_for_products(&reactants, false, false, true, state1, 1.01325, false, false, Some(&only));

        assert!(
            result.converged,
            "did not converge within {} iterations",
            result.iterations_used
        );

        assert!(
            (result.temperature_k - 3181.2589964650856).abs() < 1e-3,
            "T = {}, expected 3181.2589964650856",
            result.temperature_k
        );

        // [H, H2, H2O, O, O2, OH]
        let expected_x = [
            0.06590368364478945,
            0.06152807617606959,
            0.37638785326991414,
            0.09771754800744899,
            0.23381507062970186,
            0.16464776827208294,
        ];
        let total: f64 = result.nj.iter().sum();
        for (j, &expected_xj) in expected_x.iter().enumerate() {
            let xj = result.nj[j] / total;
            assert!(
                (xj - expected_xj).abs() < 1e-6,
                "species {j}: x = {xj}, expected {expected_xj}"
            );
        }
    }
}

#[cfg(test)]
mod e_sp_validation {
    //! Section E validation, test_sp_problem (J.4): equimolar H2/O2, SP at
    //! 1.01325 bar. Unlike HP, the guide gives the converged ln_nj directly
    //! (not just mole fractions), so this checks ln_nj to the guide's own
    //! 1e-6 relative tolerance rather than reconstructing mole fractions.
    use super::*;

    #[test]
    fn converges_h2_o2_sp_to_guide_numbers() {
        let reactants = [(1.0, Species::H2), (1.0, Species::O2)];
        let only = [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ];
        let state1 = 1.804; // S/R, given directly by the guide (J.4).

        let result =
            solve_for_products(&reactants, false, true, true, state1, 1.01325, false, false, Some(&only));

        assert!(
            result.converged,
            "did not converge within {} iterations",
            result.iterations_used
        );
        assert!(
            (result.temperature_k - 3244.7454823232674).abs() < 1e-3,
            "T = {}, expected 3244.7454823232674",
            result.temperature_k
        );

        // order: H, H2, H2O, O, O2, OH
        let expected_ln_nj = [
            -5.3971690318524281,
            -5.5822036166964706,
            -3.9841684012236138,
            -5.0642456920269865,
            -4.4084879383149573,
            -4.6538971360500030,
        ];
        for (j, &expected) in expected_ln_nj.iter().enumerate() {
            let got = result.nj[j].ln();
            assert!(
                (got - expected).abs() < 1e-6,
                "species {j}: ln_nj = {got}, expected {expected}"
            );
        }
    }
}

#[cfg(test)]
mod e_sv_validation {
    //! Section E validation, test_sv_problem (J.4): same reactants/entropy
    //! target as SP, but V fixed instead of P. The guide gives the
    //! converged ln_nj directly, so — like SP — this checks ln_nj at 1e-6
    //! rather than reconstructing mole fractions.
    use super::*;

    #[test]
    fn converges_h2_o2_sv_to_guide_numbers() {
        let reactants = [(1.0, Species::H2), (1.0, Species::O2)];
        let only = [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ];
        let state1 = 1.804; // S/R, same as SP (J.4).
        let volume = 1.0 / 0.072081; // m^3/kg

        let result = solve_for_products(
            &reactants, false, true, false, state1, volume, false, false, Some(&only),
        );

        assert!(
            result.converged,
            "did not converge within {} iterations",
            result.iterations_used
        );
        assert!(
            (result.temperature_k - 3258.5925134713739).abs() < 1e-3,
            "T = {}, expected 3258.5925134713739",
            result.temperature_k
        );

        // order: H, H2, H2O, O, O2, OH
        let expected_ln_nj = [
            -5.3856081862190726,
            -5.5747945411551276,
            -3.9891190704470452,
            -5.0532891277264751,
            -4.4106610522293108,
            -4.6460395167270923,
        ];
        for (j, &expected) in expected_ln_nj.iter().enumerate() {
            let got = result.nj[j].ln();
            assert!(
                (got - expected).abs() < 1e-6,
                "species {j}: ln_nj = {got}, expected {expected}"
            );
        }
    }
}

#[cfg(test)]
mod e_tv_validation {
    //! Section E validation, TV. The guide's own J.2/J.3 pair (test_h2_air /
    //! test_example2) is exactly a TP-at-(T,P) vs TV-at-the-equivalent-V
    //! cross-check ("TV at the equivalent volume reproduces the TP state" —
    //! an "excellent cross-check" per the guide) — but it uses the H2/air
    //! 20-species system we can't reproduce (see d1_assembly_validation's
    //! note). So: do the same cross-check ourselves on H2/O2. Solve TP,
    //! derive V from its converged n via the same ideal-gas closure
    //! solve_for_products itself uses, then confirm TV at that V reproduces
    //! the identical composition.
    use super::*;

    #[test]
    fn tv_at_equivalent_volume_reproduces_tp_state() {
        const R_CEA: f64 = 8314.51;
        let of_ratio = 15.87336;
        let mw_h2 = Species::H2.data().mw();
        let mw_o2 = Species::O2.data().mw();
        let reactants = [(1.0 / mw_h2, Species::H2), (of_ratio / mw_o2, Species::O2)];
        let only = [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ];
        let temperature_k = 3800.0;
        let pressure_bar = 1.01325;

        let tp = solve_for_products(
            &reactants,
            true,
            false,
            true,
            temperature_k,
            pressure_bar,
            false,
            false,
            Some(&only),
        );
        assert!(tp.converged);

        // P = 1e-5 * R * n * T / V  =>  V = 1e-5 * R * n * T / P
        let equivalent_volume = 1e-5 * R_CEA * tp.n * temperature_k / pressure_bar;

        let tv = solve_for_products(
            &reactants,
            true,
            false,
            false,
            temperature_k,
            equivalent_volume,
            false,
            false,
            Some(&only),
        );
        assert!(tv.converged, "TV did not converge in {} iters", tv.iterations_used);

        assert!(
            (tv.n - tp.n).abs() / tp.n < 1e-6,
            "n: TV = {}, TP = {}",
            tv.n,
            tp.n
        );
        for j in 0..6 {
            let tv_ln = tv.nj[j].ln();
            let tp_ln = tp.nj[j].ln();
            assert!(
                (tv_ln - tp_ln).abs() < 1e-6,
                "species {j}: TV ln_nj = {tv_ln}, TP ln_nj = {tp_ln}"
            );
        }
    }
}

#[cfg(test)]
mod e_uv_validation {
    //! Section E validation, UV. J.6 (test_example4) is the guide's only UV
    //! case, but it's a 92-species trace-mode stress test with condensed
    //! species — not reproducible here any more than J.2's 20-species TV
    //! case was. So: the same self-derived cross-check as TV, but for the
    //! other closure. Solve TP, derive UV's (state1 = U⁰/R, V) from that
    //! converged state, then confirm UV reproduces it.
    use super::*;

    #[test]
    fn uv_at_equivalent_state_reproduces_tp_state() {
        const R_CEA: f64 = 8314.51;
        let of_ratio = 15.87336;
        let mw_h2 = Species::H2.data().mw();
        let mw_o2 = Species::O2.data().mw();
        let reactants = [(1.0 / mw_h2, Species::H2), (of_ratio / mw_o2, Species::O2)];
        let species = [
            Species::H,
            Species::H2,
            Species::H2O,
            Species::O,
            Species::O2,
            Species::OH,
        ];
        let temperature_k = 3800.0;
        let pressure_bar = 1.01325;

        let tp = solve_for_products(
            &reactants,
            true,
            false,
            true,
            temperature_k,
            pressure_bar,
            false,
            false,
            Some(&species),
        );
        assert!(tp.converged);

        // At convergence, hsu_delta = state1/T - Σ nj*u_j == 0, so
        // state1 = T * Σ nj*(u_over_r(T)/T) = Σ nj*u_over_r(T) — the T
        // cancels, so this is just a plain weighted sum of u_over_r.
        let state1_uv: f64 = (0..6)
            .map(|j| tp.nj[j] * species[j].data().u_over_r(tp.temperature_k))
            .sum();
        // Same P/V closure as TV.
        let equivalent_volume = 1e-5 * R_CEA * tp.n * temperature_k / pressure_bar;

        let uv = solve_for_products(
            &reactants,
            false,
            false,
            false,
            state1_uv,
            equivalent_volume,
            false,
            false,
            Some(&species),
        );
        assert!(uv.converged, "UV did not converge in {} iters", uv.iterations_used);

        assert!(
            (uv.temperature_k - tp.temperature_k).abs() < 1e-3,
            "T: UV = {}, TP = {}",
            uv.temperature_k,
            tp.temperature_k
        );
        for j in 0..6 {
            let uv_ln = uv.nj[j].ln();
            let tp_ln = tp.nj[j].ln();
            assert!(
                (uv_ln - tp_ln).abs() < 1e-6,
                "species {j}: UV ln_nj = {uv_ln}, TP ln_nj = {tp_ln}"
            );
        }
    }
}

#[cfg(test)]
mod f_condensed_validation {
    //! Section F validation. The guide's only small-system condensed case
    //! (J.5, test_condensed_problem) needs F.4's solid/liquid phase-transition
    //! logic to pick H2O(L) over H2O(cr) near 273K — not implemented here
    //! (see the module note on Section F's scope) — so there's no
    //! guide-provided ground truth usable as-is.
    //!
    //! Instead: a fuel-rich carbon/oxygen TP case, deliberately with more
    //! carbon than the available O2 can fully oxidize even to CO, so some
    //! C(gr) should remain solid at equilibrium (F.3.3 should activate it).
    //! Validated the same way as D.5's gas-only case: element conservation
    //! and Lagrangian stationarity, computed independently from the public
    //! thermo_species API — but now checked for the condensed species too,
    //! which is exactly what F.2's condensed row represents.
    use super::*;

    #[test]
    fn carbon_deposits_as_condensed_graphite() {
        let reactants = [(2.5, Species::C__1), (1.0, Species::O2)];
        let temperature_k = 1200.0;
        let pressure_bar = 1.0;

        let result = solve_for_products(
            &reactants,
            true,
            false,
            true,
            temperature_k,
            pressure_bar,
            true, // include_condensed
            false,
            None, // full auto-generated C/O species pool — see note below
        );

        assert!(
            result.converged,
            "did not converge within {} iterations",
            result.iterations_used
        );

        // graphite is the only condensed candidate in the C/O pool — it
        // should have been activated by F.3.3, since 2.5 mol C is more than
        // 1 mol O2 (= 2 mol O) can fully consume even as CO (max 2 mol C).
        let graphite_idx = result
            .condensed_species
            .iter()
            .position(|&s| s == Species::C__1)
            .expect("C(gr) should be a condensed candidate for a C/O system");
        assert!(result.condensed_active[graphite_idx], "graphite never got activated");
        let nj_graphite = result.condensed_nj[graphite_idx];
        assert!(nj_graphite > 0.0, "nj_graphite = {nj_graphite}, expected positive");

        // (c_count, o_count) for any C/O species, straight from the public API.
        let co_counts = |sp: Species| -> (f64, f64) {
            let mut c = 0.0;
            let mut o = 0.0;
            for (count, constituent) in sp.data().constituents() {
                match constituent {
                    Constituent::C => c = *count,
                    Constituent::O => o = *count,
                    other => panic!("unexpected constituent {other:?} in a C/O-only system"),
                }
            }
            (c, o)
        };

        // --- Check 1: element conservation, gas + condensed ---
        let b0 = element_amounts(&reactants, &[Constituent::C, Constituent::O]);
        let gas_b: [f64; 2] = [0, 1].map(|k| {
            (0..result.gas_species.len())
                .map(|j| {
                    let (c, o) = co_counts(result.gas_species[j]);
                    (if k == 0 { c } else { o }) * result.nj[j]
                })
                .sum()
        });
        let (gr_c, gr_o) = co_counts(Species::C__1);
        let cond_b = [gr_c * nj_graphite, gr_o * nj_graphite];
        for k in 0..2 {
            assert!(
                (b0[k] - gas_b[k] - cond_b[k]).abs() < 1e-8,
                "element {k}: b0 = {}, gas+cond = {}",
                b0[k],
                gas_b[k] + cond_b[k]
            );
        }

        // --- Check 2: Lagrangian stationarity, gas AND condensed ---
        // mu_j/RT for every gas species actually present in the solution.
        // Species reported as exactly 0 (below the D.7 LOG_MIN cutoff) have
        // no defined mu and are skipped below.
        let mu_gas: Vec<f64> = (0..result.gas_species.len())
            .map(|j| {
                let d = result.gas_species[j].data();
                d.h_over_r(temperature_k) / temperature_k - d.s_over_r(temperature_k)
                    + result.nj[j].ln()
                    + (pressure_bar / result.n).ln()
            })
            .collect();
        // Solve (pi_C, pi_O) from the first two species with independent
        // (c, o) stoichiometry, then confirm every OTHER gas species (and
        // graphite) satisfies mu_j/RT == a_Cj*pi_C + a_Oj*pi_O with that
        // same pair — exactly the condensed row's own equation (F.2).
        let stoich_gas: Vec<(f64, f64)> =
            (0..result.gas_species.len()).map(|j| co_counts(result.gas_species[j])).collect();
        let (i1, i2) = (0..stoich_gas.len())
            .flat_map(|i| (0..stoich_gas.len()).map(move |j| (i, j)))
            .find(|&(i, j)| {
                let ((c1, o1), (c2, o2)) = (stoich_gas[i], stoich_gas[j]);
                result.nj[i] > 0.0 && result.nj[j] > 0.0 && (c1 * o2 - c2 * o1).abs() > 1e-9
            })
            .expect("need two gas species with independent C/O stoichiometry");
        let ((c1, o1), (c2, o2)) = (stoich_gas[i1], stoich_gas[i2]);
        let det = c1 * o2 - c2 * o1;
        let pi_c = (mu_gas[i1] * o2 - mu_gas[i2] * o1) / det;
        let pi_o = (c1 * mu_gas[i2] - c2 * mu_gas[i1]) / det;

        for j in 0..result.gas_species.len() {
            if result.nj[j] == 0.0 {
                continue; // below LOG_MIN reporting cutoff — no defined mu
            }
            let (c, o) = stoich_gas[j];
            let predicted = c * pi_c + o * pi_o;
            assert!(
                (mu_gas[j] - predicted).abs() < 1e-6,
                "gas species {j} ({:?}): mu/RT = {}, a.pi = {predicted}",
                result.gas_species[j],
                mu_gas[j]
            );
        }

        // Condensed: mu_c/RT = h_c - s_c (no ln(n) or ln(P) term, D.1/F.2).
        let graphite = Species::C__1.data();
        let mu_graphite = graphite.h_over_r(temperature_k) / temperature_k
            - graphite.s_over_r(temperature_k);
        let predicted_graphite = gr_c * pi_c + gr_o * pi_o;
        assert!(
            (mu_graphite - predicted_graphite).abs() < 1e-6,
            "C(gr): mu/RT = {mu_graphite}, a.pi = {predicted_graphite}"
        );
    }
}

#[cfg(test)]
mod g_singular_recovery_validation {
    //! Section G validation.
    //!
    //! `recovers_gauss_singular_via_element_reduction` directly exercises
    //! the recovery ladder on a case it's actually built to fix: two gas
    //! species with IDENTICAL (C, O) stoichiometry (CO2 twice, effectively)
    //! give a genuinely rank-deficient element block from the very first
    //! iteration — `gauss` reports a real Err(ierr), G.2 kicks in, and (since
    //! `iteration <= 1` on the very first hit) branch (e)'s element
    //! reduction fires immediately and the solver still reaches a valid
    //! equilibrium.
    //!
    //! `narrow_species_list_fails_gracefully_instead_of_panicking` is the
    //! scenario that used to panic before Section G existed (restricting
    //! products to only [CO, CO2, O2, C(gr)] drives CO2/O2 toward zero,
    //! leaving only CO — a 2-element system pinned by one populated species,
    //! which is rank-1). It remains genuinely unsolvable with what's
    //! implemented here: G.2(e) firing immediately would drop the O-element
    //! constraint entirely with nothing left to backstop it (divergence),
    //! and (f) trace-reseed alone can't fix a structural degeneracy that
    //! condensed insertion (only available post-convergence) is what's
    //! really needed to break. The guide's own branches (a) and (d) — which
    //! need F.4's rank/j_switch bookkeeping and Section H's ions,
    //! respectively — may be what actually rescues cases like this in real
    //! CEA. What Section G *does* guarantee here, and is worth testing on
    //! its own: no panic, and a clean `converged: false` instead.
    use super::*;

    #[test]
    fn recovers_gauss_singular_via_element_reduction() {
        // Two chemically distinct species that both happen to be pure CO2
        // in (C, O) terms doesn't exist in the crate, so instead: force the
        // degeneracy directly by using `only` with just CO2 and O2 at a
        // temperature/pressure where O2 gets truncated away almost
        // immediately, leaving CO2 as the only populated species — the
        // same "one species, two elements" rank-1 block as the narrow
        // carbon case, but recoverable because there's no third element
        // (C) whose row still needs O's constraint to stay meaningful: once
        // O is excluded, CO2's own amount is still fully pinned by the C row
        // alone (CO2 is the ONLY carbon-bearing species in this list).
        let reactants = [(1.0, Species::C__1), (5.0, Species::O2)];
        let only = [Species::CO2, Species::O2];
        let temperature_k = 2000.0;
        let pressure_bar = 1.0;

        let result = solve_for_products(
            &reactants,
            true,
            false,
            true,
            temperature_k,
            pressure_bar,
            false,
            false,
            Some(&only),
        );

        assert!(
            result.converged,
            "did not converge within {} iterations",
            result.iterations_used
        );
        for (j, &nj) in result.nj.iter().enumerate() {
            assert!(nj.is_finite() && nj >= 0.0, "species {j}: nj = {nj}");
        }

        // Element conservation must still hold even though (if the
        // reduction fired) one element's row was dropped from the matrix —
        // the guide's own contract for this fallback.
        let b0 = element_amounts(&reactants, &[Constituent::C, Constituent::O]);
        let stoich = [[1.0, 2.0], [0.0, 2.0]]; // CO2, O2
        for (k, &b0_k) in b0.iter().enumerate() {
            let b_k: f64 = (0..2).map(|j| stoich[j][k] * result.nj[j]).sum();
            assert!((b0_k - b_k).abs() < 1e-6, "element {k}: b0 = {b0_k}, Σa*nj = {b_k}");
        }
    }

    #[test]
    fn narrow_species_list_fails_gracefully_instead_of_panicking() {
        let reactants = [(2.5, Species::C__1), (1.0, Species::O2)];
        let only = [Species::CO, Species::CO2, Species::O2, Species::C__1];
        let temperature_k = 1200.0;
        let pressure_bar = 1.0;

        let result = solve_for_products(
            &reactants,
            true,
            false,
            true,
            temperature_k,
            pressure_bar,
            true, // include_condensed
            false,
            Some(&only),
        );

        // The real assertion: this returns (doesn't panic) with an honest
        // "didn't converge" rather than a singular-matrix crash.
        assert!(
            !result.converged,
            "expected this known-hard case to fail gracefully, but it unexpectedly converged \
             (nice — the recovery ladder is now doing more than this test assumed)"
        );
    }
}

#[cfg(test)]
mod h_ions_validation {
    //! Section H validation. The guide's own suggestion (a hot-air plasma
    //! with NO+/e-/N+/O+/O2-) needs species/validation data we don't have
    //! access to, so instead: the minimal possible ionization system —
    //! atomic hydrogen at a temperature high enough (15000K) that H <=> H+
    //! + e- proceeds significantly at 1 bar. Validated the same way as
    //! d5_convergence_validation: charge neutrality, element conservation,
    //! and Lagrangian stationarity (now with the E pseudo-element folded
    //! in), all computed independently from the public API.
    use super::*;

    #[test]
    fn hydrogen_ionizes_and_stays_charge_neutral() {
        let reactants = [(1.0, Species::H)];
        let only = [Species::H, Species::HPlus, Species::E];
        let temperature_k = 15000.0;
        let pressure_bar = 1.0;

        let result = solve_for_products(
            &reactants,
            true,
            false,
            true,
            temperature_k,
            pressure_bar,
            false,
            true, // include_ions
            Some(&only),
        );

        assert!(
            result.converged,
            "did not converge within {} iterations",
            result.iterations_used
        );

        let h_idx = result.gas_species.iter().position(|&s| s == Species::H).unwrap();
        let hplus_idx = result.gas_species.iter().position(|&s| s == Species::HPlus).unwrap();
        let e_idx = result.gas_species.iter().position(|&s| s == Species::E).unwrap();

        // Sanity: meaningful ionization actually occurred (not just H).
        assert!(
            result.nj[hplus_idx] > 1e-6,
            "expected significant ionization at 15000K, nj_H+ = {}",
            result.nj[hplus_idx]
        );

        // --- Check 1: charge neutrality — H+ has charge +1 (a_E = -1), e-
        // has charge -1 (a_E = +1), H is neutral (a_E = 0): nj_H+ == nj_e-.
        assert!(
            (result.nj[hplus_idx] - result.nj[e_idx]).abs() < 1e-9,
            "nj_H+ = {}, nj_e- = {}",
            result.nj[hplus_idx],
            result.nj[e_idx]
        );

        // --- Check 2: element (H) conservation, counting H+ as hydrogen too.
        let b0_h = element_amounts(&reactants, &[Constituent::H])[0];
        let total_h = result.nj[h_idx] + result.nj[hplus_idx];
        assert!((b0_h - total_h).abs() < 1e-8, "b0_H = {b0_h}, nj_H + nj_H+ = {total_h}");

        // --- Check 3: Lagrangian stationarity, now over (pi_H, pi_E).
        // mu_j/RT = h_j - s_j + ln(nj) + ln(P/n), same formula regardless of
        // species charge (D.1) — the charge only shows up in the element
        // (E) stoichiometry used to relate mu_j across species.
        let mu = |sp: Species, nj: f64| -> f64 {
            let d = sp.data();
            d.h_over_r(temperature_k) / temperature_k - d.s_over_r(temperature_k)
                + nj.ln()
                + (pressure_bar / result.n).ln()
        };
        let mu_h = mu(Species::H, result.nj[h_idx]);
        let mu_hplus = mu(Species::HPlus, result.nj[hplus_idx]);
        let mu_e = mu(Species::E, result.nj[e_idx]);

        // H: a_H=1, a_E=0  => mu_H = pi_H
        // H+: a_H=1, a_E=-1 => mu_H+ = pi_H - pi_E
        let pi_h = mu_h;
        let pi_e_multiplier = pi_h - mu_hplus;
        // e-: a_H=0, a_E=1 => mu_e- should equal pi_E
        assert!(
            (mu_e - pi_e_multiplier).abs() < 1e-6,
            "e-: mu/RT = {mu_e}, pi_E = {pi_e_multiplier}"
        );
    }
}

#[cfg(test)]
mod i_postprocess_validation {
    //! Section I validation. No guide-provided ground truth exists for I.1's
    //! reported properties or I.2's equilibrium partials, so instead this
    //! checks physical invariants and cross-validates against independent
    //! finite-difference derivatives computed from already-validated solves
    //! (Section E's HP/SP/SV cases):
    //!   1. Mole/mass fractions sum to 1.
    //!   2. HP's own conserved quantity: I.1's reported enthalpy must equal
    //!      the reactants' enthalpy (state1*R) — this is exactly what HP's
    //!      Newton solve enforces, checked here via the independent I.1
    //!      post-processing path rather than any of the solver's internal
    //!      state.
    //!   3. cp_eq (I.2): central finite difference of I.1's H over two
    //!      nearby TP solves (same P) must match the analytic cp_eq at the
    //!      midpoint solve.
    //!   4. gamma_s (I.2): central finite difference of ln P over ln V from
    //!      two nearby SV solves (same S) must match the analytic gamma_s at
    //!      the midpoint solve.
    use super::*;

    fn h2_o2_reactants() -> [(f64, Species); 2] {
        let of_ratio = 15.87336;
        let mw_h2 = Species::H2.data().mw();
        let mw_o2 = Species::O2.data().mw();
        [(1.0 / mw_h2, Species::H2), (of_ratio / mw_o2, Species::O2)]
    }

    fn only_h2_o2_products() -> [Species; 6] {
        [Species::H, Species::H2, Species::H2O, Species::O, Species::O2, Species::OH]
    }

    #[test]
    fn fractions_sum_to_one_and_hp_conserves_enthalpy() {
        let reactants = h2_o2_reactants();
        let only = only_h2_o2_products();
        let state1 = reactant_enthalpy_over_r(&reactants, 2000.0);

        let result =
            solve_for_products(&reactants, false, false, true, state1, 1.01325, false, false, Some(&only));
        assert!(result.converged);

        let mole_sum: f64 = result.gas_mole_fractions.iter().sum::<f64>()
            + result.condensed_mole_fractions.iter().sum::<f64>();
        assert!((mole_sum - 1.0).abs() < 1e-9, "mole fractions sum to {mole_sum}, expected 1");

        let mass_sum: f64 = result.gas_mass_fractions.iter().sum::<f64>()
            + result.condensed_mass_fractions.iter().sum::<f64>();
        assert!((mass_sum - 1.0).abs() < 1e-9, "mass fractions sum to {mass_sum}, expected 1");

        const R_CEA: f64 = 8314.51;
        let expected_h_kj_per_kg = state1 * R_CEA / 1000.0;
        assert!(
            (result.enthalpy_kj_per_kg - expected_h_kj_per_kg).abs()
                < 1e-4 * expected_h_kj_per_kg.abs(),
            "I.1 H = {} kJ/kg, expected (HP-conserved) {} kJ/kg",
            result.enthalpy_kj_per_kg,
            expected_h_kj_per_kg
        );
    }

    #[test]
    fn cp_eq_matches_finite_difference_of_enthalpy_over_t() {
        let reactants = h2_o2_reactants();
        let only = only_h2_o2_products();
        let t0 = 3181.2589964650856; // J.1's converged HP temperature
        let p = 1.01325;
        let dt = 0.5;

        let solve_at = |t: f64| {
            solve_for_products(&reactants, true, false, true, t, p, false, false, Some(&only))
        };
        let mid = solve_at(t0);
        let minus = solve_at(t0 - dt);
        let plus = solve_at(t0 + dt);
        assert!(mid.converged && minus.converged && plus.converged);

        let cp_findiff = (plus.enthalpy_kj_per_kg - minus.enthalpy_kj_per_kg) / (2.0 * dt);
        assert!(
            (mid.cp_eq_kj_per_kg_k - cp_findiff).abs() < 1e-2 * cp_findiff.abs(),
            "cp_eq = {}, finite-difference dH/dT = {}",
            mid.cp_eq_kj_per_kg_k,
            cp_findiff
        );
    }

    #[test]
    fn gamma_s_matches_finite_difference_of_lnp_over_lnv() {
        // Equimolar H2/O2, matching J.4's SV problem (state1 = 1.804 is the
        // entropy target for THIS mixture; equal-mass amounts here would
        // drive T below CEA's 160 K validity floor and fail the G.1 bounds
        // check).
        let reactants = [(1.0, Species::H2), (1.0, Species::O2)];
        let only = only_h2_o2_products();
        let s_over_r = 1.804; // J.4's SV state1
        let v0 = 1.0 / 0.072081; // J.4's SV converged volume
        let eps = 1e-4;

        let solve_at = |v: f64| {
            solve_for_products(&reactants, false, true, false, s_over_r, v, false, false, Some(&only))
        };
        let mid = solve_at(v0);
        let minus = solve_at(v0 * (1.0 - eps));
        let plus = solve_at(v0 * (1.0 + eps));
        assert!(
            mid.converged && minus.converged && plus.converged,
            "converged: mid={} (T={}, iters={}), minus={} (T={}, iters={}), plus={} (T={}, iters={})",
            mid.converged, mid.temperature_k, mid.iterations_used,
            minus.converged, minus.temperature_k, minus.iterations_used,
            plus.converged, plus.temperature_k, plus.iterations_used,
        );

        let gamma_findiff = -(plus.pressure_bar.ln() - minus.pressure_bar.ln())
            / (plus.volume_m3_per_kg.ln() - minus.volume_m3_per_kg.ln());
        assert!(
            (mid.gamma_s - gamma_findiff).abs() < 1e-3 * gamma_findiff.abs(),
            "gamma_s = {}, finite-difference -dlnP/dlnV = {}",
            mid.gamma_s,
            gamma_findiff
        );
    }
}

#[cfg(test)]
mod cea_validation {
    //! test_h2_o2 (CEA_RUST_PORT_GUIDE.md Section J.1): H2/O2 at
    //! O/F mass ratio 15.87336. b0 computed independently in Python from
    //! the same mass-conservation formula (B.3) using this crate's MW values.
    use super::*;

    #[test]
    fn b0_matches_h2_o2_of_ratio() {
        // O/F is a MASS ratio: pick moles so mass_O2 / mass_H2 == of_ratio.
        let of_ratio = 15.87336;
        let mw_h2 = Species::H2.data().mw();
        let mw_o2 = Species::O2.data().mw();
        let reactants = [(1.0 / mw_h2, Species::H2), (of_ratio / mw_o2, Species::O2)];
        let constituent_pool = [Constituent::H, Constituent::O];

        let b0 = element_amounts(&reactants, &constituent_pool);

        let expected_b0_h = 0.058798161538484495;
        let expected_b0_o = 0.058798141246477996;

        assert!(
            (b0[0] - expected_b0_h).abs() < 1e-9,
            "b0_H = {}, expected {expected_b0_h}",
            b0[0]
        );
        assert!(
            (b0[1] - expected_b0_o).abs() < 1e-9,
            "b0_O = {}, expected {expected_b0_o}",
            b0[1]
        );
    }
}

#[cfg(test)]
mod tp_stability_suite {
    //! Same stability sweep as thermo::mixture::tests::solve_for_products_stability
    //! (16 reactant mixtures x 15 (T, P) conditions, all solved as TP), ported
    //! over so the two solvers' convergence behavior can be compared directly
    //! on identical inputs. Reports the same per-condition/per-scenario grid;
    //! unlike mixture.rs's version this doesn't hard-fail on <100% — it's
    //! diagnostic (run with `-- --nocapture` to see the grid), since we don't
    //! yet know what pass rate to expect from the new solver.
    use super::*;

    fn reactant_scenarios() -> Vec<Vec<(f64, Species)>> {
        vec![
            vec![(1.0, Species::H2O)],
            vec![(1.0, Species::CH4)],
            vec![(1.0, Species::CH4), (2.0, Species::O2)],
            vec![(1.0, Species::CH4), (1.0, Species::O)],
            vec![(1.0, Species::CH4), (1.0, Species::HPlus), (1.0, Species::OHNeg)],
            vec![(1.0, Species::N2)],
            vec![(1.0, Species::O2)],
            vec![(1.0, Species::CO2)],
            vec![(1.0, Species::NH3)],
            vec![(2.0, Species::H2), (1.0, Species::O2)],
            vec![(1.0, Species::CO), (1.0, Species::H2O)],
            vec![(1.0, Species::N2), (3.0, Species::H2)],
            vec![(1.0, Species::CO2), (1.0, Species::H2)],
            vec![(2.0, Species::NH3), (1.5, Species::O2)],
            vec![(1.0, Species::N2), (1.0, Species::O2), (1.0, Species::H2O)],
            vec![(1.0, Species::E), (1.0, Species::HPlus)],
        ]
    }

    fn conditions() -> Vec<(f64, f64)> {
        vec![
            (300.0, 1.0),
            (300.0, 50.0),
            (500.0, 50.0),
            (1_000.0, 1.0),
            (1_000.0, 50.0),
            (2_000.0, 50.0),
            (3_000.0, 50.0),
            (4_000.0, 50.0),
            (5_000.0, 1.0),
            (5_000.0, 50.0),
            (5_000.0, 200.0),
            (7_000.0, 50.0),
            (10_000.0, 50.0),
            (15_000.0, 50.0),
            (20_000.0, 50.0),
        ]
    }

    #[test]
    fn solve_for_products_stability() {
        let reactants_list = reactant_scenarios();
        let conditions = conditions();
        let mut pass_fail = vec![vec![false; reactants_list.len()]; conditions.len()];

        for (i, &(temperature_k, pressure_bar)) in conditions.iter().enumerate() {
            for (j, reactants) in reactants_list.iter().enumerate() {
                let result = solve_for_products(
                    reactants,
                    true,  // const_t: TP
                    false, // const_s: unused for TP
                    true,  // const_p: TP
                    temperature_k,
                    pressure_bar,
                    true, // include_condensed, matching mixture.rs's always-built condensed pool
                    true, // include_ions, matching mixture.rs's always-included E in element_pool
                    None,
                );
                pass_fail[i][j] = result.converged;
                if !result.converged {
                    println!(
                        "[debug] FAILED reactants={:?} T={:.1} P={:.1}",
                        reactants
                            .iter()
                            .map(|(n, s)| format!("{:.1}*{}", n, s.data().symbol()))
                            .collect::<Vec<_>>(),
                        temperature_k,
                        pressure_bar
                    );
                }
            }
        }

        println!("========== Summary (equilibrium.rs TP solver) ==========");
        for (i, &(temperature_k, pressure_bar)) in conditions.iter().enumerate() {
            let total = pass_fail[i].len();
            let passed = pass_fail[i].iter().filter(|&&x| x).count();
            let row_header = format!("({:.1} K; {:.1} bar) Results:", temperature_k, pressure_bar);
            println!(
                "{:>30} [{}/{}] {:.2}%",
                row_header,
                passed,
                total,
                100.0 * (passed as f64 / total as f64)
            );
        }
        let total_passed: usize = pass_fail.iter().map(|row| row.iter().filter(|&&x| x).count()).sum();
        let total = conditions.len() * reactants_list.len();
        println!(
            "========== Overall: {}/{} ({:.2}%) ==========",
            total_passed,
            total,
            100.0 * total_passed as f64 / total as f64
        );
    }
}




