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
            out.enthalpy[j] = species[j].data().h_over_r(temperature_k);
            out.entropy[j] = species[j].data().s_over_r(temperature_k);
            out.energy[j] = species[j].data().u_over_r(temperature_k);
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

fn solve_for_products(
    reactants: &[(f64, Species)],
    temperature_k: f64,
    pressure_bar: f64,
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
    // contributes 0 to every sum, though its ln_nj keeps evolving.
    const SIZE: f64 = 18.420681;
    // D.5 tolerances (guide's Global conventions).
    const NJ_TOL: f64 = 0.5e-5;
    const B_TOL: f64 = 1e-6;

    let mut ln_nj = vec![(0.1 / ng as f64).ln(); ns];
    let mut n: f64 = 0.1; // D.6: total gas moles, tracked as its own state var.
    // D.3/D.5's "stored" nj: truncated/effective gas amounts, updated once
    // per iteration and read by check_convergence's tests 1 and 3. Distinct
    // from the RAW exp(ln_nj[j]) that test 4 (element balance) uses instead.
    let mut nj_stored = vec![0.1 / ng as f64; ng];

    // D.7 default budget. D.5 below now returns early on real convergence;
    // this budget only bounds the worst case (never converges).
    const MAX_ITERATIONS: usize = 50;
    let mut assembled_first_iteration: Option<DMatrix<f64>> = None;
    let mut converged = false;
    let mut iterations_used = 0;

    for iteration in 0..MAX_ITERATIONS {
        iterations_used = iteration + 1;
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

        // mu_g[j] (D.1): uses raw ln_nj, not the truncated nj_eff.
        let ln_p_over_n = (pressure_bar / n).ln();
        let mu_g: Vec<f64> = (0..ng)
            .map(|j| {
                mixture_thermo.enthalpy[j] - mixture_thermo.entropy[j] + ln_nj[j] + ln_p_over_n
            })
            .collect();

        // nj_eff[j] (D.4): truncated/effective gas amount used in every sum.
        let nj_eff: Vec<f64> = (0..ng)
            .map(|j| {
                let threshold = n.ln() - SIZE;
                if ln_nj[j] > threshold {
                    ln_nj[j].exp()
                } else {
                    0.0
                }
            })
            .collect();

        // b_delta[k] = b0[k] - Σ_j a_kj * n_j   (gas-only for now; condensed
        // contributions get added in when Section F activates a species).
        let b_delta: Vec<f64> = (0..ne)
            .map(|k| b0[k] - (0..ng).map(|j| stoich[j][k] * nj_eff[j]).sum::<f64>())
            .collect();
        let n_delta = n - nj_eff.iter().sum::<f64>();

        let neq = ne + 1;
        let mut g: DMatrix<f64> = DMatrix::zeros(neq, neq + 1);

        for k in 0..ne {
            // tmp[j] = a_kj * n_j, reused for every column of this row.
            let tmp: Vec<f64> = (0..ng).map(|j| nj_eff[j] * stoich[j][k]).collect();
            for i in 0..ne {
                // Σ_j a_kj * a_ij * n_j
                g[(k, i)] = (0..ng).map(|j| tmp[j] * stoich[j][i]).sum::<f64>();
            }
            g[(k, ne)] = tmp.iter().sum::<f64>(); // Δln n column
            g[(ne, k)] = g[(k, ne)]; // symmetric fill of the moles row
            g[(k, neq)] = b_delta[k] + (0..ng).map(|j| tmp[j] * mu_g[j]).sum::<f64>();
        }
        // Moles row (Eq. 2.26); pi coefficients already filled symmetrically above.
        g[(ne, ne)] = -n_delta;
        g[(ne, neq)] = n_delta + (0..ng).map(|j| nj_eff[j] * mu_g[j]).sum::<f64>();

        // Snapshot the assembled (pre-solve) matrix from the very first pass
        // only — that's the canonical-guess matrix Section D validation
        // checks against. gauss() overwrites g in place with the solution.
        if iteration == 0 {
            assembled_first_iteration = Some(g.clone());
        }

        // Solve for [pi_0..pi_{ne-1}, dln_n], then recover dln_nj (Eq. 2.18).
        gauss(&mut g).expect("singular matrix recovery is Section G, not implemented yet");
        let pi: Vec<f64> = (0..ne).map(|i| g[(i, neq)]).collect();
        let dln_n = g[(ne, neq)];
        let dln_nj: Vec<f64> = (0..ng)
            .map(|j| -mu_g[j] + dln_n + (0..ne).map(|i| stoich[j][i] * pi[i]).sum::<f64>())
            .collect();

        // D.2: cap the step before it's applied. dln_t = 0.0 — no ΔlnT
        // unknown yet for this TP-only driver (Section E).
        let lambda = damping_lambda(&ln_nj[..ng], &dln_nj, dln_n, 0.0, n, SIZE);

        // D.3: apply the damped update. nj_stored's truncation threshold
        // uses the OLD n (not yet updated below) — matches the guide's
        // exact pseudocode order: nj[j] update happens inside the same loop
        // that bumps ln_nj[j], before n itself is reassigned.
        let applied_dln_nj: Vec<f64> = (0..ng).map(|j| lambda * dln_nj[j]).collect();
        for j in 0..ng {
            ln_nj[j] += applied_dln_nj[j];
            let threshold = n.ln() - SIZE;
            nj_stored[j] = if ln_nj[j] > threshold { ln_nj[j].exp() } else { 0.0 };
        }
        let applied_dln_n = lambda * dln_n;
        // const_p is always true for TP — the const_v branch (n = Σ nj)
        // arrives with Section E's TV/UV/SV problem types.
        n = (n.ln() + applied_dln_n).exp();

        // D.5: convergence check. Only the tests that apply to a gas-only,
        // fixed-T, fixed-P problem: no condensed (test 2), always const_p
        // (test 3 applies), const_t (test 5 skipped), no entropy/trace/ions
        // (tests 6-7 skipped).
        let sum_nj: f64 = nj_stored.iter().sum(); // no condensed yet, so gas-only sum

        // Test 1: gas species.
        let test1 = (0..ng).all(|j| nj_stored[j] * applied_dln_nj[j].abs() / sum_nj <= NJ_TOL);

        // Test 3: total moles (const_p only, which TP always is).
        let test3 = n * applied_dln_n.abs() / sum_nj <= NJ_TOL;

        // Test 4: elements, using RAW nj (not the stored/truncated values).
        let b_max = b0.iter().cloned().fold(0.0_f64, f64::max);
        let test4 = (0..ne).all(|k| {
            if b0[k] <= B_TOL {
                return true;
            }
            let b_k: f64 = (0..ng).map(|j| stoich[j][k] * ln_nj[j].exp()).sum();
            (b0[k] - b_k).abs() <= B_TOL * b_max
        });

        if test1 && test3 && test4 {
            converged = true;
            break;
        }
    }

    let nj: Vec<f64> = (0..ng).map(|j| ln_nj[j].exp()).collect();
    TpSolveDebug {
        assembled_first_iteration: assembled_first_iteration
            .expect("MAX_ITERATIONS > 0, so iteration 0 always runs"),
        converged,
        iterations_used,
        nj,
        n,
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

        let result = solve_for_products(&reactants, 3800.0, 1.01325, false, false, Some(&only));
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
