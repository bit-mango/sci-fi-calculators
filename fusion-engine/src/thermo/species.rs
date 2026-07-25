use crate::constants::*;
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::fluid_properties::{TemperatureDependentProperty, ThermoReference};
use crate::thermo::reactions::{get_rxn_enthalpy, get_rxn_entropy};
use crate::thermo::species;
use nalgebra::{DMatrix, DVector, dmatrix, dvector};
use std::collections::HashMap;
use std::f64;

// TODO rename file to mixture.rs and move species out of it into its own file.
#[derive(Default)]
pub struct MixtureState {
    pub products: Mixture,
    pub h_total: f64, // Total enthalpy.
    pub s_total: f64, // Total entropy.
    pub n_total: f64, // Total moles.
    pub avg_mw: f64,
    pub avg_cp: f64,
}

#[derive(Clone, Default)]
pub struct Mixture {
    pub products: Vec<(f64, Species)>,
    pub temperature_k: f64,
    pub pressure_bar: f64,
    pub h_total: f64, // Total enthalpy.
    pub s_total: f64, // Total entropy.
    pub n_total: f64,
    pub avg_mw: f64,
    pub avg_cp: f64,
}

impl Mixture {
    // TODO add a new method that accepts some initial pi and ln_n guess for carry voer work?
    pub fn new(
        tr: &ThermoReference,
        reactants: &Vec<(f64, Species)>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> Self {
        let (species_pool, element_pool) = Self::reaction_pool(&reactants);
        // Precompute mu_not.
        let mu_not = species_pool
            .iter()
            .map(|species| {
                let tdp = tr.get_tdp(&species.symbol());
                let mu_not_i = tdp.h(temperature_k) - temperature_k * tdp.s(temperature_k);
                mu_not_i
            })
            .collect::<Vec<f64>>();

        // Each row is a species, and each column is how many moles of an element are in that species.
        // Example:
        //      H   O
        // H2O  2   1
        // OH   1   1
        // H2   2   0
        // H    1   0
        // O2   0   2
        // O    0   1
        let a = species_pool
            .iter()
            .map(|species| {
                let children = species.constituents();
                element_pool
                    .iter()
                    .map(|(_, elem)| {
                        children
                            .iter()
                            .find(|(_, child)| child == elem)
                            .unwrap_or(&(0.0, *elem))
                            .0
                    })
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<Vec<f64>>>();

        let j_len = a[0].len();
        // b_j is element_pool[j].0
        let b = element_pool
            .iter()
            .map(|(moles, _)| *moles)
            .collect::<Vec<f64>>();
        // Compute products given reactants and temperature_k.
        // Compute guesses.
        let (pi, ln_n) = Self::guess_pi_and_ln_n(
            &reactants,
            temperature_k,
            pressure_bar,
            j_len,
            &a,
            &mu_not,
            &species_pool,
        );
        // Find the anchor, try target temperature first, if that fails iterate up and down from it.
        let products = if let Some(ap) = Self::solve_for_products(
            tr,
            temperature_k,
            pressure_bar,
            &pi,
            ln_n,
            &species_pool,
            &a,
            &b,
        ) {
            ap.0
        } else {
            // Couldn't find anchor products so now recursively iterate to find the anchor, the walk it toward temperature_k
            let step = 100.0;
            let mut t_low = (temperature_k - (temperature_k % step)).max(200.0);
            let mut t_high = t_low + step;
            let mut anchor = None;
            let mut anchor_temperature_k = 0.0;

            while t_low > 200.0 || t_high < 20_000.0 {
                // Try t_low.
                let low_anchor = Self::solve_for_products(
                    tr,
                    t_low,
                    pressure_bar,
                    &pi,
                    ln_n,
                    &species_pool,
                    &a,
                    &b,
                );
                if low_anchor.is_some() {
                    anchor = low_anchor;
                    anchor_temperature_k = t_low;
                    break;
                } else {
                    // Try t_high.
                    let high_anchor = Self::solve_for_products(
                        tr,
                        t_high,
                        pressure_bar,
                        &pi,
                        ln_n,
                        &species_pool,
                        &a,
                        &b,
                    );
                    if high_anchor.is_some() {
                        anchor = high_anchor;
                        anchor_temperature_k = t_high;
                        break;
                    } else {
                        // Neither worked, adjust t_low, and t_high.
                        t_low = (t_low - step).max(200.0);
                        t_high = (t_high + step).min(20_000.0);
                    }
                }
            }
            let mut products = None;
            let mut anchor_iterations = 0;
            if let Some(anc) = anchor {
                // Now walk the anchor some max iterations.
                let mut anc = anc;
                while anchor_iterations < 60 {
                    let next_anchor = Self::walk_anchor(
                        tr,
                        anchor_temperature_k,
                        &anc.1,
                        anc.2,
                        &species_pool,
                        &a,
                        &b,
                        temperature_k,
                        pressure_bar,
                    );
                    if let Some(next) = next_anchor {
                        if next.3 == temperature_k {
                            // We are done!
                            products = Some(next.0);
                        } else {
                            anc.0 = next.0;
                            anc.1 = next.1;
                            anc.2 = next.2;
                            anchor_temperature_k = next.3;
                        }
                    } else {
                        panic!(
                            "Could not converge to answer!
                                Temperature wanted: {:.3} K
                                Temperature walked to: {:.3} K",
                            temperature_k, anchor_temperature_k
                        )
                    }

                    anchor_iterations += 1;
                }
            } else {
                println!(
                    "Failed to find anchor. t_low={:.3} t_high={:.3}\nreactants={:?}",
                    t_low, t_high, reactants
                );
                panic!("Unable to find anchor point!");
            }
            products.expect(&format!(
                "Anchor walk did not converge to {:.3} K within {} iterations (last anchor at {:.3} K)",
                temperature_k, anchor_iterations, anchor_temperature_k
            ))
        };
        // every 1_000 K until we find a solution. This is our anchor
        // Move the anchor closer to temperature_k until we get it backing off half the distance on fail.
        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(tr, &products, temperature_k, pressure_bar);

        Self {
            products,
            temperature_k,
            pressure_bar,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        }
    }

    fn walk_anchor(
        tr: &ThermoReference,
        anchor_temperature_k: f64,
        anchor_pi: &Vec<f64>,
        anchor_ln_n: f64,
        species_pool: &Vec<Species>,
        a: &Vec<Vec<f64>>,
        b: &Vec<f64>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> Option<(Vec<(f64, Species)>, Vec<f64>, f64, f64)> {
        let direction = if temperature_k > anchor_temperature_k {
            1.0
        } else {
            -1.0
        };
        let mut t = anchor_temperature_k;
        let mut pi = anchor_pi.clone();
        let mut ln_n = anchor_ln_n;
        let mut step = 1_000.0;

        while (temperature_k - t).abs() > 1.0e-6 {
            let remaining = temperature_k - t;
            let next_t = if remaining.abs() <= step {
                temperature_k
            } else {
                t + direction * step
            };
            match Self::solve_for_products(tr, next_t, pressure_bar, &pi, ln_n, species_pool, a, b)
            {
                Some((products, next_pi, next_ln_n)) => {
                    t = next_t;
                    pi = next_pi;
                    ln_n = next_ln_n;
                    if t == temperature_k {
                        return Some((products, pi, ln_n, t));
                    }
                }
                None => {
                    step /= 2.0;
                    if step < 1.0e-3 {
                        // Truly stuck
                        return None; // let the outer loops could not converge panic report where we got stuck
                    }
                }
            }
        }
        None
    }

    // Assumes reactants == products
    pub fn new_with_frozen_reactants(
        tr: &ThermoReference,
        reactants: &Vec<(f64, Species)>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> Self {
        let products = reactants.clone();

        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(tr, &products, temperature_k, pressure_bar);

        Self {
            products,
            temperature_k,
            pressure_bar,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        }
    }

    fn state(
        tr: &ThermoReference,
        products: &Vec<(f64, Species)>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> (f64, f64, f64, f64, f64) {
        let n_total: f64 = products.iter().map(|(moles, _)| *moles).sum();
        let products: Vec<(f64, Species, &TemperatureDependentProperty, f64)> = products
            .iter()
            .map(|(moles, species)| {
                (
                    *moles, // number of moles, n
                    *species,
                    tr.get_tdp(&species.symbol()),
                    moles / n_total, // mole fraction, x
                )
            })
            .collect();
        let h_total: f64 = products
            .iter()
            .map(|(n_i, _, tdp_i, _)| n_i * tdp_i.h(temperature_k))
            .sum();
        let s_total: f64 = products
            .iter()
            .map(|(n_i, _, tdp_i, x_i)| {
                n_i * (tdp_i.s(temperature_k)
                    - R * (x_i * pressure_bar / STD_REFERENCE_PRESSURE).ln())
            })
            .sum();
        let avg_mw: f64 = products
            .iter()
            .map(|(_, species_i, _, x_i)| x_i * species_i.mw())
            .sum();
        let avg_cp: f64 = products
            .iter()
            .map(|(_, _, tdp_i, x_i)| x_i * tdp_i.cp(temperature_k))
            .sum();

        (h_total, s_total, n_total, avg_mw, avg_cp)
    }

    pub fn print_products(&self) {
        println!("Products Temperature: {:.3} K", self.temperature_k);
        for (mole, species) in self.products.iter() {
            if *mole < 1.0e-6 {
                continue;
            }
            println!("{:.6} {}", mole, species.symbol());
        }
    }

    pub fn solve_for_target_enthalpy(&self, tr: &ThermoReference, target_enthalpy: f64) -> Mixture {
        // Find a bracket where f(T) = h(T) - starting_h changes sign
        let step = 500.0;
        // Walk down from the hot temperature.
        let mut t_low = self.temperature_k;
        let mut low_scratch = Mixture::new(tr, &self.products, t_low, self.pressure_bar);
        let mut f_low = low_scratch.h_total - target_enthalpy;
        while f_low > 0.0 && t_low > 200.0 {
            let next_t = (t_low - step).max(200.0);
            low_scratch = Mixture::new(tr, &low_scratch.products, next_t, self.pressure_bar);
            t_low = next_t;
            f_low = low_scratch.h_total - target_enthalpy;
        }
        // Walk up from the hot temperature.
        let mut t_high = self.temperature_k;
        let mut high_scratch = Mixture::new(tr, &self.products, t_high, self.pressure_bar);
        let mut f_high = high_scratch.h_total - target_enthalpy;
        while f_high < 0.0 && t_high < 20_000.0 {
            let next_t = (t_high + step).min(20_000.0);
            high_scratch = Mixture::new(tr, &high_scratch.products, next_t, self.pressure_bar);
            t_high = next_t;
            f_high = high_scratch.h_total - target_enthalpy;
        }

        if f_low > 0.0 || f_high < 0.0 {
            println!("Temperature Range: {:.3} <-> {:.3}", t_low, t_high);
            panic!("Could not bracket adiabatic mix temperature");
        }

        // Now determine actual new mix
        let mut t = t_low + (t_high - t_low) / 2.0;
        let mut iterations = 0;
        let mut final_products = None;
        while iterations < 100 {
            if (t_high - t_low).abs() < 1.0e-6 {
                final_products = Some(low_scratch);
                break;
            }
            let seed = if (t - t_low).abs() <= (t_high - t).abs() {
                &low_scratch.products
            } else {
                &high_scratch.products
            };
            let mid_scratch = Mixture::new(tr, seed, t, self.pressure_bar);
            let f_mid = mid_scratch.h_total - target_enthalpy;

            if f_mid < 0.0 {
                // Raise temperature
                t_low = t;
                low_scratch = mid_scratch;
            } else {
                // Lower temperature
                t_high = t;
                high_scratch = mid_scratch;
            }
            t = t_low + (t_high - t_low) / 2.0;
            iterations += 1;
        }

        final_products.expect("Failed to mix products")
    }

    pub fn mix(&self, tr: &ThermoReference, other: &Self) -> Self {
        let starting_h = self.h_total + other.h_total;
        if self.pressure_bar != other.pressure_bar {
            panic!("Must be Isobaric to mix!");
        }
        let pressure_bar = self.pressure_bar;

        // Before mixing, bring the colder mix up to the same temperature as the hotter mix
        // for better initial guesses.
        let (cold, hot) = if self.temperature_k <= other.temperature_k {
            (self, other)
        } else {
            (other, self)
        };

        let cold_adj = if self.temperature_k == other.temperature_k {
            // Already same temperature nothing to do.
            cold
        } else {
            // Calculate new cold mixture products
            &Mixture::new(tr, &cold.products, hot.temperature_k, pressure_bar)
        };

        // Create reactant mix using adjusted other.
        let mut reatants_mix: HashMap<String, (f64, Species)> = HashMap::new();
        // Add all of original reatants to mix.
        for s in hot.products.iter() {
            let key = s.1.symbol();
            reatants_mix.insert(key, s.clone());
        }
        for o in cold_adj.products.iter() {
            let key = o.1.symbol();
            if let Some(entry) = reatants_mix.get_mut(&key) {
                // reatants already exists! Increment moles.
                entry.0 += o.0;
            } else {
                // reatants is new, add them.
                reatants_mix.insert(key, o.clone());
            }
        }
        let mut reactants = reatants_mix
            .drain()
            .map(|(_, v)| v)
            .collect::<Vec<(f64, Species)>>();
        reactants.sort_by(|a, b| a.1.symbol().cmp(&b.1.symbol()));

        let combined = Mixture::new(tr, &reactants, hot.temperature_k, pressure_bar);
        combined.solve_for_target_enthalpy(tr, starting_h)
    }

    pub fn scale(&self, tr: &ThermoReference, factor: f64) -> Self {
        let scaled_products = self
            .products
            .iter()
            .map(|(moles, species)| (*moles * factor, *species))
            .collect();
        Self::new_with_frozen_reactants(tr, &scaled_products, self.temperature_k, self.pressure_bar)
    }

    fn reaction_pool(products: &Vec<(f64, Species)>) -> (Vec<Species>, Vec<(f64, Species)>) {
        let mut element_pool: HashMap<Species, f64> = HashMap::new();
        products.iter().for_each(|(parent_moles, parent_species)| {
            let children = parent_species.constituents();
            children.iter().for_each(|(child_moles, child_species)| {
                if let Some(existing) = element_pool.get_mut(child_species) {
                    *existing += parent_moles * child_moles;
                } else {
                    element_pool.insert(*child_species, parent_moles * child_moles);
                }
            });
        });
        // Add electrons to the element pool.
        element_pool.entry(Species::E).or_insert(0.0);
        // Filter out any species we don't have constituent elements for.
        let species_pool = Species::all()
            .iter()
            .filter(|species| {
                let children = species.constituents();
                for child in children {
                    if !element_pool.contains_key(&child.1) {
                        return false;
                    }
                }
                return true;
            })
            .map(|s| *s)
            .collect::<Vec<Species>>();
        let mut element_pool = element_pool
            .drain()
            .map(|(species, moles)| (moles, species))
            .collect::<Vec<(f64, Species)>>();
        element_pool.sort_by(|a, b| a.1.symbol().cmp(&b.1.symbol()));
        (species_pool, element_pool)
    }

    pub fn solve_for_products(
        tr: &ThermoReference,
        temperature_k: f64,
        pressure_bar: f64,
        pi: &Vec<f64>,
        ln_n: f64,
        species_pool: &Vec<Species>,
        a: &Vec<Vec<f64>>,
        b: &Vec<f64>,
    ) -> Option<(Vec<(f64, Species)>, Vec<f64>, f64)> {
        let mut pi = pi.clone();
        let mut ln_n = ln_n;
        // Precompute mu_not.
        let mu_not = species_pool
            .iter()
            .map(|species| {
                let tdp = tr.get_tdp(&species.symbol());
                let mu_not_i = tdp.h(temperature_k) - temperature_k * tdp.s(temperature_k);
                mu_not_i
            })
            .collect::<Vec<f64>>();

        // n_i = n * c_i * exp(Σⱼ πⱼ·a[i][j])
        // i is in range form 0..species_pool.length
        // j is in range form 0..element_pool.length
        let n = |i: usize, pi: &[f64], ln_n: f64| {
            let sum = pi
                .iter()
                .zip(a[i].iter())
                .map(|(pi_i, a_i_j)| pi_i * a_i_j)
                .sum::<f64>();
            ln_n.exp()
                * (STD_REFERENCE_PRESSURE / pressure_bar)
                * (sum - mu_not[i] / (R * temperature_k)).exp()
        };
        let j_len = a[0].len();

        let mut iterations = 0;
        let mut final_n_vec = vec![];
        let mut lambda = 1.0e-3;

        loop {
            let (n_vec, f, j) = residual_and_jacobian(&pi, ln_n, &a, &b, &n);
            if f.norm() < 1.0e-6 {
                final_n_vec = n_vec;
                break;
            }

            if iterations == 4999 {
                return None;
            }

            let current_norm = f.norm();
            // Solve via normal equations (JᵀJ + λ·diag(JᵀJ))·δ = -Jᵀf rather than
            // damping J directly. J here is symmetric but not necessarily positive
            // definite away from the solution (e.g. the last diagonal entry
            // n_sum - exp(ln_n) can go negative), so a diagonal shift on J itself
            // is not guaranteed to produce a descent direction for ‖f‖ — that's
            // why the line search could exhaust all 50 lambda-doublings and find
            // nothing. JᵀJ + λD is always positive definite for λ>0, so -Jᵀf-based
            // steps are always a descent direction of the merit function ½‖f‖².
            let jt = j.transpose();
            let jtj = &jt * &j;
            let neg_jtf = -(&jt * &f);
            let mut accepted = false;
            for _ in 0..50 {
                let mut jtj_damped = jtj.clone();
                for k in 0..(j_len + 1) {
                    jtj_damped[(k, k)] += lambda * jtj[(k, k)].max(1.0e-10);
                }
                if let Some(delta) = jtj_damped.clone().lu().solve(&neg_jtf) {
                    // Trust-region step limiter (NASA CEA-style): since pi/ln_n
                    // feed an exp(), an unbounded Newton step can send some
                    // species' moles to ~0 and others to overflow in one shot,
                    // which then singularizes J on the next iteration. Cap the
                    // largest per-step change so no n_i can move by more than
                    // ~e^2 in a single iteration.
                    let max_step = delta.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                    let alpha = if max_step > 25.0 {
                        25.0 / max_step
                    } else {
                        1.0
                    };

                    let trial_pi: Vec<f64> = (0..j_len).map(|k| pi[k] + alpha * delta[k]).collect();
                    let trial_ln_n = ln_n + alpha * delta[j_len];
                    let (_, trial_f, _) = residual_and_jacobian(&trial_pi, trial_ln_n, &a, &b, &n);
                    if trial_f.iter().all(|v| v.is_finite()) && trial_f.norm() < current_norm {
                        pi = trial_pi;
                        ln_n = trial_ln_n;
                        lambda *= 0.5;
                        accepted = true;
                        break;
                    }
                    lambda *= 2.0;
                }
            }
            if !accepted {
                return None;
            }

            iterations += 1;
        }

        let products = final_n_vec
            .iter()
            .zip(species_pool.iter())
            .map(|(&moles, &species)| (moles, species))
            .collect::<Vec<(f64, Species)>>();

        Some((products, pi, ln_n))
    }

    fn guess_pi_and_ln_n(
        products: &Vec<(f64, Species)>,
        temperature_k: f64,
        pressure_bar: f64,
        j_len: usize,
        a: &Vec<Vec<f64>>,
        mu_not: &Vec<f64>,
        species_pool: &Vec<Species>,
    ) -> (Vec<f64>, f64) {
        // Base ln_n_guess off of starting moles of the system.
        let ln_n_guess = products.iter().map(|(moles, _)| *moles).sum::<f64>().ln();

        let total_moles: f64 = products.iter().map(|(moles, _)| *moles).sum();
        let threshold = 1.0e-9; // relative to total system moles.

        // Identify anchor species.
        let anchors: Vec<(usize, f64)> = species_pool
            .iter()
            .enumerate() // TODO convert to filter_map?
            .map(|(i, species)| {
                if species.is_charged() {
                    return None;
                }
                let mut anchors_i = None;
                for (moles, reactant) in products.iter() {
                    if reactant == species && *moles / total_moles > threshold {
                        anchors_i = Some((i, *moles));
                    }
                }
                anchors_i
            })
            .filter_map(|anchors_i| anchors_i)
            .collect();

        // Build A (anchors.len() × j_len) and rhs (anchors.len()) from the anchor list:
        //   A[row] = a[i]                                    for (i, feed_moles) in anchors
        //   rhs[row] = ln(feed_moles) - ln_n_guess + mu_not[i]/(R*T) - ln(P_ref/P)

        // Solve (AᵀA + eps·I) · pi = Aᵀ·rhs   instead of a plain square solve.
        let mut a_ = vec![vec![0.0; j_len]; anchors.len()];
        let mut rhs = vec![0.0; anchors.len()];
        for (row, (i, feed_moles)) in anchors.iter().enumerate() {
            for j in 0..j_len {
                a_[row][j] = a[*i][j];
                rhs[row] = feed_moles.ln() - ln_n_guess + mu_not[*i] / (R * temperature_k)
                    - (STD_REFERENCE_PRESSURE / pressure_bar).ln();
            }
        }

        let rows = anchors.len();
        let flattened: Vec<f64> = a_.into_iter().flatten().collect();
        let a_matrix = DMatrix::from_row_slice(rows, j_len, &flattened);
        let rhs_vector = DVector::from_vec(rhs);

        let ata = a_matrix.transpose() * &a_matrix; // j_len * j_len, symmeric
        let atb = a_matrix.transpose() * &rhs_vector; // j_len

        // ridge term, scaled to A transpose A's magnitude rather than a fixed constant
        let eps = 1.0e-8 * ata.amax().max(1.0e-10);
        let mut ata_reg = ata.clone();
        for k in 0..j_len {
            ata_reg[(k, k)] += eps;
        }

        let pi_guess: Vec<f64> = match ata_reg.lu().solve(&atb) {
            Some(pi) => pi.iter().copied().collect(),
            None => vec![0.0; j_len], // degenerate anchors (e.g. duplicate rows) - fall back
        };

        (pi_guess, ln_n_guess)
    }

    pub fn feed_mass(&self) -> f64 {
        self.products
            .iter()
            .map(|(moles, specie)| moles * specie.mw())
            .sum()
    }
}

fn residual_and_jacobian(
    pi: &[f64],
    ln_n: f64,
    a: &Vec<Vec<f64>>,
    b: &Vec<f64>,
    n_i: &impl Fn(usize, &[f64], f64) -> f64,
) -> (Vec<f64>, DVector<f64>, DMatrix<f64>) {
    let i_len = a.len();
    let j_len = a[0].len();
    let mut n_vec = vec![0.0; i_len];
    for i in 0..i_len {
        n_vec[i] = n_i(i, pi, ln_n);
    }

    let mut f_vec = vec![0.0; j_len + 1];
    for j in 0..j_len {
        let mut sum = 0.0;
        for i in 0..i_len {
            sum += n_vec[i] * a[i][j];
        }
        sum -= b[j];
        f_vec[j] = sum;
    }
    let n_sum = n_vec.iter().sum::<f64>();
    f_vec[j_len] = n_sum - ln_n.exp();
    let mut j_vec = vec![vec![0.0; j_len + 1]; j_len + 1];
    for j in 0..j_len {
        for k in 0..j_len {
            let mut sum = 0.0;
            for i in 0..i_len {
                sum += n_vec[i] * a[i][j] * a[i][k];
            }
            j_vec[j][k] = sum;
        }
    }

    for j in 0..j_len {
        let mut sum = 0.0;
        for i in 0..i_len {
            sum += n_vec[i] * a[i][j];
        }
        j_vec[j][j_len] = sum;
        j_vec[j_len][j] = sum;
    }

    j_vec[j_len][j_len] = n_sum - ln_n.exp();
    let f_vec = DVector::from_vec(f_vec);
    let flattened_j = j_vec.iter().flatten().map(|e| *e).collect();
    // TODO from_vec expects column major data, but flatten gives row major data. It is fine because j_vec is symmetrical ie row count == col count
    // but thats messy lets fix it.
    let j_matrix = DMatrix::from_vec(j_len + 1, j_len + 1, flattened_j);
    (n_vec, f_vec, j_matrix)
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub enum Species {
    H,
    H2,
    C,
    O,
    CO,
    N2,
    N,
    H2O,
    OH,
    CO2,
    O2,
    CH4,
    NH3,
    E,
    HPlus,
    CPlus,
    OPlus,
    NPlus,
}

impl Species {
    pub fn all() -> Vec<Self> {
        vec![
            Self::H,
            Self::H2,
            Self::C,
            Self::O,
            Self::CO,
            Self::N2,
            Self::N,
            Self::H2O,
            Self::OH,
            Self::CO2,
            Self::O2,
            Self::CH4,
            Self::NH3,
            Self::E,
            Self::HPlus,
            Self::CPlus,
            Self::OPlus,
            Self::NPlus,
        ]
    }
    pub fn symbol(&self) -> String {
        match self {
            Species::H => "H".to_string(),
            Species::H2 => "H2".to_string(),
            Species::C => "C".to_string(),
            Species::O => "O".to_string(),
            Species::O2 => "O2".to_string(),
            Species::CO => "CO".to_string(),
            Species::CO2 => "CO2".to_string(),
            Species::N2 => "N2".to_string(),
            Species::N => "N".to_string(),
            Species::H2O => "H2O".to_string(),
            Species::OH => "OH".to_string(),
            Species::CH4 => "CH4".to_string(),
            Species::NH3 => "NH3".to_string(),
            Species::E => "e-".to_string(),
            Species::HPlus => "H+".to_string(),
            Species::CPlus => "C+".to_string(),
            Species::OPlus => "O+".to_string(),
            Species::NPlus => "N+".to_string(),
        }
    }

    pub fn mw(&self) -> f64 {
        match self {
            Species::H => H_MW,
            Species::H2 => H2_MW,
            Species::C => C_MW,
            Species::O => O_MW,
            Species::O2 => O2_MW,
            Species::CO => CO_MW,
            Species::CO2 => CO2_MW,
            Species::N2 => N2_MW,
            Species::N => N_MW,
            Species::H2O => H2O_MW,
            Species::OH => OH_MW,
            Species::CH4 => CH4_MW,
            Species::NH3 => NH3_MW,
            Species::E => 0.0,
            Species::HPlus => H_MW,
            Species::CPlus => C_MW,
            Species::OPlus => O_MW,
            Species::NPlus => N_MW,
        }
    }

    pub fn constituents(&self) -> Vec<(f64, Self)> {
        match self {
            Species::H => vec![(1.0, Species::H)],
            Species::H2 => vec![(2.0, Species::H)],
            Species::C => vec![(1.0, Species::C)],
            Species::O => vec![(1.0, Species::O)],
            Species::O2 => vec![(2.0, Species::O)],
            Species::CO => vec![(1.0, Species::C), (1.0, Species::O)],
            Species::CO2 => vec![(1.0, Species::C), (2.0, Species::O)],
            Species::N2 => vec![(2.0, Species::N)],
            Species::N => vec![(1.0, Species::N)],
            Species::H2O => vec![(2.0, Species::H), (1.0, Species::O)],
            Species::OH => vec![(1.0, Species::O), (1.0, Species::H)],
            Species::CH4 => vec![(1.0, Species::C), (4.0, Species::H)],
            Species::NH3 => vec![(1.0, Species::N), (3.0, Species::H)],
            Species::E => vec![(1.0, Species::E)],
            Species::HPlus => vec![(1.0, Species::H), (-1.0, Species::E)],
            Species::CPlus => vec![(1.0, Species::C), (-1.0, Species::E)],
            Species::OPlus => vec![(1.0, Species::O), (-1.0, Species::E)],
            Species::NPlus => vec![(1.0, Species::N), (-1.0, Species::E)],
        }
    }
    pub fn is_charged(&self) -> bool {
        matches!(
            self,
            Species::E | Species::HPlus | Species::CPlus | Species::OPlus | Species::NPlus
        )
    }
}
