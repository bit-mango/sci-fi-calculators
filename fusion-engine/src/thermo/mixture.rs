use crate::constants::*;
use crate::thermo::fluid_properties::{TemperatureDependentProperty, ThermoReference};
use crate::thermo::species;
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;
use std::f64;
use thermo_species::{AnySpeciesData, Constituent, Species, species};

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
    pub fn new(reactants: &Vec<(f64, Species)>, temperature_k: f64, pressure_bar: f64) -> Self {
        let (gas_species_pool, condensed_species_pool, element_pool) =
            Mixture::reaction_pool(reactants);
        let mu_not = gas_species_pool
            .iter()
            .map(|species| {
                let mu_not_i = species.data().h(temperature_k)
                    - temperature_k * species.data().s(temperature_k);
                mu_not_i
            })
            .collect::<Vec<f64>>();
        let a = element_pool
            .iter()
            .map(|(_, element)| {
                gas_species_pool
                    .iter()
                    .map(|species| {
                        let children = species.data().constituents();
                        children
                            .iter()
                            .find(|(_, child)| child == element)
                            .unwrap_or(&(0.0, *element))
                            .0
                    })
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<Vec<f64>>>();
        let b = element_pool
            .iter()
            .map(|(moles, _)| *moles)
            .collect::<Vec<f64>>();

        let result = Mixture::solve_for_products(
            temperature_k,
            pressure_bar,
            &gas_species_pool,
            &condensed_species_pool,
            &element_pool,
            &a,
            &b,
            &mu_not,
        );

        let products = result.expect("Failed to converge.");

        // every 1_000 K until we find a solution. This is our anchor
        // Move the anchor closer to temperature_k until we get it backing off half the distance on fail.
        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(&products, temperature_k, pressure_bar);

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

    // Assumes reactants == products
    pub fn new_with_frozen_reactants(
        tr: &ThermoReference,
        reactants: &Vec<(f64, Species)>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> Self {
        let products = reactants.clone();

        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(&products, temperature_k, pressure_bar);

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
        products: &Vec<(f64, Species)>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> (f64, f64, f64, f64, f64) {
        let n_total: f64 = products.iter().map(|(moles, _)| *moles).sum();
        let products: Vec<(f64, Species, f64)> = products
            .iter()
            .map(|(moles, species)| {
                (
                    *moles, // number of moles, n
                    *species,
                    moles / n_total, // mole fraction, x
                )
            })
            .collect();
        let h_total: f64 = products
            .iter()
            .map(|(n_i, species_i, _)| n_i * species_i.data().h(temperature_k))
            .sum();
        let s_total: f64 = products
            .iter()
            .map(|(n_i, species_i, x_i)| {
                n_i * (species_i.data().s(temperature_k)
                    - R * (x_i * pressure_bar / STD_REFERENCE_PRESSURE).ln())
            })
            .sum();
        let avg_mw: f64 = products
            .iter()
            .map(|(_, species_i, x_i)| x_i * species_i.data().mw())
            .sum();
        let avg_cp: f64 = products
            .iter()
            .map(|(_, species_i, x_i)| x_i * species_i.data().cp(temperature_k))
            .sum();

        (h_total, s_total, n_total, avg_mw, avg_cp)
    }

    pub fn print_products(&self) {
        println!("Products Temperature: {:.3} K", self.temperature_k);
        for (mole, species) in self.products.iter() {
            if *mole < 1.0e-6 {
                continue;
            }
            println!("{:.6} {}", mole, species.data().symbol());
        }
    }

    pub fn solve_for_target_enthalpy(&self, tr: &ThermoReference, target_enthalpy: f64) -> Mixture {
        // Find a bracket where f(T) = h(T) - starting_h changes sign
        let step = 500.0;
        // Walk down from the hot temperature.
        let mut t_low = self.temperature_k;
        let mut low_scratch = Mixture::new(&self.products, t_low, self.pressure_bar);
        let mut f_low = low_scratch.h_total - target_enthalpy;
        while f_low > 0.0 && t_low > 200.0 {
            let next_t = (t_low - step).max(200.0);
            low_scratch = Mixture::new(&low_scratch.products, next_t, self.pressure_bar);
            t_low = next_t;
            f_low = low_scratch.h_total - target_enthalpy;
        }
        // Walk up from the hot temperature.
        let mut t_high = self.temperature_k;
        let mut high_scratch = Mixture::new(&self.products, t_high, self.pressure_bar);
        let mut f_high = high_scratch.h_total - target_enthalpy;
        while f_high < 0.0 && t_high < 20_000.0 {
            let next_t = (t_high + step).min(20_000.0);
            high_scratch = Mixture::new(&high_scratch.products, next_t, self.pressure_bar);
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
            let mid_scratch = Mixture::new(seed, t, self.pressure_bar);
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
            &Mixture::new(&cold.products, hot.temperature_k, pressure_bar)
        };

        // Create reactant mix using adjusted other.
        let mut reatants_mix: HashMap<String, (f64, Species)> = HashMap::new();
        // Add all of original reatants to mix.
        for s in hot.products.iter() {
            let key = s.1.data().symbol();
            reatants_mix.insert(key.to_string(), s.clone());
        }
        for o in cold_adj.products.iter() {
            let key = o.1.data().symbol();
            if let Some(entry) = reatants_mix.get_mut(key) {
                // reatants already exists! Increment moles.
                entry.0 += o.0;
            } else {
                // reatants is new, add them.
                reatants_mix.insert(key.to_string(), o.clone());
            }
        }
        let mut reactants = reatants_mix
            .drain()
            .map(|(_, v)| v)
            .collect::<Vec<(f64, Species)>>();
        reactants.sort_by(|a, b| a.1.cmp(&b.1));

        let combined = Mixture::new(&reactants, hot.temperature_k, pressure_bar);
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

    fn reaction_pool(
        products: &Vec<(f64, Species)>,
    ) -> (Vec<Species>, Vec<Species>, Vec<(f64, Constituent)>) {
        let mut element_pool: HashMap<Constituent, f64> = HashMap::new();
        products.iter().for_each(|(parent_moles, parent_species)| {
            let children = parent_species.data().constituents();
            children.iter().for_each(|(child_moles, child_species)| {
                if let Some(existing) = element_pool.get_mut(child_species) {
                    *existing += parent_moles * child_moles;
                } else {
                    element_pool.insert(*child_species, parent_moles * child_moles);
                }
            });
        });
        // Add electrons to the element pool.
        element_pool.entry(Constituent::E).or_insert(0.0);
        // Filter out any species we don't have constituent elements for,
        // and condensed phase species.
        let gas_species_pool = Species::all()
            .iter()
            .filter(|species| {
                if species.data().phase() > 0 {
                    return false;
                }
                let children = species.data().constituents();
                for child in children {
                    if !element_pool.contains_key(&child.1) {
                        return false;
                    }
                }
                return true;
            })
            .map(|s| *s)
            .collect::<Vec<Species>>();
        let condensed_species_pool = Species::all()
            .iter()
            .filter(|species| {
                // Skip gas species.
                if species.data().phase() == 0 {
                    return false;
                }
                let children = species.data().constituents();
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
            .collect::<Vec<(f64, Constituent)>>();
        element_pool.sort_by(|a, b| a.1.cmp(&b.1));
        (gas_species_pool, condensed_species_pool, element_pool)
    }

    fn solve_for_products(
        temperature_k: f64,
        pressure_bar: f64,
        gas_species_pool: &Vec<Species>,       // gas species pool
        condensed_species_pool: &Vec<Species>, // condensed species pool
        n_lm: &Vec<(f64, Constituent)>,        // element pool
        a: &Vec<Vec<f64>>,                     // moles of element i in species j
        b: &Vec<f64>,                          // element totals, n_lm.0
        mu_not: &Vec<f64>,
    ) -> Option<Vec<(f64, Species)>> {
        let mut active_condensed_species: Vec<Species> = vec![];
        let mut rejected_condensed_species: Vec<Species> = vec![];
        let mut n_c: Vec<f64> = vec![]; // Linear moles for condensed phases
        let mut iterations = 0;
        let mut outer_iterations = 0;
        'outer: loop {
            outer_iterations += 1;
            println!(
                "[debug] outer_iterations={} active_condensed={:?} n_c={:?}",
                outer_iterations,
                active_condensed_species
                    .iter()
                    .map(|s| s.data().symbol())
                    .collect::<Vec<_>>(),
                n_c
            );
            if outer_iterations > 50 {
                println!("[debug] bailing after 50 outer iterations");
                break 'outer None;
            }
            // Initialize
            let mut enn: f64 = 0.1;
            let mut ln_n = enn.ln();
            let mut enln = vec![(enn / gas_species_pool.len() as f64).ln(); gas_species_pool.len()];
            let m = n_lm.len();
            let mut converged_pi: Vec<f64> = vec![0.0; m];
            let mut trace_damped: Vec<bool> = vec![false; gas_species_pool.len()];
            'newton: loop {
                iterations += 1;

                let c = active_condensed_species.len();
                let matrix_size = m + 1 + c;

                // Build dynamically sized G, RHS.
                let mut g = DMatrix::zeros(matrix_size, matrix_size);
                let mut rhs = DVector::zeros(matrix_size);

                let tm = (pressure_bar / enn).ln();
                let en: Vec<f64> = enln.iter().map(|x| x.exp()).collect();
                let mu: Vec<f64> = (0..gas_species_pool.len())
                    .map(|j| mu_not[j] / (R * temperature_k) + enln[j] + tm)
                    .collect();
                let sumn: f64 = en.iter().sum();

                (0..m).for_each(|i| {
                    (0..m).for_each(|k| {
                        g[(i, k)] = (0..gas_species_pool.len())
                            .map(|j| a[i][j] * a[k][j] * en[j])
                            .sum::<f64>();
                    });
                    let sum = (0..gas_species_pool.len()).map(|j| a[i][j] * en[j]).sum();
                    g[(i, m)] = sum;
                    g[(m, i)] = sum;
                });
                g[(m, m)] = sumn - enn;

                (0..m).for_each(|i| {
                    let left: f64 = (0..gas_species_pool.len())
                        .map(|j| a[i][j] * en[j] * mu[j])
                        .sum();
                    let right: f64 = (0..gas_species_pool.len()).map(|j| a[i][j] * en[j]).sum();
                    rhs[i] = left + b[i] - right;
                });
                rhs[m] = (0..gas_species_pool.len())
                    .map(|j| en[j] * mu[j])
                    .sum::<f64>()
                    + (enn - sumn);

                // Build condensed species ros/cols.
                for (idx, cond_species) in active_condensed_species.iter().enumerate() {
                    let row = m + 1 + idx;
                    let constituents = cond_species.data().constituents();
                    for i in 0..m {
                        // Atoms of element i in condensed.
                        let a_ic = constituents
                            .iter()
                            .find(|cond| cond.1 == n_lm[i].1)
                            .map(|c| c.0)
                            .unwrap_or(0.0);
                        g[(row, i)] = a_ic;
                        g[(i, row)] = a_ic;
                        rhs[i] -= a_ic * n_c[idx];
                    }
                    // RHS for condensed is mu_c_not / RT
                    let mu_not_cond = cond_species.data().h(temperature_k)
                        - temperature_k * cond_species.data().s(temperature_k);
                    rhs[row] = mu_not_cond / (R * temperature_k);
                }

                // Solve G*X = RHS using LU decomposition.
                let x = g.lu().solve(&rhs)?;
                for i in 0..m {
                    converged_pi[i] = x[i];
                }

                // Calculate Deln.
                let dln_n = x[m];
                let mut deln_gas: Vec<f64> = (0..gas_species_pool.len())
                    .map(|j| {
                        let sum_a_pi: f64 = (0..m).map(|i| a[i][j] * x[i]).sum();
                        -mu[j] + sum_a_pi + dln_n
                    })
                    .collect();

                let mut dn_cond: Vec<f64> = ((m + 1)..matrix_size).map(|i| x[i]).collect();

                let mut ambda: f64 = 1.0;
                let mut threshold = 5.0 * dln_n.abs();

                // -9.21 corresponds to a mole fraction of ~1e-4.
                const TRACE_LN_THRESHOLD: f64 = -9.21;
                const TRACE_LN_THRESHOLD_MARGIN: f64 = 1.0; // hysteresis band
                const FLOOR_LN_THRESHOLD: f64 = -40.0;

                // Step-size control for trace species.
                for j in 0..gas_species_pool.len() {
                    if deln_gas[j] > 0.0 {
                        let should_damp = if trace_damped[j] {
                            enln[j] < ln_n + TRACE_LN_THRESHOLD + TRACE_LN_THRESHOLD_MARGIN
                        } else {
                            enln[j] < ln_n + TRACE_LN_THRESHOLD
                        };
                        trace_damped[j] = should_damp;

                        if should_damp {
                            let gap = (deln_gas[j] - dln_n).abs();
                            if gap > 1.0e-8 {
                                let candidate = (TRACE_LN_THRESHOLD - enln[j] + ln_n).abs() / gap;
                                ambda = ambda.min(candidate);
                            }
                        } else if deln_gas[j] > threshold {
                            threshold = deln_gas[j];
                        }
                    } else if deln_gas[j] < 0.0 {
                        trace_damped[j] = false;
                        // Species is shrinking. If it's already far below the floor,
                        // stop letting Newton push it any further down. Clamp the
                        // per species step instead of touching the global ambda.
                        if enln[j] < ln_n + FLOOR_LN_THRESHOLD {
                            deln_gas[j] = 0.0; // freeze this species update this iteration
                        }
                    }
                }

                // Step-size control for condensed species: don't let ambda push any
                // active condensed species' moles below zero in a single step.
                for idx in 0..c {
                    if dn_cond[idx] < 0.0 {
                        let candidate = -n_c[idx] / dn_cond[idx];
                        ambda = ambda.min(candidate * 0.99);
                    }
                }

                if threshold > 2.0 {
                    ambda = ambda.min(2.0 / threshold);
                }

                // Update state variables using damped step.
                ln_n += ambda * dln_n;
                enn = ln_n.exp();
                for idx in 0..c {
                    n_c[idx] += ambda * dn_cond[idx];
                }
                if c > 0 {
                    println!(
                        "[debug]   newton_iter={} ambda={:.4e} dn_cond={:?} n_c={:?}",
                        iterations, ambda, dn_cond, n_c
                    );
                }

                // Phase removal check, now checking the post-damping amount, so a
                // species is only removed once it's actually converged to ~0, not
                // because a single raw Newton step overshot past zero.
                let mut phase_removed = false;
                for idx in 0..c {
                    if n_c[idx] < 1.0e-8 {
                        println!(
                            "[debug] removing {} at newton_iter={} n_c={:.4e}",
                            active_condensed_species[idx].data().symbol(),
                            iterations,
                            n_c[idx]
                        );
                        let rejected = active_condensed_species.remove(idx);
                        rejected_condensed_species.push(rejected);
                        n_c.remove(idx);
                        phase_removed = true;
                        break;
                    }
                }
                if phase_removed {
                    iterations = 0;
                    continue 'newton;
                }

                let mut max_gas_update: f64 = 0.0;
                for j in 0..gas_species_pool.len() {
                    let step = ambda * deln_gas[j];
                    enln[j] += step;
                    max_gas_update = max_gas_update.max(step.abs());
                }

                // TODO do I need to do this mid flight check??

                // Convergence check.
                if max_gas_update < 1.0e-5 && (ambda * dln_n).abs() < 1.0e-5 {
                    break 'newton;
                }
                if c > 0 && iterations > 4990 {
                    println!(
                        "[debug]   near-cap iter={} max_gas_update={:.6e} ambda*dln_n={:.6e}",
                        iterations,
                        max_gas_update,
                        (ambda * dln_n).abs()
                    );
                }

                if iterations == 5_000 {
                    break 'outer None;
                }
            }
            // Check all inactive condensed species to guarantee global minimum gibbs.
            let mut most_eager_species: Option<Species> = None;
            let mut max_supersaturation = 0.0;
            for candidate in condensed_species_pool.iter() {
                if active_condensed_species.contains(candidate) {
                    // Skip candidates that are already active.
                    continue;
                }
                if rejected_condensed_species.contains(candidate) {
                    // Skip candidates that were already rejected.
                    continue;
                }
                let valid_temperature_range = candidate.data().valid_temperature_range();
                if temperature_k < valid_temperature_range.0
                    || temperature_k > valid_temperature_range.1
                {
                    // Skip candidate if outside its valid temperature range.
                    continue;
                }
                let constituents = candidate.data().constituents();
                let sum_a_pi: f64 = (0..m)
                    .map(|i| {
                        let ac_i = constituents
                            .iter()
                            .find(|c| c.1 == n_lm[i].1)
                            .map(|c| c.0)
                            .unwrap_or(0.0);
                        ac_i * converged_pi[i]
                    })
                    .sum();
                let mu_not_cond = candidate.data().h(temperature_k)
                    - temperature_k * candidate.data().s(temperature_k);
                let mu_c_rt = mu_not_cond / (R * temperature_k);
                let supersaturation = sum_a_pi - mu_c_rt;
                if supersaturation > 0.0 {
                    println!(
                        "[debug] candidate={} mu_c_rt={:.3e} sum_a_pi={:.3e} supersaturation={:.3e}",
                        candidate.data().symbol(),
                        mu_c_rt,
                        sum_a_pi,
                        supersaturation
                    );
                }
                if supersaturation > max_supersaturation {
                    max_supersaturation = supersaturation;
                    most_eager_species = Some(*candidate);
                }
            }
            if let Some(species) = most_eager_species {
                println!(
                    "[debug] inserting {} (T={:.1}K) supersaturation={:.3e}",
                    species.data().symbol(),
                    temperature_k,
                    max_supersaturation
                );
                active_condensed_species.push(species);
                n_c.push(0.001);
                iterations = 0;
                continue 'outer;
            }
            // Return the converged product moles
            let result: Vec<(f64, Species)> = enln
                .iter()
                .zip(gas_species_pool.iter())
                .map(|(log_n, species)| (log_n.exp(), species.clone()))
                .collect();
            break 'outer Some(result); // Converged!
        }
    }

    pub fn feed_mass(&self) -> f64 {
        self.products
            .iter()
            .map(|(moles, specie)| moles * specie.data().mw())
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_stability_suite() {
        let reactants = vec![
            vec![(1.0, Species::H2O)],
            vec![(1.0, Species::CH4)],
            vec![(1.0, Species::CH4), (2.0, Species::O2)],
            vec![(1.0, Species::CH4), (1.0, Species::O)],
            vec![
                (1.0, Species::CH4),
                (1.0, Species::HPlus),
                (1.0, Species::OHNeg),
            ],
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
        ];
        let conditions = vec![
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
        ];
        let mut pass_fail = vec![vec![false; reactants.len()]; conditions.len()];
        for (i, (temperature_k, pressure_bar)) in conditions.iter().enumerate() {
            for (j, species) in reactants.iter().enumerate() {
                let (gas_species_pool, condensed_species_pool, element_pool) =
                    Mixture::reaction_pool(&species);
                let mu_not = gas_species_pool
                    .iter()
                    .map(|species| {
                        let mu_not_i = species.data().h(*temperature_k)
                            - temperature_k * species.data().s(*temperature_k);
                        mu_not_i
                    })
                    .collect::<Vec<f64>>();
                let a = element_pool
                    .iter()
                    .map(|(_, element)| {
                        gas_species_pool
                            .iter()
                            .map(|species| {
                                let children = species.data().constituents();
                                children
                                    .iter()
                                    .find(|(_, child)| child == element)
                                    .unwrap_or(&(0.0, *element))
                                    .0
                            })
                            .collect::<Vec<f64>>()
                    })
                    .collect::<Vec<Vec<f64>>>();
                let b = element_pool
                    .iter()
                    .map(|(moles, _)| *moles)
                    .collect::<Vec<f64>>();
                let result = Mixture::solve_for_products(
                    *temperature_k,
                    *pressure_bar,
                    &gas_species_pool,
                    &condensed_species_pool,
                    &element_pool,
                    &a,
                    &b,
                    &mu_not,
                );
                pass_fail[i][j] = result.is_some();
                if result.is_none() {
                    println!(
                        "[debug] FAILED reactants={:?} T={:.1} P={:.1}",
                        species
                            .iter()
                            .map(|(n, s)| format!("{:.1}*{}", n, s.data().symbol()))
                            .collect::<Vec<_>>(),
                        temperature_k,
                        pressure_bar
                    );
                }
            }
        }
        println!("========== Summary ==========");
        for (i, (temperature_k, pressure_bar)) in conditions.iter().enumerate() {
            let total = pass_fail[i].len();
            let test_passed = pass_fail[i].iter().filter(|&x| *x).count();
            let row_header = &format!("({:.1} K; {:.1} bar) Results:", temperature_k, pressure_bar);
            println!(
                "{:>30} [{}/{}] {:.2}%",
                row_header,
                test_passed,
                total,
                100.0 * (test_passed as f64 / total as f64)
            )
        }
        let test_passed = !pass_fail.iter().any(|inner| inner.iter().any(|x| !x));
        let percentage: f64 = pass_fail
            .iter()
            .map(|inner| inner.iter().filter(|&x| *x).count())
            .sum::<usize>() as f64
            / (conditions.len() * reactants.len()) as f64;
        assert!(
            test_passed,
            "Only {:.2}% of scenarios passed!",
            percentage * 100.0
        );
    }

    #[test]
    fn solve_for_products_stability() {
        run_stability_suite();
    }

    #[test]
    fn solve_methane_combustion() {
        let reactants = vec![(1.0, Species::CH4), (1.0, Species::O2)];
        let mixture = Mixture::new(&reactants, 3_000.0, 50.0);
        mixture.print_products();
    }

    #[test]
    fn solve_ammonia_combustion() {
        let reactants = vec![(2.0, Species::NH3), (1.5, Species::O2)];
        let mixture = Mixture::new(&reactants, 300.0, 50.0);
        mixture.print_products();
    }

    #[test]
    fn solve_ions() {
        let reactants = vec![
            (1.0, Species::CH4),
            (1.0, Species::HPlus),
            (1.0, Species::OHNeg),
        ];
        let mixture = Mixture::new(&reactants, 300.0, 1.0);
        mixture.print_products();
    }

    #[test]
    fn solve_smr() {
        // Steam Reforming Methane
        let reactants = vec![(1.0, Species::CH4), (1.0, Species::H2O)];
        let mixture = Mixture::new(&reactants, 1_273.0, 1.0);
        mixture.print_products();
    }

    #[test]
    fn solve_water() {
        let reactants = vec![(1.0, Species::H2O)];
        let mixture = Mixture::new(&reactants, 400.0, 50.0);
        mixture.print_products();
    }
}
