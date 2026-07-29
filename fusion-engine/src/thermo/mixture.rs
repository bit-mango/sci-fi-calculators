use super::species::Species;
use crate::constants::*;
use crate::thermo::fluid_properties::{TemperatureDependentProperty, ThermoReference};
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;
use std::f64;

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
    pub fn new(
        tr: &ThermoReference,
        reactants: &Vec<(f64, Species)>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> Self {
        let (species_pool, element_pool) = Mixture::reaction_pool(reactants);
        let mu_not = species_pool
            .iter()
            .map(|species| {
                let tdp = tr.get_tdp(&species.symbol());
                let mu_not_i = tdp.h(temperature_k) - temperature_k * tdp.s(temperature_k);
                mu_not_i
            })
            .collect::<Vec<f64>>();
        let a = element_pool
            .iter()
            .map(|(_, element)| {
                species_pool
                    .iter()
                    .map(|species| {
                        let children = species.constituents();
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
            &species_pool,
            &element_pool,
            &a,
            &b,
            &mu_not,
        );

        let products = result.expect("Failed to converge.");

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

    fn solve_for_products(
        temperature_k: f64,
        pressure_bar: f64,
        n_g: &Vec<Species>,         // species pool
        n_lm: &Vec<(f64, Species)>, // element pool
        a: &Vec<Vec<f64>>,          // moles of element i in species j
        b: &Vec<f64>,               // element totals, n_lm.0
        mu_not: &Vec<f64>,
    ) -> Option<Vec<(f64, Species)>> {
        // Initialize
        let mut enn: f64 = 0.1;
        let mut ln_n = enn.ln();
        let mut enln = vec![(enn / n_g.len() as f64).ln(); n_g.len()];

        let mut iterations = 0;
        loop {
            let tm = (pressure_bar / enn).ln();
            let en: Vec<f64> = enln.iter().map(|x| x.exp()).collect();
            let mu: Vec<f64> = (0..n_g.len())
                .map(|j| mu_not[j] / (R * temperature_k) + enln[j] + tm)
                .collect();
            let sumn: f64 = en.iter().sum();

            let m = n_lm.len();

            let mut g = DMatrix::zeros(m + 1, m + 1);
            let mut rhs = DVector::zeros(m + 1);
            (0..m).for_each(|i| {
                (0..m).for_each(|k| {
                    g[(i, k)] = (0..n_g.len())
                        .map(|j| a[i][j] * a[k][j] * en[j])
                        .sum::<f64>();
                });
                let sum = (0..n_g.len()).map(|j| a[i][j] * en[j]).sum();
                g[(i, m)] = sum;
                g[(m, i)] = sum;
            });
            g[(m, m)] = sumn - enn;

            (0..m).for_each(|i| {
                let left: f64 = (0..n_g.len()).map(|j| a[i][j] * en[j] * mu[j]).sum();
                let right: f64 = (0..n_g.len()).map(|j| a[i][j] * en[j]).sum();
                rhs[i] = left + b[i] - right;
            });
            rhs[m] = (0..n_g.len()).map(|j| en[j] * mu[j]).sum::<f64>() + (enn - sumn);

            // Solve G*X = RHS using LU decomposition.
            let x = g.lu().solve(&rhs)?;

            // Calculate Deln.
            let dln_n = x[m];
            let deln: Vec<f64> = (0..n_g.len())
                .map(|j| {
                    let sum_a_pi: f64 = (0..m).map(|i| a[i][j] * x[i]).sum();
                    -mu[j] + sum_a_pi + dln_n
                })
                .collect();

            let mut ambda: f64 = 1.0;
            let mut threshold = 5.0 * dln_n.abs();

            // Step-size control for trace species.
            for j in 0..n_g.len() {
                if deln[j] > 0.0 {
                    // Check if species is currently trace.
                    // -9.21 corresponds to a mole fraction of ~1e-4.
                    // Enln[j] far below Ennl -> enln[j] < ln_n - 9.21
                    if enln[j] < ln_n - 9.21 {
                        let gap = (deln[j] - dln_n).abs();
                        if gap > 1.0e-8 {
                            let candidate = (-9.21 - enln[j] + ln_n).abs() / gap;
                            ambda = ambda.min(candidate);
                        }
                    } else if deln[j] > threshold {
                        threshold = deln[j];
                    }
                }
            }
            if threshold > 2.0 {
                ambda = ambda.min(2.0 / threshold);
            }

            // Update state variables using damped step.
            ln_n += ambda * dln_n;
            enn = ln_n.exp();

            let mut max_update: f64 = 0.0;
            for j in 0..n_g.len() {
                let step = ambda * deln[j];
                enln[j] += step;
                max_update = max_update.max(step.abs());
            }

            // Convergence check.
            if max_update < 1.0e-5 && (ambda * dln_n).abs() < 1.0e-5 {
                // Return the converged product moles
                let result: Vec<(f64, Species)> = enln
                    .iter()
                    .zip(n_g.iter())
                    .map(|(log_n, species)| (log_n.exp(), species.clone()))
                    .collect();
                break Some(result); // Converged!
            }

            iterations += 1;

            if iterations == 5_000 {
                break None;
            }
        }
    }

    pub fn feed_mass(&self) -> f64 {
        self.products
            .iter()
            .map(|(moles, specie)| moles * specie.mw())
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_stability_suite(tr: &ThermoReference) {
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
                let (species_pool, element_pool) = Mixture::reaction_pool(&species);
                let mu_not = species_pool
                    .iter()
                    .map(|species| {
                        let tdp = tr.get_tdp(&species.symbol());
                        let mu_not_i =
                            tdp.h(*temperature_k) - temperature_k * tdp.s(*temperature_k);
                        mu_not_i
                    })
                    .collect::<Vec<f64>>();
                let a = element_pool
                    .iter()
                    .map(|(_, element)| {
                        species_pool
                            .iter()
                            .map(|species| {
                                let children = species.constituents();
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
                    &species_pool,
                    &element_pool,
                    &a,
                    &b,
                    &mu_not,
                );
                pass_fail[i][j] = result.is_some();
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
        let tr = &ThermoReference::new();
        run_stability_suite(tr);
    }

    #[test]
    fn solve_methane_combustion() {
        let tr = &ThermoReference::new();
        let reactants = vec![(1.0, Species::CH4), (1.0, Species::O2)];
        let mixture = Mixture::new(tr, &reactants, 3_000.0, 50.0);
        mixture.print_products();
    }

    #[test]
    fn solve_ammonia_combustion() {
        let tr = &ThermoReference::new();
        let reactants = vec![(2.0, Species::NH3), (1.5, Species::O2)];
        let mixture = Mixture::new(tr, &reactants, 300.0, 50.0);
        mixture.print_products();
    }

    #[test]
    fn solve_ions() {
        let tr = &ThermoReference::new();
        let reactants = vec![
            (1.0, Species::CH4),
            (1.0, Species::HPlus),
            (1.0, Species::OHNeg),
        ];
        let mixture = Mixture::new(tr, &reactants, 300.0, 1.0);
        mixture.print_products();
    }

    #[test]
    fn solve_water() {
        let tr = &ThermoReference::new();
        let reactants = vec![(1.0, Species::H2O)];
        let mixture = Mixture::new(tr, &reactants, 400.0, 50.0);
        mixture.print_products();
    }
}
