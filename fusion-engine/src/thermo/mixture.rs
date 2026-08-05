use super::equilibrium::{EquilibriumError, EquilibriumMode, solve_for_products};
use crate::constants::*;
use std::collections::HashMap;
use std::f64;
use thermo_species::Species;

#[derive(Clone, Default)]
pub struct Mixture {
    pub products: Vec<(f64, Species)>,
    pub temperature_k: f64,
    pub pressure_bar: f64,
    pub h_total: f64, // Total enthalpy.
    #[allow(dead_code)]
    pub s_total: f64, // Total entropy.
    pub n_total: f64,
    pub avg_mw: f64,
    pub avg_cp: f64,
}

impl Mixture {
    pub fn new(reactants: &[(f64, Species)], temperature_k: f64, pressure_bar: f64) -> Self {
        let result = solve_for_products(
            reactants,
            EquilibriumMode::TP {
                temperature_k,
                pressure_bar,
            },
            true,
            true,
            None,
        );
        match result {
            Ok(state) => {
                let (h_total, s_total, n_total, avg_mw, avg_cp) =
                    Self::state(&state.products, temperature_k, pressure_bar);

                Self {
                    products: state.products,
                    // Temperature and pressure won't change because mode is TP
                    // but writing it here for future modes where mode could change.
                    temperature_k: state.temperature_k,
                    pressure_bar: state.pressure_bar,
                    h_total,
                    s_total,
                    n_total,
                    avg_mw,
                    avg_cp,
                }
            }
            Err(err) => match err {
                EquilibriumError::FailedToConverge { iterations_used } => {
                    panic!("Failed to converge in {} iterations!", iterations_used);
                }
            },
        }
    }

    // Assumes reactants == products
    pub fn new_with_frozen_reactants(
        reactants: &[(f64, Species)],
        temperature_k: f64,
        pressure_bar: f64,
    ) -> Self {
        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(reactants, temperature_k, pressure_bar);

        Self {
            products: reactants.to_vec(),
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
        products: &[(f64, Species)],
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

    pub fn solve_for_target_enthalpy(&self, target_enthalpy: f64) -> Mixture {
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

    pub fn mix(&self, other: &Self) -> Self {
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
            reatants_mix.insert(key.to_string(), *s);
        }
        for o in cold_adj.products.iter() {
            let key = o.1.data().symbol();
            if let Some(entry) = reatants_mix.get_mut(key) {
                // reatants already exists! Increment moles.
                entry.0 += o.0;
            } else {
                // reatants is new, add them.
                reatants_mix.insert(key.to_string(), *o);
            }
        }
        let mut reactants = reatants_mix
            .drain()
            .map(|(_, v)| v)
            .collect::<Vec<(f64, Species)>>();
        reactants.sort_by_key(|a| a.1);

        let combined = Mixture::new(&reactants, hot.temperature_k, pressure_bar);
        combined.solve_for_target_enthalpy(starting_h)
    }

    pub fn scale(&self, factor: f64) -> Self {
        let scaled_products: Vec<(f64, Species)> = self
            .products
            .iter()
            .map(|(moles, species)| (*moles * factor, *species))
            .collect();
        Self::new_with_frozen_reactants(&scaled_products, self.temperature_k, self.pressure_bar)
    }

    pub fn feed_mass(&self) -> f64 {
        self.products
            .iter()
            .map(|(moles, specie)| moles * specie.data().mw())
            .sum()
    }
}
