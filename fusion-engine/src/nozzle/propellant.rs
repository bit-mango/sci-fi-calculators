use crate::constants::*;
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::fluid_properties::{TemperatureDependentProperty, ThermoReference};
use crate::thermo::reactions::{get_rxn_enthalpy, get_rxn_entropy};
use std::collections::HashMap;

#[derive(Default)]
pub struct PropellantState {
    pub alphas: Vec<f64>,
    pub h_total: f64, // Total enthalpy.
    pub s_total: f64, // Total entropy.
    pub n_total: f64, // Total moles.
    pub avg_mw: f64,
    pub avg_cp: f64,
}

#[derive(Clone)]
pub struct Propellant<'a> {
    pub species: Vec<(
        f64,
        Species,
        &'a TemperatureDependentProperty,
        Vec<(f64, Species, &'a TemperatureDependentProperty)>,
    )>, // moles, species, vec(disassociated (moles,species)
}

impl<'a> Propellant<'a> {
    pub fn new(
        thermo_reference: &'a ThermoReference,
        species: Vec<(f64, Species, Vec<(f64, Species)>)>,
    ) -> Self {
        for specie in species.iter() {
            if specie.2.len() > 2 {
                panic!("Can only disassociate into at most two species!");
            }
        }
        let species_with_tdp = species
            .iter()
            .map(|(mol, specie, disassociated)| {
                let disassociated_with_tdp = disassociated
                    .iter()
                    .map(|(d_mol, d_specie)| {
                        (
                            *d_mol,
                            *d_specie,
                            thermo_reference.get_tdp(&d_specie.symbol()),
                        )
                    })
                    .collect();
                (
                    *mol,
                    *specie,
                    thermo_reference.get_tdp(&specie.symbol()),
                    disassociated_with_tdp,
                )
            })
            .collect();
        Self {
            species: species_with_tdp,
        }
    }

    pub fn mix(&self, other: &Self) -> Self {
        let mut species_mix: HashMap<
            String,
            (
                f64,
                Species,
                &'a TemperatureDependentProperty,
                Vec<(f64, Species, &'a TemperatureDependentProperty)>,
            ),
        > = HashMap::new();
        // Add all of original species to mix.
        for s in self.species.iter() {
            let key = s.1.symbol();
            species_mix.insert(key, s.clone());
        }
        for o in other.species.iter() {
            let key = o.1.symbol();
            if let Some(entry) = species_mix.get_mut(&key) {
                // Species already exists! Increment moles.
                entry.0 += o.0;
            } else {
                // Species is new, add them.
                species_mix.insert(key, o.clone());
            }
        }
        let species = species_mix.drain().map(|(_, v)| v).collect();
        Self { species }
    }

    pub fn alphas(&self, temperature_k: f64, pressure_bar: f64) -> Vec<f64> {
        let alphas: Vec<f64> = self
            .species
            .iter()
            .map(|feed_stock_species_i| {
                // First get disassociation rxn enthalpy and entropy.
                let reactants = vec![feed_stock_species_i.2];

                let products = if feed_stock_species_i.3.len() == 2 {
                    vec![feed_stock_species_i.3[0].2, feed_stock_species_i.3[1].2]
                } else {
                    vec![feed_stock_species_i.3[0].2, feed_stock_species_i.3[0].2]
                };
                let product_factor = if feed_stock_species_i.3.len() == 2 {
                    1.0
                } else {
                    4.0
                };
                let rxn_enthalpy = get_rxn_enthalpy(temperature_k, &reactants, &products);
                let rxn_entropy = get_rxn_entropy(temperature_k, &reactants, &products);

                let alpha = calculate_disassociation_fraction(
                    temperature_k,
                    pressure_bar,
                    rxn_enthalpy,
                    rxn_entropy,
                    product_factor,
                );
                alpha
            })
            .collect();
        alphas
    }

    pub fn n(&self, alphas: &Vec<f64>) -> Vec<f64> {
        if alphas.len() != self.species.len() {
            panic!("Species and alphas must be same length!");
        }
        let mut n: Vec<f64> = vec![];
        for i in 0..self.species.len() {
            let alpha = alphas[i];
            n.push(self.species[i].0 * (1.0 - alpha));
            for disassociated_specie in self.species[i].3.iter() {
                n.push(self.species[i].0 * disassociated_specie.0 * alpha);
            }
        }

        n
    }

    pub fn h_total(&self, n: &Vec<f64>, temperature_k: f64) -> f64 {
        let mut h_total = 0.0;
        let mut i = 0;
        for specie in self.species.iter() {
            h_total += n[i] * specie.2.h(temperature_k);
            i += 1;
            for disassociated_specie in specie.3.iter() {
                h_total += n[i] * disassociated_specie.2.h(temperature_k);
                i += 1;
            }
        }

        h_total
    }

    pub fn x(&self, n: &Vec<f64>) -> Vec<f64> {
        let n_sum: f64 = n.iter().sum();
        n.iter().map(|n_i| n_i / n_sum).collect()
    }

    pub fn s_total(
        &self,
        x: &Vec<f64>,
        n: &Vec<f64>,
        temperature_k: f64,
        pressure_bar: f64,
    ) -> f64 {
        let mut s_total = 0.0;
        let mut i = 0;
        for specie in self.species.iter() {
            s_total += n[i]
                * (specie.2.s(temperature_k)
                    - R * (x[i] * pressure_bar / STD_REFERENCE_PRESSURE).ln());
            i += 1;
            for disassociated_specie in specie.3.iter() {
                s_total += n[i]
                    * (disassociated_specie.2.s(temperature_k)
                        - R * (x[i] * pressure_bar / STD_REFERENCE_PRESSURE).ln());
                i += 1;
            }
        }

        s_total
    }

    pub fn avg_mw(&self, x: &Vec<f64>) -> f64 {
        let mut mw_total = 0.0;
        let mut i = 0;
        for specie in self.species.iter() {
            mw_total += x[i] * specie.1.mw();
            i += 1;
            for disassociated_specie in specie.3.iter() {
                mw_total += x[i] * disassociated_specie.1.mw();
                i += 1;
            }
        }

        mw_total
    }

    pub fn avg_cp(&self, x: &Vec<f64>, temperature_k: f64) -> f64 {
        let mut cp_total = 0.0;
        let mut i = 0;
        for specie in self.species.iter() {
            cp_total += x[i] * specie.2.cp(temperature_k);
            i += 1;
            for disassociated_specie in specie.3.iter() {
                cp_total += x[i] * disassociated_specie.2.cp(temperature_k);
                i += 1;
            }
        }

        cp_total
    }

    pub fn feed_mass(&self) -> f64 {
        self.species
            .iter()
            .map(|specie| specie.0 * specie.1.mw())
            .sum()
    }

    pub fn state(&self, temperature_k: f64, pressure_bar: f64) -> PropellantState {
        let alphas = self.alphas(temperature_k, pressure_bar);

        let n = self.n(&alphas);
        let h_total = self.h_total(&n, temperature_k);
        let x = self.x(&n);
        let s_total = self.s_total(&x, &n, temperature_k, pressure_bar);

        let n_total = n.iter().sum();
        let avg_mw: f64 = self.avg_mw(&x);
        let avg_cp: f64 = self.avg_cp(&x, temperature_k);

        PropellantState {
            alphas,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        }
    }
}

#[derive(Copy, Clone, PartialEq)]
pub enum Species {
    H,
    H2,
    C,
    O,
    CO,
    N2,
    N,
}

impl Species {
    pub fn symbol(&self) -> String {
        match self {
            Species::H => "H".to_string(),
            Species::H2 => "H2".to_string(),
            Species::C => "C".to_string(),
            Species::O => "O".to_string(),
            Species::CO => "CO".to_string(),
            Species::N2 => "N2".to_string(),
            Species::N => "N".to_string(),
        }
    }

    pub fn mw(&self) -> f64 {
        match self {
            Species::H => H_MW,
            Species::H2 => H2_MW,
            Species::C => C_MW,
            Species::O => O_MW,
            Species::CO => CO_MW,
            Species::N2 => N2_MW,
            Species::N => N_MW,
        }
    }
}
