use crate::constants::{
    C_MW, CO_MW, ENTHALPY_CARBON_MONOXIDE, ENTHALPY_HYDROGEN, G_0, H_MW, H2_MW, O_MW, R,
    STD_REFERENCE_PRESSURE,
};
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::fluid_properties::{TemperatureDependentProperty, ThermoReference};
use crate::thermo::reactions::{get_rxn_enthalpy, get_rxn_entropy};

pub struct Propellant<'a> {
    species: Vec<(
        f64,
        Species,
        &'a TemperatureDependentProperty,
        [(f64, Species, &'a TemperatureDependentProperty); 2],
    )>, // moles, species, vec(disassociated (moles,species)
    starting_temperature_k: f64,
    chamber_temperature_k: f64,
    chamber_pressure_bar: f64,
    exit_pressure_bar: f64,
    m_dot_kg_s: f64,
}

impl<'a> Propellant<'a> {
    pub fn new(
        species: Vec<(
            f64,
            Species,
            &'a TemperatureDependentProperty,
            [(f64, Species, &'a TemperatureDependentProperty); 2],
        )>,
        starting_temperature_k: f64,
        chamber_temperature_k: f64,
        chamber_pressure_bar: f64,
        exit_pressure_bar: f64,
        m_dot_kg_s: f64,
    ) -> Self {
        Self {
            species,
            starting_temperature_k,
            chamber_temperature_k,
            chamber_pressure_bar,
            exit_pressure_bar,
            m_dot_kg_s,
        }
    }

    pub fn alphas(&self, temperature_k: f64, pressure_bar: f64) -> Vec<f64> {
        let alphas: Vec<f64> = self
            .species
            .iter()
            .map(|feed_stock_species_i| {
                // First get disassociation rxn enthalpy and entropy.
                let rxn_enthalpy = get_rxn_enthalpy(
                    temperature_k,
                    vec![feed_stock_species_i.2],
                    vec![feed_stock_species_i.3[0].2, feed_stock_species_i.3[1].2],
                );
                let rxn_entropy = get_rxn_entropy(
                    temperature_k,
                    vec![feed_stock_species_i.2],
                    vec![feed_stock_species_i.3[0].2, feed_stock_species_i.3[1].2],
                );

                let alpha = calculate_disassociation_fraction(
                    temperature_k,
                    pressure_bar,
                    rxn_enthalpy,
                    rxn_entropy,
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
            n.push(self.species[i].0 * self.species[i].3[0].0 * alpha);
            n.push(self.species[i].0 * self.species[i].3[1].0 * alpha);
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
            mw_total += x[i] * specie.3[0].1.mw();
            i += 1;
            mw_total += x[i] * specie.3[1].1.mw();
            i += 1;
        }

        mw_total
    }

    pub fn avg_cp(&self, x: &Vec<f64>, temperature_k: f64) -> f64 {
        let mut cp_total = 0.0;
        let mut i = 0;
        for specie in self.species.iter() {
            cp_total += x[i] * specie.2.cp(temperature_k);
            i += 1;
            cp_total += x[i] * specie.3[0].2.cp(temperature_k);
            i += 1;
            cp_total += x[i] * specie.3[1].2.cp(temperature_k);
            i += 1;
        }

        cp_total
    }

    pub fn feed_mass(&self) -> f64 {
        self.species
            .iter()
            .map(|specie| specie.0 * specie.1.mw())
            .sum()
    }
}

pub enum Species {
    H,
    H2,
    C,
    O,
    CO,
}

impl Species {
    pub fn symbol(&self) -> String {
        match self {
            Species::H => "H".to_string(),
            Species::H2 => "H2".to_string(),
            Species::C => "C".to_string(),
            Species::O => "O".to_string(),
            Species::CO => "CO".to_string(),
        }
    }

    pub fn mw(&self) -> f64 {
        match self {
            Species::H => H_MW,
            Species::H2 => H2_MW,
            Species::C => C_MW,
            Species::O => O_MW,
            Species::CO => CO_MW,
        }
    }
}
