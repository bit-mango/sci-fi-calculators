use crate::constants::{
    C_MW, CO_MW, ENTHALPY_CARBON_MONOXIDE, ENTHALPY_HYDROGEN, G_0, H_MW, H2_MW, O_MW, R,
    STD_REFERENCE_PRESSURE,
};

pub struct Propellant {
    pub species_feed_stock: Vec<(f64, Species)>, // moles, species
    pub species_with_disassociation: Vec<(f64, Species)>, // moles, species
    pub starting_temperature_k: f64,
    pub chamber_temperature_k: f64,
    pub chamber_pressure_bar: f64,
    pub exit_pressure_bar: f64,
    pub m_dot_kg_s: f64,
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
