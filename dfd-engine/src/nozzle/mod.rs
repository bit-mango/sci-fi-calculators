pub mod frozen_flow;
pub mod full_equilibrium_flow;
pub mod propellant;

use super::Propellant;
use crate::constants::{R, STD_REFERENCE_PRESSURE};
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::fluid_properties::TemperatureDependentProperty;
use crate::thermo::reactions::{get_rxn_enthalpy, get_rxn_entropy};

#[derive(Default)]
pub struct State {
    pub alpha_h2: f64,
    pub alpha_co: f64,
    pub n: Vec<f64>,  // Number of moles for each species.
    pub x: Vec<f64>,  // Mol % for each species.
    pub h_total: f64, // Total enthalpy.
    pub s_total: f64, // Total entropy.
}

// TODO this seems like it should be in some propellant specific struct maybe? With defining characteristics of it?
// Chamber temp, pressure, propellant starting temp, etc. Then this type has two differnt flow logics you can use to analyze it?
// Oh and this type computes the average mw and stuff too. Maybe just propellant.rs?
pub fn calculate_state(
    propellant: &Propellant,
    temperature_k: f64,
    pressure_bar: f64,
    h_tdp: &TemperatureDependentProperty,
    h2_tdp: &TemperatureDependentProperty,
    c_tdp: &TemperatureDependentProperty,
    o_tdp: &TemperatureDependentProperty,
    co_tdp: &TemperatureDependentProperty,
) -> State {
    let alphas = propellant.alphas(temperature_k, pressure_bar);

    let n = propellant.n(&alphas);
    let h_total = propellant.h_total(&n, temperature_k);
    let x = propellant.x(&n);
    let s_total = propellant.s_total(&x, &n, temperature_k, pressure_bar);

    State {
        alpha_h2: alphas[0],
        alpha_co: alphas[1],
        n,
        x,
        h_total,
        s_total,
    }
}
