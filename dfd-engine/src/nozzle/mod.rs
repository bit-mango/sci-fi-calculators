pub mod frozen_flow;
pub mod full_equilibrium_flow;
pub mod propellant;

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
    temperature_k: f64,
    pressure_bar: f64,
    h_tdp: &TemperatureDependentProperty,
    h2_tdp: &TemperatureDependentProperty,
    c_tdp: &TemperatureDependentProperty,
    o_tdp: &TemperatureDependentProperty,
    co_tdp: &TemperatureDependentProperty,
) -> State {
    // Reaction enthalpy H2 => 2H.
    let enthalpy_hydrogen_disassociation_rxn =
        get_rxn_enthalpy(temperature_k, vec![&h2_tdp], vec![&h_tdp, &h_tdp]);
    // Reaction entropy H2 => 2H.
    let entropy_hydrogen_disassociation_rxn =
        get_rxn_entropy(temperature_k, vec![&h2_tdp], vec![&h_tdp, &h_tdp]);
    // Reaction enthalpy CO => C + O.
    let enthalpy_carbon_monoxide_disassociation_rxn =
        get_rxn_enthalpy(temperature_k, vec![&co_tdp], vec![&c_tdp, &o_tdp]);
    // Reaction entropy CO => C + O.
    let entropy_carbon_monoxide_disassociation_rxn =
        get_rxn_entropy(temperature_k, vec![&co_tdp], vec![&c_tdp, &o_tdp]);

    let alpha_h2 = calculate_disassociation_fraction(
        temperature_k,
        pressure_bar,
        enthalpy_hydrogen_disassociation_rxn,
        entropy_hydrogen_disassociation_rxn,
    );

    let alpha_co = calculate_disassociation_fraction(
        temperature_k,
        pressure_bar,
        enthalpy_carbon_monoxide_disassociation_rxn,
        entropy_carbon_monoxide_disassociation_rxn,
    );

    // println!("⍺_H2: {:.4}%", 100.0 * alpha_h2);
    // println!("⍺_CO: {:.4}%", 100.0 * alpha_co);

    // The propellant mixture is.
    // W * CO + X * C + X * O + Y*H2 + Z * 2H
    // Where:
    //  W: 1-⍺_CO
    //  X: ⍺_CO
    //  Y: 3*(1-⍺_H2)
    //  Z: 3*⍺_H2
    // 3 comes from the original composition with no disassociation, CO + 3H2.
    let n = vec![
        1.0 - alpha_co,
        alpha_co,
        alpha_co,
        3.0 * (1.0 - alpha_h2),
        3.0 * alpha_h2 * 2.0,
    ]; // 2 hydrogen per 1 H2 disassociated.

    // Total enthalpy at chamber per 1 CO + 3H2 feed unit (34g) NOT mole fraction normalized.
    let h_sepcies = vec![
        co_tdp.h(temperature_k),
        c_tdp.h(temperature_k),
        o_tdp.h(temperature_k),
        h2_tdp.h(temperature_k),
        h_tdp.h(temperature_k),
    ];
    let h_total: f64 = n
        .iter()
        .zip(h_sepcies.iter())
        .map(|(n_i, h_i)| n_i * h_i)
        .sum();

    let n_sum: f64 = n.iter().sum();
    let x: Vec<f64> = n.iter().map(|n_i| n_i / n_sum).collect();

    let s_species: Vec<f64> = vec![
        co_tdp.s(temperature_k),
        c_tdp.s(temperature_k),
        o_tdp.s(temperature_k),
        h2_tdp.s(temperature_k),
        h_tdp.s(temperature_k),
    ]
    .iter()
    .zip(x.iter())
    .map(|(s_i, x_i)| s_i - R * (x_i * pressure_bar / STD_REFERENCE_PRESSURE).ln()) // Account for pressure change to entropy.
    .collect();

    let s_total: f64 = n
        .iter()
        .zip(s_species.iter())
        .map(|(n_i, s_i)| n_i * s_i)
        .sum();

    State {
        alpha_h2,
        alpha_co,
        n,
        x,
        h_total,
        s_total,
    }
}
