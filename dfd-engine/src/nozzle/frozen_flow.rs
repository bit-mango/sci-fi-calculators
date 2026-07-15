use crate::constants::{
    C_MW, CO_MW, ENTHALPY_CARBON_MONOXIDE, ENTHALPY_HYDROGEN, G_0, H_MW, H2_MW, O_MW, R,
    STD_REFERENCE_PRESSURE,
};

use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::material_properties::ThermoReference;
use crate::thermo::reactions::{get_rxn_enthalpy, get_rxn_entropy};
use std::fmt;

use integrate::prelude::*;

pub struct FrozenFlowResult {
    pub chamber_temperature_k: f64,
    pub chamber_pressure_bar: f64,
    pub propellant_m_dot: f64,
    pub propellant_mean_mw: f64,
    pub exit_temperature_k: f64,
    pub exit_pressure_bar: f64,
    pub engine_isp: f64,
    pub engine_thrust: f64,
    pub engine_power_use: f64,

    // Used for other calculations but not displayed.
    pub h_total_chamber: f64,
    pub s_total_chamber: f64,
}
impl fmt::Display for FrozenFlowResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "
            ======= Frozen Flow Results =======
            Chamber Temperature:    {:.0} K
            Chamber Pressure:       {:.3} bar
            Propellant m_dot:       {:.3} kg/s
            Propeallnt Mean Mol Wt: {:.3} g
            Exit Temperature:       {:.0} K
            Exit Pressure:          {:.3} bar
            Engine Isp:             {:.0} s
            Engine Thrust:          {:.3} kN
            Engine Power Draw:      {:.3} MW
            ",
            self.chamber_temperature_k,
            self.chamber_pressure_bar,
            self.propellant_m_dot,
            self.propellant_mean_mw,
            self.exit_temperature_k,
            self.exit_pressure_bar,
            self.engine_isp,
            self.engine_thrust / 1.0e3,
            self.engine_power_use / 1.0e6
        )
    }
}

pub fn calculate_frozen_flow_results(
    chamber_temperature: f64,
    chamber_pressure: f64,
    propellant_m_dot: f64,
) -> FrozenFlowResult {
    let thermo_reference = ThermoReference::new();

    let h_tdp = thermo_reference.get_tdp("H");
    let h2_tdp = thermo_reference.get_tdp("H2");
    let c_tdp = thermo_reference.get_tdp("C");
    let o_tdp = thermo_reference.get_tdp("O");
    let co_tdp = thermo_reference.get_tdp("CO");

    // Reaction enthalpy H2 => 2H.
    let enthalpy_hydrogen_disassociation_rxn =
        get_rxn_enthalpy(chamber_temperature, vec![&h2_tdp], vec![&h_tdp, &h_tdp]);
    // Reaction entropy H2 => 2H.
    let entropy_hydrogen_disassociation_rxn =
        get_rxn_entropy(chamber_temperature, vec![&h2_tdp], vec![&h_tdp, &h_tdp]);
    // Reaction enthalpy CO => C + O.
    let enthalpy_carbon_monoxide_disassociation_rxn =
        get_rxn_enthalpy(chamber_temperature, vec![&co_tdp], vec![&c_tdp, &o_tdp]);
    // Reaction entropy CO => C + O.
    let entropy_carbon_monoxide_disassociation_rxn =
        get_rxn_entropy(chamber_temperature, vec![&co_tdp], vec![&c_tdp, &o_tdp]);

    let alpha_h2_chamber = calculate_disassociation_fraction(
        chamber_temperature,
        chamber_pressure,
        enthalpy_hydrogen_disassociation_rxn,
        entropy_hydrogen_disassociation_rxn,
    );

    let alpha_co_chamber = calculate_disassociation_fraction(
        chamber_temperature,
        chamber_pressure,
        enthalpy_carbon_monoxide_disassociation_rxn,
        entropy_carbon_monoxide_disassociation_rxn,
    );

    // println!("⍺_H2: {:.4}%", 100.0 * alpha_h2_chamber);
    // println!("⍺_CO: {:.4}%", 100.0 * alpha_co_chamber);

    // The propellant mixture is.
    // W * CO + X * C + X * O + Y*H2 + Z * 2H
    // Where:
    //  W: 1-⍺_CO
    //  X: ⍺_CO
    //  Y: 3*(1-⍺_H2)
    //  Z: 3*⍺_H2
    // 3 comes from the original composition with no disassociation, CO + 3H2.
    let n_chamber = vec![
        1.0 - alpha_co_chamber,
        alpha_co_chamber,
        alpha_co_chamber,
        3.0 * (1.0 - alpha_h2_chamber),
        3.0 * alpha_h2_chamber * 2.0,
    ]; // 2 hydrogen per 1 H2 disassociated.

    // Total enthalpy at chamber per 1 CO + 3H2 feed unit (34g) NOT mole fraction normalized.
    let h_sepcies_chamber = vec![
        co_tdp.h(chamber_temperature),
        c_tdp.h(chamber_temperature),
        o_tdp.h(chamber_temperature),
        h2_tdp.h(chamber_temperature),
        h_tdp.h(chamber_temperature),
    ];
    let h_total_chamber: f64 = n_chamber
        .iter()
        .zip(h_sepcies_chamber.iter())
        .map(|(n_i, h_i)| n_i * h_i)
        .sum();

    let n_sum: f64 = n_chamber.iter().sum();
    let x: Vec<f64> = n_chamber.iter().map(|n_i| n_i / n_sum).collect();

    let s_sepcies_chamber = vec![
        co_tdp.s(chamber_temperature) - R * (x[0] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        c_tdp.s(chamber_temperature) - R * (x[1] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        o_tdp.s(chamber_temperature) - R * (x[2] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        h2_tdp.s(chamber_temperature) - R * (x[3] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        h_tdp.s(chamber_temperature) - R * (x[4] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
    ];

    let s_total_chamber: f64 = n_chamber
        .iter()
        .zip(s_sepcies_chamber.iter())
        .map(|(n_i, s_i)| n_i * s_i)
        .sum();

    // Mixture properties.
    let mw = vec![CO_MW, C_MW, O_MW, H2_MW, H_MW];
    let mixture_mean_molecular_weight: f64 =
        x.iter().zip(mw.iter()).map(|(x_i, mw_i)| x_i * mw_i).sum();

    let cp = vec![
        co_tdp.cp(chamber_temperature),
        c_tdp.cp(chamber_temperature),
        o_tdp.cp(chamber_temperature),
        h2_tdp.cp(chamber_temperature),
        h_tdp.cp(chamber_temperature),
    ];
    let mixture_cp: f64 = x.iter().zip(cp.iter()).map(|(x_i, cp_i)| x_i * cp_i).sum();
    let mixture_mean_molecular_weight_kg = mixture_mean_molecular_weight / 1.0e3;
    let mixture_cp_mass_basis = mixture_cp / mixture_mean_molecular_weight_kg;
    let mixture_specific_gas_constant = R / mixture_mean_molecular_weight_kg;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    let mixture_starting_temperature = 1_000.0;

    // Cp is temperature dependent, so to find the total sensible heat we need to integrate.
    let cp_t_co = |temperature: f64| x[0] * co_tdp.cp(temperature);
    let cp_t_c = |temperature: f64| x[1] * c_tdp.cp(temperature);
    let cp_t_o = |temperature: f64| x[2] * o_tdp.cp(temperature);
    let cp_t_h2 = |temperature: f64| x[3] * h2_tdp.cp(temperature);
    let cp_t_h = |temperature: f64| x[4] * h_tdp.cp(temperature);
    let integrate_steps: u32 = 10_000;

    let sensible_co = trapezoidal_rule(
        cp_t_co,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_c = trapezoidal_rule(
        cp_t_c,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_o = trapezoidal_rule(
        cp_t_o,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_h2 = trapezoidal_rule(
        cp_t_h2,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_h = trapezoidal_rule(
        cp_t_h,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );

    let disassociation_h2 = ENTHALPY_HYDROGEN * 3.0 * alpha_h2_chamber / n_sum; // Enthalpy of hydrogen times the number of H2 moles actually split.
    let disassociation_co = ENTHALPY_CARBON_MONOXIDE * alpha_co_chamber / n_sum; // Enthalpy of carbon monoxide times the number of CO moles actually split.
    // Divide by n_sum so it is based off per mol basis.

    // Only sensible heat, no phase changes for any species.

    let engine_power = propellant_m_dot
        * (sensible_co
            + sensible_c
            + sensible_o
            + sensible_h2
            + sensible_h
            + disassociation_h2
            + disassociation_co)
        / mixture_mean_molecular_weight_kg;

    let exit_pressure = chamber_pressure * 5.0e-5;

    // Nozzle expansion. Assume frozen flow for lower bound Isp.
    // Guesstimate exit pressure for simplicity. TODO use area ratio later so we are bound by some nozzle size.
    let exit_temperature = chamber_temperature
        * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
    let exit_velocity =
        (2.0 * mixture_cp_mass_basis * (chamber_temperature - exit_temperature)).sqrt();
    let isp = exit_velocity / G_0;

    let engine_thrust = exit_velocity * propellant_m_dot;

    FrozenFlowResult {
        chamber_temperature_k: chamber_temperature,
        chamber_pressure_bar: chamber_pressure,
        propellant_m_dot: propellant_m_dot,
        propellant_mean_mw: mixture_mean_molecular_weight,
        exit_temperature_k: exit_temperature,
        exit_pressure_bar: exit_pressure,
        engine_isp: isp,
        engine_thrust: engine_thrust,
        engine_power_use: engine_power,
        h_total_chamber: h_total_chamber,
        s_total_chamber: s_total_chamber,
    }
}
