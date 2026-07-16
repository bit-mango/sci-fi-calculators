use crate::constants::{
    C_MW, CO_MW, ENTHALPY_CARBON_MONOXIDE, ENTHALPY_HYDROGEN, G_0, H_MW, H2_MW, O_MW, R,
    STD_REFERENCE_PRESSURE,
};

use super::{calculate_state, propellant::Propellant};
use crate::thermo::fluid_properties::ThermoReference;
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
    pub h_total: f64,
    pub s_total: f64,
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
            Propellant Mean Mol Wt: {:.3} g
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
    propellant: &Propellant,
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

    let state_chamber = calculate_state(
        chamber_temperature,
        chamber_pressure,
        h_tdp,
        h2_tdp,
        c_tdp,
        o_tdp,
        co_tdp,
    );

    // Mixture properties.
    let mw = vec![CO_MW, C_MW, O_MW, H2_MW, H_MW];
    let mixture_mean_molecular_weight: f64 = state_chamber
        .x
        .iter()
        .zip(mw.iter())
        .map(|(x_i, mw_i)| x_i * mw_i)
        .sum();

    let cp = vec![
        co_tdp.cp(chamber_temperature),
        c_tdp.cp(chamber_temperature),
        o_tdp.cp(chamber_temperature),
        h2_tdp.cp(chamber_temperature),
        h_tdp.cp(chamber_temperature),
    ];
    let mixture_cp: f64 = state_chamber
        .x
        .iter()
        .zip(cp.iter())
        .map(|(x_i, cp_i)| x_i * cp_i)
        .sum();
    let mixture_mean_molecular_weight_kg = mixture_mean_molecular_weight / 1.0e3;
    let mixture_cp_mass_basis = mixture_cp / mixture_mean_molecular_weight_kg;
    let mixture_specific_gas_constant = R / mixture_mean_molecular_weight_kg;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    let mixture_starting_temperature = 1_000.0;

    let state_start = calculate_state(
        mixture_starting_temperature,
        chamber_pressure,
        h_tdp,
        h2_tdp,
        c_tdp,
        o_tdp,
        co_tdp,
    );

    let feed_mass_kg = 0.034;
    let engine_power =
        propellant_m_dot * (state_chamber.h_total - state_start.h_total) / feed_mass_kg;

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
        h_total: state_chamber.h_total,
        s_total: state_chamber.s_total,
    }
}
