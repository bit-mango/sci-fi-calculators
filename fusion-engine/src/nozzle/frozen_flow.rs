// use crate::constants::{G_0, R};

// use super::area_ratio::exit_pressure_from_area_ratio;
// use super::propellant::Propellant;
// use std::fmt;

// #[derive(Default)]
// pub struct FrozenFlowResult {
//     pub chamber_temperature_k: f64,
//     pub chamber_pressure_bar: f64,
//     pub chamber_alphas: Vec<(String, f64)>,
//     pub propellant_m_dot: f64,
//     pub propellant_mean_mw: f64,
//     pub exit_temperature_k: f64,
//     pub exit_pressure_bar: f64,
//     pub engine_isp: f64,
//     pub engine_thrust: f64,
//     pub engine_power_use: f64,

//     // Used for other calculations but not displayed.
//     pub h_total: f64,
//     pub s_total: f64,
// }

// impl fmt::Display for FrozenFlowResult {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         let mut alphas = "Chamber Alphas: \n".to_string();
//         for alpha in self.chamber_alphas.iter() {
//             alphas += &format!("              ⍺_{}: {:.4}% \n", alpha.0, alpha.1);
//         }
//         write!(
//             f,
//             "
//             ======= Frozen Flow Results =======
//             Chamber Temperature:    {:.0} K
//             Chamber Pressure:       {:.3} bar
//             {}
//             Propellant m_dot:       {:.3} kg/s
//             Propellant Mean Mol Wt: {:.3} g
//             Exit Temperature:       {:.0} K
//             Exit Pressure:          {:.3} mbar
//             Engine Isp:             {:.0} s
//             Engine Thrust:          {:.3} kN
//             Engine Power Draw:      {:.3} MW
//             ",
//             self.chamber_temperature_k,
//             self.chamber_pressure_bar,
//             alphas,
//             self.propellant_m_dot,
//             self.propellant_mean_mw * 1.0e3,
//             self.exit_temperature_k,
//             self.exit_pressure_bar * 1.0e3,
//             self.engine_isp,
//             self.engine_thrust / 1.0e3,
//             self.engine_power_use / 1.0e6
//         )
//     }
// }

// pub fn calculate_frozen_flow_results(
//     propellant: &Propellant,
//     starting_temperature: f64,
//     chamber_temperature: f64,
//     chamber_pressure: f64,
//     propellant_m_dot: f64,
//     target_area_ratio: f64,
// ) -> FrozenFlowResult {
//     let state_chamber = propellant.state(chamber_temperature, chamber_pressure);

//     let mut chamber_alphas = vec![];
//     for i in 0..state_chamber.alphas.len() {
//         chamber_alphas.push((
//             propellant.species[i].1.symbol(),
//             state_chamber.alphas[i] * 100.0,
//         ));
//     }

//     // Mixture properties.
//     let mixture_cp_mass_basis = state_chamber.avg_cp / state_chamber.avg_mw;
//     let mixture_specific_gas_constant = R / state_chamber.avg_mw;
//     let mixture_gamma =
//         mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

//     let state_start = propellant.state(starting_temperature, chamber_pressure);

//     let feed_mass_kg = propellant.feed_mass();
//     let engine_power =
//         propellant_m_dot * (state_chamber.h_total - state_start.h_total) / feed_mass_kg;

//     let exit_pressure =
//         exit_pressure_from_area_ratio(chamber_pressure, target_area_ratio, mixture_gamma);

//     // Nozzle expansion. Assume frozen flow for lower bound Isp.
//     let exit_temperature = chamber_temperature
//         * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
//     let exit_velocity =
//         (2.0 * mixture_cp_mass_basis * (chamber_temperature - exit_temperature)).sqrt();
//     let isp = exit_velocity / G_0;

//     let engine_thrust = exit_velocity * propellant_m_dot;

//     FrozenFlowResult {
//         chamber_temperature_k: chamber_temperature,
//         chamber_pressure_bar: chamber_pressure,
//         chamber_alphas: chamber_alphas,
//         propellant_m_dot: propellant_m_dot,
//         propellant_mean_mw: state_chamber.avg_mw,
//         exit_temperature_k: exit_temperature,
//         exit_pressure_bar: exit_pressure,
//         engine_isp: isp,
//         engine_thrust: engine_thrust,
//         engine_power_use: engine_power,
//         h_total: state_chamber.h_total,
//         s_total: state_chamber.s_total,
//     }
// }
