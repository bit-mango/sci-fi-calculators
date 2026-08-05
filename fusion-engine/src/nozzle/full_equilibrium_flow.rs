// use crate::{constants::G_0, nozzle::frozen_flow::FrozenFlowResult};

// use super::{
//     frozen_flow::calculate_frozen_flow_results,
//     propellant::{Propellant, PropellantState},
// };
// use std::fmt;

// #[derive(Default)]
// pub struct FullEquilibriumFlowResult {
//     pub exit_temperature_k: f64,
//     pub exit_pressure_bar: f64,
//     pub engine_isp: f64,
//     pub engine_thrust: f64,
//     pub frozen_flow_results: FrozenFlowResult,
// }

// impl fmt::Display for FullEquilibriumFlowResult {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         let mut alphas = "Chamber Alphas: \n".to_string();
//         for alpha in self.frozen_flow_results.chamber_alphas.iter() {
//             alphas += &format!("              ⍺_{}: {:.4}% \n", alpha.0, alpha.1);
//         }
//         let estimated_isp = 0.85 * (self.engine_isp - self.frozen_flow_results.engine_isp)
//             + self.frozen_flow_results.engine_isp;
//         let estimated_thrust = estimated_isp * G_0 * self.frozen_flow_results.propellant_m_dot;
//         write!(
//             f,
//             "
//             ========= Frozen Flow Results =========
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

//             ==== Full Equilibrium Flow Results ====
//             Exit Temperature:       {:.0} K
//             Exit Pressure:          {:.3} mbar
//             Engine Isp:             {:.0} s
//             Engine Thrust:          {:.3} kN

//             =============== Overall ===============
//             Engine Isp Range(s):     {:.0} <-> {:.0}
//             Engine Thrust Range(kN): {:.3} <-> {:.3}
//             Estimated Isp(s):        {:.0}
//             Estimated Thrust(kN):    {:.3}
//             ",
//             self.frozen_flow_results.chamber_temperature_k,
//             self.frozen_flow_results.chamber_pressure_bar,
//             alphas,
//             self.frozen_flow_results.propellant_m_dot,
//             self.frozen_flow_results.propellant_mean_mw * 1.0e3,
//             self.frozen_flow_results.exit_temperature_k,
//             self.frozen_flow_results.exit_pressure_bar * 1.0e3,
//             self.frozen_flow_results.engine_isp,
//             self.frozen_flow_results.engine_thrust / 1.0e3,
//             self.frozen_flow_results.engine_power_use / 1.0e6,
//             self.exit_temperature_k,
//             self.exit_pressure_bar * 1.0e3,
//             self.engine_isp,
//             self.engine_thrust / 1.0e3,
//             self.frozen_flow_results.engine_isp,
//             self.engine_isp,
//             self.frozen_flow_results.engine_thrust / 1.0e3,
//             self.engine_thrust / 1.0e3,
//             estimated_isp,
//             estimated_thrust / 1.0e3
//         )
//     }
// }
// pub fn calculate_full_quilibrium_flow_results(
//     propellant: &Propellant,
//     starting_temperature: f64,
//     chamber_temperature: f64,
//     chamber_pressure: f64,
//     propellant_m_dot: f64,
//     target_area_ratio: f64,
// ) -> FullEquilibriumFlowResult {
//     let frozen_flow_results = calculate_frozen_flow_results(
//         propellant,
//         starting_temperature,
//         chamber_temperature,
//         chamber_pressure,
//         propellant_m_dot,
//         target_area_ratio,
//     );

//     // Assume full equilibrium flow for upper bound Isp. This is the exit velocity assuming we get full recombination.
//     // Exit state recompute equilibrium at exit, T, P.
//     // Iterate exit temperature until s_total_chamber ~ s_total_exit.
//     let mut exit_temperature_low = 300.0;
//     let mut exit_temperature_mid = 0.0;
//     let mut exit_temperature_high = chamber_temperature;

//     let mut state_low;
//     let mut state_mid = PropellantState::default();
//     let mut state_high;

//     for i in 0..100 {
//         state_low = propellant.state(exit_temperature_low, frozen_flow_results.exit_pressure_bar);
//         state_high = propellant.state(exit_temperature_high, frozen_flow_results.exit_pressure_bar);

//         if frozen_flow_results.s_total < state_low.s_total
//             || frozen_flow_results.s_total > state_high.s_total
//         {
//             panic!(
//                 "Entropy outside bracket. Chamber Entropy: {:.2}, Bracket: {:.2} <-> {:.2}",
//                 frozen_flow_results.s_total, state_low.s_total, state_high.s_total
//             );
//         } else {
//             // Compute middle entropy
//             exit_temperature_mid =
//                 exit_temperature_low + (exit_temperature_high - exit_temperature_low) / 2.0;
//             state_mid =
//                 propellant.state(exit_temperature_mid, frozen_flow_results.exit_pressure_bar);

//             if (state_mid.s_total - frozen_flow_results.s_total).abs() < 0.001 {
//                 break;
//             } else {
//                 if frozen_flow_results.s_total > state_mid.s_total {
//                     // Low is set to mid.
//                     exit_temperature_low = exit_temperature_mid;
//                 } else {
//                     // hihg is set to mid.
//                     exit_temperature_high = exit_temperature_mid;
//                 }
//             }
//         }

//         if i == 99 {
//             panic!("Failed to find exit temperature.")
//         }
//     }

//     let feed_mass_kg = propellant.feed_mass(); // 1 mol CO (28g) + 3 mol H2 (6g) = 34g, fixed regardless of dissociation state
//     let delta_h_per_kg = (frozen_flow_results.h_total - state_mid.h_total) / feed_mass_kg;
//     let exit_velocity = (2.0 * delta_h_per_kg).sqrt();
//     let isp = exit_velocity / G_0;
//     let engine_thrust = exit_velocity * propellant_m_dot;

//     FullEquilibriumFlowResult {
//         exit_temperature_k: exit_temperature_mid,
//         exit_pressure_bar: frozen_flow_results.exit_pressure_bar,
//         engine_isp: isp,
//         engine_thrust: engine_thrust,
//         frozen_flow_results,
//     }
// }
