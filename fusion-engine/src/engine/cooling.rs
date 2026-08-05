// use crate::constants::{CH4_MW, H2O_MW, NH3_MW, UREA_LATENT_HEAT_OF_FUSION, UREA_MW};
// use crate::thermo::fluid_properties::ThermoReference;
// use crate::{Propellant, Species, nozzle::full_equilibrium_flow::FullEquilibriumFlowResult};

// pub struct HeatState {
//     low_grade_in: f64,
//     high_grade_in: f64,
//     low_grade_out: f64,
//     high_grade_out: f64,
// }

// pub fn calculate_heat_state(engine_power_draw: f64, engine_efficiency: f64, propellant_m_dot: f64) {
//     let thermo_reference = ThermoReference::new();
//     // Just handle case for CH4 + H2O propellant.
//     // 1:1 stoich so n dots are the same.
//     let n_dot_ch4 = propellant_m_dot / (CH4_MW + H2O_MW);
//     let n_dot_h2o = propellant_m_dot / (CH4_MW + H2O_MW);

//     // Assume CH4 starts as a liquid at 100 K.
//     let enthalpy_start = -0.64438e3; // J/mol
//     // Room temp is 20 C.
//     let enthalpy_room = 14.414e3; // J/mol
//     let ch4_low_grade_in = (enthalpy_room - enthalpy_start) * n_dot_ch4;

//     let enthalpy_room = thermo_reference.get_tdp("CH4").h(293.0);
//     let enthalpy_end = thermo_reference.get_tdp("CH4").h(1_000.0);
//     let ch4_high_grade_in = (enthalpy_end - enthalpy_room) * n_dot_ch4;
//     let enthalpy_start = 0.0018350e3; // J/mol
//     let enthalpy_room = 1.5141e3; // J/mol
//     let h2o_low_grade_in = (enthalpy_room - enthalpy_start) * n_dot_h2o;
//     let enthalpy_room = thermo_reference.get_tdp("H2O").h(293.0);
//     let enthalpy_end = thermo_reference.get_tdp("H2O").h(1_000.0);
//     let h2o_high_grade_in = (enthalpy_end - enthalpy_room) * n_dot_h2o;
//     // Finally we have energy needed in SMR.
//     let smr_heat_req = 206.0e3; // J/mol
//     let smr_high_grade_in = smr_heat_req * n_dot_ch4;
//     // Find waste heat from engine, for now assume all power to the engines!
//     let total_engine_power = engine_power_draw / engine_efficiency;
//     let engine_high_grade_out = total_engine_power * (1.0 - engine_efficiency);

//     println!("CH4 Low grade:        {:.3} kW", ch4_low_grade_in / 1.0e3);
//     println!("CH4 High grade:       {:.3} kW", ch4_high_grade_in / 1.0e3);
//     println!("H2O Low grade:        {:.3} kW", h2o_low_grade_in / 1.0e3);
//     println!("H2O High grade:       {:.3} kW", h2o_high_grade_in / 1.0e3);
//     println!("SMR High grade:       {:.3} kW", smr_high_grade_in / 1.0e3);
//     println!(
//         "Total Low Grade In:   {:.3} MW",
//         (ch4_low_grade_in + h2o_low_grade_in) / 1.0e6
//     );
//     println!(
//         "Total High Grade In:  {:.3} MW",
//         (ch4_high_grade_in + h2o_high_grade_in + smr_high_grade_in) / 1.0e6
//     );
//     println!(
//         "Total High Grade Out: {:.3} MW",
//         engine_high_grade_out / 1.0e6
//     );
// }

// pub fn calculate_heat_state_urea(
//     engine_power_draw: f64,
//     engine_efficiency: f64,
//     propellant_m_dot: f64,
// ) {
//     let thermo_reference = ThermoReference::new();
//     // Just handle case for CH4 + H2O + Urea propellant.
//     let n_dot = propellant_m_dot / (CH4_MW + H2O_MW + UREA_MW);
//     let n_dot_urea = n_dot;

//     // Assume CH4 starts as a liquid at 100 K.
//     let enthalpy_start = -0.64438e3; // J/mol
//     // Room temp is 20 C.
//     let enthalpy_room = 14.414e3; // J/mol
//     let ch4_low_grade_in = (enthalpy_room - enthalpy_start) * n_dot;

//     let enthalpy_start = 0.0018350e3; // J/mol
//     let enthalpy_room = 1.5141e3; // J/mol
//     let h2o_low_grade_in = (enthalpy_room - enthalpy_start) * n_dot;

//     let urea_low_grade_in = 1.49e3 * n_dot; // Estimated 0C to 20C.

//     // Assume all products start at room temp then raise to 1,000 K for
//     // sensible heat approx.
//     let enthalpy_room = thermo_reference.get_tdp("CO").h(293.0);
//     let enthalpy_end = thermo_reference.get_tdp("CO").h(1_000.0);
//     let mut sensible_products = (enthalpy_end - enthalpy_room) * 2.0 * n_dot; // 2 CO in products
//     let enthalpy_room = thermo_reference.get_tdp("N2").h(293.0);
//     let enthalpy_end = thermo_reference.get_tdp("N2").h(1_000.0);
//     sensible_products += (enthalpy_end - enthalpy_room) * n_dot; // 1 N2 in products
//     let enthalpy_room = thermo_reference.get_tdp("H2").h(293.0);
//     let enthalpy_end = thermo_reference.get_tdp("H2").h(1_000.0);
//     sensible_products += (enthalpy_end - enthalpy_room) * 5.0 * n_dot; // 5 H2 in products

//     // Not facturing in sensible heat of Urea as I dont know it.
//     let urea_latent_heat = UREA_LATENT_HEAT_OF_FUSION * n_dot_urea;
//     // Urea also decomposes.
//     let urea_decomposition = 85.0e3 * n_dot_urea; // 185 - 100 from HNCO + H2O exothermic reaction.
//     let methane_cracking = 74.8e3 * n_dot_urea;
//     let ammonia_cracking = 46.0e3 * 2.0 * n_dot_urea; // Get two moles of ammonia from 1 urea.
//     let boudard = 173.0e3 * n_dot_urea; // CO2 + C -> 2 CO.

//     // Find waste heat from engine, for now assume all power to the engines!
//     let total_engine_power = engine_power_draw / engine_efficiency;
//     let engine_high_grade_out = total_engine_power * (1.0 - engine_efficiency);

//     println!(
//         "[Urea] Total Low Grade In:    {:.3} MW",
//         (ch4_low_grade_in + h2o_low_grade_in + urea_low_grade_in) / 1.0e6
//     );

//     println!(
//         "[Urea] Total High Grade In:   {:.3} MW",
//         (sensible_products
//             + urea_latent_heat
//             + urea_decomposition
//             + methane_cracking
//             + ammonia_cracking
//             + boudard)
//             / 1.0e6
//     );
//     println!(
//         "[Urea] Total High Grade Out: {:.3} MW",
//         engine_high_grade_out / 1.0e6
//     );
// }

// pub fn calculate_heat_state_ammonia(
//     engine_power_draw: f64,
//     engine_efficiency: f64,
//     propellant_m_dot: f64,
// ) {
//     let thermo_reference = ThermoReference::new();
//     let n_dot = propellant_m_dot / NH3_MW;

//     // Assums Ammonia starts at 200 K.
//     let enthalpy_start = 0.33278e3;
//     let enthalpy_room = 28.632e3;
//     let low_grade = (enthalpy_room - enthalpy_start) * n_dot;

//     // Assume all products start at room temp then raise to 1,000 K for
//     // sensible heat approx.
//     let enthalpy_room = thermo_reference.get_tdp("N2").h(293.0);
//     let enthalpy_end = thermo_reference.get_tdp("N2").h(1_000.0);
//     let mut sensible_products = (enthalpy_end - enthalpy_room) * n_dot * 0.5; // 0.5 N2 in products
//     let enthalpy_room = thermo_reference.get_tdp("H2").h(293.0);
//     let enthalpy_end = thermo_reference.get_tdp("H2").h(1_000.0);
//     sensible_products += (enthalpy_end - enthalpy_room) * 1.5 * n_dot; // 1.5 H2 in products

//     let ammonia_cracking = 46.0e3 * n_dot;

//     // Find waste heat from engine, for now assume all power to the engines!
//     let total_engine_power = engine_power_draw / engine_efficiency;
//     let engine_high_grade_out = total_engine_power * (1.0 - engine_efficiency);

//     println!(
//         "[Ammonia] Total Low Grade In:    {:.3} MW",
//         low_grade / 1.0e6
//     );

//     println!(
//         "[Ammonia] Total High Grade In:   {:.3} MW",
//         (sensible_products + ammonia_cracking) / 1.0e6
//     );
//     println!(
//         "[Ammonia] Total High Grade Out: {:.3} MW",
//         engine_high_grade_out / 1.0e6
//     );
// }
