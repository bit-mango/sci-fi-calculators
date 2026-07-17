use crate::constants::{CH4_MW, H2O_MW};
use crate::thermo::fluid_properties::ThermoReference;
use crate::{Propellant, Species, nozzle::full_equilibrium_flow::FullEquilibriumFlowResult};

pub struct HeatState {
    low_grade_in: f64,
    high_grade_in: f64,
    low_grade_out: f64,
    high_grade_out: f64,
}

pub fn calculate_heat_state(
    engine_power_draw: f64,
    engine_efficiency: f64,
    propellant_m_dot: f64,
) -> HeatState {
    let thermo_reference = ThermoReference::new();
    // Just handle case for CH4 + H2O propellant.
    // 1:1 stoich so n dots are the same.
    let n_dot_ch4 = propellant_m_dot / (CH4_MW + H2O_MW);
    let n_dot_h2o = propellant_m_dot / (CH4_MW + H2O_MW);

    // Assume CH4 starts as a liquid at 100 K.
    let enthalpy_start = -0.64438e3; // J/mol
    // Room temp is 20 C.
    let enthalpy_room = 14.414e3; // J/mol
    let ch4_low_grade_in = (enthalpy_room - enthalpy_start) * n_dot_ch4;

    let enthalpy_room = thermo_reference.get_tdp("CH4").h(293.0);
    let enthalpy_end = thermo_reference.get_tdp("CH4").h(1_000.0);
    let ch4_high_grade_in = (enthalpy_end - enthalpy_room) * n_dot_ch4;
    let enthalpy_start = 0.0018350e3; // J/mol
    let enthalpy_room = 1.5141e3; // J/mol
    let h2o_low_grade_in = (enthalpy_room - enthalpy_start) * n_dot_h2o;
    let enthalpy_room = thermo_reference.get_tdp("H2O").h(293.0);
    let enthalpy_end = thermo_reference.get_tdp("H2O").h(1_000.0);
    let h2o_high_grade_in = (enthalpy_end - enthalpy_room) * n_dot_h2o;
    // Finally we have energy needed in SMR.
    let smr_heat_req = 206.0e3; // J/mol
    let smr_high_grade_in = smr_heat_req * n_dot_ch4;
    // Find waste heat from engine, for now assume all power to the engines!
    let total_engine_power = engine_power_draw / engine_efficiency;
    let engine_high_grade_out = total_engine_power * (1.0 - engine_efficiency);

    println!("CH4 Low grade:        {:.3} kW", ch4_low_grade_in / 1.0e3);
    println!("CH4 High grade:       {:.3} kW", ch4_high_grade_in / 1.0e3);
    println!("H2O Low grade:        {:.3} kW", h2o_low_grade_in / 1.0e3);
    println!("H2O High grade:       {:.3} kW", h2o_high_grade_in / 1.0e3);
    println!("SMR High grade:       {:.3} kW", smr_high_grade_in / 1.0e3);
    println!(
        "Total Low Grade In:   {:.3} MW",
        (ch4_low_grade_in + h2o_low_grade_in) / 1.0e6
    );
    println!(
        "Total High Grade In:  {:.3} MW",
        (ch4_high_grade_in + h2o_high_grade_in + smr_high_grade_in) / 1.0e6
    );
    println!(
        "Total High Grade Out: {:.3} MW",
        engine_high_grade_out / 1.0e6
    );

    todo!()
}
