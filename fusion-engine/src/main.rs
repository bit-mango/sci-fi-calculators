// dfd-engine.rs
// A DFD engine, or Direct Fusion Drive, utilizes a fusion reactor, in this case a
// Tritium, He4 variant, to super heat some propellant, and expel it out a magnetic nozzle.
// Notably the DFD engine can also serve to generate power as well.
// Normally DFD engines use Hydrogen as the propellant, and that is the ideal case, but
// Hydrogen is extremely difficult to store, especially for long space voyages.
// This calcualtor will explore using Methane(CH4), and water(H2O) as the propellant,
// by first using Methane Steam reforming to convert CH4 + H20 => CO + 3H2, then feeding that gas
// into the DFD engine to heat and expel it. Importantly, CO is extremely thermodynamically stable
// so Hydrogen will disassociate much more readily than CO, leading to less coking(Carbon build up in the nozzle),
// and a lower overall average molecular weight of the propellant(increases efficiency).
//
// This calculate will compare the reference design(Hydrogen only), to the proposed Methane + Water design.
//
// • -> opt + 8
// ∆ -> opt + j
// ⍺ -> ctrl + cmd + space, then alpha
//

mod constants;
mod engine;
mod hybrid_engine;
mod nozzle;
mod thermo;
// use crate::engine::cooling::{
//     calculate_heat_state, calculate_heat_state_ammonia, calculate_heat_state_urea,
// };
use crate::hybrid_engine::sweep_engine;
// use crate::nozzle::{
// full_equilibrium_flow::calculate_full_quilibrium_flow_results,
// };

// use crate::thermo::equilibrium::*;

fn main() {
    sweep_engine();
    // let target_area_ratio = 100.0;
    // let thermo_reference = ThermoReference::new();
    // println!("Scenario A: As planned");
    // let propellant = Propellant::new(
    //     &thermo_reference,
    //     vec![
    //         (1.0, Species::CO, vec![(1.0, Species::C), (1.0, Species::O)]),
    //         (3.0, Species::H2, vec![(2.0, Species::H)]),
    //     ],
    // );

    // let chamber_temperature = 8_000.0; // 4_800 K lowers engine waste == endothermic heat + dimished radiator. was 5600 K
    // let chamber_pressure = 200.0;
    // // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    // let starting_temperature = 1_000.0;
    // let propellant_m_dot = 1.0;

    // let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
    //     &propellant,
    //     starting_temperature,
    //     chamber_temperature,
    //     chamber_pressure,
    //     propellant_m_dot,
    //     target_area_ratio,
    // );

    // println!("{}", full_equilibrium_flow_results);
    // calculate_heat_state(
    //     full_equilibrium_flow_results
    //         .frozen_flow_results
    //         .engine_power_use,
    //     0.8,
    //     propellant_m_dot,
    // );

    // println!("Scenario B: Reduce Engine Power");
    // let propellant = Propellant::new(
    //     &thermo_reference,
    //     vec![
    //         (1.0, Species::CO, vec![(1.0, Species::C), (1.0, Species::O)]),
    //         (3.0, Species::H2, vec![(2.0, Species::H)]),
    //     ],
    // );

    // let chamber_temperature = 4_790.0; // 4_800 K lowers engine waste == endothermic heat + dimished radiator. was 5600 K
    // let chamber_pressure = 10.0;
    // // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    // let starting_temperature = 1_000.0;
    // let propellant_m_dot = 1.0;

    // let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
    //     &propellant,
    //     starting_temperature,
    //     chamber_temperature,
    //     chamber_pressure,
    //     propellant_m_dot,
    //     target_area_ratio,
    // );

    // println!("{}", full_equilibrium_flow_results);
    // calculate_heat_state(
    //     full_equilibrium_flow_results
    //         .frozen_flow_results
    //         .engine_power_use,
    //     0.8,
    //     propellant_m_dot,
    // );
    // println!("Scenario C: Urea");

    // // Use m dot of 2.264 to add urea, and drop temp to 4_600 to keep engine power the same.
    // let propellant = Propellant::new(
    //     &thermo_reference,
    //     vec![
    //         (2.0, Species::CO, vec![(1.0, Species::C), (1.0, Species::O)]),
    //         (1.0, Species::N2, vec![(2.0, Species::N)]),
    //         (5.0, Species::H2, vec![(2.0, Species::H)]),
    //     ],
    // );

    // let chamber_temperature = 5_050.0; // 4_800 K lowers engine waste == endothermic heat + dimished radiator. was 5600 K
    // let chamber_pressure = 10.0;
    // // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    // let starting_temperature = 1_000.0;
    // let propellant_m_dot = 2.0;

    // let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
    //     &propellant,
    //     starting_temperature,
    //     chamber_temperature,
    //     chamber_pressure,
    //     propellant_m_dot,
    //     target_area_ratio,
    // );

    // println!("{}", full_equilibrium_flow_results);
    // calculate_heat_state_urea(
    //     full_equilibrium_flow_results
    //         .frozen_flow_results
    //         .engine_power_use,
    //     0.8,
    //     propellant_m_dot,
    // );

    // println!("Scenario D: Ammonia Propellant");
    // let propellant = Propellant::new(
    //     &thermo_reference,
    //     vec![
    //         (1.0, Species::N2, vec![(2.0, Species::N)]),
    //         (3.0, Species::H2, vec![(2.0, Species::H)]),
    //     ],
    // );

    // let chamber_temperature = 8_000.0; // 4_800 K lowers engine waste == endothermic heat + dimished radiator. was 5600 K
    // let chamber_pressure = 10.0;
    // // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    // let starting_temperature = 1_000.0;
    // let propellant_m_dot = 1.0;

    // let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
    //     &propellant,
    //     starting_temperature,
    //     chamber_temperature,
    //     chamber_pressure,
    //     propellant_m_dot,
    //     target_area_ratio,
    // );

    // println!("{}", full_equilibrium_flow_results);
    // calculate_heat_state_ammonia(
    //     full_equilibrium_flow_results
    //         .frozen_flow_results
    //         .engine_power_use,
    //     0.8,
    //     propellant_m_dot,
    // );
    // // TODO ammonia propellant
}
