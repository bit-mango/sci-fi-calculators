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
mod nozzle;
mod thermo;

use crate::nozzle::{
    full_equilibrium_flow::calculate_full_quilibrium_flow_results,
    propellant::{Propellant, Species},
};
use crate::thermo::fluid_properties::ThermoReference;

fn main() {
    let target_area_ratio = 100.0;
    let thermo_reference = ThermoReference::new();
    let propellant = Propellant::new(
        &thermo_reference,
        vec![
            (1.0, Species::CO, vec![(1.0, Species::C), (1.0, Species::O)]),
            (3.0, Species::H2, vec![(2.0, Species::H)]),
        ],
    );

    let chamber_temperature = 5_600.0;
    let chamber_pressure = 10.0;
    // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    let starting_temperature = 1_273.15;

    let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
        &propellant,
        starting_temperature,
        chamber_temperature,
        chamber_pressure,
        1.0,
        target_area_ratio,
    );

    println!("{}", full_equilibrium_flow_results);

    let propellant = Propellant::new(
        &thermo_reference,
        vec![
            (1.0, Species::N2, vec![(2.0, Species::N)]),
            (3.0, Species::H2, vec![(2.0, Species::H)]),
        ],
    );

    let chamber_temperature = 7_000.0;
    let chamber_pressure = 50.0;
    // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    let starting_temperature = 1_273.15;

    let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
        &propellant,
        starting_temperature,
        chamber_temperature,
        chamber_pressure,
        0.85,
        target_area_ratio,
    );

    println!("{}", full_equilibrium_flow_results);

    let propellant = Propellant::new(
        &thermo_reference,
        vec![(3.0, Species::H2, vec![(2.0, Species::H)])],
    );

    let chamber_temperature = 20_000.0;
    let chamber_pressure = 10.0;
    // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    let starting_temperature = 1_273.15;

    let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
        &propellant,
        starting_temperature,
        chamber_temperature,
        chamber_pressure,
        0.09,
        target_area_ratio,
    );

    println!("{}", full_equilibrium_flow_results);
}
