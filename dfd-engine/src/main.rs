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

use crate::constants::{G_0, R, STD_REFERENCE_PRESSURE};
use crate::nozzle::{
    frozen_flow::calculate_frozen_flow_results,
    full_equilibrium_flow::calculate_full_quilibrium_flow_results,
    propellant::{Propellant, Species},
};
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::fluid_properties::ThermoReference;

// Use 10-100 bar for chamber pressure for now.

fn main() {
    let thermo_reference = ThermoReference::new();
    let propellant = Propellant::new(
        vec![
            (
                1.0,
                Species::CO,
                thermo_reference.get_tdp("CO"),
                [
                    (1.0, Species::C, thermo_reference.get_tdp("C")),
                    (1.0, Species::O, thermo_reference.get_tdp("O")),
                ],
            ),
            (
                3.0,
                Species::H2,
                thermo_reference.get_tdp("H2"),
                [
                    (1.0, Species::H, thermo_reference.get_tdp("H")),
                    (1.0, Species::H, thermo_reference.get_tdp("H")),
                ],
            ),
        ],
        1_000.0,
        5_750.0,
        1.0,
        5.0e-5,
        1.0,
    );

    let chamber_temperature = 5_750.0;
    let chamber_pressure = 10.0;

    let frozen_flow_results =
        calculate_frozen_flow_results(&propellant, chamber_temperature, chamber_pressure, 1.0);

    println!("{}", frozen_flow_results);

    let full_equilibrium_flow_results = calculate_full_quilibrium_flow_results(
        &propellant,
        chamber_temperature,
        chamber_pressure,
        1.0,
    );

    println!("{}", full_equilibrium_flow_results);
}
