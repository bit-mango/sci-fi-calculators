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
    let species = Species::H;
    println!("Species: {}, mw: {}", species.symbol(), species.mw());
    let propellant = Propellant {
        species_feed_stock: vec![(1.0, Species::CO), (3.0, Species::H2)],
        species_with_disassociation: vec![
            (1.0, Species::CO),
            (1.0, Species::C),
            (1.0, Species::O),
            (3.0, Species::H2),
            (2.0, Species::H),
        ],
        starting_temperature_k: 1_000.0,
        chamber_temperature_k: 5_750.0,
        chamber_pressure_bar: 1.0,
        exit_pressure_bar: 5.0e-5,
        m_dot_kg_s: 1.0,
    };

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
