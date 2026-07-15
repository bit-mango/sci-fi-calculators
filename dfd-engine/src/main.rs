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
use crate::nozzle::frozen_flow::calculate_frozen_flow_results;
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::material_properties::ThermoReference;

// Use 10-100 bar for chamber pressure for now.

fn main() {
    let thermo_reference = ThermoReference::new();
    let chamber_temperature = 5_750.0;
    let chamber_pressure = 10.0;

    let frozen_flow_results =
        calculate_frozen_flow_results(chamber_temperature, chamber_pressure, 1.0);

    println!("{}", frozen_flow_results);

    let h_tdp = thermo_reference.get_tdp("H");
    let h2_tdp = thermo_reference.get_tdp("H2");
    let c_tdp = thermo_reference.get_tdp("C");
    let o_tdp = thermo_reference.get_tdp("O");
    let co_tdp = thermo_reference.get_tdp("CO");

    // Assume full equilibrium flow for upper bound Isp. This is the exit velocity assuming we get full recombination.
    // Exit state recompute equilibrium at exit, T, P.
    // Iterate exit temperature until s_total_chamber ~ s_total_exit.
    let mut exit_temperature_low = 300.0;
    let mut exit_temperature_high = chamber_temperature;
    let mut real_exit_temperaute = 0.0;
    let mut enthalpy_rxn = (0.0, 0.0);
    let mut entropy_rxn = (0.0, 0.0);
    for i in 0..100 {
        // Reaction enthalpy H2 => 2H.
        let enthalpy_hydrogen_disassociation_rxn =
            2.0 * h_tdp.h(exit_temperature_low) - h2_tdp.h(exit_temperature_low);
        // Reaction entropy H2 => 2H.
        let entropy_hydrogen_disassociation_rxn =
            2.0 * h_tdp.s(exit_temperature_low) - h2_tdp.s(exit_temperature_low);
        // Reaction enthalpy CO => C + O.
        let enthalpy_carbon_monoxide_disassociation_rxn = c_tdp.h(exit_temperature_low)
            + o_tdp.h(exit_temperature_low)
            - co_tdp.h(exit_temperature_low);
        // Reaction entropy CO => C + O.
        let entropy_carbon_monoxide_disassociation_rxn = c_tdp.s(exit_temperature_low)
            + o_tdp.s(exit_temperature_low)
            - co_tdp.s(exit_temperature_low);

        let alpha_h2_exit = calculate_disassociation_fraction(
            exit_temperature_low,
            frozen_flow_results.exit_pressure_bar,
            enthalpy_hydrogen_disassociation_rxn,
            entropy_hydrogen_disassociation_rxn,
        );

        let alpha_co_exit = calculate_disassociation_fraction(
            exit_temperature_low,
            frozen_flow_results.exit_pressure_bar,
            enthalpy_carbon_monoxide_disassociation_rxn,
            entropy_carbon_monoxide_disassociation_rxn,
        );
        let n_exit = vec![
            1.0 - alpha_co_exit,
            alpha_co_exit,
            alpha_co_exit,
            3.0 * (1.0 - alpha_h2_exit),
            3.0 * alpha_h2_exit * 2.0,
        ]; // 2 hydrogen per 1 H2 disassociated.

        let n_sum_exit: f64 = n_exit.iter().sum();
        let x_exit: Vec<f64> = n_exit.iter().map(|n_i| n_i / n_sum_exit).collect();

        let h_sepcies_exit = vec![
            co_tdp.h(exit_temperature_low),
            c_tdp.h(exit_temperature_low),
            o_tdp.h(exit_temperature_low),
            h2_tdp.h(exit_temperature_low),
            h_tdp.h(exit_temperature_low),
        ];
        let h_total_exit: f64 = n_exit
            .iter()
            .zip(h_sepcies_exit.iter())
            .map(|(n_i, h_i)| n_i * h_i)
            .sum();

        let s_sepcies_exit = vec![
            co_tdp.s(exit_temperature_low)
                - R * (x_exit[0] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            c_tdp.s(exit_temperature_low)
                - R * (x_exit[1] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            o_tdp.s(exit_temperature_low)
                - R * (x_exit[2] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            h2_tdp.s(exit_temperature_low)
                - R * (x_exit[3] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            h_tdp.s(exit_temperature_low)
                - R * (x_exit[4] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
        ];

        let s_total_exit_low: f64 = n_exit
            .iter()
            .zip(s_sepcies_exit.iter())
            .map(|(n_i, s_i)| n_i * s_i)
            .sum();

        // Reaction enthalpy H2 => 2H.
        let enthalpy_hydrogen_disassociation_rxn =
            2.0 * h_tdp.h(exit_temperature_high) - h2_tdp.h(exit_temperature_high);
        // Reaction entropy H2 => 2H.
        let entropy_hydrogen_disassociation_rxn =
            2.0 * h_tdp.s(exit_temperature_high) - h2_tdp.s(exit_temperature_high);
        // Reaction enthalpy CO => C + O.
        let enthalpy_carbon_monoxide_disassociation_rxn = c_tdp.h(exit_temperature_high)
            + o_tdp.h(exit_temperature_high)
            - co_tdp.h(exit_temperature_high);
        // Reaction entropy CO => C + O.
        let entropy_carbon_monoxide_disassociation_rxn = c_tdp.s(exit_temperature_high)
            + o_tdp.s(exit_temperature_high)
            - co_tdp.s(exit_temperature_high);

        let alpha_h2_exit = calculate_disassociation_fraction(
            exit_temperature_high,
            frozen_flow_results.exit_pressure_bar,
            enthalpy_hydrogen_disassociation_rxn,
            entropy_hydrogen_disassociation_rxn,
        );

        let alpha_co_exit = calculate_disassociation_fraction(
            exit_temperature_high,
            frozen_flow_results.exit_pressure_bar,
            enthalpy_carbon_monoxide_disassociation_rxn,
            entropy_carbon_monoxide_disassociation_rxn,
        );
        let n_exit = vec![
            1.0 - alpha_co_exit,
            alpha_co_exit,
            alpha_co_exit,
            3.0 * (1.0 - alpha_h2_exit),
            3.0 * alpha_h2_exit * 2.0,
        ]; // 2 hydrogen per 1 H2 disassociated.

        let n_sum_exit: f64 = n_exit.iter().sum();
        let x_exit: Vec<f64> = n_exit.iter().map(|n_i| n_i / n_sum_exit).collect();

        let h_sepcies_exit = vec![
            co_tdp.h(exit_temperature_high),
            c_tdp.h(exit_temperature_high),
            o_tdp.h(exit_temperature_high),
            h2_tdp.h(exit_temperature_high),
            h_tdp.h(exit_temperature_high),
        ];
        let h_total_exit: f64 = n_exit
            .iter()
            .zip(h_sepcies_exit.iter())
            .map(|(n_i, h_i)| n_i * h_i)
            .sum();

        let s_sepcies_exit = vec![
            co_tdp.s(exit_temperature_high)
                - R * (x_exit[0] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            c_tdp.s(exit_temperature_high)
                - R * (x_exit[1] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            o_tdp.s(exit_temperature_high)
                - R * (x_exit[2] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            h2_tdp.s(exit_temperature_high)
                - R * (x_exit[3] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
            h_tdp.s(exit_temperature_high)
                - R * (x_exit[4] * frozen_flow_results.exit_pressure_bar / STD_REFERENCE_PRESSURE)
                    .ln(),
        ];

        let s_total_exit_high: f64 = n_exit
            .iter()
            .zip(s_sepcies_exit.iter())
            .map(|(n_i, s_i)| n_i * s_i)
            .sum();

        if frozen_flow_results.s_total_chamber < s_total_exit_low
            || frozen_flow_results.s_total_chamber > s_total_exit_high
        {
            panic!(
                "Entropy outside bracket. Chamber Entropy: {:.2}, Bracket: {:.2} <-> {:.2}",
                frozen_flow_results.s_total_chamber, s_total_exit_low, s_total_exit_high
            );
        } else {
            // Compute middle entropy
            let exit_temperature_mid =
                exit_temperature_low + (exit_temperature_high - exit_temperature_low) / 2.0;
            // Reaction enthalpy H2 => 2H.
            let enthalpy_hydrogen_disassociation_rxn =
                2.0 * h_tdp.h(exit_temperature_mid) - h2_tdp.h(exit_temperature_mid);
            // Reaction entropy H2 => 2H.
            let entropy_hydrogen_disassociation_rxn =
                2.0 * h_tdp.s(exit_temperature_mid) - h2_tdp.s(exit_temperature_mid);
            // Reaction enthalpy CO => C + O.
            let enthalpy_carbon_monoxide_disassociation_rxn = c_tdp.h(exit_temperature_mid)
                + o_tdp.h(exit_temperature_mid)
                - co_tdp.h(exit_temperature_mid);
            // Reaction entropy CO => C + O.
            let entropy_carbon_monoxide_disassociation_rxn = c_tdp.s(exit_temperature_mid)
                + o_tdp.s(exit_temperature_mid)
                - co_tdp.s(exit_temperature_mid);

            let alpha_h2_exit = calculate_disassociation_fraction(
                exit_temperature_mid,
                frozen_flow_results.exit_pressure_bar,
                enthalpy_hydrogen_disassociation_rxn,
                entropy_hydrogen_disassociation_rxn,
            );

            let alpha_co_exit = calculate_disassociation_fraction(
                exit_temperature_mid,
                frozen_flow_results.exit_pressure_bar,
                enthalpy_carbon_monoxide_disassociation_rxn,
                entropy_carbon_monoxide_disassociation_rxn,
            );
            let n_exit = vec![
                1.0 - alpha_co_exit,
                alpha_co_exit,
                alpha_co_exit,
                3.0 * (1.0 - alpha_h2_exit),
                3.0 * alpha_h2_exit * 2.0,
            ]; // 2 hydrogen per 1 H2 disassociated.

            let n_sum_exit: f64 = n_exit.iter().sum();
            let x_exit: Vec<f64> = n_exit.iter().map(|n_i| n_i / n_sum_exit).collect();

            let h_sepcies_exit = vec![
                co_tdp.h(exit_temperature_mid),
                c_tdp.h(exit_temperature_mid),
                o_tdp.h(exit_temperature_mid),
                h2_tdp.h(exit_temperature_mid),
                h_tdp.h(exit_temperature_mid),
            ];
            let h_total_exit: f64 = n_exit
                .iter()
                .zip(h_sepcies_exit.iter())
                .map(|(n_i, h_i)| n_i * h_i)
                .sum();

            let s_sepcies_exit = vec![
                co_tdp.s(exit_temperature_mid)
                    - R * (x_exit[0] * frozen_flow_results.exit_pressure_bar
                        / STD_REFERENCE_PRESSURE)
                        .ln(),
                c_tdp.s(exit_temperature_mid)
                    - R * (x_exit[1] * frozen_flow_results.exit_pressure_bar
                        / STD_REFERENCE_PRESSURE)
                        .ln(),
                o_tdp.s(exit_temperature_mid)
                    - R * (x_exit[2] * frozen_flow_results.exit_pressure_bar
                        / STD_REFERENCE_PRESSURE)
                        .ln(),
                h2_tdp.s(exit_temperature_mid)
                    - R * (x_exit[3] * frozen_flow_results.exit_pressure_bar
                        / STD_REFERENCE_PRESSURE)
                        .ln(),
                h_tdp.s(exit_temperature_mid)
                    - R * (x_exit[4] * frozen_flow_results.exit_pressure_bar
                        / STD_REFERENCE_PRESSURE)
                        .ln(),
            ];

            let s_total_exit_mid: f64 = n_exit
                .iter()
                .zip(s_sepcies_exit.iter())
                .map(|(n_i, s_i)| n_i * s_i)
                .sum();

            if (s_total_exit_mid - frozen_flow_results.s_total_chamber).abs() < 0.001 {
                real_exit_temperaute = exit_temperature_mid;
                enthalpy_rxn.0 = enthalpy_carbon_monoxide_disassociation_rxn;
                entropy_rxn.0 = entropy_carbon_monoxide_disassociation_rxn;
                enthalpy_rxn.1 = enthalpy_hydrogen_disassociation_rxn;
                entropy_rxn.1 = entropy_hydrogen_disassociation_rxn;
                break;
            } else {
                if frozen_flow_results.s_total_chamber > s_total_exit_mid {
                    // Low is set to mid.
                    exit_temperature_low = exit_temperature_mid;
                } else {
                    // hihg is set to mid.
                    exit_temperature_high = exit_temperature_mid;
                }
            }
        }

        if i == 99 {
            panic!("Failed to find exit temperature.")
        }
    }

    let alpha_co_exit = calculate_disassociation_fraction(
        real_exit_temperaute,
        frozen_flow_results.exit_pressure_bar,
        enthalpy_rxn.0,
        entropy_rxn.0,
    );

    let alpha_h2_exit = calculate_disassociation_fraction(
        real_exit_temperaute,
        frozen_flow_results.exit_pressure_bar,
        enthalpy_rxn.1,
        entropy_rxn.1,
    );

    let n_exit = vec![
        1.0 - alpha_co_exit,
        alpha_co_exit,
        alpha_co_exit,
        3.0 * (1.0 - alpha_h2_exit),
        3.0 * alpha_h2_exit * 2.0,
    ]; // 2 hydrogen per 1 H2 disassociated.

    let h_sepcies_exit = vec![
        co_tdp.h(real_exit_temperaute),
        c_tdp.h(real_exit_temperaute),
        o_tdp.h(real_exit_temperaute),
        h2_tdp.h(real_exit_temperaute),
        h_tdp.h(real_exit_temperaute),
    ];
    let h_total_exit: f64 = n_exit
        .iter()
        .zip(h_sepcies_exit.iter())
        .map(|(n_i, h_i)| n_i * h_i)
        .sum();

    let feed_mass_kg = 0.034; // 1 mol CO (28g) + 3 mol H2 (6g) = 34g, fixed regardless of dissociation state
    let delta_h_per_kg = (frozen_flow_results.h_total_chamber - h_total_exit) / feed_mass_kg;
    let exit_velocity = (2.0 * delta_h_per_kg).sqrt();
    let isp = exit_velocity / G_0;

    println!("Upper Isp: {:.0}", isp);
    println!("Exit Temperature: {:.0} K", real_exit_temperaute);
    println!(
        "Exit Pressure: {:.3} mbar",
        frozen_flow_results.exit_pressure_bar * 1.0e3
    );
    println!("⍺_H2: {:.4}%", 100.0 * alpha_h2_exit);
    println!("⍺_CO: {:.4}%", 100.0 * alpha_co_exit);
}
