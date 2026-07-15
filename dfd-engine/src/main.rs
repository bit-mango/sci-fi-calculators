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

mod thermo;

use crate::thermo::{
    constants::{C_MW, CO_MW, G_0, H_MW, H2_MW, O_MW, R},
    material_properties::{PolynomialCoefficients, TemperatureDependentProperty},
};

use integrate::prelude::*;

const STD_REFERENCE_PRESSURE: f64 = 1.0; // [bar]
const ENTHALPY_HYDROGEN: f64 = 435.998e3; // ∆H [J/mol]
const ENTROPY_HYDROGEN: f64 = 98.753; // ∆S [J/K•mol]

const ENTHALPY_CARBON_MONOXIDE: f64 = 1076.375e3; // ∆H [J/mol]
const ENTROPY_CARBON_MONOXIDE: f64 = 121.498; // ∆S [J/K•mol]

// Use 10-100 bar for chamber pressure for now.

// Step 1: Determine chemical equilibrium of Propellants(H2 ⇌ 2H, CO ⇌ C + O  disassociation fractions).
// Need to maximize Hydrogen disassociation but minimize Carbon Monoxide disassociation(to minimize coking).
// The disassociation fraction is given by:
// Kp(T) = exp(-∆G(T)/RT)
// Kp = [4⍺^2 / (1-⍺^2)]* (P/P_0)
// where ⍺ is the disassociation fraction.
// Solve for ⍺.
// ⍺ = sqrt(N/D)
// where N = Kp * P_0 / P, D = 4 + Kp * P_0 / P
// where ∆G(T) = ∆H-T∆S
fn calculate_disassociation_fraction(
    chamber_temperature_k: f64,
    chamber_pressure_bar: f64,
    enthalpy: f64,
    entropy: f64,
) -> f64 {
    let gibbs = enthalpy - chamber_temperature_k * entropy;
    let kp = (-gibbs / (R * chamber_temperature_k)).exp();

    let numerator = kp * STD_REFERENCE_PRESSURE / chamber_pressure_bar;
    let denominator = 4.0 + numerator;
    let alpha = (numerator / denominator).sqrt();

    alpha
}

fn main() {
    let chamber_temperature = 5_750.0;
    let chamber_pressure = 10.0;

    let h_tdp = TemperatureDependentProperty::new().with_hydrogen_coefficients();
    let h2_tdp = TemperatureDependentProperty::new().with_hydrogen2_coefficients();
    let c_tdp = TemperatureDependentProperty::new().with_carbon_coefficients();
    let o_tdp = TemperatureDependentProperty::new().with_oxygen_coefficients();
    let co_tdp = TemperatureDependentProperty::new().with_carbon_monoxide_coefficients();

    // Reaction enthalpy H2 => 2H.
    let enthalpy_hydrogen_disassociation_rxn =
        2.0 * h_tdp.h(chamber_temperature) - h2_tdp.h(chamber_temperature);
    // Reaction entropy H2 => 2H.
    let entropy_hydrogen_disassociation_rxn =
        2.0 * h_tdp.s(chamber_temperature) - h2_tdp.s(chamber_temperature);
    // Reaction enthalpy CO => C + O.
    let enthalpy_carbon_monoxide_disassociation_rxn =
        c_tdp.h(chamber_temperature) + o_tdp.h(chamber_temperature) - co_tdp.h(chamber_temperature);
    // Reaction entropy CO => C + O.
    let entropy_carbon_monoxide_disassociation_rxn =
        c_tdp.s(chamber_temperature) + o_tdp.s(chamber_temperature) - co_tdp.s(chamber_temperature);

    let alpha_h2_chamber = calculate_disassociation_fraction(
        chamber_temperature,
        chamber_pressure,
        enthalpy_hydrogen_disassociation_rxn,
        entropy_hydrogen_disassociation_rxn,
    );

    let alpha_co_chamber = calculate_disassociation_fraction(
        chamber_temperature,
        chamber_pressure,
        enthalpy_carbon_monoxide_disassociation_rxn,
        entropy_carbon_monoxide_disassociation_rxn,
    );

    println!("Chamber Temperature: {:.0} K", chamber_temperature);
    println!("Chamber Pressure: {:.0} bar", chamber_pressure);
    println!("⍺_H2: {:.4}%", 100.0 * alpha_h2_chamber);
    println!("⍺_CO: {:.4}%", 100.0 * alpha_co_chamber);

    // The propellant mixture is.
    // W * CO + X * C + X * O + Y*H2 + Z * 2H
    // Where:
    //  W: 1-⍺_CO
    //  X: ⍺_CO
    //  Y: 3*(1-⍺_H2)
    //  Z: 3*⍺_H2
    // 3 comes from the original composition with no disassociation, CO + 3H2.
    let n_chamber = vec![
        1.0 - alpha_co_chamber,
        alpha_co_chamber,
        alpha_co_chamber,
        3.0 * (1.0 - alpha_h2_chamber),
        3.0 * alpha_h2_chamber * 2.0,
    ]; // 2 hydrogen per 1 H2 disassociated.

    // Total enthalpy at chamber per 1 CO + 3H2 feed unit (34g) NOT mole fraction normalized.
    let h_sepcies_chamber = vec![
        co_tdp.h(chamber_temperature),
        c_tdp.h(chamber_temperature),
        o_tdp.h(chamber_temperature),
        h2_tdp.h(chamber_temperature),
        h_tdp.h(chamber_temperature),
    ];
    let h_total_chamber: f64 = n_chamber
        .iter()
        .zip(h_sepcies_chamber.iter())
        .map(|(n_i, h_i)| n_i * h_i)
        .sum();

    let n_sum: f64 = n_chamber.iter().sum();
    let x: Vec<f64> = n_chamber.iter().map(|n_i| n_i / n_sum).collect();

    let s_sepcies_chamber = vec![
        co_tdp.s(chamber_temperature) - R * (x[0] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        c_tdp.s(chamber_temperature) - R * (x[1] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        o_tdp.s(chamber_temperature) - R * (x[2] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        h2_tdp.s(chamber_temperature) - R * (x[3] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
        h_tdp.s(chamber_temperature) - R * (x[4] * chamber_pressure / STD_REFERENCE_PRESSURE).ln(),
    ];

    let s_total_chamber: f64 = n_chamber
        .iter()
        .zip(s_sepcies_chamber.iter())
        .map(|(n_i, s_i)| n_i * s_i)
        .sum();

    // Mixture properties.
    let mw = vec![CO_MW, C_MW, O_MW, H2_MW, H_MW];
    let mixture_mean_molecular_weight: f64 =
        x.iter().zip(mw.iter()).map(|(x_i, mw_i)| x_i * mw_i).sum();
    println!(
        "Mixture Mean Molecular Weight: {:.3} g",
        mixture_mean_molecular_weight
    );

    let cp = vec![
        co_tdp.cp(chamber_temperature),
        c_tdp.cp(chamber_temperature),
        o_tdp.cp(chamber_temperature),
        h2_tdp.cp(chamber_temperature),
        h_tdp.cp(chamber_temperature),
    ];
    let mixture_cp: f64 = x.iter().zip(cp.iter()).map(|(x_i, cp_i)| x_i * cp_i).sum();
    let mixture_mean_molecular_weight_kg = mixture_mean_molecular_weight / 1.0e3;
    let mixture_cp_mass_basis = mixture_cp / mixture_mean_molecular_weight_kg;
    let mixture_specific_gas_constant = R / mixture_mean_molecular_weight_kg;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    // Assume the propellant mixture starts at 1,000 K as it was preheated using waste heat.
    let mixture_starting_temperature = 1_000.0;

    // Cp is temperature dependent, so to find the total sensible heat we need to integrate.
    let cp_t_co = |temperature: f64| x[0] * co_tdp.cp(temperature);
    let cp_t_c = |temperature: f64| x[1] * c_tdp.cp(temperature);
    let cp_t_o = |temperature: f64| x[2] * o_tdp.cp(temperature);
    let cp_t_h2 = |temperature: f64| x[3] * h2_tdp.cp(temperature);
    let cp_t_h = |temperature: f64| x[4] * h_tdp.cp(temperature);
    let integrate_steps: u32 = 10_000;

    let sensible_co = trapezoidal_rule(
        cp_t_co,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_c = trapezoidal_rule(
        cp_t_c,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_o = trapezoidal_rule(
        cp_t_o,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_h2 = trapezoidal_rule(
        cp_t_h2,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );
    let sensible_h = trapezoidal_rule(
        cp_t_h,
        mixture_starting_temperature,
        chamber_temperature,
        integrate_steps,
    );

    let disassociation_h2 = ENTHALPY_HYDROGEN * 3.0 * alpha_h2_chamber / n_sum; // Enthalpy of hydrogen times the number of H2 moles actually split.
    let disassociation_co = ENTHALPY_CARBON_MONOXIDE * alpha_co_chamber / n_sum; // Enthalpy of carbon monoxide times the number of CO moles actually split.
    // Divide by n_sum so it is based off per mol basis.

    // Only sensible heat, no phase changes for any species.
    let mixture_m_dot = 1.0;
    println!("Propellant Mixture m_dot: {:.3} kg/s", mixture_m_dot);

    let engine_power = mixture_m_dot
        * (sensible_co
            + sensible_c
            + sensible_o
            + sensible_h2
            + sensible_h
            + disassociation_h2
            + disassociation_co)
        / mixture_mean_molecular_weight_kg;

    println!("Engine Power: {:.3} MW", engine_power / 1.0e6);

    let exit_pressure = chamber_pressure * 5.0e-5;

    // Nozzle expansion. Assume frozen flow for lower bound Isp.
    // Guesstimate exit pressure for simplicity. TODO use area ratio later so we are bound by some nozzle size.
    let exit_temperature = chamber_temperature
        * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
    println!("Exit temperature: {} K", exit_temperature);
    let exit_velocity =
        (2.0 * mixture_cp_mass_basis * (chamber_temperature - exit_temperature)).sqrt();
    let isp = exit_velocity / G_0;

    println!("Calculated Isp: {:.0} s", isp);

    let engine_thrust = exit_velocity * mixture_m_dot;
    println!("Engine Thrust: {:.3} kN", engine_thrust / 1.0e3);

    // Assume full equilibrium flow for upper bound Isp. This is the exit velocity assuming we get full recombination.
    // Exit state recompute equilibrium at exit, T, P.
    // Iterate exit temperature until s_total_chamber ~ s_total_exit.
    let mut exit_temperature_low = exit_temperature;
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
            exit_pressure,
            enthalpy_hydrogen_disassociation_rxn,
            entropy_hydrogen_disassociation_rxn,
        );

        let alpha_co_exit = calculate_disassociation_fraction(
            exit_temperature_low,
            exit_pressure,
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
                - R * (x_exit[0] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            c_tdp.s(exit_temperature_low)
                - R * (x_exit[1] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            o_tdp.s(exit_temperature_low)
                - R * (x_exit[2] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            h2_tdp.s(exit_temperature_low)
                - R * (x_exit[3] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            h_tdp.s(exit_temperature_low)
                - R * (x_exit[4] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
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
            exit_pressure,
            enthalpy_hydrogen_disassociation_rxn,
            entropy_hydrogen_disassociation_rxn,
        );

        let alpha_co_exit = calculate_disassociation_fraction(
            exit_temperature_high,
            exit_pressure,
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
                - R * (x_exit[0] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            c_tdp.s(exit_temperature_high)
                - R * (x_exit[1] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            o_tdp.s(exit_temperature_high)
                - R * (x_exit[2] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            h2_tdp.s(exit_temperature_high)
                - R * (x_exit[3] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            h_tdp.s(exit_temperature_high)
                - R * (x_exit[4] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
        ];

        let s_total_exit_high: f64 = n_exit
            .iter()
            .zip(s_sepcies_exit.iter())
            .map(|(n_i, s_i)| n_i * s_i)
            .sum();

        if s_total_chamber < s_total_exit_low || s_total_chamber > s_total_exit_high {
            panic!(
                "Entropy outside bracket. Chamber Entropy: {:.2}, Bracket: {:.2} <-> {:.2}",
                s_total_chamber, s_total_exit_low, s_total_exit_high
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
                exit_pressure,
                enthalpy_hydrogen_disassociation_rxn,
                entropy_hydrogen_disassociation_rxn,
            );

            let alpha_co_exit = calculate_disassociation_fraction(
                exit_temperature_mid,
                exit_pressure,
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
                    - R * (x_exit[0] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
                c_tdp.s(exit_temperature_mid)
                    - R * (x_exit[1] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
                o_tdp.s(exit_temperature_mid)
                    - R * (x_exit[2] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
                h2_tdp.s(exit_temperature_mid)
                    - R * (x_exit[3] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
                h_tdp.s(exit_temperature_mid)
                    - R * (x_exit[4] * exit_pressure / STD_REFERENCE_PRESSURE).ln(),
            ];

            let s_total_exit_mid: f64 = n_exit
                .iter()
                .zip(s_sepcies_exit.iter())
                .map(|(n_i, s_i)| n_i * s_i)
                .sum();

            if (s_total_exit_mid - s_total_chamber).abs() < 0.001 {
                real_exit_temperaute = exit_temperature_mid;
                enthalpy_rxn.0 = enthalpy_carbon_monoxide_disassociation_rxn;
                entropy_rxn.0 = entropy_carbon_monoxide_disassociation_rxn;
                enthalpy_rxn.1 = enthalpy_hydrogen_disassociation_rxn;
                entropy_rxn.1 = entropy_hydrogen_disassociation_rxn;
                break;
            } else {
                if s_total_chamber > s_total_exit_mid {
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
        exit_pressure,
        enthalpy_rxn.0,
        entropy_rxn.0,
    );

    let alpha_h2_exit = calculate_disassociation_fraction(
        real_exit_temperaute,
        exit_pressure,
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
    let delta_h_per_kg = (h_total_chamber - h_total_exit) / feed_mass_kg;
    let exit_velocity = (2.0 * delta_h_per_kg).sqrt();
    let isp = exit_velocity / G_0;

    println!("Upper Isp: {:.0}", isp);
    println!("Exit Temperature: {:.0} K", real_exit_temperaute);
    println!("Exit Pressure: {:.3} mbar", exit_pressure * 1.0e3);
    println!("⍺_H2: {:.4}%", 100.0 * alpha_h2_exit);
    println!("⍺_CO: {:.4}%", 100.0 * alpha_co_exit);
}
