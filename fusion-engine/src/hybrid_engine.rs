use crate::constants::*;
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::fluid_properties::{TemperatureDependentProperty, ThermoReference};
use crate::thermo::reactions::{get_rxn_enthalpy, get_rxn_entropy};
use std::f64::consts::PI;

use crate::nozzle::{
    full_equilibrium_flow::calculate_full_quilibrium_flow_results,
    propellant::{Propellant, Species},
};

fn o_sourcing_cost_j_per_mol() -> f64 {
    O2_DISSOCIATION + H2O_ELECTROLYSIS
}

// Where:
// v_o is the voltage used to accelerate the oxygen
// v_ch4 is the voltage used to accelerate the methane
fn energy_in(v_o: f64, v_ch4: f64) -> f64 {
    let o_sourcing_cost = o_sourcing_cost_j_per_mol();
    o_sourcing_cost + CH4_FIRST_IONIZATION + F * (v_o + v_ch4) - O_ELECTRON_AFFINITY
}

fn energy_out(v_o: f64, v_ch4: f64) -> f64 {
    // We get none of the oxygen power back regardless of the source
    // O electron affinity is accounted for in partial oxidation energy.
    CH4_FIRST_IONIZATION + CH4_PARTIAL_OXIDATION_ENERGY + F * (v_o + v_ch4)
}

fn calculate_engine_output(v_o: f64, v_ch4: f64, collision_theta_deg: f64, engine_power: f64) {
    // Start with CH4 + H2O propellant feeds. 1 mole propellant is CH4 + H2O
    let mut total_energy_in_per_propellant_mole = 0.0;

    // First we need to split the water, and disassociate the resulting O2.
    total_energy_in_per_propellant_mole += H2O_ELECTROLYSIS + O2_DISSOCIATION;

    // Now we need to ionize the CH4 and give the electron to the O.
    total_energy_in_per_propellant_mole += CH4_FIRST_IONIZATION - O_ELECTRON_AFFINITY;

    // Now each species is accelerated using an electric field.
    let electrostatic_energy = F * (v_o + v_ch4);
    total_energy_in_per_propellant_mole += electrostatic_energy;

    let propellant_mw = CH4_MW + H2O_MW;
    let propellant_n_dot = engine_power / total_energy_in_per_propellant_mole;
    let propellant_m_dot = propellant_n_dot * propellant_mw;

    // Species are now at the reaction chamber, where they collide and combust.
    // collision_theta_deg controls how much kinetic energy is used to smash the species together,
    // versus how much is carried through out the nozzle.
    // A theta of 90 deg is a head on collision, great reaction rate, but kinetic energy turns into heat.
    // A theta of 0 degrees means the species run parallel and have a poor reaction rate, but kinetic energy is retained.

    let collision_theta_rad = collision_theta_deg * PI / 180.0;
    let species_mw = CH4_MW + O_MW;
    let v_y = (2.0 * electrostatic_energy / species_mw).sqrt() * collision_theta_rad.cos();

    // Any heat from collision plus reaction plus ionization. will go to raising the temperature of the products.
    // Products are CO + 3H2 (H2 injected from water electrolysis).
    // Assume all products start at room temperature.
    let thermo_reference = ThermoReference::new();
    let q_chamber = electrostatic_energy * collision_theta_rad.sin().powi(2)
        + CH4_PARTIAL_OXIDATION_ENERGY
        + CH4_FIRST_IONIZATION
        - O_ELECTRON_AFFINITY; // J/mol
    // Use bisection method to estimate
    let t_start = 300.0;
    let chamber_pressure = 100.0;
    let mut t_low = 300.0;
    let mut t_high = 20_000.0;
    let mut t_chamber = t_low + (t_high - t_low) / 2.0;
    let co_tdp = thermo_reference.get_tdp("CO");
    let c_tdp = thermo_reference.get_tdp("C");
    let o_tdp = thermo_reference.get_tdp("O");
    let h2_tdp = thermo_reference.get_tdp("H2");
    let h_tdp = thermo_reference.get_tdp("H");

    let rxn_enthalpy = get_rxn_enthalpy(t_start, &vec![co_tdp], &vec![c_tdp, o_tdp]);
    let rxn_entropy = get_rxn_entropy(t_start, &vec![co_tdp], &vec![c_tdp, o_tdp]);

    let mut alpha_co = calculate_disassociation_fraction(
        t_start,
        chamber_pressure,
        rxn_enthalpy,
        rxn_entropy,
        4.0,
    );

    let rxn_enthalpy = get_rxn_enthalpy(t_start, &vec![h2_tdp], &vec![h_tdp, h_tdp]);
    let rxn_entropy = get_rxn_entropy(t_start, &vec![h2_tdp], &vec![h_tdp, h_tdp]);

    let mut alpha_h2 = calculate_disassociation_fraction(
        t_start,
        chamber_pressure,
        rxn_enthalpy,
        rxn_entropy,
        1.0,
    );

    // Weighted by mole fraction.
    let mut enthalpy_start = (1.0 - alpha_co) * co_tdp.h(t_start);
    enthalpy_start += alpha_co * c_tdp.h(t_start);
    enthalpy_start += alpha_co * o_tdp.h(t_start);
    enthalpy_start += 3.0 * (1.0 - alpha_h2) * h2_tdp.h(t_start);
    enthalpy_start += 3.0 * 2.0 * alpha_h2 * h_tdp.h(t_start);

    for i in 0..200 {
        // println!(
        //     "Trying Temperature Range: {:.3} K <-> {:.3} K",
        //     t_low, t_high
        // );
        let rxn_enthalpy = get_rxn_enthalpy(t_chamber, &vec![co_tdp], &vec![c_tdp, o_tdp]);
        let rxn_entropy = get_rxn_entropy(t_chamber, &vec![co_tdp], &vec![c_tdp, o_tdp]);

        alpha_co = calculate_disassociation_fraction(
            t_chamber,
            chamber_pressure,
            rxn_enthalpy,
            rxn_entropy,
            4.0,
        );

        let rxn_enthalpy = get_rxn_enthalpy(t_chamber, &vec![h2_tdp], &vec![h_tdp, h_tdp]);
        let rxn_entropy = get_rxn_entropy(t_chamber, &vec![h2_tdp], &vec![h_tdp, h_tdp]);

        alpha_h2 = calculate_disassociation_fraction(
            t_chamber,
            chamber_pressure,
            rxn_enthalpy,
            rxn_entropy,
            1.0,
        );

        let mut enthalpy_end = (1.0 - alpha_co) * co_tdp.h(t_chamber);
        enthalpy_end += alpha_co * c_tdp.h(t_chamber);
        enthalpy_end += alpha_co * o_tdp.h(t_chamber);
        enthalpy_end += 3.0 * (1.0 - alpha_h2) * h2_tdp.h(t_chamber);
        enthalpy_end += 3.0 * 2.0 * alpha_h2 * h_tdp.h(t_chamber);

        let delta_enthalpy = enthalpy_end - enthalpy_start;

        // Change in enthalpy must equal heat energy in.
        if (delta_enthalpy - q_chamber).abs() < 1.0e-6 {
            // We found it!
            break;
        } else {
            if delta_enthalpy < q_chamber {
                // Need to raise enthalpy, making it hotter.
                t_low = t_chamber;
            } else {
                // Need to lower enthalpy, making it colder.
                t_high = t_chamber;
            }
            t_chamber = t_low + (t_high - t_low) / 2.0;
        }
        if i == 99 {
            panic!("Failed to find chamber temperature.")
        }
    }

    // Find average mw of propellant with disassociation.
    let n = vec![
        (1.0 - alpha_co),
        alpha_co,
        alpha_co,
        3.0 * (1.0 - alpha_h2),
        3.0 * 2.0 * alpha_h2,
    ];
    let propellant_mw: f64 = vec![CO_MW, C_MW, O_MW, H2_MW, H_MW]
        .iter()
        .zip(n.iter())
        .map(|(mw_i, n_i)| mw_i * n_i)
        .sum::<f64>()
        / n.iter().sum::<f64>();

    let v_heat = (2.0 * q_chamber / propellant_mw).sqrt();

    let v_e = v_y + v_heat;

    let isp = v_e / G_0;
    let thrust = propellant_m_dot * v_e;

    println!("Chamber Temperature: {:.3} K", t_chamber);
    println!("Alpha CO: {:.2} %", alpha_co * 100.0);
    println!("Alpha H2: {:.2} %", alpha_h2 * 100.0);
    println!("Isp: {:.0} s", isp);
    println!("Thrust: {:.3} kN", thrust / 1.0e3);
    println!("Propellant m dot: {:.3} g/s", propellant_m_dot * 1.0e3);
}

pub fn sweep_engine() {
    let voltages = vec![100.0, 250.0, 500.0, 1_000.0];
    let collision_theta_deg = 10.0;
    let engine_power = 50.0e6;
    for voltage in voltages {
        calculate_engine_output(voltage, voltage, collision_theta_deg, engine_power);
    }
}
