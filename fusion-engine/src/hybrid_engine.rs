use crate::constants::*;
use crate::nozzle::area_ratio::exit_pressure_from_area_ratio;
use std::f64::consts::PI;
use thermo_species::{Species, species};

use crate::thermo::equilibrium::{EquilibriumError, EquilibriumMode, solve_for_products};
use crate::thermo::mixture::Mixture;

fn required_max_channel_area(
    n_dot_reaction: f64,
    species_mw: f64,
    voltage: f64,
    chamber_pressure_bar: f64,
) -> f64 {
    let m_dot_channel = n_dot_reaction * species_mw; // kg/s through this one channel.
    let v_beam = (2.0 * F * voltage / species_mw).sqrt(); // m/s Full beam speed, not the axial-projected v_y
    let momentum_flux = m_dot_channel * v_beam; // N
    let chamber_pressure_pa = chamber_pressure_bar * 1.0e5; // Convert to Pa
    momentum_flux / chamber_pressure_pa
}

fn min_channel_area_child_langmuir(
    total_current: f64,
    species_mw: f64,
    voltage: f64,
    gap_spacing: f64, // meters, accelerator gap, chosen design parameter
) -> f64 {
    let q_over_m = F / species_mw; // C/kg
    let j_max =
        (4.0 / 9.0) * EPSILON_0 * (2.0 * q_over_m).sqrt() * voltage.powf(1.5) / gap_spacing.powi(2);
    total_current / j_max // m^2 Minimum area needed to pass this current.
}

fn apertures_needed(
    total_current: f64,
    per_aperture_area: f64,
    gap_spacing: f64,
    species_mw: f64,
    voltage: f64,
) -> f64 {
    let q_over_m = F / species_mw;
    let j_max =
        (4.0 / 9.0) * EPSILON_0 * (2.0 * q_over_m).sqrt() * voltage.powf(1.5) / gap_spacing.powi(2);
    total_current / (j_max * per_aperture_area)
}

fn funnel_convergence_half_angle_deg(
    accelerator_area: f64,
    combustion_inlet_area: f64,
    funnel_length: f64,
) -> f64 {
    let accel_radius = (accelerator_area / PI).sqrt();
    let inlet_radius = (combustion_inlet_area / PI).sqrt();
    ((accel_radius - inlet_radius) / funnel_length)
        .atan()
        .to_degrees()
}

fn solve_for_chamber_state(q_chamber: f64, ions: &Mixture, diluent: Option<&Mixture>) -> Mixture {
    // Mix ions and diluent.
    let mixture = if let Some(diluent) = diluent {
        ions.mix(diluent)
    } else {
        ions.clone()
    };
    let target_enthalpy = mixture.h_total + q_chamber; // J

    // Use HP mode to reequilibrate at the new enthalpy and pressure.
    // This allows species to dissociate/recombine as needed with the added heat.
    // HP mode wants H/R per unit feed mass [K·kmol/kg = K·mol/g], not J.
    let result = solve_for_products(
        &mixture.products,
        EquilibriumMode::HP {
            h_over_r: target_enthalpy / (R * mixture.feed_mass()),
            pressure_bar: mixture.pressure_bar,
        },
        true,
        true,
        None,
    );

    match result {
        Ok(state) => {
            let (h_total, s_total, n_total, avg_mw, avg_cp) =
                Mixture::state(&state.products, state.temperature_k, state.pressure_bar);

            Mixture {
                products: state.products,
                temperature_k: state.temperature_k,
                pressure_bar: state.pressure_bar,
                h_total,
                s_total,
                n_total,
                avg_mw,
                avg_cp,
            }
        }
        Err(err) => match err {
            EquilibriumError::FailedToConverge { iterations_used } => {
                panic!(
                    "solve_for_chamber_state failed to converge after {} iterations",
                    iterations_used
                )
            }
        },
    }
}

// Isentropic full-equilibrium expansion (CEA-style SP solve) from a chamber
// state to the nozzle exit pressure. The frozen-gamma temperature estimate
// breaks down for near-gamma=1 plasma: T barely drops across the nozzle and
// low-pressure dissociation can push exit enthalpy above the chamber's,
// making v_heat imaginary.
fn expand_to_exit(chamber: &Mixture, exit_pressure_bar: f64) -> Mixture {
    let result = solve_for_products(
        &chamber.products,
        EquilibriumMode::SP {
            // s_total is J/K on the chamber's mole basis; SP mode wants
            // S/R per unit feed mass [K*kmol/kg = K*mol/g].
            s_over_r: chamber.s_total / (R * chamber.feed_mass()),
            pressure_bar: exit_pressure_bar,
        },
        true,
        true,
        None,
    );
    match result {
        Ok(state) => {
            let (h_total, s_total, n_total, avg_mw, avg_cp) =
                Mixture::state(&state.products, state.temperature_k, state.pressure_bar);
            Mixture {
                products: state.products,
                temperature_k: state.temperature_k,
                pressure_bar: state.pressure_bar,
                h_total,
                s_total,
                n_total,
                avg_mw,
                avg_cp,
            }
        }
        Err(EquilibriumError::FailedToConverge { iterations_used }) => {
            panic!("expand_to_exit failed to converge after {iterations_used} iterations")
        }
    }
}

fn choked_throat_area(
    m_dot: f64,
    pressure_bar: f64,
    temperature: f64,
    gamma: f64,
    avg_mw: f64,
) -> f64 {
    let r_specific = R / avg_mw;
    let pressure_pa = pressure_bar * 1.0e5;
    let choke_factor = (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    m_dot / (pressure_pa * (gamma / (r_specific * temperature)).sqrt() * choke_factor)
}

fn choked_chamber_pressure_bar(
    m_dot: f64,
    area: f64,
    temperature: f64,
    gamma: f64,
    avg_mw: f64,
) -> f64 {
    let r_specific = R / avg_mw;
    let choke_factor = (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    let pressure_pa = m_dot / (area * (gamma / (r_specific * temperature)).sqrt() * choke_factor);
    pressure_pa / 1.0e5
}

#[allow(clippy::too_many_arguments)]
fn calculate_engine_output(
    field_voltage: f64,
    collision_theta_deg: f64,
    engine_power: f64,
    chamber_pressure: f64,
    target_area_ratio: f64,
    ions: &Mixture,
    // Moles of CH5+/OH- pair units in the `ions` feed basis (the CH4 amount
    // in its reactant list). Used to normalize `ions` to a per-pair basis so
    // the J/mol energy terms below survive changes to the feed mole ratios.
    ion_pair_moles: f64,
    // `diluent` amounts are interpreted per mole of ion pairs.
    diluent: &Mixture,
    gap_spacing: f64,
    aperture_diameter: f64,
    funnel_length: f64,
    coupling_efficiency: f64,
    fixed_total_throat_area: Option<f64>,
) -> (Option<f64>, Option<f64>) {
    let ions = ions.scale(1.0 / ion_pair_moles);

    // Start with CH4 + H20 ion feeds.
    let mut total_energy_in_per_ions_mole = 0.0;

    // CH4 + H2O => CH5+ + OH-
    total_energy_in_per_ions_mole += CH4_PARTIAL_OXIDATION_ENERGY;

    // Now each species is accelerated using an electric field.
    let electrostatic_energy = F * (2.0 * field_voltage);
    total_energy_in_per_ions_mole += electrostatic_energy;

    // This needs to account for the Hydrogen in the water, but this isn't very clean...
    // TODO DOES THIS need H2 MW added?
    let ions_mw = ions.feed_mass() * 1.0e-3; // feed_mass is grams -> kg per mole of ion pairs
    let ions_n_dot = engine_power / total_energy_in_per_ions_mole;
    let ions_m_dot = ions_n_dot * ions_mw;

    // Species are now at the reaction chamber, where they collide and combust.
    // collision_theta_deg controls how much kinetic energy is used to smash the species together,
    // versus how much is carried through out the nozzle.
    // A theta of 90 deg is a head on collision, great reaction rate, but kinetic energy turns into heat.
    // A theta of 0 degrees means the species run parallel and have a poor reaction rate, but kinetic energy is retained.

    let collision_theta_rad = collision_theta_deg * PI / 180.0;
    let v_y = (2.0 * electrostatic_energy / ions_mw).sqrt() * collision_theta_rad.cos();

    // Accelerated species enter chamber, where some of them can entrain with diluent.
    let diluent_m_dot = (diluent.n_total * diluent.avg_mw * 1.0e-3 * ions_n_dot).abs(); // kg/s (avg_mw is g/mol)
    let entrained_ion_m_dot = coupling_efficiency * ions_m_dot;
    let p_ion_entrained = entrained_ion_m_dot * v_y;
    let v_common = p_ion_entrained / (entrained_ion_m_dot + diluent_m_dot);
    // coupling_efficiency
    // 0 -> No coupling, accelerated species pass through without hitting anything.
    // 1 -> Full coupling, accelerated species share all ke with chamber species.

    // Dissipated kinetic energy becomes chamber heat.
    let ke_before = 0.5 * entrained_ion_m_dot * v_y.powi(2);
    let ke_after_full = 0.5 * (entrained_ion_m_dot + diluent_m_dot) * v_common.powi(2);
    let ke_dissipated_full = ke_before - ke_after_full;
    let ke_dissipated = coupling_efficiency * ke_dissipated_full;
    let q_mixing = ke_dissipated / (coupling_efficiency * ions_n_dot); // J/mol, adds to q_chamber

    // Any heat from collision plus reaction plus ionization. will go to raising the temperature of the products.
    // Products are CO + 3H2 (H2 injected from water electrolysis).
    // Assume all products start at room temperature.
    // Heat for the accelerated species that hit eachother and react, but do not
    // hit any other chamber species.
    let q_chamber = electrostatic_energy * collision_theta_rad.sin().powi(2)
        + CH4_PARTIAL_OXIDATION_ENERGY
        + CH4_FIRST_IONIZATION;
    let q_chamber_plasma = q_chamber * (1.0 - coupling_efficiency); // J/mol
    //  Heat for the accelerated species that hit eachother and react, then
    // fully mix with other chamber species.
    let q_chamber_diluent = q_chamber * coupling_efficiency + q_mixing; // J/mol

    let plasma = &ions;
    let plasma_species = plasma.scale(1.0 - coupling_efficiency);
    let chamber_plasma_state = solve_for_chamber_state(q_chamber_plasma, &plasma_species, None);
    let diluent_species = plasma.scale(coupling_efficiency);
    let chamber_diluent_state =
        solve_for_chamber_state(q_chamber_diluent, &diluent_species, Some(diluent));

    // Handle fast flow.
    // Mixture properties.
    let mixture_cp_mass_basis = chamber_plasma_state.avg_cp / chamber_plasma_state.avg_mw;
    let mixture_specific_gas_constant = R / chamber_plasma_state.avg_mw;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    let exit_pressure =
        exit_pressure_from_area_ratio(chamber_pressure, target_area_ratio, mixture_gamma);

    // Nozzle expansion. full equilibrium flow for upper bound Isp.
    let exit_temperature = chamber_plasma_state.temperature_k
        * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
    let plasma_species = expand_to_exit(&chamber_plasma_state, exit_pressure);
    let plasma_species_low = Mixture::new_with_frozen_reactants(
        &chamber_plasma_state.products,
        exit_temperature,
        exit_pressure,
    );
    let total_pool_mass = ions_mw * (1.0 - coupling_efficiency);
    // let total_pool_mass = propellant_mw + diluent.feed_mass();
    let v_heat =
        (2.0 * (chamber_plasma_state.h_total - plasma_species.h_total) / total_pool_mass).sqrt();
    let v_e = v_y + v_heat;
    let plasma_electrostatic_isp = v_y / G_0;
    let plasma_combustion_isp = v_heat / G_0;
    let plasma_isp = plasma_electrostatic_isp + plasma_combustion_isp;
    let plasma_engine_thrust = v_e * (ions_m_dot * (1.0 - coupling_efficiency));

    let v_heat_low = (2.0 * (chamber_plasma_state.h_total - plasma_species_low.h_total)
        / total_pool_mass)
        .sqrt();
    let v_e_low = v_y + v_heat_low;
    let plasma_combustion_isp_low = v_heat_low / G_0;
    let plasma_isp_low = plasma_electrostatic_isp + plasma_combustion_isp_low;
    let plasma_engine_thrust_low = v_e_low * (ions_m_dot * (1.0 - coupling_efficiency));

    println!(
        "[Fast] Chamber Plasma Temperature: {:.3} K",
        chamber_plasma_state.temperature_k
    );
    chamber_plasma_state.print_products();
    println!(
        "[Fast] Electrostatic Isp: {:.0} s",
        plasma_electrostatic_isp
    );
    println!("[Fast] Combustion Isp: {:.0} s", plasma_combustion_isp);
    println!("[Fast] Isp: {:.0} s", plasma_isp);
    println!("[Fast] Thrust: {:.3} kN", plasma_engine_thrust / 1.0e3);
    println!("[Fast] Lower Isp: {:.0} s", plasma_isp_low);
    println!(
        "[Fast] Lower Thrust: {:.3} kN",
        plasma_engine_thrust_low / 1.0e3
    );

    // Handle slow flow.
    // Mixture properties.
    let mixture_cp_mass_basis = chamber_diluent_state.avg_cp / chamber_diluent_state.avg_mw;
    let mixture_specific_gas_constant = R / chamber_diluent_state.avg_mw;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    let exit_pressure =
        exit_pressure_from_area_ratio(chamber_pressure, target_area_ratio, mixture_gamma);

    // Nozzle expansion. full equilibrium flow for upper bound Isp.
    let exit_temperature = chamber_diluent_state.temperature_k
        * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
    let diluent_species = expand_to_exit(&chamber_diluent_state, exit_pressure);
    let diluent_species_low = Mixture::new_with_frozen_reactants(
        &chamber_diluent_state.products,
        exit_temperature,
        exit_pressure,
    );
    let total_pool_mass = ions_mw * coupling_efficiency + diluent.feed_mass() * 1.0e-3; // kg
    let v_heat =
        (2.0 * (chamber_diluent_state.h_total - diluent_species.h_total) / total_pool_mass).sqrt();
    let v_e = v_common + v_heat;
    let slow_electrostatic_isp = v_common / G_0;
    let slow_combustion_isp = v_heat / G_0;
    let slow_isp = slow_electrostatic_isp + slow_combustion_isp;
    let slow_engine_thrust = v_e * (ions_m_dot * coupling_efficiency + diluent_m_dot);

    let v_heat_low = (2.0 * (chamber_diluent_state.h_total - diluent_species_low.h_total)
        / total_pool_mass)
        .sqrt();
    let v_e_low = v_common + v_heat_low;
    let slow_ecombustion_isp_low = v_heat_low / G_0;
    let slow_isp_low = slow_electrostatic_isp + slow_ecombustion_isp_low;
    let slow_engine_thrust_low = v_e_low * (ions_m_dot * coupling_efficiency + diluent_m_dot);

    println!(
        "[Slow] Chamber Temperature: {:.3} K",
        chamber_diluent_state.temperature_k
    );
    chamber_diluent_state.print_products();
    println!("[Slow] Electrostatic Isp: {:.0} s", slow_electrostatic_isp);
    println!("[Slow] Combustion Isp: {:.0} s", slow_combustion_isp);
    println!("[Slow] Isp: {:.0} s", slow_isp);
    println!("[Slow] Thrust: {:.3} kN", slow_engine_thrust / 1.0e3);
    println!("[Slow] Lower Isp: {:.0} s", slow_isp_low);
    println!(
        "[Slow] Lower Thrust: {:.3} kN",
        slow_engine_thrust_low / 1.0e3
    );

    println!("Ions m dot: {:.3} g/s", ions_m_dot * 1.0e3);
    println!("Diluent m dot: {:.3} g/s", diluent_m_dot * 1.0e3);
    let thrust_combined = slow_engine_thrust + plasma_engine_thrust;
    let v_combined = thrust_combined / (ions_m_dot + diluent_m_dot);
    let isp_combined = v_combined / G_0;
    println!("[Combined] Isp: {:.0} s", isp_combined);
    println!("[Combined] Thrust: {:.3} kN", thrust_combined / 1.0e3);

    // This can just use CH4_MW because O_MW == CH4_MW.
    let aperture_area = PI * aperture_diameter.powi(2) / 4.0;
    let total_current = ions_n_dot * F;
    let apertures_needed = apertures_needed(
        total_current,
        aperture_area,
        gap_spacing,
        CH4_MW,
        field_voltage,
    );
    let max_channel_area = required_max_channel_area(
        ions_n_dot, // After acceleration grids, ions are funneled
        // down into a tighter package, could be magnetic funnel,
        // or could use a concave accelerator plate to naturally squeeze them together.
        CH4_MW,
        field_voltage,
        chamber_pressure,
    );
    let max_channel_diameter = (4.0 * max_channel_area / PI).sqrt();

    let min_channel_area =
        min_channel_area_child_langmuir(total_current, CH4_MW, field_voltage, gap_spacing);

    let convergence_angle_deg =
        funnel_convergence_half_angle_deg(min_channel_area, max_channel_area, funnel_length);
    println!(
        "Max Channel Diameter: {:.3} mm",
        max_channel_diameter * 1.0e3
    );
    println!(
        "Accelerator Aperture Total Area: {:.3} m^2",
        min_channel_area
    );
    println!(
        "Combustion Chamber Inlet Area: {:.3} mm^2",
        max_channel_area * 1.0e6
    );
    let compression_ratio = min_channel_area / max_channel_area;
    println!(
        "Compression Ratio Aperture Area:Combustion Inlet {:.0}:{:.0}",
        compression_ratio, 1.0
    );
    println!("Funnel Length: {:.3} m", funnel_length);
    println!(
        "Funnel Convergence Half-Angle: {:.1} deg",
        convergence_angle_deg
    );
    println!("Total Channel Area: {:.3} m^2", {
        PI * (aperture_diameter / 2.0).powi(2) * apertures_needed
    });
    println!(
        "Configured Aperture Diameter: {:.3} mm",
        aperture_diameter * 1.0e3
    );
    println!("Configured Gap Spacing: {:.3} mm", gap_spacing * 1.0e3);
    println!("Apertures Needed: {:.0}", apertures_needed);

    // Nozzle Design
    let slow_m_dot = ions_m_dot * coupling_efficiency + diluent_m_dot;
    if let Some(fixed_area) = fixed_total_throat_area {
        let available_annulus_area = fixed_area - max_channel_area;
        let actual_chamber_pressure = choked_chamber_pressure_bar(
            slow_m_dot,
            available_annulus_area,
            chamber_diluent_state.temperature_k,
            mixture_gamma,
            chamber_diluent_state.avg_mw * 1.0e-3, // g/mol -> kg/mol
        );
        println!(
            "Actual Chamber Pressure: {:.3} bar",
            actual_chamber_pressure
        );
        (Some(actual_chamber_pressure), None)
    } else {
        let annulus_area_design = choked_throat_area(
            slow_m_dot,
            chamber_pressure,
            chamber_diluent_state.temperature_k,
            mixture_gamma,
            chamber_diluent_state.avg_mw * 1.0e-3, // g/mol -> kg/mol
        );
        let fixed_total_area_throat = max_channel_area + annulus_area_design;
        println!(
            "Fixed Throat Area: {:.3} mm^2",
            fixed_total_area_throat * 1.0e6
        );
        (None, Some(fixed_total_area_throat))
    }
}

// TODO RENAME OOGA BOOGA VARIABLES IN THIS like all the state stuff.
// TODO add some effiency number?
// Where the engine power is not accounting for
// losses at all
// Waste heat needs to account for this loss, and loss from fusion to electricity, and
// the heat from disassociation and water electrolysis
pub fn sweep_engine() {
    let voltage = 1_000.0; // Increase this for higher temps
    let engine_power = 50.0e6;
    let target_area_ratio = 100.0;
    let gap_spacing = 0.05e-3;
    let aperture_diameter = 0.5e-3;
    let funnel_length = 3.0; // meters
    let chamber_pressure = 10.0;
    let collision_theta_deg = 15.0;

    let ions = Mixture::new(
        &[(1.0, species!("CH4")), (1.0, species!("O"))],
        273.0,
        chamber_pressure,
    );

    let diluent = Mixture::new(
        &[(1.0, species!("H2")), (19.0, Species::H2O)],
        1_000.0,
        chamber_pressure,
    );
    // TODO seems like raising chamber pressure, increases temperature, which lets
    // me lower theta, which converses more axial velocity!
    let coupling_efficiency = 0.05;
    let fixed_total_throat_area = None; // Derive area
    println!("===== Thrust Mode =====");
    calculate_engine_output(
        voltage,
        collision_theta_deg,
        engine_power,
        chamber_pressure,
        target_area_ratio,
        &ions,
        1.0, // ion_pair_moles: the CH4 amount in `ions` above
        &diluent,
        gap_spacing,
        aperture_diameter,
        funnel_length,
        coupling_efficiency,
        fixed_total_throat_area,
    );

    let diluent = Mixture::new(&[(1.0, species!("H2"))], 1_000.0, chamber_pressure);

    let coupling_efficiency = 0.002; // Lower because the length of actual plasma in the chamber is less
    // So there is less time for the plasma to couple energy with the diluetnt.
    println!("===== Isp Mode =====");
    calculate_engine_output(
        voltage,
        collision_theta_deg,
        engine_power,
        chamber_pressure,
        target_area_ratio,
        &ions,
        1.0, // ion_pair_moles: the CH4 amount in `ions` above
        &diluent,
        gap_spacing,
        aperture_diameter,
        funnel_length,
        coupling_efficiency,
        None,
    );

    // Using methane clathrate is CH4•5.75H2O
    // where 1 H2O is split, O goes to ions, H2 is the diluent
    let diluent = Mixture::new(
        &[(1.0, species!("H2")), (4.75, species!("H2O"))],
        1_000.0,
        chamber_pressure,
    );

    let coupling_efficiency = 0.01; // Lower because the length of actual plasma in the chamber is less
    // So there is less time for the plasma to couple energy with the diluetnt.
    println!("===== Mixed Mode =====");
    calculate_engine_output(
        voltage,
        collision_theta_deg,
        engine_power,
        chamber_pressure,
        target_area_ratio,
        &ions,
        1.0, // ion_pair_moles: the CH4 amount in `ions` above
        &diluent,
        gap_spacing,
        aperture_diameter,
        funnel_length,
        coupling_efficiency,
        None,
    );
}
