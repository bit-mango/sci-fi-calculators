use crate::constants::*;
use crate::nozzle::area_ratio::exit_pressure_from_area_ratio;
use crate::nozzle::propellant::PropellantState;
use crate::thermo::fluid_properties::ThermoReference;
use std::f64::consts::PI;

use crate::nozzle::propellant::{Propellant, Species};

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

fn calculate_engine_output(
    field_voltage: f64,
    collision_theta_deg: f64,
    engine_power: f64,
    t_mixture_start: f64,
    chamber_pressure: f64,
    target_area_ratio: f64,
    propellant_start: &Propellant,
    diluent: &Propellant,
    gap_spacing: f64,
    aperture_diameter: f64,
    funnel_length: f64,
) {
    let propellant = propellant_start.mix(diluent);
    // Start with CH4 + H2O propellant feeds. 1 mole propellant is CH4 + H2O
    let mut total_energy_in_per_propellant_mole = 0.0;

    // First we need to split the water, and disassociate the resulting O2.
    total_energy_in_per_propellant_mole += H2O_ELECTROLYSIS + O2_DISSOCIATION;

    // Now we need to ionize the CH4 and give the electron to the O.
    // TODO why do we not subtract the O electron affinity here?
    total_energy_in_per_propellant_mole += CH4_FIRST_IONIZATION;

    // Now each species is accelerated using an electric field.
    let electrostatic_energy = F * (2.0 * field_voltage); // Each species accelerated through field.
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
    let collision_species_mw = CH4_MW + O_MW;
    let v_y =
        (2.0 * electrostatic_energy / collision_species_mw).sqrt() * collision_theta_rad.cos();

    // Any heat from collision plus reaction plus ionization. will go to raising the temperature of the products.
    // Products are CO + 3H2 (H2 injected from water electrolysis).
    // Assume all products start at room temperature.
    let q_chamber = electrostatic_energy * collision_theta_rad.sin().powi(2)
        + CH4_PARTIAL_OXIDATION_ENERGY
        + CH4_FIRST_IONIZATION; // J/mol
    // Use bisection method to estimate
    let mut t_low = t_mixture_start;
    let mut t_high = 20_000.0;
    let mut t_chamber = t_low + (t_high - t_low) / 2.0;
    let start_state = propellant.state(t_mixture_start, chamber_pressure);

    let mut chamber_state = PropellantState::default();

    for i in 0..200 {
        // println!(
        //     "Trying Temperature Range: {:.3} K <-> {:.3} K",
        //     t_low, t_high
        // );
        chamber_state = propellant.state(t_chamber, chamber_pressure);

        let delta_enthalpy = chamber_state.h_total - start_state.h_total;

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
        if i == 199 {
            panic!("Failed to find chamber temperature.")
        }
    }

    // Mixture properties.
    let mixture_cp_mass_basis = chamber_state.avg_cp / chamber_state.avg_mw;
    let mixture_specific_gas_constant = R / chamber_state.avg_mw;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    let exit_pressure =
        exit_pressure_from_area_ratio(chamber_pressure, target_area_ratio, mixture_gamma);

    // Nozzle expansion. Assume frozen flow for lower bound Isp.
    let exit_temperature =
        t_chamber * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
    let exit_state = propellant.state(exit_temperature, exit_pressure);
    let total_pool_mass = propellant_mw + diluent.feed_mass();
    let v_heat = (2.0 * (chamber_state.h_total - exit_state.h_total) / total_pool_mass).sqrt();
    let v_e = v_y + v_heat;
    let electrostatic_isp = v_y / G_0;
    let ecombustion_isp = v_heat / G_0;
    let isp = electrostatic_isp + ecombustion_isp;

    let diluent_state_start = diluent.state(t_mixture_start, chamber_pressure);
    let diluent_m_dot =
        (diluent_state_start.n_total * diluent_state_start.avg_mw * propellant_n_dot).abs();
    let engine_thrust = v_e * (propellant_m_dot + diluent_m_dot);

    // This can just use CH4_MW because O_MW == CH4_MW.
    let aperture_area = PI * aperture_diameter.powi(2) / 4.0;
    let total_current = propellant_n_dot * F;
    let apertures_needed = apertures_needed(
        total_current,
        aperture_area,
        gap_spacing,
        CH4_MW,
        field_voltage,
    );
    let max_channel_area = required_max_channel_area(
        propellant_n_dot, // After acceleration grids, ions are funneled
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

    println!("Chamber Temperature: {:.3} K", t_chamber);
    println!("Alpha CO: {:.2} %", chamber_state.alphas[0] * 100.0);
    println!("Alpha H2: {:.2} %", chamber_state.alphas[1] * 100.0);
    println!("Electrostatic Isp: {:.0} s", electrostatic_isp);
    println!("Combustion Isp: {:.0} s", ecombustion_isp);
    println!("Isp: {:.0} s", isp);
    println!("Thrust: {:.3} kN", engine_thrust / 1.0e3);
    println!("Propellant m dot: {:.3} g/s", propellant_m_dot * 1.0e3);
    println!("Diluent m dot: {:.3} g/s", diluent_m_dot * 1.0e3);
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
}

pub fn sweep_engine() {
    let thermo_reference = ThermoReference::new();
    let propellant = Propellant::new(
        &thermo_reference,
        vec![
            (1.0, Species::CO, vec![(1.0, Species::C), (1.0, Species::O)]),
            (3.0, Species::H2, vec![(2.0, Species::H)]),
        ],
    );
    // TODO should let you specify a diluent inlet temperature...
    // Think if I literally just split the start state into 2 h_total calls
    // one for the material coming from the accelerator and one for diluent which
    // could be much hotter than room temp.
    let diluent = Propellant::new(
        &thermo_reference,
        vec![
            // (1.0, Species::CO, vec![(1.0, Species::C), (1.0, Species::O)]),
            // (3.0, Species::H2, vec![(2.0, Species::H)]),
        ],
    );

    let voltages = vec![300.0, 400.0, 500.0];
    let collision_theta_deg = 5.0;
    let engine_power = 50.0e6;
    let t_mixture_start = 300.0;
    let chamber_pressure = 25.0;
    let target_area_ratio = 100.0;
    let gap_spacing = 0.05e-3;
    let aperture_diameter = 0.5e-3;
    let funnel_length = 7.0; // meters
    for voltage in voltages {
        println!("===== Voltage: {:.0} V =====", voltage);
        calculate_engine_output(
            voltage,
            collision_theta_deg,
            engine_power,
            t_mixture_start,
            chamber_pressure,
            target_area_ratio,
            &propellant,
            &diluent,
            gap_spacing,
            aperture_diameter,
            funnel_length,
        );
    }
}
