use crate::nozzle::area_ratio::exit_pressure_from_area_ratio;
use crate::nozzle::propellant::PropellantState;
use crate::thermo::fluid_properties::ThermoReference;
use crate::{constants::*, nozzle::propellant};
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

fn solve_for_chamber_state(
    t_ions_start: f64,
    t_diluent_start: f64,
    t_allowed_max_chamber: f64,
    chamber_pressure: f64,
    q_chamber: f64,
    ions: &Propellant,
    diluent: Option<&Propellant>,
) -> (PropellantState, f64) {
    let mut t_low = t_ions_start;
    let mut t_high = t_allowed_max_chamber;
    let mut t_chamber = t_low + (t_high - t_low) / 2.0;
    // Ions and diluent can start out at different temperatures, so
    // get initial state separately.
    let mixture_start_state = ions.state(t_ions_start, chamber_pressure);
    let diluent_start_h = if let Some(d) = diluent {
        d.state(t_diluent_start, chamber_pressure).h_total
    } else {
        0.0
    };
    let starting_enthalpy = mixture_start_state.h_total + diluent_start_h;

    let mixture = if let Some(d) = diluent {
        &ions.mix(d)
    } else {
        ions
    };

    let mut chamber_state = PropellantState::default();

    for i in 0..200 {
        chamber_state = mixture.state(t_chamber, chamber_pressure);

        let delta_enthalpy = chamber_state.h_total - starting_enthalpy;

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
            panic!(
                "Failed to find chamber temperature. Last Guess: {:.3}",
                t_chamber
            );
        }
    }
    (chamber_state, t_chamber)
}

fn calculate_engine_output(
    field_voltage: f64,
    collision_theta_deg: f64,
    engine_power: f64,
    t_mixture_start: f64,
    t_diluent_start: f64,
    chamber_pressure: f64,
    target_area_ratio: f64,
    propellant: &Propellant,
    diluent: &Propellant,
    gap_spacing: f64,
    aperture_diameter: f64,
    funnel_length: f64,
    t_allowed_max_chamber: f64,
    coupling_efficiency: f64,
) {
    // It is kind of like there are two exhaust streams.
    // 1) Accelerated stream that does not mix with the diluent. (1.0 - coupling_efficiency)
    // 2) Accelerated stream that mixes with the all diluent. coupling_efficiency
    // Accelerated species only partially mix with diluent.

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

    // This needs to account for the Hydrogen in the water, but this isn't very clean...
    // TODO DOES THIS need H2 MW added?
    let propellant_mw = propellant.feed_mass();
    let propellant_n_dot = engine_power / total_energy_in_per_propellant_mole;
    let propellant_m_dot = propellant_n_dot * propellant_mw;

    // Species are now at the reaction chamber, where they collide and combust.
    // collision_theta_deg controls how much kinetic energy is used to smash the species together,
    // versus how much is carried through out the nozzle.
    // A theta of 90 deg is a head on collision, great reaction rate, but kinetic energy turns into heat.
    // A theta of 0 degrees means the species run parallel and have a poor reaction rate, but kinetic energy is retained.

    let collision_theta_rad = collision_theta_deg * PI / 180.0;
    let collision_species_mw = propellant.feed_mass();
    let v_y =
        (2.0 * electrostatic_energy / collision_species_mw).sqrt() * collision_theta_rad.cos();

    // Accelerated species enter chamber, where some of them can slam into diluent.
    let diluent_state_start = diluent.state(t_mixture_start, chamber_pressure);
    let diluent_m_dot =
        (diluent_state_start.n_total * diluent_state_start.avg_mw * propellant_n_dot).abs();
    let entrained_ion_m_dot = coupling_efficiency * propellant_m_dot;
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
    let q_chamber_mixing = ke_dissipated / (coupling_efficiency * propellant_n_dot); // J/mol, adds to q_chamber

    // Any heat from collision plus reaction plus ionization. will go to raising the temperature of the products.
    // Products are CO + 3H2 (H2 injected from water electrolysis).
    // Assume all products start at room temperature.

    // Heat for the accelerated species that hit eachother and react, but do not
    // hit any other chamber species.
    let q_chamber_fast = (electrostatic_energy * collision_theta_rad.sin().powi(2)
        + CH4_PARTIAL_OXIDATION_ENERGY
        + CH4_FIRST_IONIZATION)
        * (1.0 - coupling_efficiency); // J/mol
    //  Heat for the accelerated species that hit eachother and react, then
    // fully mix with other chamber species.
    let q_chamber_slow = (electrostatic_energy * collision_theta_rad.sin().powi(2)
        + CH4_PARTIAL_OXIDATION_ENERGY
        + CH4_FIRST_IONIZATION)
        * coupling_efficiency
        + q_chamber_mixing; // J/mol

    let fast_species = propellant.scale(1.0 - coupling_efficiency);
    let (chamber_fast_state, t_chamber_fast) = solve_for_chamber_state(
        t_mixture_start,
        t_diluent_start,
        20_000.0,
        chamber_pressure,
        q_chamber_fast,
        &fast_species,
        None,
    );
    let slow_species = propellant.scale(coupling_efficiency);
    let (chamber_slow_state, t_chamber_slow) = solve_for_chamber_state(
        t_mixture_start,
        t_diluent_start,
        t_allowed_max_chamber,
        chamber_pressure,
        q_chamber_slow,
        &slow_species,
        Some(diluent),
    );

    // Handle fast flow.
    // Mixture properties.
    let mixture_cp_mass_basis = chamber_fast_state.avg_cp / chamber_fast_state.avg_mw;
    let mixture_specific_gas_constant = R / chamber_fast_state.avg_mw;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    let exit_pressure =
        exit_pressure_from_area_ratio(chamber_pressure, target_area_ratio, mixture_gamma);

    // Nozzle expansion. Assume frozen flow for lower bound Isp.
    let exit_temperature = t_chamber_fast
        * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
    let exit_state = fast_species.state(exit_temperature, exit_pressure);
    let total_pool_mass = propellant_mw * (1.0 - coupling_efficiency);
    // let total_pool_mass = propellant_mw + diluent.feed_mass();
    let v_heat = (2.0 * (chamber_fast_state.h_total - exit_state.h_total) / total_pool_mass).sqrt();
    let v_e = v_y + v_heat;
    let fast_electrostatic_isp = v_y / G_0;
    let fast_ecombustion_isp = v_heat / G_0;
    let fast_isp = fast_electrostatic_isp + fast_ecombustion_isp;

    let fast_engine_thrust = v_e * (propellant_m_dot * (1.0 - coupling_efficiency));

    println!("[Fast] Chamber Temperature: {:.3} K", t_chamber_fast);
    for (i, specie) in propellant.species.iter().enumerate() {
        println!(
            "[Fast] ⍺ {}: {:.2} %",
            specie.1.symbol(),
            chamber_fast_state.alphas[i] * 100.0
        );
    }
    println!("[Fast] Electrostatic Isp: {:.0} s", fast_electrostatic_isp);
    println!("[Fast] Combustion Isp: {:.0} s", fast_ecombustion_isp);
    println!("[Fast] Isp: {:.0} s", fast_isp);
    println!("[Fast] Thrust: {:.3} kN", fast_engine_thrust / 1.0e3);

    // Handle slow flow.
    // Mixture properties.
    let mixture_cp_mass_basis = chamber_slow_state.avg_cp / chamber_slow_state.avg_mw;
    let mixture_specific_gas_constant = R / chamber_slow_state.avg_mw;
    let mixture_gamma =
        mixture_cp_mass_basis / (mixture_cp_mass_basis - mixture_specific_gas_constant);

    let exit_pressure =
        exit_pressure_from_area_ratio(chamber_pressure, target_area_ratio, mixture_gamma);

    // Nozzle expansion. Assume frozen flow for lower bound Isp.
    let exit_temperature = t_chamber_slow
        * (exit_pressure / chamber_pressure).powf((mixture_gamma - 1.0) / mixture_gamma);
    let slow_species = slow_species.mix(diluent);
    let exit_state = slow_species.state(exit_temperature, exit_pressure);
    let total_pool_mass = propellant_mw * coupling_efficiency + diluent.feed_mass();
    let v_heat = (2.0 * (chamber_slow_state.h_total - exit_state.h_total) / total_pool_mass).sqrt();
    let v_e = v_common + v_heat;
    let slow_electrostatic_isp = v_common / G_0;
    let slow_combustion_isp = v_heat / G_0;
    let slow_isp = slow_electrostatic_isp + slow_combustion_isp;

    let slow_engine_thrust = v_e * (propellant_m_dot * coupling_efficiency + diluent_m_dot);

    println!("[Slow] Chamber Temperature: {:.3} K", t_chamber_slow);
    for (i, specie) in slow_species.species.iter().enumerate() {
        println!(
            "[Slow] ⍺ {}: {:.2} %",
            specie.1.symbol(),
            chamber_slow_state.alphas[i] * 100.0
        );
    }
    println!("[Slow] Electrostatic Isp: {:.0} s", slow_electrostatic_isp);
    println!("[Slow] Combustion Isp: {:.0} s", slow_combustion_isp);
    println!("[Slow] Isp: {:.0} s", slow_isp);
    println!("[Slow] Thrust: {:.3} kN", slow_engine_thrust / 1.0e3);

    println!("Propellant m dot: {:.3} g/s", propellant_m_dot * 1.0e3);
    println!("Diluent m dot: {:.3} g/s", diluent_m_dot * 1.0e3);
    let thrust_combined = slow_engine_thrust + fast_engine_thrust;
    let v_combined = thrust_combined / (propellant_m_dot + diluent_m_dot);
    let isp_combined = v_combined / G_0;
    println!("[Combined] Isp: {:.0} s", isp_combined);
    println!("[Combined] Thrust: {:.3} kN", thrust_combined / 1.0e3);

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

    println!("[Slow] Avg MW: {:.4} kg/mol", chamber_slow_state.avg_mw);
    println!("[Slow] Gamma: {:.4}", mixture_gamma);
    println!(
        "[Slow] Mass flow: {:.3} g/s",
        (propellant_m_dot * coupling_efficiency + diluent_m_dot) * 1.0e3
    );
}

// TODO thermo reference for water at higher temperatures?
// TODO add some effiency number?
// Where the engine power is not accounting for
// losses at all
// Waste heat needs to account for this loss, and loss from fusion to electricity, and
// the heat from disassociation and water electrolysis
pub fn sweep_engine() {
    let thermo_reference = ThermoReference::new();
    let propellant = Propellant::new(
        &thermo_reference,
        vec![
            (1.0, Species::CO, vec![(1.0, Species::C), (1.0, Species::O)]),
            (2.0, Species::H2, vec![(2.0, Species::H)]),
        ],
    );
    let diluent = Propellant::new(
        &thermo_reference,
        vec![(1.0, Species::H2, vec![(2.0, Species::H)])],
    );

    let voltage = 400.0;
    let collision_theta_deg = 5.0;
    let engine_power = 50.0e6;
    let t_mixture_start = 300.0;
    let t_diluent_start = 1_000.0;
    let chamber_pressure = 25.0;
    let target_area_ratio = 100.0;
    let gap_spacing = 0.05e-3;
    let aperture_diameter = 0.5e-3;
    let funnel_length = 7.0; // meters
    let t_allowed_max_chamber = 6_000.0; // When using water diluent
    let coupling_efficiency = 0.005; // Lower because chamber only has 1 mole H2
    println!("===== Isp Mode =====");
    calculate_engine_output(
        voltage,
        collision_theta_deg,
        engine_power,
        t_mixture_start,
        t_diluent_start,
        chamber_pressure,
        target_area_ratio,
        &propellant,
        &diluent,
        gap_spacing,
        aperture_diameter,
        funnel_length,
        t_allowed_max_chamber,
        coupling_efficiency,
    );

    let diluent = Propellant::new(
        &thermo_reference,
        vec![
            (1.0, Species::H2, vec![(2.0, Species::H)]),
            (
                5.0,
                Species::H2O,
                vec![(1.0, Species::H2), (1.0, Species::O)],
            ),
        ],
    );

    let voltage = 400.0;
    let collision_theta_deg = 5.0;
    let engine_power = 50.0e6;
    let t_mixture_start = 300.0;
    let t_diluent_start = 1_000.0;
    let chamber_pressure = 25.0;
    let target_area_ratio = 100.0;
    let gap_spacing = 0.05e-3;
    let aperture_diameter = 0.5e-3;
    let funnel_length = 7.0; // meters
    let t_allowed_max_chamber = 6_000.0; // When using water diluent
    let coupling_efficiency = 0.03; // Higher because chamber has 6 moles of diluent plus disassociation
    println!("===== Thrust Mode =====");
    calculate_engine_output(
        voltage,
        collision_theta_deg,
        engine_power,
        t_mixture_start,
        t_diluent_start,
        chamber_pressure,
        target_area_ratio,
        &propellant,
        &diluent,
        gap_spacing,
        aperture_diameter,
        funnel_length,
        t_allowed_max_chamber,
        coupling_efficiency,
    );
}
