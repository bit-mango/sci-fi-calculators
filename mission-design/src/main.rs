mod lambert_solver;

use std::f64::consts::PI;

use crate::lambert_solver::solver::solve;

const AU: f64 = 1.495978707e11; // m
const R_MARS: f64 = 1.524 * AU;
const R_EARTH: f64 = 1.0 * AU;

const ENGINE_ISP: f64 = 908.0;
const ENGINE_THRUST: f64 = 8.912e3; // 8.8912 kN
const G_0: f64 = 9.81;
const DRY_MASS: f64 = 300.0e3; // 500 mt
const CARGO_MASS: f64 = 500.0e3; // 500 mt

fn main() {
    let tof = 90.0 * 86_400.0;
    println!("Mission Length: {:.0} days", tof / 86_400.0);
    let mut best: Option<(f64, f64, f64, f64)> = None; // total, dv_mars, dv_earth, theta_deg
    let mut worst: Option<(f64, f64, f64, f64)> = None; // total, dv_mars, dv_earth, theta_deg

    let mut theta_deg = 20.0;
    while theta_deg < 340.0 {
        let theta = theta_deg * PI / 180.0;
        let (dv_mars, dv_earth) = solve(tof, R_MARS, R_EARTH, theta);
        let total = dv_mars + dv_earth;
        if total.is_finite() {
            if best.is_none() || total < best.unwrap().0 {
                best = Some((total, dv_mars, dv_earth, theta_deg));
            }
            if worst.is_none() || total > worst.unwrap().0 {
                worst = Some((total, dv_mars, dv_earth, theta_deg));
            }
        }
        theta_deg += 1.0;
    }

    // let (total, dv_mars, dv_earth, theta_deg) = worst.unwrap();
    // println!("Worst theta = {:.0} deg", theta_deg);
    // println!("Delta V Mars: {:.3} km/s", dv_mars / 1.0e3);
    // println!("Delta V Earth: {:.3} km/s", dv_earth / 1.0e3);
    // println!("Total: {:.3} km/s", total / 1.0e3);

    let (_, _, _, theta_deg) = best.unwrap();
    // Choose a theta that is close to the best but not perfect.
    let theta_deg = theta_deg + 10.0;
    let theta = theta_deg * PI / 180.0;
    let (dv_mars, dv_earth) = solve(tof, R_MARS, R_EARTH, theta);
    let total = dv_mars + dv_earth;

    println!("Mission theta = {:.0} deg", theta_deg);
    println!("Delta V Mars: {:.3} km/s", dv_mars / 1.0e3);
    println!("Delta V Earth: {:.3} km/s", dv_earth / 1.0e3);
    println!("Total: {:.3} km/s", total / 1.0e3);

    let exhaust_velocity = ENGINE_ISP * G_0;

    // Work backwards between every burn stage.
    let m_final = DRY_MASS + CARGO_MASS;
    let stage_1_m_initial = calculate_stage_initial_mass(exhaust_velocity, dv_earth, m_final);
    let stage_0_m_initial =
        calculate_stage_initial_mass(exhaust_velocity, dv_mars, stage_1_m_initial);

    let m_dot = ENGINE_THRUST / exhaust_velocity;

    let stage_0_propellant_used = stage_0_m_initial - stage_1_m_initial;
    let stage_1_propellant_used = stage_1_m_initial - m_final;
    let m_propellant = stage_0_m_initial - m_final;
    // As mission planned.
    println!("Propellant m_dot: {:.3} kg/s", m_dot);
    println!("Total Propellant Mass: {:.3} mt", m_propellant / 1.0e3);
    println!(
        "Stage 0 Propellant Used: {:.3} mt",
        stage_0_propellant_used / 1.0e3
    );
    println!(
        "Stage 1 Propellant Used: {:.3} mt",
        stage_1_propellant_used / 1.0e3
    );
    println!(
        "Stage 0 Burn Time: {:.3} days",
        (stage_0_propellant_used / m_dot) / 86_400.0
    );
    println!(
        "Stage 1 Burn Time: {:.3} days",
        (stage_1_propellant_used / m_dot) / 86_400.0
    );

    // Lose 75% of radiators.
    let exhaust_velocity = 809.0 * G_0;

    // Work backwards between every burn stage.
    let m_final = DRY_MASS + CARGO_MASS;
    let stage_1_m_initial = calculate_stage_initial_mass(exhaust_velocity, dv_earth, m_final);
    let stage_0_m_initial =
        calculate_stage_initial_mass(exhaust_velocity, dv_mars, stage_1_m_initial);

    let m_dot = 7_936.0 / exhaust_velocity;

    let stage_0_propellant_used = stage_0_m_initial - stage_1_m_initial;
    let stage_1_propellant_used = stage_1_m_initial - m_final;
    let m_propellant = stage_0_m_initial - m_final;
    // As mission planned.
    println!("Propellant m_dot: {:.3} kg/s", m_dot);
    println!("Total Propellant Mass: {:.3} mt", m_propellant / 1.0e3);
    println!(
        "Stage 0 Propellant Used: {:.3} mt",
        stage_0_propellant_used / 1.0e3
    );
    println!(
        "Stage 1 Propellant Used: {:.3} mt",
        stage_1_propellant_used / 1.0e3
    );
    println!(
        "Stage 0 Burn Time: {:.3} days",
        (stage_0_propellant_used / m_dot) / 86_400.0
    );
    println!(
        "Stage 1 Burn Time: {:.3} days",
        (stage_1_propellant_used / m_dot) / 86_400.0
    );

    // Back up using urea.
    let exhaust_velocity = 697.0 * G_0;

    // Work backwards between every burn stage.
    // Assume 300 mt of Urea is used as propellant.
    let m_final = DRY_MASS + CARGO_MASS - 300.0e3;
    let stage_1_m_initial = calculate_stage_initial_mass(exhaust_velocity, dv_earth, m_final);
    let stage_0_m_initial =
        calculate_stage_initial_mass(exhaust_velocity, dv_mars, stage_1_m_initial);

    let m_dot = 13_666.0 / exhaust_velocity;

    let stage_0_propellant_used = stage_0_m_initial - stage_1_m_initial;
    let stage_1_propellant_used = stage_1_m_initial - m_final;
    let m_propellant = stage_0_m_initial - m_final;
    // As mission planned.
    println!("Propellant m_dot: {:.3} kg/s", m_dot);
    println!("Total Propellant Mass: {:.3} mt", m_propellant / 1.0e3);
    println!(
        "Stage 0 Propellant Used: {:.3} mt",
        stage_0_propellant_used / 1.0e3
    );
    println!(
        "Stage 1 Propellant Used: {:.3} mt",
        stage_1_propellant_used / 1.0e3
    );
    println!(
        "Stage 0 Burn Time: {:.3} days",
        (stage_0_propellant_used / m_dot) / 86_400.0
    );
    println!(
        "Stage 1 Burn Time: {:.3} days",
        (stage_1_propellant_used / m_dot) / 86_400.0
    );
}

fn calculate_stage_initial_mass(
    exhaust_velocity: f64,
    stage_delta_v: f64,
    stage_m_final: f64,
) -> f64 {
    let mass_ratio = (stage_delta_v / exhaust_velocity).exp();
    let stage_m_initial = stage_m_final * mass_ratio;
    stage_m_initial
}
