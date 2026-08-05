fn area_ratio_from_mach(mach: f64, gamma: f64) -> f64 {
    let term = 2.0 / (gamma + 1.0) * (1.0 + mach.powi(2) * (gamma - 1.0) / 2.0);
    let exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    (1.0 / mach) * term.powf(exponent)
}

fn solve_exit_mach(target_area_ratio: f64, gamma: f64) -> f64 {
    // Need to find a supersonic mach number: M > 1 bisect between a number just above 1, and a larger upper bound.
    let mut mach_low = 1.0001;
    let mut mach_high = 100.0;

    let mut mach_mid = 0.0;
    for i in 0..100 {
        mach_mid = mach_low + (mach_high - mach_low) / 2.0;
        let ratio_mid = area_ratio_from_mach(mach_mid, gamma);

        if (ratio_mid - target_area_ratio).abs() < 1.0e-6 {
            break;
        }

        if ratio_mid < target_area_ratio {
            mach_low = mach_mid;
        } else {
            mach_high = mach_mid;
        }

        if i == 99 {
            panic!("Failed to find mach");
        }
    }

    mach_mid
}

pub fn exit_pressure_from_area_ratio(
    chamber_pressure: f64,
    target_area_ratio: f64,
    gamma: f64,
) -> f64 {
    let mach_exit = solve_exit_mach(target_area_ratio, gamma);
    let pressure_ratio =
        (1.0 + mach_exit.powi(2) * (gamma - 1.0) / 2.0).powf(-gamma / (gamma - 1.0));
    chamber_pressure * pressure_ratio
}
