use std::f64::{self, consts::PI};

const MICRO_SUN: f64 = 1.32712440018e20; // m^3 /  s^2
const MICRO_MARS: f64 = 4.282837e13; // m^3 /  s^2
const MICRO_EARTH: f64 = 3.986004e14; // m^3 /  s^2
const R_PARKING_MARS: f64 = 3_789.5e3; // m, 400 km altitude.
const R_PARKING_EARTH: f64 = 6_678.0e3; // m, 300 km altitude.
const AU: f64 = 1.495978707e11; // m
const R_MARS: f64 = 1.524 * AU;
const R_EARTH: f64 = 1.0 * AU;

fn calculate_circular_orbit_speed(r: f64) -> f64 {
    (MICRO_SUN / r).sqrt()
}

#[derive(Clone, Copy)]
struct Vec2 {
    x: f64,
    y: f64,
}
impl Vec2 {
    fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    fn sub(self, o: Vec2) -> Vec2 {
        Vec2::new(self.x - o.x, self.y - o.y)
    }
    fn scale(self, s: f64) -> Vec2 {
        Vec2::new(self.x * s, self.y * s)
    }
    fn mag(self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}

pub fn solve(target_tof: f64, starting_r: f64, ending_r: f64, delta_theta: f64) -> (f64, f64) {
    let c_z = |z: f64| {
        if z > 0.0 {
            (1.0 - z.sqrt().cos()) / z
        } else if z < 0.0 {
            ((-z).sqrt().cosh() - 1.0) / (-z)
        } else {
            0.5
        }
    };
    let s_z = |z: f64| {
        if z > 0.0 {
            let sz = z.sqrt();
            (sz - sz.sin()) / sz.powi(3)
        } else if z < 0.0 {
            let sz = (-z).sqrt();
            (sz.sinh() - sz) / sz.powi(3)
        } else {
            1.0 / 6.0
        }
    };

    let a = delta_theta.sin() * (starting_r * ending_r / (1.0 - delta_theta.cos())).sqrt();
    let mut z_low = -4.0 * PI * PI;
    let mut z_high = 4.0 * PI * PI;
    let mut z = 0.0;

    let (mut f, mut g, mut g_dot) = (0.0, 0.0, 0.0);

    for i in 0..200 {
        let y = starting_r + ending_r + a * (z * s_z(z) - 1.0) / (c_z(z).sqrt());

        if a > 0.0 && y < 0.0 {
            // Bump z up.
            z_low = z;
            z = z_low + (z_high - z_low) / 2.0;
            continue;
        }

        let x = (y / c_z(z)).sqrt();
        let t = (x.powi(3) * s_z(z) + a * y.sqrt()) / MICRO_SUN.sqrt();
        if t < target_tof {
            z_low = z;
        } else {
            z_high = z;
        }

        if (t - target_tof).abs() < 1.0e-3 {
            f = 1.0 - y / starting_r;
            g = a * (y / MICRO_SUN).sqrt();
            g_dot = 1.0 - y / ending_r;
            break;
        }
        z = z_low + (z_high - z_low) / 2.0;
        if i == 199 {
            panic!("Failed to find solution!");
        }
    }

    let r1_vec = Vec2::new(starting_r, 0.0);
    let r2_vec = Vec2::new(ending_r * delta_theta.cos(), ending_r * delta_theta.sin());

    let v1_vec = r2_vec.sub(r1_vec.scale(f)).scale(1.0 / g);
    let v2_vec = r2_vec.scale(g_dot).sub(r1_vec).scale(1.0 / g);

    let v_mars_circ_vec = Vec2::new(0.0, calculate_circular_orbit_speed(R_MARS));
    let v_earth_circ_vec = Vec2::new(
        -delta_theta.sin() * calculate_circular_orbit_speed(R_EARTH),
        delta_theta.cos() * calculate_circular_orbit_speed(R_EARTH),
    );

    let vinf_depart = v1_vec.sub(v_mars_circ_vec).mag();
    let vinf_arrive = v2_vec.sub(v_earth_circ_vec).mag();

    let v_hyp_mars = (vinf_depart.powi(2) + 2.0 * MICRO_MARS / R_PARKING_MARS).sqrt();
    let dv_mars = v_hyp_mars - (MICRO_MARS / R_PARKING_MARS).sqrt();
    let v_hyp_earth = (vinf_arrive.powi(2) + 2.0 * MICRO_EARTH / R_PARKING_EARTH).sqrt();
    let dv_earth = v_hyp_earth - (MICRO_EARTH / R_PARKING_EARTH).sqrt();

    (dv_mars, dv_earth)
}
