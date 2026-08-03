use crate::{Constituent, Species};
use std::fmt;

const R: f64 = 8.314; // Ideal Gas Constant [J/mol•K]

impl fmt::Display for Species {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let d = self.data();
        writeln!(f, "{} (phase {})", d.symbol(), d.phase())?;
        writeln!(f, "  mw: {} g/mol, ΔHf: {} J/mol", d.mw(), d.h_formation())?;
        writeln!(f, "  constituents: {:?}", d.constituents())?;
        writeln!(f, "  temperature ranges: {:?}", d.temperature_data())?;
        write!(f, "  sibling phases: {:?}", self.phases())
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseSpeciesError;

#[derive(Debug)]
pub struct SpeciesData<const C: usize, const T: usize> {
    pub symbol: &'static str,
    pub constituents: [(f64, Constituent); C],
    pub temperature_data: [(f64, f64, [f64; 9]); T],
    pub mw: f64,
    pub h_formation: f64,
    /// 0 = gas, non-zero = condensed. The specific non-zero value is not
    /// semantically meaningful (NASA's raw record ordinal, not a phase ID)
    /// and is arbitrary after lambda-transition merging.
    pub phase: u8,
}

fn find_coefficients(ranges: &[(f64, f64, [f64; 9])], temperature_k: f64) -> &(f64, f64, [f64; 9]) {
    // Check if temperature_k is greater than the maximum range we have.
    if let Some(last) = ranges.last()
        && temperature_k >= last.1
    {
        return last;
    }
    // Check if temperature_k is less than the minimum range we have.
    if let Some(first) = ranges.first()
        && temperature_k <= first.0
    {
        return first;
    }
    // Falls within a range we have.
    for coeff in ranges.iter() {
        if temperature_k >= coeff.0 && temperature_k <= coeff.1 {
            return coeff;
        }
    }

    panic!("Not sure how we got here");
}

pub trait AnySpeciesData {
    fn symbol(&self) -> &str;
    fn constituents(&self) -> &[(f64, Constituent)];
    fn temperature_data(&self) -> &[(f64, f64, [f64; 9])];
    fn valid_temperature_range(&self) -> (f64, f64);
    fn mw(&self) -> f64;
    fn h_formation(&self) -> f64;
    fn phase(&self) -> u8;
    fn cp(&self, temperature_k: f64) -> f64 {
        let coeff = find_coefficients(self.temperature_data(), temperature_k);
        let t_1 = temperature_k;
        let t_2 = t_1 * t_1;
        let t_3 = t_2 * t_1;
        let t_4 = t_2 * t_2;

        let res = coeff.2[0] / t_2
            + coeff.2[1] / t_1
            + coeff.2[2]
            + coeff.2[3] * t_1
            + coeff.2[4] * t_2
            + coeff.2[5] * t_3
            + coeff.2[6] * t_4;

        res * R
    }
    fn h(&self, temperature_k: f64) -> f64 {
        let coeff = find_coefficients(self.temperature_data(), temperature_k);
        let t_1 = temperature_k;
        let t_2 = t_1 * t_1;
        let t_3 = t_2 * t_1;
        let t_4 = t_2 * t_2;

        let res = -coeff.2[0] / t_2
            + coeff.2[1] * t_1.ln() / t_1
            + coeff.2[2]
            + coeff.2[3] * t_1 / 2.0
            + coeff.2[4] * t_2 / 3.0
            + coeff.2[5] * t_3 / 4.0
            + coeff.2[6] * t_4 / 5.0
            + coeff.2[7] / t_1;

        res * R * t_1
    }
    fn s(&self, temperature_k: f64) -> f64 {
        let coeff = find_coefficients(self.temperature_data(), temperature_k);
        let t_1 = temperature_k;
        let t_2 = t_1 * t_1;
        let t_3 = t_2 * t_1;
        let t_4 = t_2 * t_2;

        let res = -coeff.2[0] / (2.0 * t_2) - coeff.2[1] / t_1
            + coeff.2[2] * t_1.ln()
            + coeff.2[3] * t_1
            + coeff.2[4] * t_2 / 2.0
            + coeff.2[5] * t_3 / 3.0
            + coeff.2[6] * t_4 / 4.0
            + coeff.2[8];

        res * R
    }

    fn cp_over_r(&self, temperature_k: f64) -> f64 {
        let (_, _, coeff) = find_coefficients(self.temperature_data(), temperature_k);
        let mut cp = coeff[6];
        cp = temperature_k * cp + coeff[5];
        cp = temperature_k * cp + coeff[4];
        cp = temperature_k * cp + coeff[3];
        cp = temperature_k * cp + coeff[2];
        cp = temperature_k * cp + coeff[1];
        cp = temperature_k * cp + coeff[0];

        cp / (temperature_k * temperature_k)
    }

    fn h_over_r(&self, temperature_k: f64) -> f64 {
        let (_, _, coeff) = find_coefficients(self.temperature_data(), temperature_k);
        let mut h = coeff[6] / 5.0;
        h = temperature_k * h + coeff[5] / 4.0;
        h = temperature_k * h + coeff[4] / 3.0;
        h = temperature_k * h + coeff[3] / 2.0;
        h = temperature_k * h + coeff[2];
        h = temperature_k * h + coeff[7];
        h = temperature_k * h - coeff[0];

        h / temperature_k + coeff[1] * temperature_k.ln()
    }

    fn s_over_r(&self, temperature_k: f64) -> f64 {
        let (_, _, coeff) = find_coefficients(self.temperature_data(), temperature_k);
        let mut s = coeff[6] / 4.0;
        s = temperature_k * s + coeff[5] / 3.0;
        s = temperature_k * s + coeff[4] / 2.0;
        s = temperature_k * s + coeff[3];
        s = temperature_k * s + coeff[8];
        s = temperature_k * s - coeff[1];
        s = temperature_k * s - coeff[0] / 2.0;

        s / (temperature_k * temperature_k) + coeff[2] * temperature_k.ln()
    }

    fn u_over_r(&self, temperature_k: f64) -> f64 {
        self.h_over_r(temperature_k) - temperature_k
    }
}

impl<const C: usize, const T: usize> AnySpeciesData for SpeciesData<C, T> {
    fn symbol(&self) -> &str {
        self.symbol
    }
    fn constituents(&self) -> &[(f64, Constituent)] {
        &self.constituents
    }
    fn temperature_data(&self) -> &[(f64, f64, [f64; 9])] {
        &self.temperature_data
    }
    fn valid_temperature_range(&self) -> (f64, f64) {
        let first = self.temperature_data.first().unwrap();
        let last = self.temperature_data.last().unwrap();
        (first.0, last.1)
    }
    fn mw(&self) -> f64 {
        self.mw
    }
    fn h_formation(&self) -> f64 {
        self.h_formation
    }
    fn phase(&self) -> u8 {
        self.phase
    }
}

#[cfg(test)]
mod cea_validation {
    //! Ground truth pulled from NASA CEA2 (cea2.f, compiled from the
    //! `rocketcea` distribution) via `output debug 1`, which prints each
    //! species' H0j/RT and S0j/R directly from the thermo.inp fits.
    use crate::Species;
    use std::str::FromStr;

    // (species, T_k, h_over_rt, s_over_r)
    const CASES: &[(&str, f64, f64, f64)] = &[
        ("H2", 300.0, 0.0213920, 15.7387),
        ("H2", 1000.0, 2.48710, 19.9912),
        ("H2", 3000.0, 3.55726, 24.4017),
        ("O2", 300.0, 0.0217926, 24.6955),
        ("O2", 1000.0, 2.73103, 29.2967),
        ("O2", 3000.0, 3.93358, 34.2198),
        ("H2O", 300.0, -96.9245, 22.7358),
        ("H2O", 1000.0, -25.9574, 27.9916),
        ("H2O", 3000.0, -4.57705, 34.5172),
    ];

    #[test]
    fn h_and_s_over_r_match_cea() {
        for &(symbol, t, expected_h_over_rt, expected_s_over_r) in CASES {
            let species = Species::from_str(symbol).unwrap();
            let data = species.data();
            let h_over_rt = data.h_over_r(t) / t;
            let s_over_r = data.s_over_r(t);

            assert!(
                (h_over_rt - expected_h_over_rt).abs() < 1e-3,
                "{symbol} at {t}K: h_over_rt = {h_over_rt}, expected {expected_h_over_rt}"
            );
            assert!(
                (s_over_r - expected_s_over_r).abs() < 1e-3,
                "{symbol} at {t}K: s_over_r = {s_over_r}, expected {expected_s_over_r}"
            );
        }
    }
}
