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
