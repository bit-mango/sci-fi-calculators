pub const R: f64 = 8.314; // Ideal Gas Constant [J/mol•K]
pub const G_0: f64 = 9.81;

pub const H_MW: f64 = 1.0e-3;
pub const H2_MW: f64 = 2.0 * H_MW;
pub const C_MW: f64 = 12.0e-3;
pub const O_MW: f64 = 16.0e-3;
pub const CO_MW: f64 = C_MW + O_MW;

pub const STD_REFERENCE_PRESSURE: f64 = 1.0; // [bar]
pub const ENTHALPY_HYDROGEN: f64 = 435.998e3; // ∆H [J/mol]
pub const ENTROPY_HYDROGEN: f64 = 98.753; // ∆S [J/K•mol]

pub const ENTHALPY_CARBON_MONOXIDE: f64 = 1076.375e3; // ∆H [J/mol]
pub const ENTROPY_CARBON_MONOXIDE: f64 = 121.498; // ∆S [J/K•mol]
