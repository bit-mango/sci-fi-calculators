pub const R: f64 = 8.314; // Ideal Gas Constant [J/mol•K]
pub const G_0: f64 = 9.80665;
pub const F: f64 = 96_485.332; // Faraday constant, C/mol
pub const EPSILON_0: f64 = 8.854e-12; // F/m

pub const O2_DISSOCIATION: f64 = 249.2e3; // J/mol O, so double for O2
pub const O_ELECTRON_AFFINITY: f64 = 141.0e3; // O + e- => O-
pub const O_FIRST_IONIZATION: f64 = 1_313.9e3; // O =? O+ + e-
pub const CH4_FIRST_IONIZATION: f64 = 1_216.0e3; // CH4 =? CH4+ + e-
pub const H2O_ELECTROLYSIS: f64 = 285.8e3; // J/mol H2O => H2 + 1/2 O2
pub const CH4_PARTIAL_OXIDATION_ENERGY: f64 = 284.7e3; // CH4 + O => CO + 2H2

pub const H_MW: f64 = 1.0e-3;
pub const H2_MW: f64 = 2.0 * H_MW;
pub const C_MW: f64 = 12.0e-3;
pub const O_MW: f64 = 16.0e-3;
pub const CO_MW: f64 = C_MW + O_MW;
pub const N_MW: f64 = 14.0e-3;
pub const N2_MW: f64 = 2.0 * N_MW;
pub const CH4_MW: f64 = 16.0e-3;
pub const H2O_MW: f64 = 18.0e-3;
pub const UREA_MW: f64 = 60.0e-3;
pub const NH3_MW: f64 = 17.0e-3;
pub const UREA_LATENT_HEAT_OF_FUSION: f64 = 14.5e3;

pub const STD_REFERENCE_PRESSURE: f64 = 1.0; // [bar]
