use crate::constants::*;
use std::f64;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub enum Species {
    H,
    H2,
    C,
    O,
    CO,
    N2,
    N,
    H2O,
    OH,
    CO2,
    O2,
    CH4,
    NH3,
    E,
    HPlus,
    CPlus,
    OPlus,
    NPlus,
    OHNeg,
}

impl Species {
    pub fn all() -> Vec<Self> {
        vec![
            Self::H,
            Self::H2,
            Self::C,
            Self::O,
            Self::CO,
            Self::N2,
            Self::N,
            Self::H2O,
            Self::OH,
            Self::CO2,
            Self::O2,
            Self::CH4,
            Self::NH3,
            Self::E,
            Self::HPlus,
            Self::CPlus,
            Self::OPlus,
            Self::NPlus,
            Self::OHNeg,
        ]
    }
    pub fn symbol(&self) -> String {
        match self {
            Species::H => "H".to_string(),
            Species::H2 => "H2".to_string(),
            Species::C => "C".to_string(),
            Species::O => "O".to_string(),
            Species::O2 => "O2".to_string(),
            Species::CO => "CO".to_string(),
            Species::CO2 => "CO2".to_string(),
            Species::N2 => "N2".to_string(),
            Species::N => "N".to_string(),
            Species::H2O => "H2O".to_string(),
            Species::OH => "OH".to_string(),
            Species::CH4 => "CH4".to_string(),
            Species::NH3 => "NH3".to_string(),
            Species::E => "e-".to_string(),
            Species::HPlus => "H+".to_string(),
            Species::CPlus => "C+".to_string(),
            Species::OPlus => "O+".to_string(),
            Species::NPlus => "N+".to_string(),
            Species::OHNeg => "OH-".to_string(),
        }
    }

    pub fn mw(&self) -> f64 {
        match self {
            Species::H => H_MW,
            Species::H2 => H2_MW,
            Species::C => C_MW,
            Species::O => O_MW,
            Species::O2 => O2_MW,
            Species::CO => CO_MW,
            Species::CO2 => CO2_MW,
            Species::N2 => N2_MW,
            Species::N => N_MW,
            Species::H2O => H2O_MW,
            Species::OH => OH_MW,
            Species::CH4 => CH4_MW,
            Species::NH3 => NH3_MW,
            Species::E => 0.0,
            Species::HPlus => H_MW,
            Species::CPlus => C_MW,
            Species::OPlus => O_MW,
            Species::NPlus => N_MW,
            Species::OHNeg => OH_MW,
        }
    }

    pub fn constituents(&self) -> Vec<(f64, Self)> {
        match self {
            Species::H => vec![(1.0, Species::H)],
            Species::H2 => vec![(2.0, Species::H)],
            Species::C => vec![(1.0, Species::C)],
            Species::O => vec![(1.0, Species::O)],
            Species::O2 => vec![(2.0, Species::O)],
            Species::CO => vec![(1.0, Species::C), (1.0, Species::O)],
            Species::CO2 => vec![(1.0, Species::C), (2.0, Species::O)],
            Species::N2 => vec![(2.0, Species::N)],
            Species::N => vec![(1.0, Species::N)],
            Species::H2O => vec![(2.0, Species::H), (1.0, Species::O)],
            Species::OH => vec![(1.0, Species::O), (1.0, Species::H)],
            Species::CH4 => vec![(1.0, Species::C), (4.0, Species::H)],
            Species::NH3 => vec![(1.0, Species::N), (3.0, Species::H)],
            Species::E => vec![(1.0, Species::E)],
            Species::HPlus => vec![(1.0, Species::H), (-1.0, Species::E)],
            Species::CPlus => vec![(1.0, Species::C), (-1.0, Species::E)],
            Species::OPlus => vec![(1.0, Species::O), (-1.0, Species::E)],
            Species::NPlus => vec![(1.0, Species::N), (-1.0, Species::E)],
            Species::OHNeg => vec![(1.0, Species::O), (1.0, Species::H), (1.0, Species::E)],
        }
    }
    pub fn is_charged(&self) -> bool {
        matches!(
            self,
            Species::E
                | Species::HPlus
                | Species::CPlus
                | Species::OPlus
                | Species::NPlus
                | Species::OHNeg
        )
    }
}
