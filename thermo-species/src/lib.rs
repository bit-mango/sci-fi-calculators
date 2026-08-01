mod parsing;
mod types;
#[allow(non_upper_case_globals, non_camel_case_types)]
mod generated {
    include!(concat!(env!("OUT_DIR"), "/thermo.rs"));
}
pub use generated::{AnySpeciesData, Constituent, Species, SpeciesData};

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_species_macro() {
        // Compile time Species
        assert_eq!(species!("H2O(cr)"), Species::H2O__1);
    }
    #[test]
    fn test_from_str() {
        // Run time dynamic Species
        assert_eq!(Species::from_str("H2O").unwrap(), Species::H2O);

        println!("{}", species!("CH4"));
    }
}
