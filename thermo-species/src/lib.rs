mod parsing;
mod types;
#[allow(non_upper_case_globals, non_camel_case_types)]
mod generated {
    include!(concat!(env!("OUT_DIR"), "/thermo.rs"));
}
pub use generated::{AnySpeciesData, Element, Species, SpeciesData};
