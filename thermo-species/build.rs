use std::collections::HashMap;
use std::fmt::Write;
use std::{env, fs, path::Path};

#[path = "src/parsing.rs"]
mod parsing;
#[path = "src/types.rs"]
mod types;

fn main() {
    println!("cargo::rerun-if-changed=build.rs");
    println!("cargo::rerun-if-changed=data/thermo.inp");
    println!("cargo::rerun-if-changed=src/parsing.rs");
    println!("cargo::rerun-if-changed=src/types.rs");
    let custom_thermo_inp = env::var_os("THERMO_INP_PATH");
    let thermo_inp_contents = if let Some(custom) = custom_thermo_inp {
        println!("cargo::rerun-if-changed={:?}", &custom);
        let contents =
            fs::read_to_string(&custom).expect(&format!("Failed to read file: {:?}", custom));
        contents
    } else {
        let default = "data/thermo.inp";
        let contents =
            fs::read_to_string(&default).expect(&format!("Failed to read file: {:?}", default));
        contents
    };
    let thermo_inp_cleaned: Vec<String> = thermo_inp_contents
        .split("\n") // Convert file into vec of full lines
        .filter(|row| !row.starts_with("!")) // Remove comments
        .map(|line| line.to_string())
        .collect();

    let mut iter = thermo_inp_cleaned.into_iter().peekable();

    loop {
        if iter.peek().unwrap().split(" ").collect::<Vec<&str>>()[0] == "e-" {
            break;
        } else {
            iter.next(); // Throw out every line until we hit the first specimen.
        }
    }
    // types::SpeciesData

    let mut buffer = String::new();
    // writeln!(&mut buffer, "#![allow(non_snake_case)]").unwrap();
    writeln!(&mut buffer, "use crate::types::{{Element, SpeciesData}};").unwrap();
    let mut condensed_phases: HashMap<String, Vec<String>> = HashMap::new();
    // Key is the identifier without '_phase' at end.
    // Value is vec of identifiers.
    // Example:
    // Key: H2O
    // Value: vec!["H2O_cr", "H2O_L",...]
    // First capture all products/reactants.

    // TODO store all values in a hashmap, concating temperature data
    // together if exact same species name, and other metadata matches
    // Also maybe add some loud warning about when an entries temperature data
    // is not in ascending order.
    // TODO rather than different constants for each phase, what if I had 1 constant
    // where phase 0 was gas, then the next entry is the next phase and so on?
    while iter.peek().unwrap() != "END PRODUCTS" {
        let (
            (cleaned_species_identifier, dirty_species_identifier),
            (t_intervals, constituents, phase, mw, h_formation),
            cleaned_temperature_data,
        ) = parsing::process_entry(&mut iter);

        // Update condensed phase.
        let key = if phase == 0 {
            // Key is just cleaned identifier.
            cleaned_species_identifier.clone()
        } else {
            let idx_phase_start = cleaned_species_identifier.rfind("_").unwrap();
            let identifier_no_phase = cleaned_species_identifier
                .clone()
                .split_off(idx_phase_start);
            identifier_no_phase
        };
        if let Some(entry) = condensed_phases.get_mut(&key) {
            entry.push(cleaned_species_identifier.clone());
        } else {
            condensed_phases.insert(key, vec![cleaned_species_identifier.clone()]);
        }

        let c = constituents.len();
        // Add beginning of static species.
        write!(
            &mut buffer,
            "pub static {}: SpeciesData<{}, {}> = SpeciesData {{",
            cleaned_species_identifier, c, t_intervals
        )
        .unwrap();
        // Add symbol.
        write!(&mut buffer, "symbol: \"{}\",", dirty_species_identifier).unwrap();
        // Add constituents.
        let mut constituent_str = "".to_string();
        for (moles, element) in constituents {
            constituent_str += &format!("({:?}, Element::{}),", moles, element);
        }
        write!(&mut buffer, "constituents: [{}],", constituent_str).unwrap();
        // Add temperature data.
        let mut temperature_data_str = "".to_string();
        for (low_temperature, high_temperature, coefficients) in cleaned_temperature_data {
            let mut coefficients_str = "".to_string();
            for coefficient in coefficients {
                coefficients_str += &format!("{:?},", coefficient);
            }
            temperature_data_str += &format!(
                "({:?}, {:?}, [{}],),",
                low_temperature, high_temperature, coefficients_str
            );
        }
        write!(&mut buffer, "temperature_data: [{}],", temperature_data_str).unwrap();
        // Add mw.
        write!(&mut buffer, "mw: {:?},", mw).unwrap();
        // Add h_formation.
        write!(&mut buffer, "h_formation: {:?},", h_formation).unwrap();
        // Add phase.
        write!(&mut buffer, "phase: {:?},", phase).unwrap();
        // Close out struct definition.
        writeln!(&mut buffer, "}};").unwrap();
    }

    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("generated_sample.rs");
    fs::write(&dest_path, buffer).unwrap();

    // use write! and writeln! macros for building string.
    // format parsed floats using {:?} to keep precision.
}

// TODO need to setup grouping, where we enforce that when I strip the _phase at the end of the
// identifier, all matching species get grouped together.

// TODO build Species enum.
// TODO All const array of Species enum. Note all speces data is in the binary, but only copied once!
// impl Species {
// pub fn data(&self) -> &'static SpeciesData {
// match self => statics NOT a copy
// }
// }
