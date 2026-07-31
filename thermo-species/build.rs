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
    // writeln!(&mut buffer, "#[allow(non_upper_case)]").unwrap();
    writeln!(
        &mut buffer,
        "pub use crate::types::{{Element, SpeciesData, AnySpeciesData}};"
    )
    .unwrap();
    // Key is the identifier without '_phase' at end.
    // Value is vec of identifiers.
    // Example:
    // Key: H2O
    // Value: vec!["H2O_cr", "H2O_L",...]
    // First capture all products/reactants.

    let mut species: HashMap<
        String,
        (
            String,
            (usize, Vec<(f64, String)>, u8, f64, f64),
            Vec<(f64, f64, Vec<f64>)>,
        ),
    > = HashMap::new();
    // TODO store all values in a hashmap, concating temperature data
    // together if exact same species name, and other metadata matches
    // Also maybe add some loud warning about when an entries temperature data
    // is not in ascending order.
    while iter.peek().unwrap() != "END PRODUCTS" {
        let (
            (cleaned_species_identifier, dirty_species_identifier),
            (t_intervals, constituents, phase, mw, h_formation),
            cleaned_temperature_data,
        ) = parsing::process_entry(&mut iter);

        if cleaned_species_identifier.starts_with("Inert") {
            continue;
        }

        // Consolidate species.
        if let Some(entry) = species.get_mut(&cleaned_species_identifier) {
            // Exact species already exists, check they are the same, and append temperature data together.
            if entry.0 != dirty_species_identifier || entry.1.3 != mw || entry.1.4 != h_formation {
                panic!("Species have the exact same name, but different data!");
            }
            let current_temperature_data = entry.2.clone();
            let new_temperature_data = cleaned_temperature_data;
            let mut c_iter = current_temperature_data.iter().peekable();
            let mut n_iter = new_temperature_data.iter().peekable();
            let mut combined_temperature_data = vec![];
            while let (Some(c_peek), Some(n_peek)) = (c_iter.peek(), n_iter.peek()) {
                // Compare them, and order them.
                if c_peek.0 < n_peek.0 {
                    let c_high = c_peek.1;
                    // Add current temperature data.
                    combined_temperature_data.push(c_iter.next().unwrap().clone());
                    if c_high == n_peek.0 {
                        // Add new temperature data.
                        combined_temperature_data.push(n_iter.next().unwrap().clone());
                    }
                } else if n_peek.0 < c_peek.0 {
                    let n_high = n_peek.1;
                    // Add new temperature data.
                    combined_temperature_data.push(n_iter.next().unwrap().clone());
                    if n_high == c_peek.0 {
                        // Add current temperature data.
                        combined_temperature_data.push(c_iter.next().unwrap().clone());
                    }
                } else {
                    // We have conflicting temperature ranges, ie two different
                    // entries start with the same lower temperature.
                    panic!(
                        "Conflicting temperature ranges for {}",
                        dirty_species_identifier
                    );
                }
            }
            // Append remaining data into combined.
            c_iter.for_each(|c| combined_temperature_data.push(c.clone()));
            n_iter.for_each(|n| combined_temperature_data.push(n.clone()));

            entry.2 = combined_temperature_data;
        } else {
            species.insert(
                cleaned_species_identifier,
                (
                    dirty_species_identifier,
                    (t_intervals, constituents, phase, mw, h_formation),
                    cleaned_temperature_data,
                ),
            );
        }
    }

    let mut all = vec![];
    let mut condensed_phases: HashMap<String, Vec<String>> = HashMap::new();

    let mut consolidated_species = species.drain().collect::<Vec<(
        String,
        (
            String,
            (usize, Vec<(f64, String)>, u8, f64, f64),
            Vec<(f64, f64, Vec<f64>)>,
        ),
    )>>();
    consolidated_species.sort_by(|a, b| a.0.cmp(&b.0));
    consolidated_species.iter().for_each(
        |(
            cleaned_species_identifier,
            (
                dirty_species_identifier,
                (_, constituents, phase, mw, h_formation),
                cleaned_temperature_data,
            ),
        )| {
            // Update condensed phase.
            let key = if *phase == 0 {
                // Key is just cleaned identifier.
                cleaned_species_identifier.clone()
            } else {
                let idx_phase_start = cleaned_species_identifier.rfind("_").unwrap();
                let mut identifier_no_phase = cleaned_species_identifier.clone();
                let _ = identifier_no_phase.split_off(idx_phase_start);
                identifier_no_phase
            };
            // Update condensed_phase
            if let Some(entry) = condensed_phases.get_mut(&key) {
                if *phase != 0 {
                    entry.push(cleaned_species_identifier.clone());
                }
            } else {
                let mut phases = vec![];
                if *phase != 0 {
                    phases.push(cleaned_species_identifier.clone());
                }
                condensed_phases.insert(key, phases);
            }

            // Add each species to the all vec.
            all.push(cleaned_species_identifier.clone());

            let c = constituents.len();
            let t = cleaned_temperature_data.len();
            // Add beginning of static species.
            write!(
                &mut buffer,
                "pub static {}: SpeciesData<{}, {}> = SpeciesData {{",
                cleaned_species_identifier, c, t
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
        },
    );

    // Add Species enum.
    write!(&mut buffer, "pub enum Species {{").unwrap();
    for species in all.iter() {
        write!(&mut buffer, "{},", species).unwrap();
    }
    writeln!(&mut buffer, "}}").unwrap();

    write!(&mut buffer, "impl Species {{").unwrap();
    // Add phases
    write!(&mut buffer, "pub fn phases(&self) -> &'static [Species] {{").unwrap();
    write!(&mut buffer, "match self {{").unwrap();
    condensed_phases.iter().for_each(|(base, phases)| {
        let base_is_real = all.contains(base);
        let mut left = if base_is_real {
            format!("Species::{}", base)
        } else {
            String::new()
        };
        let mut right = "&[".to_string();
        phases.iter().for_each(|phase| {
            left += &format!(" | Species::{}", phase);
            right += &format!("Species::{},", phase);
        });
        right += "],";
        write!(&mut buffer, "{} => {}", left, right).unwrap();
    });
    // Close match arm
    write!(&mut buffer, "}}").unwrap();
    // Close fn definition
    writeln!(&mut buffer, "}}").unwrap();

    // Add data
    write!(
        &mut buffer,
        "pub fn data(&self) -> &'static dyn AnySpeciesData {{"
    )
    .unwrap();
    write!(&mut buffer, "match self {{").unwrap();
    all.iter().for_each(|species| {
        write!(&mut buffer, "Species::{} => &{},", species, species).unwrap();
    });
    // Close match arm
    write!(&mut buffer, "}}").unwrap();
    // Close fn definition
    write!(&mut buffer, "}}").unwrap();

    // Close the implementation.
    write!(&mut buffer, "}}").unwrap();

    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("thermo.rs");
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
