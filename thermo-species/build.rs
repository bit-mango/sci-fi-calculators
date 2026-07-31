use std::collections::HashMap;
use std::fmt::Write;
use std::{env, fs, path::Path};

#[path = "src/parsing.rs"]
mod parsing;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("cargo::rerun-if-changed=build.rs");
    println!("cargo::rerun-if-changed=data/thermo.inp");
    println!("cargo::rerun-if-changed=src/parsing.rs");
    println!("cargo::rerun-if-changed=src/types.rs");
    let custom_thermo_inp = env::var_os("THERMO_INP_PATH");
    let thermo_inp_contents = if let Some(custom) = custom_thermo_inp {
        println!("cargo::rerun-if-changed={:?}", custom.to_string_lossy());
        fs::read_to_string(&custom).unwrap_or_else(|_| panic!("Failed to read file: {:?}", custom))
    } else {
        let default = "data/thermo.inp";
        fs::read_to_string(default).unwrap_or_else(|_| panic!("Failed to read file: {:?}", default))
    };
    let thermo_inp_cleaned: Vec<String> = thermo_inp_contents
        .split("\n") // Convert file into vec of full lines
        .filter(|row| !row.starts_with("!")) // Remove comments
        .map(|line| line.to_string())
        .collect();

    let mut iter = thermo_inp_cleaned.into_iter().peekable();

    loop {
        if iter
            .peek()
            .ok_or("reached end of input before finding first species (\"e-\")")?
            .split(" ")
            .collect::<Vec<&str>>()[0]
            == "e-"
        {
            break;
        } else {
            iter.next(); // Throw out every line until we hit the first specimen.
        }
    }

    let mut buffer = String::new();
    writeln!(
        &mut buffer,
        "pub use crate::types::{{Element, SpeciesData, AnySpeciesData, ParseSpeciesError}};"
    )?;
    let mut species: HashMap<String, (String, parsing::Metadata, parsing::TemperatureData)> =
        HashMap::new();
    while iter
        .peek()
        .ok_or("reached end of input before \"END PRODUCTS\"")?
        != "END PRODUCTS"
    {
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
            if entry.0 != dirty_species_identifier
                || entry.1.3 != mw
                || entry.1.4 != h_formation
                || entry.1.1 != constituents
            {
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
    for (
        cleaned_species_identifier,
        (
            dirty_species_identifier,
            (_, constituents, phase, mw, h_formation),
            cleaned_temperature_data,
        ),
    ) in consolidated_species.iter()
    {
        // Update condensed phase.
        let key = if *phase == 0 {
            // Key is just cleaned identifier.
            cleaned_species_identifier.clone()
        } else {
            let idx_phase_start = cleaned_species_identifier.rfind("__").unwrap_or_else(|| {
                panic!("Failed to find phase in: {}", cleaned_species_identifier)
            });
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
        )?;
        // Add symbol.
        write!(&mut buffer, "symbol: \"{}\",", dirty_species_identifier)?;
        // Add constituents.
        let mut constituent_str = "".to_string();
        for (moles, element) in constituents {
            constituent_str += &format!("({:?}, Element::{}),", moles, element);
        }
        write!(&mut buffer, "constituents: [{}],", constituent_str)?;
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
        write!(&mut buffer, "temperature_data: [{}],", temperature_data_str)?;
        // Add mw.
        write!(&mut buffer, "mw: {:?},", mw)?;
        // Add h_formation.
        write!(&mut buffer, "h_formation: {:?},", h_formation)?;
        // Add phase.
        write!(&mut buffer, "phase: {:?},", phase)?;
        // Close out struct definition.
        writeln!(&mut buffer, "}};")?;
    }

    // Add Species enum.
    write!(
        &mut buffer,
        "#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]"
    )?;
    write!(&mut buffer, "pub enum Species {{")?;
    for species in all.iter() {
        write!(&mut buffer, "{},", species)?;
    }
    writeln!(&mut buffer, "}}")?;

    write!(&mut buffer, "pub static ALL: [Species; {}] = [", all.len())?;
    for species in all.iter() {
        write!(&mut buffer, "Species::{},", species)?;
    }
    writeln!(&mut buffer, "];")?;

    // Start Species impl
    write!(&mut buffer, "impl Species {{")?;

    // Add all
    writeln!(&mut buffer, "pub fn all() -> &'static [Species] {{&ALL}}")?;

    // Add phases
    write!(&mut buffer, "pub fn phases(&self) -> &'static [Species] {{")?;
    write!(&mut buffer, "match self {{")?;
    for (base, phases) in condensed_phases.iter() {
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
        write!(&mut buffer, "{} => {}", left, right)?;
    }
    // Close match arm
    write!(&mut buffer, "}}")?;
    // Close fn definition
    writeln!(&mut buffer, "}}")?;

    // Add data
    write!(
        &mut buffer,
        "pub fn data(&self) -> &'static dyn AnySpeciesData {{"
    )?;
    write!(&mut buffer, "match self {{")?;
    for species in all.iter() {
        write!(&mut buffer, "Species::{} => &{},", species, species)?;
    }
    // Close match arm
    write!(&mut buffer, "}}")?;
    // Close fn definition
    write!(&mut buffer, "}}")?;

    // Close the implementation.
    writeln!(&mut buffer, "}}")?;

    // Implement from_str for Species.
    writeln!(&mut buffer, "impl std::str::FromStr for Species {{")?;
    writeln!(&mut buffer, "type Err = ParseSpeciesError;")?;
    writeln!(
        &mut buffer,
        "fn from_str(s: &str) -> Result<Self, Self::Err> {{"
    )?;
    write!(&mut buffer, "match s {{")?;
    for (cleaned_species_identifier, (dirty_species_identifier, (_, _, _, _, _), _)) in
        consolidated_species.iter()
    {
        write!(
            &mut buffer,
            "\"{}\" => Ok(Species::{}),",
            dirty_species_identifier, cleaned_species_identifier
        )?;
    }
    write!(&mut buffer, "_ => Err(ParseSpeciesError), }}",)?;

    // Close functin definition.
    writeln!(&mut buffer, "}}")?;

    // Close the implementation
    writeln!(&mut buffer, "}}")?;

    // Define species!() macro.
    writeln!(&mut buffer, "#[macro_export]")?;
    writeln!(&mut buffer, "macro_rules! species {{")?;
    for (cleaned_species_identifier, (dirty_species_identifier, (_, _, _, _, _), _)) in
        consolidated_species.iter()
    {
        write!(
            &mut buffer,
            "(\"{}\") => {{ $crate::Species::{}}};",
            dirty_species_identifier, cleaned_species_identifier
        )?;
    }
    write!(
        &mut buffer,
        "($other:literal) => {{ compile_error!(concat!(\"Unknown species: \", $other))}}",
    )?;
    // Close the macro
    writeln!(&mut buffer, "}}")?;

    let out_dir =
        env::var_os("OUT_DIR").ok_or("OUT_DIR not set — build scripts must be run via cargo")?;
    let dest_path = Path::new(&out_dir).join("thermo.rs");
    fs::write(&dest_path, buffer)?;

    Ok(())
}
