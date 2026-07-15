use crate::constants::R;
use std::{fmt, fs};

pub struct PolynomialCoefficients {
    pub a_1: f64,
    pub a_2: f64,
    pub a_3: f64,
    pub a_4: f64,
    pub a_5: f64,
    pub a_6: f64,
    pub a_7: f64,
    pub b_1: f64,
    pub b_2: f64,
    pub min_temperature: f64,
    pub max_temperature: f64,
}

impl PolynomialCoefficients {
    pub fn new(
        a_1: f64,
        a_2: f64,
        a_3: f64,
        a_4: f64,
        a_5: f64,
        a_6: f64,
        a_7: f64,
        b_1: f64,
        b_2: f64,
        min_temperature: f64,
        max_temperature: f64,
    ) -> Self {
        Self {
            a_1,
            a_2,
            a_3,
            a_4,
            a_5,
            a_6,
            a_7,
            b_1,
            b_2,
            min_temperature,
            max_temperature,
        }
    }
}

pub struct TemperatureDependentProperty {
    coefficients: Vec<PolynomialCoefficients>,
}

impl fmt::Display for TemperatureDependentProperty {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Use the write! macro to send formatted text to the formatter buffer
        let mut buffer = String::new();
        for coefficients in &self.coefficients {
            buffer += &format!(
                "
                Temperature Range(K): {} <=> {}
                Coefficients:
                    a_1: {}
                    a_2: {}
                    a_3: {}
                    a_4: {}
                    a_5: {}
                    a_6: {}
                    a_7: {}
                    b_1: {}
                    b_2: {}
            ",
                coefficients.min_temperature,
                coefficients.max_temperature,
                coefficients.a_1,
                coefficients.a_2,
                coefficients.a_3,
                coefficients.a_4,
                coefficients.a_5,
                coefficients.a_6,
                coefficients.a_7,
                coefficients.b_1,
                coefficients.b_2,
            );
        }
        write!(f, "{}", buffer)
    }
}

impl TemperatureDependentProperty {
    pub fn new() -> Self {
        Self {
            coefficients: vec![],
        }
    }
    // TODO could add more checks to make sure we dont have overlaps, and to order them.
    pub fn add_coefficients(&mut self, pc: PolynomialCoefficients) {
        self.coefficients.push(pc);
    }

    fn find_coefficients(&self, temperature_k: f64) -> Option<&PolynomialCoefficients> {
        self.coefficients.iter().find(|&coefficients| {
            coefficients.min_temperature <= temperature_k
                && coefficients.max_temperature >= temperature_k
        })
    }
    pub fn cp(&self, temperature_k: f64) -> f64 {
        let c = self.find_coefficients(temperature_k).expect(&format!(
            "Could not find coefficients for temperature_k {temperature_k}"
        ));
        let t_1 = temperature_k;
        let t_2 = t_1 * t_1;
        let t_3 = t_2 * t_1;
        let t_4 = t_2 * t_2;

        let res = c.a_1 / t_2
            + c.a_2 / t_1
            + c.a_3
            + c.a_4 * t_1
            + c.a_5 * t_2
            + c.a_6 * t_3
            + c.a_7 * t_4;

        res * R
    }

    pub fn h(&self, temperature_k: f64) -> f64 {
        let c = self.find_coefficients(temperature_k).expect(&format!(
            "Could not find coefficients for temperature_k {temperature_k}"
        ));
        let t_1 = temperature_k;
        let t_2 = t_1 * t_1;
        let t_3 = t_2 * t_1;
        let t_4 = t_2 * t_2;

        let res = -c.a_1 / t_2
            + c.a_2 * t_1.ln() / t_1
            + c.a_3
            + c.a_4 * t_1 / 2.0
            + c.a_5 * t_2 / 3.0
            + c.a_6 * t_3 / 4.0
            + c.a_7 * t_4 / 5.0
            + c.b_1 / t_1;

        res * R * t_1
    }

    pub fn s(&self, temperature_k: f64) -> f64 {
        let c = self.find_coefficients(temperature_k).expect(&format!(
            "Could not find coefficients for temperature_k {temperature_k}"
        ));
        let t_1 = temperature_k;
        let t_2 = t_1 * t_1;
        let t_3 = t_2 * t_1;
        let t_4 = t_2 * t_2;

        let res = -c.a_1 / (2.0 * t_2) - c.a_2 / t_1
            + c.a_3 * t_1.ln()
            + c.a_4 * t_1
            + c.a_5 * t_2 / 2.0
            + c.a_6 * t_3 / 3.0
            + c.a_7 * t_4 / 4.0
            + c.b_2;

        res * R
    }
}

use std::collections::HashMap;
pub struct ThermoReference {
    reference: HashMap<String, TemperatureDependentProperty>,
}

impl ThermoReference {
    pub fn new() -> Self {
        let mut path_to_file =
            std::env::current_dir().expect(&format!("Failed to read current dir."));
        path_to_file.push("thermo.inp");
        let path_to_file = path_to_file.to_str().unwrap();
        let reference = Self::create_reference_from_file(path_to_file)
            .expect(&format!("Failed to read file: {}", path_to_file));

        Self { reference }
    }

    pub fn get_tdp(&self, key: &str) -> &TemperatureDependentProperty {
        self.reference
            .get(key)
            .expect(&format!("Unsupported Species: {}", key))
    }

    fn create_reference_from_file(
        path_to_file: &str,
    ) -> std::io::Result<HashMap<String, TemperatureDependentProperty>> {
        let contents = fs::read_to_string(path_to_file)?;

        let mut content: Vec<&str> = contents
            .split("\n")
            .filter(|row| !row.starts_with("!"))
            .collect();

        let no_header_content = content.split_off(2);

        let mut content_idx = 0;
        let mut reference = HashMap::new();
        while content_idx < no_header_content.len() {
            // content_idx starts by pointing at a new species we need to add to the hashmap.
            let key: &str = no_header_content[content_idx]
                .split(" ")
                .collect::<Vec<&str>>()[0];

            if key == "END" {
                break;
            }
            let number_of_temperature_intervals = no_header_content[content_idx + 1]
                .split(" ")
                .collect::<Vec<&str>>()[1];

            let number_of_temperature_intervals = number_of_temperature_intervals
                .parse::<usize>()
                .expect(&format!(
                    "Failed to parse temperature intervals from {}",
                    no_header_content[content_idx + 1]
                ));

            // Entries can be repeated in thermo.inp, so check if it exists,
            // if not create it.
            if !reference.contains_key(key) {
                let new_tdp = TemperatureDependentProperty::new();
                reference.insert(key.to_string(), new_tdp);
            }
            let tdp = reference.get_mut(key).unwrap();
            for _ in 0..number_of_temperature_intervals {
                // Grab temperature range first.
                let temperature_row: Vec<&str> = no_header_content[content_idx + 2]
                    .split(" ")
                    .filter(|x| !x.is_empty())
                    .collect();
                let min_temperature = temperature_row[0]
                    .parse::<f64>()
                    .expect(&format!("Failed to parse: {}", temperature_row[0]));
                let max_temperature = temperature_row[1]
                    .parse::<f64>()
                    .expect(&format!("Failed to parse: {}", temperature_row[1]));

                // Grab coefficients.
                let mut coefficient_row: Vec<f64> = [
                    no_header_content[content_idx + 3],
                    no_header_content[content_idx + 4],
                ]
                .concat()
                .split(" ")
                .filter(|x| !x.is_empty())
                .map(|coeff| {
                    let mut coeff = coeff.replace("D", "e");
                    let res = if coeff.len() > 17 {
                        let mut coeffs = vec![];
                        while coeff.len() > 17 {
                            // Need to split malformed strings. Happens when consecutive number(s) are negative.
                            // Example: "1.000000000D+00-1.000000000D+00-1.000000000D+00"
                            let e_idx = coeff.find("e").unwrap();
                            let coefficients_string = coeff.split_at(e_idx + 4);

                            let first_coefficient = coefficients_string
                                .0
                                .parse::<f64>()
                                .expect(&format!("Failed to parse: {}", coefficients_string.0));
                            coeffs.push(first_coefficient);
                            coeff = coefficients_string.1.to_string();
                        }
                        // Parse final coefficient.
                        coeffs.push(
                            coeff
                                .parse::<f64>()
                                .expect(&format!("Failed to parse: {}", coeff)),
                        );
                        coeffs
                    } else {
                        // Properly formatted number
                        vec![
                            coeff
                                .parse::<f64>()
                                .expect(&format!("Failed to parse: {}", coeff)),
                        ]
                    };
                    res
                })
                .flatten()
                .collect();

                if coefficient_row.len() != 9 {
                    // Some entries have 10 coefficients and I have no idea why, but the 7th idx
                    // appears to always be zero, so if it is remove it.
                    if coefficient_row.len() == 10 && coefficient_row[7] == 0.0 {
                        let mut row_end = coefficient_row.split_off(8);
                        coefficient_row.pop(); // Removes 7th element.
                        coefficient_row.append(&mut row_end);
                    } else {
                        // Wrong number of coefficients, print row for debugging.
                        panic!("coefficients_row wrong length: {:?}", coefficient_row);
                    }
                }

                let pc = PolynomialCoefficients::new(
                    coefficient_row[0],
                    coefficient_row[1],
                    coefficient_row[2],
                    coefficient_row[3],
                    coefficient_row[4],
                    coefficient_row[5],
                    coefficient_row[6],
                    coefficient_row[7],
                    coefficient_row[8],
                    min_temperature,
                    max_temperature,
                );

                tdp.add_coefficients(pc);
                content_idx += 3;
            }
            content_idx += 2;
        }

        Ok(reference)
    }
}
