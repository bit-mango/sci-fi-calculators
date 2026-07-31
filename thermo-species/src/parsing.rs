#![allow(dead_code)]
fn clean_identifier(raw_identifier: &str) -> String {
    let mut cleaned_identifier = raw_identifier.to_string();
    // Special naming for electrons.
    if cleaned_identifier == "e-" {
        cleaned_identifier = "E".to_string();
    }
    if cleaned_identifier.contains("+") {
        cleaned_identifier = cleaned_identifier.replace("+", "Plus");
    }
    if cleaned_identifier.contains("-") {
        cleaned_identifier = cleaned_identifier.replace("-", "Neg");
    }
    if cleaned_identifier.contains(",") {
        cleaned_identifier = cleaned_identifier.replace(",", "_");
    }
    if cleaned_identifier.contains(".") {
        cleaned_identifier = cleaned_identifier.replace(".", "_");
    }
    if cleaned_identifier.contains("'") {
        cleaned_identifier = cleaned_identifier.replace("'", "Prime");
    }
    if cleaned_identifier.contains("(") {
        cleaned_identifier = cleaned_identifier.replace("(", "_");
    }
    if cleaned_identifier.contains(")") {
        cleaned_identifier = cleaned_identifier.replace(")", "_");
    }
    // Remove trailing '_' if present.
    if cleaned_identifier.ends_with("_") {
        let mut identifier = cleaned_identifier.chars();
        identifier.next_back(); // Remove trailing '_'.
        cleaned_identifier = identifier.as_str().to_string();
    }
    // Remove leading '_' if present.
    if cleaned_identifier.starts_with("_") {
        let mut identifier = cleaned_identifier.chars();
        identifier.next(); // Remove leading '_'.
        cleaned_identifier = identifier.as_str().to_string();
    }
    cleaned_identifier
}

pub fn process_entry(
    iter: &mut impl Iterator<Item = String>,
) -> (
    (String, String),
    (usize, Vec<(f64, String)>, u8, f64, f64),
    Vec<(f64, f64, Vec<f64>)>,
) {
    // Get entry header.
    let header = iter.next().unwrap();
    let species_identifier = header.split(" ").collect::<Vec<&str>>()[0];
    let cleaned_species_identifier = clean_identifier(species_identifier);

    // Get entry metadata.
    let metadata = iter.next().unwrap();
    let cleaned_metadata = parse_metadata(&metadata);

    // Get temperature data.
    let cleaned_temperature_data = (0..cleaned_metadata.0)
        .map(|_| parse_temperature_data(iter))
        .collect::<Vec<(f64, f64, Vec<f64>)>>();

    (
        (cleaned_species_identifier, species_identifier.to_string()),
        cleaned_metadata,
        cleaned_temperature_data,
    )
}

fn parse_temperature_data(iter: &mut impl Iterator<Item = String>) -> (f64, f64, Vec<f64>) {
    let temperature_header = iter.next().unwrap();
    let temperature_low = &temperature_header[0..11]
        .trim()
        .parse::<f64>()
        .expect(&format!(
            "Failed to read lower temperature range: {}",
            &temperature_header[0..11]
        ));
    let temperature_high = &temperature_header[11..22]
        .trim()
        .parse::<f64>()
        .expect(&format!(
            "Failed to read upper temperature range: {}",
            &temperature_header[11..22]
        ));

    // First 5 coefficients.
    let mut coefficients = vec![0.0; 9];
    let coefficients_row = iter.next().unwrap();
    let mut i = 0;
    while i < 5 {
        let idx = i * 16;
        coefficients[i] = *&coefficients_row[idx..(idx + 16)]
            .trim()
            .replace("D", "e")
            .parse::<f64>()
            .expect(&format!(
                "Failed to parse coefficient[{}]: {}",
                i,
                &coefficients_row[0..16]
            ));

        i += 1;
    }
    // Last 4 coefficents.
    let offset = i; // Capture last i to properly index coefficients.
    i = 0; // Reset for indexing into next row.
    let coefficients_row = iter.next().unwrap();
    while i < 2 {
        let idx = i * 16;
        coefficients[i + offset] = *&coefficients_row[idx..(idx + 16)]
            .trim()
            .replace("D", "e")
            .parse::<f64>()
            .expect(&format!(
                "Failed to parse coefficient[{}]: {}",
                i,
                &coefficients_row[idx..(idx + 16)]
            ));

        i += 1;
    }
    while i < 4 {
        // Need to skip white space separating a and b coefficients.
        // So add 16.
        let idx = i * 16 + 16;
        coefficients[i + offset] = *&coefficients_row[idx..(idx + 16)]
            .trim()
            .replace("D", "e")
            .parse::<f64>()
            .expect(&format!(
                "Failed to parse coefficient[{}]: {}",
                i,
                &coefficients_row[idx..(idx + 16)]
            ));

        i += 1;
    }
    (*temperature_low, *temperature_high, coefficients)
}

fn parse_metadata(raw_metadata: &str) -> (usize, Vec<(f64, String)>, u8, f64, f64) {
    let temperature_intervals = &raw_metadata[0..2]
        .trim() // Remove white space
        .parse::<usize>()
        .expect(&format!(
            "Failed to parse temperature intervals: {}",
            &raw_metadata[0..2]
        ));

    let mut constituents = vec![];
    let raw_constituents = &raw_metadata[10..50];
    for i in 0..5 {
        let idx = i * 8;
        let constituent = &raw_constituents[idx..(idx + 8)];
        let name = constituent[0..2].trim().to_string();
        if name == "" {
            // Once we find this we have processed all the constiuents.
            break;
        }
        let moles = constituent[2..8].trim().parse::<f64>().expect(&format!(
            "Failed to parse moles from constituent: {}",
            constituent
        ));
        let element = if name.len() >= 2 {
            // Make sure the second letter is lower case.
            name.chars()
                .enumerate()
                .map(|(i, c)| if i == 1 { c.to_ascii_lowercase() } else { c })
                .collect::<String>()
        } else {
            name
        };
        constituents.push((moles, element));
    }

    let phase = *&raw_metadata[50..52]
        .trim()
        .parse::<i32>()
        .expect(&format!("Failed to parse phase: {}", &raw_metadata[50..52])) as u8;

    let mw = &raw_metadata[52..65].trim().parse::<f64>().expect(&format!(
        "Failed to parse molecular weight: {}",
        &raw_metadata[52..65]
    ));
    let std_heat_of_formation = &raw_metadata[65..80].trim().parse::<f64>().expect(&format!(
        "Failed to parse std heat of formation: {}",
        &raw_metadata[65..80]
    ));

    (
        *temperature_intervals,
        constituents,
        phase,
        *mw,
        *std_heat_of_formation,
    )
}

// TODO some identifier sanitizer function
// input is the raw string Ag+, H2O(L), etc and the output is the sanitized name
// AgPlus, H2O
//
// TODO Add All constant which returns a fixed array of species. Or could it return an array of enums?
// enums are way cheaper in memory compared to a full struct.

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_parse_metadata() {
        // Start with a simple single element species.
        let sample =
            " 3 g 6/97 H   1.00    0.00    0.00    0.00    0.00 0    1.0079400     217998.828";
        let (temperature_intervals, constituents, is_gas, mw, std_heat_of_formation) =
            parse_metadata(sample);
        assert_eq!(
            temperature_intervals, 3,
            "Temperature intervals should equal 3"
        );
        assert_eq!(constituents.len(), 1, "Should have 1 constituent");
        let expected_constituents = vec![(1.0, "H".to_string())];
        assert_eq!(constituents, expected_constituents, "Wrong constituents",);
        assert_eq!(is_gas, true, "is_gas should be true");
        assert_eq!(mw, 1.0079400, "mw should equal expected");
        assert_eq!(
            std_heat_of_formation, 217998.828,
            "std_heat_of_formation should equal expected"
        );

        // Try a multi element species.
        let sample =
            " 2 tpis91 C   1.00H   1.00CL  2.00BR  1.00    0.00 0  163.8286400     -45000.000";
        let (temperature_intervals, constituents, is_gas, mw, std_heat_of_formation) =
            parse_metadata(sample);
        assert_eq!(
            temperature_intervals, 2,
            "Temperature intervals should equal 2"
        );
        assert_eq!(constituents.len(), 4, "Should have 4 constituents");
        let expected_constituents = vec![
            (1.0, "C".to_string()),
            (1.0, "H".to_string()),
            (2.0, "CL".to_string()),
            (1.0, "BR".to_string()),
        ];
        assert_eq!(constituents, expected_constituents, "Wrong constituents",);
        assert_eq!(is_gas, true, "is_gas should be true");
        assert_eq!(mw, 163.8286400, "mw should equal expected");
        assert_eq!(
            std_heat_of_formation, -45000.000,
            "std_heat_of_formation should equal expected"
        );

        // Try a different phase.
        let sample =
            " 1 coda89 AG  1.00    0.00    0.00    0.00    0.00 1  107.8682000          0.000";
        let (temperature_intervals, constituents, is_gas, mw, std_heat_of_formation) =
            parse_metadata(sample);
        assert_eq!(
            temperature_intervals, 1,
            "Temperature intervals should equal 1"
        );
        assert_eq!(constituents.len(), 1, "Should have 1 constituent");
        let expected_constituents = vec![(1.0, "AG".to_string())];
        assert_eq!(constituents, expected_constituents, "Wrong constituents",);
        assert_eq!(is_gas, false, "is_gas should be false");
        assert_eq!(mw, 107.8682000, "mw should equal expected");
        assert_eq!(
            std_heat_of_formation, 0.000,
            "std_heat_of_formation should equal expected"
        );
    }

    #[test]
    fn test_parse_temperature_data() {
        // Try scenario where there is a gap between a and b coefficients.
        let sample = vec![
            "    298.150   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428"
                .to_string(),
            " 0.000000000D+00 0.000000000D+00 2.500000000D+00 0.000000000D+00 0.000000000D+00"
                .to_string(),
            " 0.000000000D+00 0.000000000D+00                -7.453750000D+02-1.172081224D+01"
                .to_string(),
        ];
        let (temperature_low, temperature_high, coefficients) =
            parse_temperature_data(&mut sample.into_iter());

        assert_eq!(temperature_low, 298.150, "Wrong temperature low");
        assert_eq!(temperature_high, 1000.000, "Wrong temperature high");
        let expected_coefficients = vec![
            0.000000000e+00,
            0.000000000e+00,
            2.500000000e+00,
            0.000000000e+00,
            0.000000000e+00,
            0.000000000e+00,
            0.000000000e+00,
            -7.453750000e+02,
            -1.172081224e+01,
        ];
        assert_eq!(coefficients, expected_coefficients, "Wrong coefficients");

        // Try scenario where there is a zero number between a and b coefficients.
        let sample = vec![
            "    300.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10027.417"
                .to_string(),
            " 3.218921730D+04-2.877601815D+02 4.203583820D+00 3.455405960D-03-6.746193340D-06"
                .to_string(),
            " 7.654571640D-09-2.870328419D-12 0.000000000D+00 4.733624710D+04-2.143628603D+00"
                .to_string(),
        ];
        let (temperature_low, temperature_high, coefficients) =
            parse_temperature_data(&mut sample.into_iter());

        assert_eq!(temperature_low, 300.000, "Wrong temperature low");
        assert_eq!(temperature_high, 1000.000, "Wrong temperature high");
        let expected_coefficients = vec![
            3.218921730e+04,
            -2.877601815e+02,
            4.203583820e+00,
            3.455405960e-03,
            -6.746193340e-06,
            7.654571640e-09,
            -2.870328419e-12,
            4.733624710e+04,
            -2.143628603e+00,
        ];
        assert_eq!(coefficients, expected_coefficients, "Wrong coefficients");
    }
}
