use crate::thermo::fluid_properties::TemperatureDependentProperty;

pub fn get_rxn_enthalpy(
    rxn_temperature: f64,
    rxn_reactants: &Vec<&TemperatureDependentProperty>,
    rxn_products: &Vec<&TemperatureDependentProperty>,
) -> f64 {
    let enthalpy_products: f64 = rxn_products
        .iter()
        .map(|product_tdp| product_tdp.h(rxn_temperature))
        .sum();
    let enthalpy_reactants: f64 = rxn_reactants
        .iter()
        .map(|reactant_tdp| reactant_tdp.h(rxn_temperature))
        .sum();

    // Reaction enthalpy
    enthalpy_products - enthalpy_reactants
}

pub fn get_rxn_entropy(
    rxn_temperature: f64,
    rxn_reactants: &Vec<&TemperatureDependentProperty>,
    rxn_products: &Vec<&TemperatureDependentProperty>,
) -> f64 {
    let entropy_products: f64 = rxn_products
        .iter()
        .map(|product_tdp| product_tdp.s(rxn_temperature))
        .sum();
    let entropy_reactants: f64 = rxn_reactants
        .iter()
        .map(|reactant_tdp| reactant_tdp.s(rxn_temperature))
        .sum();

    // Reaction entropy
    entropy_products - entropy_reactants
}
