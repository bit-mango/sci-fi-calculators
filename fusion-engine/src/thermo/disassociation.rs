use crate::constants::{R, STD_REFERENCE_PRESSURE};

// Step 1: Determine chemical equilibrium of Propellants(H2 ⇌ 2H, CO ⇌ C + O  disassociation fractions).
// Need to maximize Hydrogen disassociation but minimize Carbon Monoxide disassociation(to minimize coking).
// The disassociation fraction is given by:
// Kp(T) = exp(-∆G(T)/RT)
// Kp = [4⍺^2 / (1-⍺^2)]* (P/P_0)
// where ⍺ is the disassociation fraction.
// Solve for ⍺.
// ⍺ = sqrt(N/D)
// where N = Kp * P_0 / P, D = 4 + Kp * P_0 / P
// where ∆G(T) = ∆H-T∆S
pub fn calculate_disassociation_fraction(
    chamber_temperature_k: f64,
    chamber_pressure_bar: f64,
    enthalpy: f64,
    entropy: f64,
    product_factor: f64,
) -> f64 {
    let gibbs = enthalpy - chamber_temperature_k * entropy;
    let kp = (-gibbs / (R * chamber_temperature_k)).exp();

    let numerator = kp * STD_REFERENCE_PRESSURE / chamber_pressure_bar;
    let denominator = product_factor + numerator;
    let alpha = (numerator / denominator).sqrt();

    alpha
}
