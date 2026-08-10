use super::equilibrium::{EquilibriumError, EquilibriumMode, solve_for_products};
use std::collections::HashMap;
use std::f64;
use thermo_species::Species;

const R: f64 = 8.314; // Ideal Gas Constant [J/mol•K]
pub const STD_REFERENCE_PRESSURE: f64 = 1.0; // [bar]

const MAX_FROZEN_ITERS: usize = 10;
const FROZEN_EXPANSION_TOL: f64 = 1.0e-6;

#[derive(Debug, Clone)]
pub enum MixtureError {
    /// A Gibbs-minimization solve (`solve_for_products`) failed to converge.
    EquilibriumSolveFailed {
        problem: &'static str, // e.g. "with_adiabatic_flame (HP)"
        iterations_used: usize,
    },
    /// The frozen-composition entropy root-find failed to converge.
    FrozenExpansionDidNotConverge { iterations_used: usize },
    /// `mix()` requires both mixtures to be at the same pressure.
    NotIsobaric {
        self_pressure_bar: f64,
        other_pressure_bar: f64,
    },
}

impl std::fmt::Display for MixtureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MixtureError::EquilibriumSolveFailed {
                problem,
                iterations_used,
            } => {
                write!(
                    f,
                    "{problem} failed to converge after {iterations_used} iterations"
                )
            }
            MixtureError::FrozenExpansionDidNotConverge { iterations_used } => {
                write!(
                    f,
                    "frozen expansion failed to converge after {iterations_used} iterations"
                )
            }
            MixtureError::NotIsobaric {
                self_pressure_bar,
                other_pressure_bar,
            } => {
                write!(
                    f,
                    "mix() requires equal pressures, got {self_pressure_bar} bar and {other_pressure_bar} bar"
                )
            }
        }
    }
}

impl std::error::Error for MixtureError {}

#[derive(Clone)]
pub struct Mixture {
    pub products: Vec<(f64, Species)>,
    pub temperature_k: f64,
    pub pressure_bar: f64,
    pub h_total: f64, // Total enthalpy.
    pub s_total: f64, // Total entropy.
    pub n_total: f64,
    pub avg_mw: f64,
    pub avg_cp: f64,
}

impl Mixture {
    /// No equilibrium solve — products = reactants, evaluated at (t, p).
    pub fn new(reactants: &[(f64, Species)], temperature_k: f64, pressure_bar: f64) -> Self {
        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(reactants, temperature_k, pressure_bar);

        Self {
            products: reactants.to_vec(),
            temperature_k,
            pressure_bar,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        }
    }

    /// TP problem: reequilibrate at an arbitrary, independently specified
    /// (temperature_k, pressure_bar) — may differ from self's current state.
    pub fn with_equilibrium_at(
        self,
        temperature_k: f64,
        pressure_bar: f64,
        only: Option<&[Species]>,
        insert: Option<&[Species]>,
    ) -> Result<Self, MixtureError> {
        let state = solve_for_products(
            &self.products,
            EquilibriumMode::TP {
                temperature_k,
                pressure_bar,
            },
            true,
            true,
            only,
            insert,
        )
        .map_err(|EquilibriumError::FailedToConverge { iterations_used }| {
            MixtureError::EquilibriumSolveFailed {
                problem: "with_equilibrium_at (TP)",
                iterations_used,
            }
        })?;

        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(&state.products, state.temperature_k, state.pressure_bar);

        Ok(Self {
            products: state.products,
            // Temperature and pressure won't change because mode is TP,
            // but reading them off `state` keeps this correct if the mode
            // handling in `solve_for_products` ever changes.
            temperature_k: state.temperature_k,
            pressure_bar: state.pressure_bar,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        })
    }

    /// Reequilibrates at self's current (temperature_k, pressure_bar) —
    /// useful after directly mutating `products`.
    pub fn with_equilibrium(
        self,
        only: Option<&[Species]>,
        insert: Option<&[Species]>,
    ) -> Result<Self, MixtureError> {
        let temperature_k = self.temperature_k;
        let pressure_bar = self.pressure_bar;
        self.with_equilibrium_at(temperature_k, pressure_bar, only, insert)
    }

    /// Adiabatic flame state: HP solve at self's current enthalpy (i.e. no
    /// heat added). `insert` pre-activates condensed species (CEA's
    /// `insert` keyword) — required when a condensed product dominates,
    /// e.g. MgO(cr) for metal fuels, or the gas-only solve crashes the
    /// temperature toward 0 K.
    pub fn with_adiabatic_flame(
        self,
        pressure_bar: f64,
        only: Option<&[Species]>,
        insert: Option<&[Species]>,
    ) -> Result<Self, MixtureError> {
        let h_target = self.h_total;
        self.solve_hp(
            "with_adiabatic_flame (HP)",
            h_target,
            pressure_bar,
            only,
            insert,
        )
    }

    /// Re-equilibrates the mixture after adding `q` joules of heat at
    /// constant pressure. Uses HP mode so species can dissociate/recombine
    /// as needed with the added heat.
    pub fn with_heat_addition(
        self,
        joules: f64,
        only: Option<&[Species]>,
        insert: Option<&[Species]>,
    ) -> Result<Self, MixtureError> {
        let h_target = self.h_total + joules;
        let pressure_bar = self.pressure_bar;
        self.solve_hp(
            "with_heat_addition (HP)",
            h_target,
            pressure_bar,
            only,
            insert,
        )
    }

    /// Equilibrium (shifting) expansion: SP solve — composition is allowed
    /// to re-equilibrate as it expands to `pressure_bar`.
    pub fn with_equilibrium_expansion(
        self,
        pressure_bar: f64,
        only: Option<&[Species]>,
        insert: Option<&[Species]>,
    ) -> Result<Self, MixtureError> {
        // s_total is J/K on self's mole basis; SP mode wants S/R per unit
        // feed mass [K*kmol/kg = K*mol/g].
        let s_over_r = self.s_total / (R * self.feed_mass());
        let state = solve_for_products(
            &self.products,
            EquilibriumMode::SP {
                s_over_r,
                pressure_bar,
            },
            true,
            true,
            only,
            insert,
        )
        .map_err(|EquilibriumError::FailedToConverge { iterations_used }| {
            MixtureError::EquilibriumSolveFailed {
                problem: "with_equilibrium_expansion (SP)",
                iterations_used,
            }
        })?;

        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(&state.products, state.temperature_k, state.pressure_bar);

        Ok(Self {
            products: state.products,
            temperature_k: state.temperature_k,
            pressure_bar: state.pressure_bar,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        })
    }

    /// Frozen expansion: composition fixed (no re-equilibration), only
    /// temperature is solved, from entropy conservation at fixed mole
    /// fractions. CEA's "frozen" nozzle mode — contrast with
    /// `with_equilibrium_expansion`, which lets composition shift. No
    /// chemistry solve happens here, so there's no `only`/`insert` to take.
    pub fn with_frozen_expansion(self, pressure_bar: f64) -> Result<Self, MixtureError> {
        let s_target = self.s_total;
        let mut t = self.temperature_k;

        for _ in 0..MAX_FROZEN_ITERS {
            let (h_total, s_total, n_total, avg_mw, avg_cp) =
                Self::state(&self.products, t, pressure_bar);
            let residual = s_total - s_target;
            if residual.abs() < FROZEN_EXPANSION_TOL {
                return Ok(Self {
                    products: self.products,
                    temperature_k: t,
                    pressure_bar,
                    h_total,
                    s_total,
                    n_total,
                    avg_mw,
                    avg_cp,
                });
            }
            // dS/dT = n_total * avg_cp / T for a fixed-composition ideal-gas mixture.
            t -= residual / (n_total * avg_cp / t);
        }
        Err(MixtureError::FrozenExpansionDidNotConverge {
            iterations_used: MAX_FROZEN_ITERS,
        })
    }

    pub fn mix(self, other: Self) -> Result<Self, MixtureError> {
        if self.pressure_bar != other.pressure_bar {
            return Err(MixtureError::NotIsobaric {
                self_pressure_bar: self.pressure_bar,
                other_pressure_bar: other.pressure_bar,
            });
        }
        let starting_h = self.h_total + other.h_total;
        let pressure_bar = self.pressure_bar;

        // Before mixing, bring the colder mix up to the same temperature as the hotter mix
        // for better initial guesses.
        let (cold, hot) = if self.temperature_k <= other.temperature_k {
            (self, other)
        } else {
            (other, self)
        };

        let cold_adj = if cold.temperature_k == hot.temperature_k {
            // Already same temperature nothing to do.
            cold
        } else {
            // Calculate new cold mixture products
            Mixture::new(&cold.products, hot.temperature_k, pressure_bar)
        };

        // Create reactant mix using adjusted other.
        let mut reactants_mix: HashMap<String, (f64, Species)> = HashMap::new();
        // Add all of original reatants to mix.
        for s in hot.products.iter() {
            let key = s.1.data().symbol();
            reactants_mix.insert(key.to_string(), *s);
        }
        for o in cold_adj.products.iter() {
            let key = o.1.data().symbol();
            if let Some(entry) = reactants_mix.get_mut(key) {
                // reatants already exists! Increment moles.
                entry.0 += o.0;
            } else {
                // reatants is new, add them.
                reactants_mix.insert(key.to_string(), *o);
            }
        }
        let mut reactants = reactants_mix
            .drain()
            .map(|(_, v)| v)
            .collect::<Vec<(f64, Species)>>();
        reactants.sort_by_key(|a| a.1);

        // Use HP mode to reequilibrate at the combined enthalpy and pressure.
        // This allows species to dissociate/recombine as needed (e.g., H2O → H2 + O, etc.).
        // starting_h is in J; HP mode wants H/R per unit feed mass [K·mol/g].
        let feed_mass: f64 = reactants.iter().map(|(n, s)| n * s.data().mw()).sum();
        let state = solve_for_products(
            &reactants,
            EquilibriumMode::HP {
                h_over_r: starting_h / (R * feed_mass),
                pressure_bar,
            },
            true,
            true,
            None,
            None,
        )
        .map_err(|EquilibriumError::FailedToConverge { iterations_used }| {
            MixtureError::EquilibriumSolveFailed {
                problem: "mix() (HP)",
                iterations_used,
            }
        })?;

        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(&state.products, state.temperature_k, state.pressure_bar);

        Ok(Self {
            products: state.products,
            temperature_k: state.temperature_k,
            pressure_bar: state.pressure_bar,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        })
    }

    pub fn scale(self, factor: f64) -> Self {
        let scaled_products: Vec<(f64, Species)> = self
            .products
            .iter()
            .map(|(moles, species)| (*moles * factor, *species))
            .collect();
        Self::new(&scaled_products, self.temperature_k, self.pressure_bar)
    }

    pub fn feed_mass(&self) -> f64 {
        self.products
            .iter()
            .map(|(moles, specie)| moles * specie.data().mw())
            .sum()
    }

    pub fn print_products(&self) {
        println!("Products Temperature: {:.3} K", self.temperature_k);
        for (mole, species) in self.products.iter() {
            if *mole < 1.0e-4 {
                continue;
            }
            println!("{:.4} {}", mole, species.data().symbol());
        }
    }

    fn solve_hp(
        self,
        problem: &'static str,
        h_target: f64,
        pressure_bar: f64,
        only: Option<&[Species]>,
        insert: Option<&[Species]>,
    ) -> Result<Self, MixtureError> {
        let state = solve_for_products(
            &self.products,
            EquilibriumMode::HP {
                h_over_r: h_target / (R * self.feed_mass()),
                pressure_bar,
            },
            true,
            true,
            only,
            insert,
        )
        .map_err(|EquilibriumError::FailedToConverge { iterations_used }| {
            MixtureError::EquilibriumSolveFailed {
                problem,
                iterations_used,
            }
        })?;

        let (h_total, s_total, n_total, avg_mw, avg_cp) =
            Self::state(&state.products, state.temperature_k, state.pressure_bar);

        Ok(Self {
            products: state.products,
            temperature_k: state.temperature_k,
            pressure_bar: state.pressure_bar,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        })
    }

    fn state(
        products: &[(f64, Species)],
        temperature_k: f64,
        pressure_bar: f64,
    ) -> (f64, f64, f64, f64, f64) {
        let n_total: f64 = products.iter().map(|(moles, _)| *moles).sum();
        let products: Vec<(f64, Species, f64)> = products
            .iter()
            .map(|(moles, species)| {
                (
                    *moles, // number of moles, n
                    *species,
                    moles / n_total, // mole fraction, x
                )
            })
            .collect();
        let h_total: f64 = products
            .iter()
            .map(|(n_i, species_i, _)| n_i * species_i.data().h(temperature_k))
            .sum();
        let s_total: f64 = products
            .iter()
            .map(|(n_i, species_i, x_i)| {
                // Extinct species are reported with exactly 0 moles; skip
                // them or 0 * ln(0) poisons the sum with NaN.
                if *n_i <= 0.0 {
                    return 0.0;
                }
                n_i * (species_i.data().s(temperature_k)
                    - R * (x_i * pressure_bar / STD_REFERENCE_PRESSURE).ln())
            })
            .sum();
        let avg_mw: f64 = products
            .iter()
            .map(|(_, species_i, x_i)| x_i * species_i.data().mw())
            .sum();
        let avg_cp: f64 = products
            .iter()
            .map(|(_, species_i, x_i)| x_i * species_i.data().cp(temperature_k))
            .sum();

        (h_total, s_total, n_total, avg_mw, avg_cp)
    }
}
