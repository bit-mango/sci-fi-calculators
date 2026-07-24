use crate::constants::*;
use crate::thermo::disassociation::calculate_disassociation_fraction;
use crate::thermo::fluid_properties::{TemperatureDependentProperty, ThermoReference};
use crate::thermo::reactions::{get_rxn_enthalpy, get_rxn_entropy};
use crate::thermo::species;
use nalgebra::{DMatrix, DVector, dmatrix, dvector};
use std::collections::HashMap;

#[derive(Default)]
pub struct MixtureState {
    pub products: Mixture,
    pub h_total: f64, // Total enthalpy.
    pub s_total: f64, // Total entropy.
    pub n_total: f64, // Total moles.
    pub avg_mw: f64,
    pub avg_cp: f64,
}

#[derive(Clone, Default)]
pub struct Mixture {
    pub reactants: Vec<(f64, Species)>, // moles, species, Temperature Dependent Property
}

impl Mixture {
    pub fn new(reactants: Vec<(f64, Species)>) -> Self {
        Self { reactants }
    }

    pub fn mix(&self, other: &Self) -> Self {
        let mut reatants_mix: HashMap<String, (f64, Species)> = HashMap::new();
        // Add all of original reatants to mix.
        for s in self.reactants.iter() {
            let key = s.1.symbol();
            reatants_mix.insert(key, s.clone());
        }
        for o in other.reactants.iter() {
            let key = o.1.symbol();
            if let Some(entry) = reatants_mix.get_mut(&key) {
                // reatants already exists! Increment moles.
                entry.0 += o.0;
            } else {
                // reatants is new, add them.
                reatants_mix.insert(key, o.clone());
            }
        }
        let mut reactants = reatants_mix
            .drain()
            .map(|(_, v)| v)
            .collect::<Vec<(f64, Species)>>();
        reactants.sort_by(|a, b| a.1.symbol().cmp(&b.1.symbol()));
        Self { reactants }
    }

    pub fn scale(&self, factor: f64) -> Self {
        Self {
            reactants: self
                .reactants
                .iter()
                .map(|(mol, sp)| (*mol * factor, *sp))
                .collect(),
        }
    }

    pub fn reaction_pool(&self) -> (Vec<Species>, Vec<(f64, Species)>) {
        let mut element_pool: HashMap<Species, f64> = HashMap::new();
        self.reactants
            .iter()
            .for_each(|(parent_moles, parent_species)| {
                let children = parent_species.constituents();
                children.iter().for_each(|(child_moles, child_species)| {
                    if let Some(existing) = element_pool.get_mut(child_species) {
                        *existing += parent_moles * child_moles;
                    } else {
                        element_pool.insert(*child_species, parent_moles * child_moles);
                    }
                });
            });
        // Filter out any species we don't have constituent elements for.
        let species_pool = Species::all()
            .iter()
            .filter(|species| {
                let children = species.constituents();
                for child in children {
                    if !element_pool.contains_key(&child.1) {
                        return false;
                    }
                }
                return true;
            })
            .map(|s| *s)
            .collect::<Vec<Species>>();
        let element_pool = element_pool
            .drain()
            .map(|(species, moles)| (moles, species))
            .collect::<Vec<(f64, Species)>>();
        (species_pool, element_pool)
    }

    pub fn solve_for_products(
        &self,
        temperature_k: f64,
        pressure_bar: f64,
        tr: &ThermoReference,
    ) -> Mixture {
        let (species_pool, element_pool) = self.reaction_pool();
        // Precompute mu_not.
        let mu_not = species_pool
            .iter()
            .map(|species| {
                let tdp = tr.get_tdp(&species.symbol());
                let mu_not_i = tdp.h(temperature_k) - temperature_k * tdp.s(temperature_k);
                mu_not_i
            })
            .collect::<Vec<f64>>();

        // Each row is a species, and each column is how many moles of an element are in that species.
        // Example:
        //      H   O
        // H2O  2   1
        // OH   1   1
        // H2   2   0
        // H    1   0
        // O2   0   2
        // O    0   1
        let a = species_pool
            .iter()
            .map(|species| {
                let children = species.constituents();
                element_pool
                    .iter()
                    .map(|(_, elem)| {
                        children
                            .iter()
                            .find(|(_, child)| child == elem)
                            .unwrap_or(&(0.0, *elem))
                            .0
                    })
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<Vec<f64>>>();

        // n_i = n * c_i * exp(Σⱼ πⱼ·a[i][j])
        // i is in range form 0..species_pool.length
        // j is in range form 0..element_pool.length
        let n = |i: usize, pi: &[f64], ln_n: f64| {
            let sum = pi
                .iter()
                .zip(a[i].iter())
                .map(|(pi_i, a_i_j)| pi_i * a_i_j)
                .sum::<f64>();
            ln_n.exp()
                * (STD_REFERENCE_PRESSURE / pressure_bar)
                * (sum - mu_not[i] / (R * temperature_k)).exp()
        };
        let j_len = a[0].len();
        // b_j is element_pool[j].0
        let b = element_pool
            .iter()
            .map(|(moles, _)| *moles)
            .collect::<Vec<f64>>();

        let (mut pi, mut ln_n) = self.guess_pi_and_ln_n(
            temperature_k,
            pressure_bar,
            j_len,
            &a,
            &b,
            &mu_not,
            &species_pool,
        );

        let mut iterations = 0;
        let mut final_n_vec = vec![];
        let mut lambda = 1.0e-3;

        loop {
            let (n_vec, f, j) = residual_and_jacobian(&pi, ln_n, &a, &b, &n);
            if f.norm() < 1.0e-6 {
                final_n_vec = n_vec;
                break;
            }

            if iterations <= 30 || iterations >= 55 {
                println!("[{}] f norm: {}, lambda: {}", iterations, f.norm(), lambda);
                let species_name = species_pool
                    .iter()
                    .map(|species| species.symbol())
                    .zip(n_vec.iter())
                    .map(|(species, moles)| (*moles, species))
                    .collect::<Vec<(f64, String)>>();
                println!("[{}] n: {:?}", iterations, species_name);
                println!("[{}] pi: {:?}, ln_n: {}", iterations, pi, ln_n);
            }

            if iterations == 4999 {
                println!("[5000] f: {}", f);
                println!("[5000] n: {:?}", n_vec);
                let species_name = species_pool
                    .iter()
                    .map(|species| species.symbol())
                    .zip(n_vec.iter())
                    .map(|(species, moles)| (*moles, species))
                    .collect::<Vec<(f64, String)>>();
                println!("[5000] species: {:?}", species_name);

                panic!("Failed to find solution to mixture!");
            }

            let current_norm = f.norm();
            // Solve via normal equations (JᵀJ + λ·diag(JᵀJ))·δ = -Jᵀf rather than
            // damping J directly. J here is symmetric but not necessarily positive
            // definite away from the solution (e.g. the last diagonal entry
            // n_sum - exp(ln_n) can go negative), so a diagonal shift on J itself
            // is not guaranteed to produce a descent direction for ‖f‖ — that's
            // why the line search could exhaust all 50 lambda-doublings and find
            // nothing. JᵀJ + λD is always positive definite for λ>0, so -Jᵀf-based
            // steps are always a descent direction of the merit function ½‖f‖².
            let jt = j.transpose();
            let jtj = &jt * &j;
            let neg_jtf = -(&jt * &f);
            let mut accepted = false;
            for _ in 0..50 {
                let mut jtj_damped = jtj.clone();
                for k in 0..(j_len + 1) {
                    jtj_damped[(k, k)] += lambda * jtj[(k, k)].max(1.0e-10);
                }
                if let Some(delta) = jtj_damped.clone().lu().solve(&neg_jtf) {
                    // Trust-region step limiter (NASA CEA-style): since pi/ln_n
                    // feed an exp(), an unbounded Newton step can send some
                    // species' moles to ~0 and others to overflow in one shot,
                    // which then singularizes J on the next iteration. Cap the
                    // largest per-step change so no n_i can move by more than
                    // ~e^2 in a single iteration.
                    let max_step = delta.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                    let alpha = if max_step > 25.0 {
                        25.0 / max_step
                    } else {
                        1.0
                    };

                    let trial_pi: Vec<f64> = (0..j_len).map(|k| pi[k] + alpha * delta[k]).collect();
                    let trial_ln_n = ln_n + alpha * delta[j_len];
                    let (_, trial_f, _) = residual_and_jacobian(&trial_pi, trial_ln_n, &a, &b, &n);
                    if trial_f.iter().all(|v| v.is_finite()) && trial_f.norm() < current_norm {
                        pi = trial_pi;
                        ln_n = trial_ln_n;
                        lambda *= 0.5;
                        accepted = true;
                        break;
                    }
                    lambda *= 2.0;
                }
            }
            if !accepted {
                panic!(
                    "Line search failed to find a reducing step at iteration {}",
                    iterations
                );
            }

            iterations += 1;
        }

        let products = final_n_vec
            .iter()
            .zip(species_pool.iter())
            .map(|(&moles, &species)| (moles, species))
            .collect::<Vec<(f64, Species)>>();

        Mixture::new(products)
    }

    fn guess_pi_and_ln_n(
        &self,
        temperature_k: f64,
        pressure_bar: f64,
        j_len: usize,
        a: &Vec<Vec<f64>>,
        b: &Vec<f64>,
        mu_not: &Vec<f64>,
        species_pool: &Vec<Species>,
    ) -> (Vec<f64>, f64) {
        // Base ln_n_guess off of starting moles of the system.
        let ln_n_guess = self
            .reactants
            .iter()
            .map(|(moles, _)| *moles)
            .sum::<f64>()
            .ln();

        // Identify anchor species.
        let anchors: Vec<(usize, f64)> = species_pool
            .iter()
            .enumerate() // TODO convert to filter_map?
            .map(|(i, species)| {
                let mut anchors_i = None;
                for (moles, reactant) in self.reactants.iter() {
                    if reactant == species {
                        anchors_i = Some((i, *moles));
                    }
                }
                anchors_i
            })
            .filter(|anchors_i| anchors_i.is_some())
            .map(|anchors_i| anchors_i.unwrap())
            .collect();

        // Build A (anchors.len() × j_len) and rhs (anchors.len()) from the anchor list:
        //   A[row] = a[i]                                    for (i, feed_moles) in anchors
        //   rhs[row] = ln(feed_moles) - ln_n_guess + mu_not[i]/(R*T) - ln(P_ref/P)

        // Solve (AᵀA + eps·I) · pi = Aᵀ·rhs   instead of a plain square solve.
        let mut a_ = vec![vec![0.0; j_len]; anchors.len()];
        let mut rhs = vec![0.0; anchors.len()];
        for (row, (i, feed_moles)) in anchors.iter().enumerate() {
            for j in 0..j_len {
                a_[row][j] = a[*i][j];
                rhs[row] = feed_moles.ln() - ln_n_guess + mu_not[*i] / (R * temperature_k)
                    - (STD_REFERENCE_PRESSURE / pressure_bar).ln();
            }
        }

        let rows = anchors.len();
        let flattened: Vec<f64> = a_.into_iter().flatten().collect();
        let a_matrix = DMatrix::from_row_slice(rows, j_len, &flattened);
        let rhs_vector = DVector::from_vec(rhs);

        let ata = a_matrix.transpose() * &a_matrix; // j_len * j_len, symmeric
        let atb = a_matrix.transpose() * &rhs_vector; // j_len

        // ridge term, scaled to A transpose A's magnitude rather than a fixed constant
        let eps = 1.0e-8 * ata.amax().max(1.0e-10);
        let mut ata_reg = ata.clone();
        for k in 0..j_len {
            ata_reg[(k, k)] += eps;
        }

        let pi_guess: Vec<f64> = match ata_reg.lu().solve(&atb) {
            Some(pi) => pi.iter().copied().collect(),
            None => vec![0.0; j_len], // degenerate anchors (e.g. duplicate rows) - fall back
        };

        (pi_guess, ln_n_guess)
    }

    pub fn h_total(&self, n: &Vec<f64>, temperature_k: f64, tr: &ThermoReference) -> f64 {
        let mut h_total = 0.0;
        let mut i = 0;
        for specie in self.reactants.iter() {
            h_total += n[i] * tr.get_tdp(&specie.1.symbol()).h(temperature_k);
            i += 1;
        }

        h_total
    }

    pub fn x(&self, n: &Vec<f64>) -> Vec<f64> {
        let n_sum: f64 = n.iter().sum();
        n.iter().map(|n_i| n_i / n_sum).collect()
    }

    pub fn s_total(
        &self,
        x: &Vec<f64>,
        n: &Vec<f64>,
        temperature_k: f64,
        pressure_bar: f64,
        tr: &ThermoReference,
    ) -> f64 {
        let mut s_total = 0.0;
        let mut i = 0;
        for specie in self.reactants.iter() {
            s_total += n[i]
                * (tr.get_tdp(&specie.1.symbol()).s(temperature_k)
                    - R * (x[i] * pressure_bar / STD_REFERENCE_PRESSURE).ln());
            i += 1;
        }

        s_total
    }

    pub fn avg_mw(&self, x: &Vec<f64>) -> f64 {
        let mut mw_total = 0.0;
        let mut i = 0;
        for specie in self.reactants.iter() {
            mw_total += x[i] * specie.1.mw();
            i += 1;
        }

        mw_total
    }

    pub fn avg_cp(&self, x: &Vec<f64>, temperature_k: f64, tr: &ThermoReference) -> f64 {
        let mut cp_total = 0.0;
        let mut i = 0;
        for specie in self.reactants.iter() {
            cp_total += x[i] * tr.get_tdp(&specie.1.symbol()).cp(temperature_k);
            i += 1;
        }

        cp_total
    }

    pub fn feed_mass(&self) -> f64 {
        self.reactants
            .iter()
            .map(|specie| specie.0 * specie.1.mw())
            .sum()
    }

    pub fn state(
        &self,
        temperature_k: f64,
        pressure_bar: f64,
        tr: &ThermoReference,
    ) -> MixtureState {
        // TODO state always tries to solve this, but for low temps its not needed.
        let products = if temperature_k <= 0.0 {
            self.clone()
        } else {
            self.solve_for_products(temperature_k, pressure_bar, tr)
        };

        let n = products
            .reactants
            .iter()
            .map(|(moles, _)| *moles)
            .collect::<Vec<f64>>();
        let h_total = products.h_total(&n, temperature_k, tr);
        let x = products.x(&n);
        let s_total = products.s_total(&x, &n, temperature_k, pressure_bar, tr);

        let n_total = n.iter().sum();
        let avg_mw: f64 = products.avg_mw(&x);
        let avg_cp: f64 = products.avg_cp(&x, temperature_k, tr);

        MixtureState {
            products,
            h_total,
            s_total,
            n_total,
            avg_mw,
            avg_cp,
        }
    }
}

fn residual_and_jacobian(
    pi: &[f64],
    ln_n: f64,
    a: &Vec<Vec<f64>>,
    b: &Vec<f64>,
    n_i: &impl Fn(usize, &[f64], f64) -> f64,
) -> (Vec<f64>, DVector<f64>, DMatrix<f64>) {
    let i_len = a.len();
    let j_len = a[0].len();
    let mut n_vec = vec![0.0; i_len];
    for i in 0..i_len {
        n_vec[i] = n_i(i, pi, ln_n);
    }

    let mut f_vec = vec![0.0; j_len + 1];
    for j in 0..j_len {
        let mut sum = 0.0;
        for i in 0..i_len {
            sum += n_vec[i] * a[i][j];
        }
        sum -= b[j];
        f_vec[j] = sum;
    }
    let n_sum = n_vec.iter().sum::<f64>();
    f_vec[j_len] = n_sum - ln_n.exp();
    let mut j_vec = vec![vec![0.0; j_len + 1]; j_len + 1];
    for j in 0..j_len {
        for k in 0..j_len {
            let mut sum = 0.0;
            for i in 0..i_len {
                sum += n_vec[i] * a[i][j] * a[i][k];
            }
            j_vec[j][k] = sum;
        }
    }

    for j in 0..j_len {
        let mut sum = 0.0;
        for i in 0..i_len {
            sum += n_vec[i] * a[i][j];
        }
        j_vec[j][j_len] = sum;
        j_vec[j_len][j] = sum;
    }

    j_vec[j_len][j_len] = n_sum - ln_n.exp();
    let f_vec = DVector::from_vec(f_vec);
    let flattened_j = j_vec.iter().flatten().map(|e| *e).collect();
    // TODO from_vec expects column major data, but flatten gives row major data. It is fine because j_vec is symmetrical ie row count == col count
    // but thats messy lets fix it.
    let j_matrix = DMatrix::from_vec(j_len + 1, j_len + 1, flattened_j);
    (n_vec, f_vec, j_matrix)
}

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub enum Species {
    H,
    H2,
    C,
    O,
    CO,
    N2,
    N,
    H2O,
    OH,
    CO2,
    O2,
    CH4,
    NH3,
}

impl Species {
    pub fn all() -> Vec<Self> {
        vec![
            Self::H,
            Self::H2,
            Self::C,
            Self::O,
            Self::CO,
            Self::N2,
            Self::N,
            Self::H2O,
            Self::OH,
            Self::CO2,
            Self::O2,
            Self::CH4,
            Self::NH3,
        ]
    }
    pub fn symbol(&self) -> String {
        match self {
            Species::H => "H".to_string(),
            Species::H2 => "H2".to_string(),
            Species::C => "C".to_string(),
            Species::O => "O".to_string(),
            Species::O2 => "O2".to_string(),
            Species::CO => "CO".to_string(),
            Species::CO2 => "CO2".to_string(),
            Species::N2 => "N2".to_string(),
            Species::N => "N".to_string(),
            Species::H2O => "H2O".to_string(),
            Species::OH => "OH".to_string(),
            Species::CH4 => "CH4".to_string(),
            Species::NH3 => "NH3".to_string(),
        }
    }

    pub fn mw(&self) -> f64 {
        match self {
            Species::H => H_MW,
            Species::H2 => H2_MW,
            Species::C => C_MW,
            Species::O => O_MW,
            Species::O2 => O2_MW,
            Species::CO => CO_MW,
            Species::CO2 => CO2_MW,
            Species::N2 => N2_MW,
            Species::N => N_MW,
            Species::H2O => H2O_MW,
            Species::OH => OH_MW,
            Species::CH4 => CH4_MW,
            Species::NH3 => NH3_MW,
        }
    }

    pub fn constituents(&self) -> Vec<(f64, Self)> {
        match self {
            Species::H => vec![(1.0, Species::H)],
            Species::H2 => vec![(2.0, Species::H)],
            Species::C => vec![(1.0, Species::C)],
            Species::O => vec![(1.0, Species::O)],
            Species::O2 => vec![(2.0, Species::O)],
            Species::CO => vec![(1.0, Species::C), (1.0, Species::O)],
            Species::CO2 => vec![(1.0, Species::C), (2.0, Species::O)],
            Species::N2 => vec![(2.0, Species::N)],
            Species::N => vec![(1.0, Species::N)],
            Species::H2O => vec![(2.0, Species::H), (1.0, Species::O)],
            Species::OH => vec![(1.0, Species::O), (1.0, Species::H)],
            Species::CH4 => vec![(1.0, Species::C), (4.0, Species::H)],
            Species::NH3 => vec![(1.0, Species::N), (3.0, Species::H)],
        }
    }
}
