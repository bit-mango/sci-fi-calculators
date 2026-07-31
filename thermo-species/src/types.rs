#[rustfmt::skip] // So it doesn't convert the single line periods into multi-line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Element {
    // Electron e-
    E,
    // Period 1
    H, D, He, // D is Deuterium

    // Period 2
    Li, Be, B, C, N, O, F, Ne,

    // Period 3
    Na, Mg, Al, Si, P, S, Cl, Ar,

    // Period 4
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,

    // Period 5
    Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,

    // Period 6
    Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
    Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,

    // Period 7
    Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, FmA, Md, No, Lr,
    Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
}

#[derive(Debug)]
pub struct SpeciesData<const C: usize, const T: usize> {
    pub symbol: &'static str,
    pub constituents: [(f64, Element); C],
    pub temperature_data: [(f64, f64, [f64; 9]); T],
    pub mw: f64,
    pub h_formation: f64,
    pub phase: u8,
}

pub trait AnySpeciesData {
    fn symbol(&self) -> &str;
    fn constituents(&self) -> &[(f64, Element)];
    fn temperature_data(&self) -> &[(f64, f64, [f64; 9])];
    fn mw(&self) -> f64;
    fn h_formation(&self) -> f64;
    fn phase(&self) -> u8;
}

impl<const C: usize, const T: usize> AnySpeciesData for SpeciesData<C, T> {
    fn symbol(&self) -> &str {
        self.symbol
    }
    fn constituents(&self) -> &[(f64, Element)] {
        &self.constituents
    }
    fn temperature_data(&self) -> &[(f64, f64, [f64; 9])] {
        &self.temperature_data
    }
    fn mw(&self) -> f64 {
        self.mw
    }
    fn h_formation(&self) -> f64 {
        self.h_formation
    }
    fn phase(&self) -> u8 {
        self.phase
    }
}
