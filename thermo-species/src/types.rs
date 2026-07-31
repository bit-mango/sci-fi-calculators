use std::str::FromStr;

#[rustfmt::skip] // So it doesn't convert the single line periods into multi-line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Element {
    // Period 1
    H, He,

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

#[derive(Debug, PartialEq, Eq)]
pub struct ParseElementError;

impl FromStr for Element {
    type Err = ParseElementError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(Self::H),
            "He" => Ok(Self::He),
            "Li" => Ok(Self::Li),
            "Be" => Ok(Self::Be),
            "B" => Ok(Self::B),
            "C" => Ok(Self::C),
            "N" => Ok(Self::N),
            "O" => Ok(Self::O),
            "F" => Ok(Self::F),
            "Ne" => Ok(Self::Ne),
            "Na" => Ok(Self::Na),
            "Mg" => Ok(Self::Mg),
            "Al" => Ok(Self::Al),
            "Si" => Ok(Self::Si),
            "P" => Ok(Self::P),
            "S" => Ok(Self::S),
            "Cl" => Ok(Self::Cl),
            "Ar" => Ok(Self::Ar),
            "K" => Ok(Self::K),
            "Ca" => Ok(Self::Ca),
            "Sc" => Ok(Self::Sc),
            "Ti" => Ok(Self::Ti),
            "V" => Ok(Self::V),
            "Cr" => Ok(Self::Cr),
            "Mn" => Ok(Self::Mn),
            "Fe" => Ok(Self::Fe),
            "Co" => Ok(Self::Co),
            "Ni" => Ok(Self::Ni),
            "Cu" => Ok(Self::Cu),
            "Zn" => Ok(Self::Zn),
            "Ga" => Ok(Self::Ga),
            "Ge" => Ok(Self::Ge),
            "As" => Ok(Self::As),
            "Se" => Ok(Self::Se),
            "Br" => Ok(Self::Br),
            "Kr" => Ok(Self::Kr),
            "Rb" => Ok(Self::Rb),
            "Sr" => Ok(Self::Sr),
            "Y" => Ok(Self::Y),
            "Zr" => Ok(Self::Zr),
            "Nb" => Ok(Self::Nb),
            "Mo" => Ok(Self::Mo),
            "Tc" => Ok(Self::Tc),
            "Ru" => Ok(Self::Ru),
            "Rh" => Ok(Self::Rh),
            "Pd" => Ok(Self::Pd),
            "Ag" => Ok(Self::Ag),
            "Cd" => Ok(Self::Cd),
            "In" => Ok(Self::In),
            "Sn" => Ok(Self::Sn),
            "Sb" => Ok(Self::Sb),
            "Te" => Ok(Self::Te),
            "I" => Ok(Self::I),
            "Xe" => Ok(Self::Xe),
            "Cs" => Ok(Self::Cs),
            "Ba" => Ok(Self::Ba),
            "La" => Ok(Self::La),
            "Ce" => Ok(Self::Ce),
            "Pr" => Ok(Self::Pr),
            "Nd" => Ok(Self::Nd),
            "Pm" => Ok(Self::Pm),
            "Sm" => Ok(Self::Sm),
            "Eu" => Ok(Self::Eu),
            "Gd" => Ok(Self::Gd),
            "Tb" => Ok(Self::Tb),
            "Dy" => Ok(Self::Dy),
            "Ho" => Ok(Self::Ho),
            "Er" => Ok(Self::Er),
            "Tm" => Ok(Self::Tm),
            "Yb" => Ok(Self::Yb),
            "Lu" => Ok(Self::Lu),
            "Hf" => Ok(Self::Hf),
            "Ta" => Ok(Self::Ta),
            "W" => Ok(Self::W),
            "Re" => Ok(Self::Re),
            "Os" => Ok(Self::Os),
            "Ir" => Ok(Self::Ir),
            "Pt" => Ok(Self::Pt),
            "Au" => Ok(Self::Au),
            "Hg" => Ok(Self::Hg),
            "Tl" => Ok(Self::Tl),
            "Pb" => Ok(Self::Pb),
            "Bi" => Ok(Self::Bi),
            "Po" => Ok(Self::Po),
            "At" => Ok(Self::At),
            "Rn" => Ok(Self::Rn),
            "Fr" => Ok(Self::Fr),
            "Ra" => Ok(Self::Ra),
            "Ac" => Ok(Self::Ac),
            "Th" => Ok(Self::Th),
            "Pa" => Ok(Self::Pa),
            "U" => Ok(Self::U),
            "Np" => Ok(Self::Np),
            "Pu" => Ok(Self::Pu),
            "Am" => Ok(Self::Am),
            "Cm" => Ok(Self::Cm),
            "Bk" => Ok(Self::Bk),
            "Cf" => Ok(Self::Cf),
            "Es" => Ok(Self::Es),
            "FmA" => Ok(Self::FmA),
            "Md" => Ok(Self::Md),
            "No" => Ok(Self::No),
            "Lr" => Ok(Self::Lr),
            "Rf" => Ok(Self::Rf),
            "Db" => Ok(Self::Db),
            "Sg" => Ok(Self::Sg),
            "Bh" => Ok(Self::Bh),
            "Hs" => Ok(Self::Hs),
            "Mt" => Ok(Self::Mt),
            "Ds" => Ok(Self::Ds),
            "Rg" => Ok(Self::Rg),
            "Cn" => Ok(Self::Cn),
            "Nh" => Ok(Self::Nh),
            "Fl" => Ok(Self::Fl),
            "Mc" => Ok(Self::Mc),
            "Lv" => Ok(Self::Lv),
            "Ts" => Ok(Self::Ts),
            "Og" => Ok(Self::Og),
            _ => Err(ParseElementError),
        }
    }
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
