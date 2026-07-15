use super::constants::R;

pub struct PolynomialCoefficients {
    a_1: f64,
    a_2: f64,
    a_3: f64,
    a_4: f64,
    a_5: f64,
    a_6: f64,
    a_7: f64,
    b_1: f64,
    b_2: f64,
    min_temperture: f64,
    max_temperature: f64,
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
        min_temperture: f64,
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
            min_temperture,
            max_temperature,
        }
    }
}

pub struct TemperatureDependentProperty {
    coefficients: Vec<PolynomialCoefficients>,
}

impl TemperatureDependentProperty {
    pub fn new() -> Self {
        Self {
            coefficients: vec![],
        }
    }
    // TODO could add more checks to make sure we dont have overlaps, and to order them.
    pub fn with_coefficients(mut self, pc: PolynomialCoefficients) -> Self {
        self.coefficients.push(pc);

        self
    }
    pub fn with_hydrogen_coefficients(self) -> Self {
        // H: From 200 K  to 1,000 K
        let a_1 = 0.000000000e+00;
        let a_2 = 0.000000000e+00;
        let a_3 = 2.500000000e+00;
        let a_4 = 0.000000000e+00;
        let a_5 = 0.000000000e+00;
        let a_6 = 0.000000000e+00;
        let a_7 = 0.000000000e+00;
        let b_1 = 2.547370801e+04;
        let b_2 = -4.466828530e-01;
        let h_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 190.0, 1_000.0,
        );
        let h_tdp = self.with_coefficients(h_pc);

        // H: From 1,000 K  to 6,000 K
        let a_1 = 6.078774250e+01;
        let a_2 = -1.819354417e-01;
        let a_3 = 2.500211817e+00;
        let a_4 = -1.226512864e-07;
        let a_5 = 3.732876330e-11;
        let a_6 = -5.687744560e-15;
        let a_7 = 3.410210197e-19;
        let b_1 = 2.547486398e+04;
        let b_2 = -4.481917770e-01;
        let h_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 1_000.0, 6_000.0,
        );
        let h_tdp = h_tdp.with_coefficients(h_pc);

        // H: From 6,000 K  to 20,000 K
        let a_1 = 2.173757694e+08;
        let a_2 = -1.312035403e+05;
        let a_3 = 3.399174200e+01;
        let a_4 = -3.813999680e-03;
        let a_5 = 2.432854837e-07;
        let a_6 = -7.694275540e-12;
        let a_7 = 9.644105630e-17;
        let b_1 = 1.067638086e+06;
        let b_2 = -2.742301051e+02;
        let h_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 6_000.0, 20_000.0,
        );

        h_tdp.with_coefficients(h_pc)
    }

    pub fn with_hydrogen2_coefficients(self) -> Self {
        // H2: From 200 K to 1,000 K
        let a_1 = 4.078323210e+04;
        let a_2 = -8.009186040e+02;
        let a_3 = 8.214702010e+00;
        let a_4 = -1.269714457e-02;
        let a_5 = 1.753605076e-05;
        let a_6 = -1.202860270e-08;
        let a_7 = 3.368093490e-12;
        let b_1 = 2.682484665e+03;
        let b_2 = -3.043788844e+01;
        let h2_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 190.0, 1_000.0,
        );
        let h2_tdp = self.with_coefficients(h2_pc);

        // H2: From 1,000 K to 6,000 K
        let a_1 = 5.608128010e+05;
        let a_2 = -8.371504740e+02;
        let a_3 = 2.975364532e+00;
        let a_4 = 1.252249124e-03;
        let a_5 = -3.740716190e-07;
        let a_6 = 5.936625200e-11;
        let a_7 = -3.606994100e-15;
        let b_1 = 5.339824410e+03;
        let b_2 = -2.202774769e+00;
        let h2_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 1_000.0, 6_000.0,
        );
        let h2_tdp = h2_tdp.with_coefficients(h2_pc);

        // H2: From 6,000 K to 20,000 K
        let a_1 = 4.966884120e+08;
        let a_2 = -3.147547149e+05;
        let a_3 = 7.984121880e+01;
        let a_4 = -8.414789210e-03;
        let a_5 = 4.753248350e-07;
        let a_6 = -1.371873492e-11;
        let a_7 = 1.605461756e-16;
        let b_1 = 2.488433516e+06;
        let b_2 = -6.695728110e+02;
        let h2_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 6_000.0, 20_000.0,
        );

        h2_tdp.with_coefficients(h2_pc)
    }

    pub fn with_carbon_coefficients(self) -> Self {
        // C: From 200 K to 1,000 K.
        let a_1 = 6.495031470e+02;
        let a_2 = -9.649010860e-01;
        let a_3 = 2.504675479e+00;
        let a_4 = -1.281448025e-05;
        let a_5 = 1.980133654e-08;
        let a_6 = -1.606144025e-11;
        let a_7 = 5.314483411e-15;
        let b_1 = 8.545763110e+04;
        let b_2 = 4.747924288e+00;
        let c_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 190.0, 1_000.0,
        );
        let c_tdp = self.with_coefficients(c_pc);

        // C: From 1,000 K to 6,000 K.
        let a_1 = -1.289136472e+05;
        let a_2 = 1.719528572e+02;
        let a_3 = 2.646044387e+00;
        let a_4 = -3.353068950e-04;
        let a_5 = 1.742092740e-07;
        let a_6 = -2.902817829e-11;
        let a_7 = 1.642182385e-15;
        let b_1 = 8.410597850e+04;
        let b_2 = 4.130047418e+00;
        let c_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 1_000.0, 6_000.0,
        );
        let c_tdp = c_tdp.with_coefficients(c_pc);

        // C: From 6,000 K to 20,000 K.
        let a_1 = 4.432528010e+08;
        let a_2 = -2.886018412e+05;
        let a_3 = 7.737108320e+01;
        let a_4 = -9.715281890e-03;
        let a_5 = 6.649595330e-07;
        let a_6 = -2.230078776e-11;
        let a_7 = 2.899388702e-16;
        let b_1 = 2.355273444e+06;
        let b_2 = -6.405123160e+02;

        let c_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 6_000.0, 20_000.0,
        );

        c_tdp.with_coefficients(c_pc)
    }

    pub fn with_oxygen_coefficients(self) -> Self {
        // O: From 200 K to 1,000 K.
        let a_1 = -7.953611300e+03;
        let a_2 = 1.607177787e+02;
        let a_3 = 1.966226438e+00;
        let a_4 = 1.013670310e-03;
        let a_5 = -1.110415423e-06;
        let a_6 = 6.517507500e-10;
        let a_7 = -1.584779251e-13;
        let b_1 = 2.840362437e+04;
        let b_2 = 8.404241820e+00;
        let o_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 190.0, 1_000.0,
        );
        let o_tdp = self.with_coefficients(o_pc);

        // O: From 1,000 K to 6,000 K.
        let a_1 = 2.619020262e+05;
        let a_2 = -7.298722030e+02;
        let a_3 = 3.317177270e+00;
        let a_4 = -4.281334360e-04;
        let a_5 = 1.036104594e-07;
        let a_6 = -9.438304330e-12;
        let a_7 = 2.725038297e-16;
        let b_1 = 3.392428060e+04;
        let b_2 = -6.679585350e-01;
        let o_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 1_000.0, 6_000.0,
        );
        let o_tdp = o_tdp.with_coefficients(o_pc);

        // O: From 6,000 K to 20,000 K.
        let a_1 = 1.779004264e+08;
        let a_2 = -1.082328257e+05;
        let a_3 = 2.810778365e+01;
        let a_4 = -2.975232262e-03;
        let a_5 = 1.854997534e-07;
        let a_6 = -5.796231540e-12;
        let a_7 = 7.191720164e-17;
        let b_1 = 8.890942630e+05;
        let b_2 = -2.181728151e+02;
        let o_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 6_000.0, 20_000.0,
        );

        o_tdp.with_coefficients(o_pc)
    }

    pub fn with_carbon_monoxide_coefficients(self) -> Self {
        // CO: From 200 K to 1,000 K.
        let a_1 = 1.489045326e+04;
        let a_2 = -2.922285939e+02;
        let a_3 = 5.724527170e+00;
        let a_4 = -8.176235030e-03;
        let a_5 = 1.456903469e-05;
        let a_6 = -1.087746302e-08;
        let a_7 = 3.027941827e-12;
        let b_1 = -1.303131878e+04;
        let b_2 = -7.859241350e+00;
        let co_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 190.0, 1_000.0,
        );
        let co_tdp = self.with_coefficients(co_pc);

        // CO: From 1,000 K to 6,000 K.
        let a_1 = 4.619197250e+05;
        let a_2 = -1.944704863e+03;
        let a_3 = 5.916714180e+00;
        let a_4 = -5.664282830e-04;
        let a_5 = 1.398814540e-07;
        let a_6 = -1.787680361e-11;
        let a_7 = 9.620935570e-16;
        let b_1 = -2.466261084e+03;
        let b_2 = -1.387413108e+01;
        let co_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 1_000.0, 6_000.0,
        );
        let co_tdp = co_tdp.with_coefficients(co_pc);

        // CO: From 6,000 K to 20,000 K.
        let a_1 = 8.868662960e+08;
        let a_2 = -7.500377840e+05;
        let a_3 = 2.495474979e+02;
        let a_4 = -3.956351100e-02;
        let a_5 = 3.297772080e-06;
        let a_6 = -1.318409933e-10;
        let a_7 = 1.998937948e-15;
        let b_1 = 5.701421130e+06;
        let b_2 = -2.060704786e+03;
        let co_pc = PolynomialCoefficients::new(
            a_1, a_2, a_3, a_4, a_5, a_6, a_7, b_1, b_2, 6_000.0, 20_000.0,
        );
        co_tdp.with_coefficients(co_pc)
    }

    fn find_coefficients(&self, temperature_k: f64) -> Option<&PolynomialCoefficients> {
        self.coefficients.iter().find(|&coefficients| {
            coefficients.min_temperture <= temperature_k
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
