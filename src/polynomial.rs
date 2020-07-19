pub struct Knot {
    pub x: f64,
    pub y: f64,
}

impl Knot {
    pub fn new(x: f64, y: f64) -> Knot {
        Knot { x, y }
    }
}

pub trait Evaluate {
    fn evaluate(&self, v: f64) -> f64;
}

pub trait HasDerivative {
    type DerivativeOf;
    fn derivative(&self) -> Self::DerivativeOf;
}

pub trait HasIntegral: Evaluate {
    type IntegralOf;
    fn indefinite(&self) -> Self::IntegralOf;
    fn integral(&self, knot: Knot) -> Self::IntegralOf;
}

/// Vertically translate the polynomial.
pub trait Translate {
    fn translate(&self, v: f64) -> Self;
}

#[inline]
pub fn evaluate_coeffs(coeffs: &[f64], v: f64) -> f64 {
    // Multiplying degree via `powi` seems mildly better than what we
    // did in original Haskell code which was to multiply the `v` into
    // the accumulator as we went: this is because this allows us to
    // do the instructions in parallel without having to wait before
    // processing each elemnt. I think.
    coeffs
        .iter()
        .enumerate()
        .fold(0.0, |acc, (i, e)| acc + e * v.powi(i as i32))
}

pub struct Poly0(pub f64);
impl Evaluate for Poly0 {
    fn evaluate(&self, _: f64) -> f64 {
        self.0
    }
}

pub struct Poly1(pub [f64; 2]);
impl Evaluate for Poly1 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}
pub struct Poly2(pub [f64; 3]);
impl Evaluate for Poly2 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}
pub struct Poly3(pub [f64; 4]);
impl Evaluate for Poly3 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}
pub struct Poly4(pub [f64; 5]);
impl Evaluate for Poly4 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}
pub struct Poly5(pub [f64; 6]);
impl Evaluate for Poly5 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}
pub struct Poly6(pub [f64; 7]);
impl Evaluate for Poly6 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}
pub struct Poly7(pub [f64; 8]);
impl Evaluate for Poly7 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}
pub struct Poly8(pub [f64; 9]);
impl Evaluate for Poly8 {
    fn evaluate(&self, v: f64) -> f64 {
        evaluate_coeffs(&self.0, v)
    }
}

#[cfg(test)]
mod tests {
    use crate::polynomial::*;
    #[test]
    fn evaluate_poly0() {
        let poly = Poly0(7.0);
        let v = 17.0;
        assert_eq!(poly.evaluate(v), 7.0);
    }

    #[test]
    fn evaluate_poly1() {
        let poly = Poly1([7.0, 3.0]);
        let v = 17.0;
        assert_eq!(poly.evaluate(v), 58.0);
    }

    #[test]
    fn evaluate_poly2() {
        let poly = Poly2([7.0, 3.0, 9.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 97.0);
    }

    #[test]
    fn evaluate_poly3() {
        let poly = Poly3([7.0, 3.0, 9.0, 8.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 313.0);
    }

    #[test]
    fn evaluate_poly4() {
        let poly = Poly4([7.0, 3.0, 9.0, 8.0, 6.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 799.0);
    }

    #[test]
    fn evaluate_poly5() {
        let poly = Poly5([7.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 1163.5);
    }

    #[test]
    fn evaluate_poly6() {
        let poly = Poly6([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 3715.0);
    }

    #[test]
    fn evaluate_poly7() {
        let poly = Poly7([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 13556.5);
    }

    #[test]
    fn evaluate_poly8() {
        let poly = Poly8([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 72605.5);
    }
}
