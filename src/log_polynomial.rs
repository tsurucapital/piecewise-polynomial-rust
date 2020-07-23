use crate::polynomial::*;
use std::ops::{Add, Mul, Neg};

/// Polynomial of natural log of x
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Log<T>(pub T);

impl<Scalar, T: Mul<Scalar>> Mul<Scalar> for Log<T> {
    type Output = Log<<T as Mul<Scalar>>::Output>;
    fn mul(self, rhs: Scalar) -> Self::Output {
        let r: <T as Mul<Scalar>>::Output = self.0 * rhs;
        Log(r)
    }
}

impl<T: Evaluate> Evaluate for Log<T> {
    fn evaluate(&self, v: f64) -> f64 {
        self.0.evaluate(v.ln())
    }
}

impl<T: Translate> Translate for Log<T> {
    fn translate(&self, v: f64) -> Self {
        Log(self.0.translate(v))
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct IntOfLog<T> {
    pub k: f64,
    pub poly: T,
}

impl<T: Default> Default for IntOfLog<T> {
    fn default() -> Self {
        IntOfLog {
            k: 0.0,
            poly: Default::default(),
        }
    }
}

impl<T: Add> Add for IntOfLog<T> {
    type Output = IntOfLog<<T as Add>::Output>;
    fn add(self, other: Self) -> Self::Output {
        IntOfLog {
            k: self.k + other.k,
            poly: self.poly + other.poly,
        }
    }
}

impl<Scalar, T> Mul<Scalar> for IntOfLog<T>
where
    // Allow any scalar as long as it can copy both the f64 and the
    // underlying polynomial. Also we use it twice so it needs to be
    // Copy.
    Scalar: Mul<f64, Output = f64> + Copy,
    T: Mul<Scalar>,
{
    type Output = IntOfLog<<T as Mul<Scalar>>::Output>;
    fn mul(self, rhs: Scalar) -> Self::Output {
        IntOfLog {
            k: rhs * self.k,
            poly: self.poly * rhs,
        }
    }
}

impl<T: Neg> Neg for IntOfLog<T> {
    type Output = IntOfLog<<T as Neg>::Output>;
    fn neg(self) -> Self::Output {
        IntOfLog {
            k: self.k.neg(),
            poly: self.poly.neg(),
        }
    }
}

impl<T: Evaluate> Evaluate for IntOfLog<T> {
    fn evaluate(&self, v: f64) -> f64 {
        self.k + self.poly.evaluate(v.ln())
    }
}

impl<T: Translate + Copy> Translate for IntOfLog<T> {
    fn translate(&self, v: f64) -> Self {
        IntOfLog {
            k: self.k + v,
            poly: self.poly,
        }
    }
}

#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct IntOfLogPoly4 {
    /// Constant term
    pub k: f64,
    pub coeffs: [f64; 4],
    pub u: f64,
}

impl Add for IntOfLogPoly4 {
    type Output = Self;
    fn add(self, other: Self) -> Self::Output {
        IntOfLogPoly4 {
            k: self.k + other.k,
            coeffs: [
                self.coeffs[0] + other.coeffs[0],
                self.coeffs[1] + other.coeffs[1],
                self.coeffs[2] + other.coeffs[2],
                self.coeffs[3] + other.coeffs[3],
            ],
            u: self.u + other.u,
        }
    }
}

impl Neg for IntOfLogPoly4 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        IntOfLogPoly4 {
            k: self.k.neg(),
            coeffs: [
                self.coeffs[0].neg(),
                self.coeffs[1].neg(),
                self.coeffs[2].neg(),
                self.coeffs[3].neg(),
            ],
            u: self.u.neg(),
        }
    }
}

impl Mul<f64> for IntOfLogPoly4 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        IntOfLogPoly4 {
            k: self.k * rhs,
            coeffs: [
                self.coeffs[0] * rhs,
                self.coeffs[1] * rhs,
                self.coeffs[2] * rhs,
                self.coeffs[3] * rhs,
            ],
            u: self.u * rhs,
        }
    }
}

// Hide ugly taylor expansion stuff
mod taylor {
    fn exp_5_tail_anal(x: f64) -> f64 {
        let x = x.recip();

        let c0: f64 = 0.0;
        let c1: f64 = -1.0 / 24.0;
        let c2: f64 = -1.0 / 6.0;
        let c3: f64 = -1.0 / 2.0;
        let c4: f64 = -1.0;
        let c5: f64 = x.recip().exp() - 1.0;

        let x2 = x * x;
        let x4 = x2 * x2;

        let t0 = c1.mul_add(x, c0);
        let t1 = c3.mul_add(x, c2);
        let t2 = c5.mul_add(x, c4);

        t2.mul_add(x4, t1.mul_add(x2, t0))
    }

    fn exp_5_tail_taylor(x: f64) -> f64 {
        // (((C0+C1x) + (C2+C3x)x2) + ((C4+C5x) + (C6+C7x)x2)x4) + (((C8+C9x) + (C10+C11x)x2) + ((C12+C13x) + (C14+C15x)x2)x4)x8
        let c0: f64 = 1.0 / 120.0;
        let c1: f64 = 1.0 / 720.0;
        let c2: f64 = 1.0 / 5040.0;
        let c3: f64 = 1.0 / 40320.0;
        let c4: f64 = 1.0 / 362880.0;
        let c5: f64 = 1.0 / 3628800.0;
        let c6: f64 = 1.0 / 39916800.0;
        let c7: f64 = 1.0 / 479001600.0;
        let c8: f64 = 1.0 / 6227020800.0;
        let c9: f64 = 1.0 / 87178291200.0;
        let c10: f64 = 1.0 / 1307674368000.0;
        let c11: f64 = 1.0 / 20922789888000.0;
        let c12: f64 = 1.0 / 355687428096000.0;
        let c13: f64 = 1.0 / 6402373705728000.0;
        let c14: f64 = 1.0 / 121645100408832000.0;
        let c15: f64 = 1.0 / 2432902008176640000.0;

        // Start with
        //
        //   (((C0+C1x) + (C2+C3x)x2) + ((C4+C5x) + (C6+C7x)x2)x4)
        // + (((C8+C9x) + (C10+C11x)x2) + ((C12+C13x) + (C14+C15x)x2)x4)x8

        // Smallest individual terms
        let t0 = c1.mul_add(x, c0);
        let t1 = c3.mul_add(x, c2);
        let t2 = c5.mul_add(x, c4);
        let t3 = c7.mul_add(x, c6);
        let t4 = c9.mul_add(x, c8);
        let t5 = c11.mul_add(x, c10);
        let t6 = c13.mul_add(x, c12);
        let t7 = c15.mul_add(x, c14);

        // We now have
        //
        //   ((t0 + t1x2) + (t2 + t3x2)x4)
        // + ((t4 + t5x2) + (t6 + t7x2)x4)x8

        let x2 = x * x;
        // Reduce top and bottom row terms
        let top_l = t1.mul_add(x2, t0);
        let top_r = t3.mul_add(x2, t2);
        let bot_l = t5.mul_add(x2, t4);
        let bot_r = t7.mul_add(x2, t6);

        // We now have
        //
        //   (top_l + top_rx4)
        // + (bot_l + bot_rx4)x8

        let x4 = x2 * x2;
        // Reduce the top and bottom completely
        let top = top_r.mul_add(x4, top_l);
        let bot = bot_r.mul_add(x4, bot_l);

        // We know have
        //
        // top + botx8
        let x8 = x4 * x4;
        // Final reduce, the answer.
        bot.mul_add(x8, top)
    }

    pub fn exp_5_taylor(x: f64) -> f64 {
        const LOWER_THRES: f64 = -1.71;
        const UPPER_THRES: f64 = 1.72;
        if LOWER_THRES < x && x < UPPER_THRES {
            exp_5_tail_taylor(x)
        } else {
            exp_5_tail_anal(x)
        }
    }

    #[test]
    fn exp_tail_tailor_ref() {
        use assert_approx_eq::assert_approx_eq;
        // Check some sample value against reference implementation to
        // make sure we haven't messed it up when copying down. Allow
        // 1e-16 of floating point shenanighans: this is just close to
        // the limit of the test passing: if the test fails for you
        // due to this delta being too small, feel free to bump it a
        // little.
        let diff = 1e-16;
        assert_approx_eq!(exp_5_taylor(-2.0), 6.187439065522517e-3, diff);
        assert_approx_eq!(exp_5_taylor(-1.5), 6.624834877573014e-3, diff);
        assert_approx_eq!(exp_5_taylor(-1.0), 7.120558828557678e-3, diff);
        assert_approx_eq!(exp_5_taylor(-0.5), 7.685555862397111e-3, diff);
        assert_approx_eq!(exp_5_taylor(0.0), 8.333333333333333e-3, diff);
        assert_approx_eq!(exp_5_taylor(0.5), 9.080662404100699e-3, diff);
        assert_approx_eq!(exp_5_taylor(1.0), 9.948495125711903e-3, diff);
        assert_approx_eq!(exp_5_taylor(1.5), 1.0963169756452968e-2, diff);
        assert_approx_eq!(exp_5_taylor(2.0), 1.2158003091582829e-2, diff);
    }
}

impl Evaluate for IntOfLogPoly4 {
    fn evaluate(&self, v: f64) -> f64 {
        let x = v.ln().neg();

        let cr = {
            let c0: f64 = 0.0;
            let c1: f64 = self.coeffs[0];
            let c2: f64 = self.coeffs[1];
            let c3: f64 = self.coeffs[2];
            let c4: f64 = self.coeffs[3];
            let c5: f64 = self.u * crate::log_polynomial::taylor::exp_5_taylor(x);

            let t0 = c1.mul_add(x, c0);
            let t1 = c3.mul_add(x, c2);
            let t2 = c5.mul_add(x, c4);

            let x2 = x * x;
            let x4 = x2 * x2;

            t2.mul_add(x4, t1.mul_add(x2, t0))
        };

        v.mul_add(cr, self.k)
    }
}

impl Translate for IntOfLogPoly4 {
    fn translate(&self, v: f64) -> Self {
        IntOfLogPoly4 {
            k: self.k + v,
            coeffs: self.coeffs,
            u: self.u,
        }
    }
}

impl HasIntegral for Log<Poly0> {
    type IntegralOf = IntOfLog<Poly0>;
    fn indefinite(&self) -> Self::IntegralOf {
        IntOfLog {
            k: 0.0,
            poly: Poly0((self.0).0),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl HasIntegral for Log<Poly1> {
    type IntegralOf = IntOfLog<Poly1>;
    fn indefinite(&self) -> Self::IntegralOf {
        let b = (self.0).0[1];
        let a = (self.0).0[0] - b;
        IntOfLog {
            k: 0.0,
            poly: Poly1([a, b]),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl HasIntegral for Log<Poly2> {
    type IntegralOf = IntOfLog<Poly2>;
    fn indefinite(&self) -> Self::IntegralOf {
        let c = (self.0).0[2];
        let b = (self.0).0[1] - 2.0 * c;
        let a = (self.0).0[0] - b;
        IntOfLog {
            k: 0.0,
            poly: Poly2([a, b, c]),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl HasIntegral for Log<Poly3> {
    type IntegralOf = IntOfLog<Poly3>;
    fn indefinite(&self) -> Self::IntegralOf {
        let d = (self.0).0[3];
        let c = (self.0).0[2] - 3.0 * d;
        let b = (self.0).0[1] - 2.0 * c;
        let a = (self.0).0[0] - b;
        IntOfLog {
            k: 0.0,
            poly: Poly3([a, b, c, d]),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

// NB: special case
impl HasIntegral for Log<Poly4> {
    type IntegralOf = IntOfLogPoly4;
    fn indefinite(&self) -> Self::IntegralOf {
        let a = (self.0).0[0].neg();
        let b = (a + (self.0).0[1]) * (1.0 / 2.0);
        let c = (b - (self.0).0[2]) * (1.0 / 3.0);
        let d = (c + (self.0).0[3]) * (1.0 / 4.0);
        let u = (d - (self.0).0[4]) * 24.0;
        IntOfLogPoly4 {
            k: 0.0,
            coeffs: [a, b, c, d],
            u,
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl HasIntegral for Log<Poly5> {
    type IntegralOf = IntOfLog<Poly5>;
    fn indefinite(&self) -> Self::IntegralOf {
        let f = (self.0).0[5];
        let e = (self.0).0[4] - 5.0 * f;
        let d = (self.0).0[3] - 4.0 * e;
        let c = (self.0).0[2] - 3.0 * d;
        let b = (self.0).0[1] - 2.0 * c;
        let a = (self.0).0[0] - b;
        IntOfLog {
            k: 0.0,
            poly: Poly5([a, b, c, d, e, f]),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl HasIntegral for Log<Poly6> {
    type IntegralOf = IntOfLog<Poly6>;
    fn indefinite(&self) -> Self::IntegralOf {
        let g = (self.0).0[6];
        let f = (self.0).0[5] - 6.0 * g;
        let e = (self.0).0[4] - 5.0 * f;
        let d = (self.0).0[3] - 4.0 * e;
        let c = (self.0).0[2] - 3.0 * d;
        let b = (self.0).0[1] - 2.0 * c;
        let a = (self.0).0[0] - b;
        IntOfLog {
            k: 0.0,
            poly: Poly6([a, b, c, d, e, f, g]),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl HasIntegral for Log<Poly7> {
    type IntegralOf = IntOfLog<Poly7>;
    fn indefinite(&self) -> Self::IntegralOf {
        let h = (self.0).0[7];
        let g = (self.0).0[6] - 7.0 * h;
        let f = (self.0).0[5] - 6.0 * g;
        let e = (self.0).0[4] - 5.0 * f;
        let d = (self.0).0[3] - 4.0 * e;
        let c = (self.0).0[2] - 3.0 * d;
        let b = (self.0).0[1] - 2.0 * c;
        let a = (self.0).0[0] - b;
        IntOfLog {
            k: 0.0,
            poly: Poly7([a, b, c, d, e, f, g, h]),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl HasIntegral for Log<Poly8> {
    type IntegralOf = IntOfLog<Poly8>;
    fn indefinite(&self) -> Self::IntegralOf {
        let i = (self.0).0[8];
        let h = (self.0).0[7] - 8.0 * i;
        let g = (self.0).0[6] - 7.0 * h;
        let f = (self.0).0[5] - 6.0 * g;
        let e = (self.0).0[4] - 5.0 * f;
        let d = (self.0).0[3] - 4.0 * e;
        let c = (self.0).0[2] - 3.0 * d;
        let b = (self.0).0[1] - 2.0 * c;
        let a = (self.0).0[0] - b;
        IntOfLog {
            k: 0.0,
            poly: Poly8([a, b, c, d, e, f, g, h, i]),
        }
    }

    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn evaluate_log_poly0() {
        use crate::log_polynomial::*;
        let poly = Log(Poly0(2.0));
        assert_eq!(poly.evaluate(7.0), 2.0);
    }

    #[test]
    fn integral_log_poly0() {
        use crate::log_polynomial::*;
        let poly = Log(Poly0(2.0));
        let result = IntOfLog {
            k: 0.0,
            poly: Poly0(2.0),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly1() {
        use crate::log_polynomial::*;
        let poly = Log(Poly1([2.0, 3.0]));
        assert_eq!(poly.evaluate(7.0), 7.8377304471659395);
    }

    #[test]
    fn integral_log_poly1() {
        use crate::log_polynomial::*;
        let poly = Log(Poly1([2.0, 3.0]));
        let result = IntOfLog {
            k: 3.0,
            poly: Poly1([-1.0, 3.0]),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly2() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Log(Poly2([2.0, 3.0, 4.0]));
        assert_approx_eq!(poly.evaluate(7.0), 22.98399567995183, 1e-12);
    }

    #[test]
    fn integral_log_poly2() {
        use crate::log_polynomial::*;
        let poly = Log(Poly2([2.0, 3.0, 4.0]));
        let result = IntOfLog {
            k: -5.0,
            poly: Poly2([7.0, -5.0, 4.0]),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly3() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Log(Poly3([2.0, 3.0, 4.0, 5.0]));
        assert_approx_eq!(poly.evaluate(7.0), 59.82558472590394, 1e-12);
    }

    #[test]
    fn integral_log_poly3() {
        use crate::log_polynomial::*;
        let poly = Log(Poly3([2.0, 3.0, 4.0, 5.0]));
        let result = IntOfLog {
            k: 25.0,
            poly: Poly3([-23.0, 25.0, -11.0, 5.0]),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly4() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Log(Poly4([2.0, 3.0, 4.0, 5.0, 6.0]));
        assert_approx_eq!(poly.evaluate(7.0), 145.85409116411583, 1e-12);
    }

    #[test]
    fn evaluate_int_of_log_poly_4() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = IntOfLogPoly4 {
            k: 1.0,
            coeffs: [2.0, 3.0, 4.0, 5.0],
            u: 6.0,
        };
        assert_approx_eq!(poly.evaluate(7.0), 341.4921166923076, 1e-12);
    }

    #[test]
    fn integral_log_poly4() {
        use crate::log_polynomial::*;
        let poly = Log(Poly4([2.0, 3.0, 4.0, 5.0, 6.0]));
        let result = IntOfLogPoly4 {
            k: 2.0,
            coeffs: [-2.0, 0.5, -1.1666666666666665, 0.9583333333333334],
            u: -121.0,
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly5() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Log(Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0]));
        assert_approx_eq!(poly.evaluate(7.0), 341.1584589146673, 1e-12);
    }

    #[test]
    fn integral_log_poly5() {
        use crate::log_polynomial::*;
        let poly = Log(Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0]));
        let result = IntOfLog {
            k: 721.0,
            poly: Poly5([-719.0, 721.0, -359.0, 121.0, -29.0, 7.0]),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly6() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Log(Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]));
        assert_approx_eq!(poly.evaluate(7.0), 775.4953176125293, 1e-12);
    }

    #[test]
    fn integral_log_poly6() {
        use crate::log_polynomial::*;
        let poly = Log(Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]));
        let result = IntOfLog {
            k: -5039.0,
            poly: Poly6([5041.0, -5039.0, 2521.0, -839.0, 211.0, -41.0, 8.0]),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly7() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Log(Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]));
        assert_approx_eq!(poly.evaluate(7.0), 1726.3233817426242, 1e-12);
    }

    #[test]
    fn integral_log_poly7() {
        use crate::log_polynomial::*;
        let poly = Log(Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]));
        let result = IntOfLog {
            k: 40321.0,
            poly: Poly7([
                -40319.0, 40321.0, -20159.0, 6721.0, -1679.0, 337.0, -55.0, 9.0,
            ]),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_log_poly8() {
        use crate::log_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Log(Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]));
        assert_approx_eq!(poly.evaluate(7.0), 3782.130026184144, 1e-12);
    }

    #[test]
    fn integral_log_poly8() {
        use crate::log_polynomial::*;
        let poly = Log(Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]));
        let result = IntOfLog {
            k: -362879.0,
            poly: Poly8([
                362881.0, -362879.0, 181441.0, -60479.0, 15121.0, -3023.0, 505.0, -71.0, 10.0,
            ]),
        };
        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }
}
