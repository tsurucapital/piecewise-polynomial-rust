use std::default::Default;
use std::ops::{Add, Mul, Neg};

#[derive(Debug, PartialEq, Clone, Copy)]
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

pub trait HasIntegral
where
    Self::IntegralOf: Evaluate,
{
    type IntegralOf;
    fn indefinite(&self) -> Self::IntegralOf;
    fn integral(&self, knot: Knot) -> Self::IntegralOf;
}

/// Vertically translate the polynomial.
pub trait Translate {
    fn translate(&self, v: f64) -> Self;
}

#[inline]
/// For our polynomials, we use manually unrolled Estin's scheme
/// (manually unrolled simply because it's easier than making general
/// implementation that's going to inline well). However, if you have
/// some collection of coefficients of your own, I provide FMA-enabled
/// Horner's scheme which performs fairly closely.
pub fn evaluate_horners<'a, I>(coeffs: I, v: f64) -> f64
where
    I: IntoIterator<Item = &'a f64>,
    <I as IntoIterator>::IntoIter: DoubleEndedIterator,
{
    let mut iter = coeffs.into_iter().rev();
    match iter.next() {
        None => return 0.0,
        Some(&first) => iter.fold(first, |acc, &e| acc.mul_add(v, e)),
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly0(pub f64);
impl Evaluate for Poly0 {
    fn evaluate(&self, _: f64) -> f64 {
        self.0
    }
}
impl HasDerivative for Poly0 {
    type DerivativeOf = Poly0;
    fn derivative(&self) -> Self::DerivativeOf {
        Poly0(0.0)
    }
}
impl Translate for Poly0 {
    fn translate(&self, v: f64) -> Self {
        Poly0(self.0 + v)
    }
}
impl HasIntegral for Poly0 {
    type IntegralOf = Poly1;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [0.0, self.0];
        Poly1(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl Mul<f64> for Poly0 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly0(self.0 * rhs)
    }
}
impl Neg for Poly0 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly0 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly0(self.0 + other.0)
    }
}

impl Default for Poly0 {
    fn default() -> Self {
        Poly0(0.0)
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly1(pub [f64; 2]);
impl Evaluate for Poly1 {
    fn evaluate(&self, x: f64) -> f64 {
        let c = self.0;
        c[1].mul_add(x, c[0])
    }
}
impl HasDerivative for Poly1 {
    type DerivativeOf = Poly0;
    fn derivative(&self) -> Self::DerivativeOf {
        Poly0(self.0[1])
    }
}
impl Translate for Poly1 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly1(dst)
    }
}
impl HasIntegral for Poly1 {
    type IntegralOf = Poly2;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [0.0, self.0[0], self.0[1] / 2.0];
        Poly2(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}
impl Mul<f64> for Poly1 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly1([self.0[0] * rhs, self.0[1] * rhs])
    }
}
impl Neg for Poly1 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly1 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly1([self.0[0] + other.0[0], self.0[1] + other.0[1]])
    }
}
impl Default for Poly1 {
    fn default() -> Self {
        Poly1([0.0; 2])
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly2(pub [f64; 3]);
impl Evaluate for Poly2 {
    fn evaluate(&self, x: f64) -> f64 {
        let x2 = x * x;
        let c = self.0;
        c[2].mul_add(x2, c[1].mul_add(x, c[0]))
    }
}
impl HasDerivative for Poly2 {
    type DerivativeOf = Poly1;
    fn derivative(&self) -> Self::DerivativeOf {
        let coeffs = self.0;
        let dst = [coeffs[1], 2.0 * coeffs[2]];
        Poly1(dst)
    }
}
impl Translate for Poly2 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly2(dst)
    }
}
impl HasIntegral for Poly2 {
    type IntegralOf = Poly3;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [0.0, self.0[0], self.0[1] / 2.0, self.0[2] / 3.0];
        Poly3(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}
impl Mul<f64> for Poly2 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly2([self.0[0] * rhs, self.0[1] * rhs, self.0[2] * rhs])
    }
}
impl Neg for Poly2 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly2 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly2([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
        ])
    }
}
impl Default for Poly2 {
    fn default() -> Self {
        Poly2([0.0; 3])
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly3(pub [f64; 4]);
impl Evaluate for Poly3 {
    fn evaluate(&self, x: f64) -> f64 {
        // P3(x) = (C0 + C1x) + (C2 + C3x) x2
        let x2 = x * x;
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);

        t1.mul_add(x2, t0)
    }
}
impl HasDerivative for Poly3 {
    type DerivativeOf = Poly2;
    fn derivative(&self) -> Self::DerivativeOf {
        let coeffs = self.0;
        let dst = [coeffs[1], 2.0 * coeffs[2], 3.0 * coeffs[3]];
        Poly2(dst)
    }
}
impl Translate for Poly3 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly3(dst)
    }
}
impl HasIntegral for Poly3 {
    type IntegralOf = Poly4;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [
            0.0,
            self.0[0],
            self.0[1] / 2.0,
            self.0[2] / 3.0,
            self.0[3] / 4.0,
        ];
        Poly4(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}
impl Mul<f64> for Poly3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly3([
            self.0[0] * rhs,
            self.0[1] * rhs,
            self.0[2] * rhs,
            self.0[3] * rhs,
        ])
    }
}
impl Neg for Poly3 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly3([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
            self.0[3] + other.0[3],
        ])
    }
}
impl Default for Poly3 {
    fn default() -> Self {
        Poly3([0.0; 4])
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly4(pub [f64; 5]);
impl Evaluate for Poly4 {
    fn evaluate(&self, x: f64) -> f64 {
        // P4(x) = (C0 + C1x) + (C2 + C3x) x2 + C4x4
        let x2 = x * x;
        let x4 = x2 * x2;
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = c[4];

        t2.mul_add(x4, t1.mul_add(x2, t0))
    }
}
impl HasDerivative for Poly4 {
    type DerivativeOf = Poly3;
    fn derivative(&self) -> Self::DerivativeOf {
        let coeffs = self.0;
        let dst = [coeffs[1], 2.0 * coeffs[2], 3.0 * coeffs[3], 4.0 * coeffs[4]];
        Poly3(dst)
    }
}
impl Translate for Poly4 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly4(dst)
    }
}
impl HasIntegral for Poly4 {
    type IntegralOf = Poly5;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [
            0.0,
            self.0[0],
            self.0[1] / 2.0,
            self.0[2] / 3.0,
            self.0[3] / 4.0,
            self.0[4] / 5.0,
        ];
        Poly5(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}
impl Mul<f64> for Poly4 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly4([
            self.0[0] * rhs,
            self.0[1] * rhs,
            self.0[2] * rhs,
            self.0[3] * rhs,
            self.0[4] * rhs,
        ])
    }
}
impl Neg for Poly4 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly4 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly4([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
            self.0[3] + other.0[3],
            self.0[4] + other.0[4],
        ])
    }
}
impl Default for Poly4 {
    fn default() -> Self {
        Poly4([0.0; 5])
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly5(pub [f64; 6]);
impl Evaluate for Poly5 {
    fn evaluate(&self, x: f64) -> f64 {
        // P5(x) = (C0 + C1x) + (C2 + C3x) x2 + (C4 + C5x) x4
        let x2 = x * x;
        let x4 = x2 * x2;
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = c[5].mul_add(x, c[4]);

        t2.mul_add(x4, t1.mul_add(x2, t0))
    }
}
impl HasDerivative for Poly5 {
    type DerivativeOf = Poly4;
    fn derivative(&self) -> Self::DerivativeOf {
        let coeffs = self.0;
        let dst = [
            coeffs[1],
            2.0 * coeffs[2],
            3.0 * coeffs[3],
            4.0 * coeffs[4],
            5.0 * coeffs[5],
        ];
        Poly4(dst)
    }
}
impl Translate for Poly5 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly5(dst)
    }
}
impl HasIntegral for Poly5 {
    type IntegralOf = Poly6;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [
            0.0,
            self.0[0],
            self.0[1] / 2.0,
            self.0[2] / 3.0,
            self.0[3] / 4.0,
            self.0[4] / 5.0,
            self.0[5] / 6.0,
        ];
        Poly6(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}
impl Mul<f64> for Poly5 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly5([
            self.0[0] * rhs,
            self.0[1] * rhs,
            self.0[2] * rhs,
            self.0[3] * rhs,
            self.0[4] * rhs,
            self.0[5] * rhs,
        ])
    }
}
impl Neg for Poly5 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly5 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly5([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
            self.0[3] + other.0[3],
            self.0[4] + other.0[4],
            self.0[5] + other.0[5],
        ])
    }
}
impl Default for Poly5 {
    fn default() -> Self {
        Poly5([0.0; 6])
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly6(pub [f64; 7]);
impl Evaluate for Poly6 {
    fn evaluate(&self, x: f64) -> f64 {
        // P6(x) = (C0 + C1x) + (C2 + C3x) x2 + ((C4 + C5x) + C6x2)x4
        let x2 = x * x;
        let x4 = x2 * x2;
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = c[6].mul_add(x2, c[5].mul_add(x, c[4]));

        t2.mul_add(x4, t1.mul_add(x2, t0))
    }
}
impl HasDerivative for Poly6 {
    type DerivativeOf = Poly5;
    fn derivative(&self) -> Self::DerivativeOf {
        let coeffs = self.0;
        let dst = [
            coeffs[1],
            2.0 * coeffs[2],
            3.0 * coeffs[3],
            4.0 * coeffs[4],
            5.0 * coeffs[5],
            6.0 * coeffs[6],
        ];
        Poly5(dst)
    }
}
impl Translate for Poly6 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly6(dst)
    }
}
impl HasIntegral for Poly6 {
    type IntegralOf = Poly7;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [
            0.0,
            self.0[0],
            self.0[1] / 2.0,
            self.0[2] / 3.0,
            self.0[3] / 4.0,
            self.0[4] / 5.0,
            self.0[5] / 6.0,
            self.0[6] / 7.0,
        ];
        Poly7(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}
impl Mul<f64> for Poly6 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly6([
            self.0[0] * rhs,
            self.0[1] * rhs,
            self.0[2] * rhs,
            self.0[3] * rhs,
            self.0[4] * rhs,
            self.0[5] * rhs,
            self.0[6] * rhs,
        ])
    }
}
impl Neg for Poly6 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly6 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly6([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
            self.0[3] + other.0[3],
            self.0[4] + other.0[4],
            self.0[5] + other.0[5],
            self.0[6] + other.0[6],
        ])
    }
}
impl Default for Poly6 {
    fn default() -> Self {
        Poly6([0.0; 7])
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly7(pub [f64; 8]);
impl Evaluate for Poly7 {
    fn evaluate(&self, x: f64) -> f64 {
        // P7(x) = (C0 + C1x) + (C2 + C3x) x2 + ((C4 + C5x) + (C6 + C7x) x2)x4
        let x2 = x * x;
        let x4 = x2 * x2;
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = (c[7].mul_add(x, c[6])).mul_add(x2, c[5].mul_add(x, c[4]));

        t2.mul_add(x4, t1.mul_add(x2, t0))
    }
}
impl HasDerivative for Poly7 {
    type DerivativeOf = Poly6;
    fn derivative(&self) -> Self::DerivativeOf {
        let coeffs = self.0;
        let dst = [
            coeffs[1],
            2.0 * coeffs[2],
            3.0 * coeffs[3],
            4.0 * coeffs[4],
            5.0 * coeffs[5],
            6.0 * coeffs[6],
            7.0 * coeffs[7],
        ];
        Poly6(dst)
    }
}
impl Translate for Poly7 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly7(dst)
    }
}
impl HasIntegral for Poly7 {
    type IntegralOf = Poly8;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [
            0.0,
            self.0[0],
            self.0[1] / 2.0,
            self.0[2] / 3.0,
            self.0[3] / 4.0,
            self.0[4] / 5.0,
            self.0[5] / 6.0,
            self.0[6] / 7.0,
            self.0[7] / 8.0,
        ];
        Poly8(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}
impl Mul<f64> for Poly7 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly7([
            self.0[0] * rhs,
            self.0[1] * rhs,
            self.0[2] * rhs,
            self.0[3] * rhs,
            self.0[4] * rhs,
            self.0[5] * rhs,
            self.0[6] * rhs,
            self.0[7] * rhs,
        ])
    }
}
impl Neg for Poly7 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly7 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly7([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
            self.0[3] + other.0[3],
            self.0[4] + other.0[4],
            self.0[5] + other.0[5],
            self.0[6] + other.0[6],
            self.0[7] + other.0[7],
        ])
    }
}
impl Default for Poly7 {
    fn default() -> Self {
        Poly7([0.0; 8])
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Poly8(pub [f64; 9]);
impl Evaluate for Poly8 {
    fn evaluate(&self, x: f64) -> f64 {
        // P8(x) = (C0 + C1x) + (C2 + C3x) x2 + ((C4 + C5x) + (C6 + C7x) x2)x4 + C8x8
        let x2 = x * x;
        let x4 = x2 * x2;
        let x8 = x4 * x4;
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = (c[7].mul_add(x, c[6])).mul_add(x2, c[5].mul_add(x, c[4]));
        let t3 = c[8];

        t3.mul_add(x8, t2.mul_add(x4, t1.mul_add(x2, t0)))
    }
}

impl HasDerivative for Poly8 {
    type DerivativeOf = Poly7;
    fn derivative(&self) -> Self::DerivativeOf {
        let coeffs = self.0;
        let dst = [
            coeffs[1],
            2.0 * coeffs[2],
            3.0 * coeffs[3],
            4.0 * coeffs[4],
            5.0 * coeffs[5],
            6.0 * coeffs[6],
            7.0 * coeffs[7],
            8.0 * coeffs[8],
        ];
        Poly7(dst)
    }
}
impl Translate for Poly8 {
    fn translate(&self, v: f64) -> Self {
        let mut dst = self.0;
        dst[0] += v;
        Poly8(dst)
    }
}
impl Mul<f64> for Poly8 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly8([
            self.0[0] * rhs,
            self.0[1] * rhs,
            self.0[2] * rhs,
            self.0[3] * rhs,
            self.0[4] * rhs,
            self.0[5] * rhs,
            self.0[6] * rhs,
            self.0[7] * rhs,
            self.0[8] * rhs,
        ])
    }
}
impl Neg for Poly8 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Add for Poly8 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Poly8([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
            self.0[3] + other.0[3],
            self.0[4] + other.0[4],
            self.0[5] + other.0[5],
            self.0[6] + other.0[6],
            self.0[7] + other.0[7],
            self.0[8] + other.0[8],
        ])
    }
}
impl Default for Poly8 {
    fn default() -> Self {
        Poly8([0.0; 9])
    }
}

#[cfg(test)]
mod tests {
    // Test results taken from the reference Haskell implementation.

    use crate::polynomial::*;
    #[test]
    fn evaluate_poly0() {
        let poly = Poly0(7.0);
        let v = 17.0;
        assert_eq!(poly.evaluate(v), 7.0);
    }

    #[test]
    fn derivative_poly0() {
        let poly = Poly0(7.0);
        let result = Poly0(0.0);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly0() {
        let poly = Poly0(7.0);
        let result = Poly0(14.0);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly0() {
        let poly = Poly0(7.0);
        let result = Poly1([-9.0, 7.0]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly1() {
        let poly = Poly1([7.0, 3.0]);
        let v = 17.0;
        assert_eq!(poly.evaluate(v), 58.0);
    }

    #[test]
    fn derivative_poly1() {
        let poly = Poly1([7.0, 3.0]);
        let result = Poly0(3.0);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly1() {
        let poly = Poly1([7.0, 3.0]);
        let result = Poly1([14.0, 3.0]);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly1() {
        let poly = Poly1([7.0, 3.0]);
        let result = Poly2([-15.0, 7.0, 1.5]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly2() {
        let poly = Poly2([7.0, 3.0, 9.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 97.0);
    }

    #[test]
    fn derivative_poly2() {
        let poly = Poly2([7.0, 3.0, 9.0]);
        let result = Poly1([3.0, 18.0]);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly2() {
        let poly = Poly2([7.0, 3.0, 9.0]);
        let result = Poly2([14.0, 3.0, 9.0]);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly2() {
        let poly = Poly2([7.0, 3.0, 9.0]);
        let result = Poly3([-39.0, 7.0, 1.5, 3.0]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly3() {
        let poly = Poly3([7.0, 3.0, 9.0, 8.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 313.0);
    }

    #[test]
    fn derivative_poly3() {
        let poly = Poly3([7.0, 3.0, 9.0, 8.0]);
        let result = Poly2([3.0, 18.0, 24.0]);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly3() {
        let poly = Poly3([7.0, 3.0, 9.0, 8.0]);
        let result = Poly3([14.0, 3.0, 9.0, 8.0]);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly3() {
        let poly = Poly3([7.0, 3.0, 9.0, 8.0]);
        let result = Poly4([-71.0, 7.0, 1.5, 3.0, 2.0]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly4() {
        let poly = Poly4([7.0, 3.0, 9.0, 8.0, 6.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 799.0);
    }

    #[test]
    fn derivative_poly4() {
        let poly = Poly4([7.0, 3.0, 9.0, 8.0, 6.0]);
        let result = Poly3([3.0, 18.0, 24.0, 24.0]);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly4() {
        let poly = Poly4([7.0, 3.0, 9.0, 8.0, 6.0]);
        let result = Poly4([14.0, 3.0, 9.0, 8.0, 6.0]);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly4() {
        let poly = Poly4([7.0, 3.0, 9.0, 8.0, 6.0]);
        let result = Poly5([-109.4, 7.0, 1.5, 3.0, 2.0, 1.2]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly5() {
        let poly = Poly5([7.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 1163.5);
    }

    #[test]
    fn derivative_poly5() {
        let poly = Poly5([7.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        let result = Poly4([3.0, 18.0, 24.0, 24.0, 7.5]);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly5() {
        let poly = Poly5([7.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        let result = Poly5([14.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly5() {
        let poly = Poly5([7.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        let result = Poly6([-125.4, 7.0, 1.5, 3.0, 2.0, 1.2, 0.25]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly6() {
        let poly = Poly6([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 3715.0);
    }

    #[test]
    fn derivative_poly6() {
        let poly = Poly6([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        let result = Poly5([3.0, 18.0, 24.0, 24.0, 7.5, 21.0]);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly6() {
        let poly = Poly6([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        let result = Poly6([14.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly6() {
        let poly = Poly6([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        let result = Poly7([-189.4, 7.0, 1.5, 3.0, 2.0, 1.2, 0.25, 0.5]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly7() {
        let poly = Poly7([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 13556.5);
    }

    #[test]
    fn derivative_poly7() {
        let poly = Poly7([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        let result = Poly6([3.0, 18.0, 24.0, 24.0, 7.5, 21.0, 31.5]);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly7() {
        let poly = Poly7([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        let result = Poly7([14.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        assert_eq!(poly.translate(7.0), result);
    }

    #[test]
    fn integral_poly7() {
        let poly = Poly7([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        let result = Poly8([-333.4, 7.0, 1.5, 3.0, 2.0, 1.2, 0.25, 0.5, 0.5625]);
        assert_eq!(poly.integral(Knot { x: 2.0, y: 5.0 }), result);
    }

    #[test]
    fn evaluate_poly8() {
        let poly = Poly8([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
        let v = 3.0;
        assert_eq!(poly.evaluate(v), 72605.5);
    }

    #[test]
    fn derivative_poly8() {
        let poly = Poly8([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
        let result = Poly7([3.0, 18.0, 24.0, 24.0, 7.5, 21.0, 31.5, 72.0]);
        assert_eq!(poly.derivative(), result);
    }

    #[test]
    fn translate_poly8() {
        let poly = Poly8([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
        let result = Poly8([14.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
        assert_eq!(poly.translate(7.0), result);
    }
}
