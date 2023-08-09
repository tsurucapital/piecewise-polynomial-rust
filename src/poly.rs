use approx::{AbsDiffEq, RelativeEq};
use arbitrary::Arbitrary;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul, MulAssign, Neg};

/// A join-point between two functions
#[derive(Debug, PartialEq, Clone, Copy, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
pub struct Knot {
    pub x: f64,
    pub y: f64,
}

impl Knot {
    pub fn new(x: f64, y: f64) -> Knot {
        Knot { x, y }
    }
}

/// Functions in ℝ → ℝ
pub trait Evaluate {
    fn evaluate(&self, v: f64) -> f64;
}

/// Differentiable functions
pub trait HasDerivative {
    type DerivativeOf;
    fn derivative(&self) -> Self::DerivativeOf;
}

/// Integrable functions
pub trait HasIntegral
where
    Self::IntegralOf: Evaluate,
{
    type IntegralOf;
    fn indefinite(&self) -> Self::IntegralOf;
    fn integral(&self, knot: Knot) -> Self::IntegralOf;
}

/// Functions which can be translated vertically
pub trait Translate {
    /// Vertically translate the polynomial
    fn translate(&mut self, v: f64);
}

/// Rank-n polynomial
///
/// Evaluated using FMA-enabled Horner's scheme.
#[derive(Debug, PartialEq, Clone, Arbitrary)]
pub struct PolyN(pub Vec<f64>);
impl Evaluate for PolyN {
    fn evaluate(&self, x: f64) -> f64 {
        let mut iter = self.0.iter().rev();
        match iter.next() {
            None => 0.0,
            Some(&first) => iter.fold(first, |acc, &e| acc.mul_add(x, e)),
        }
    }
}
impl Translate for PolyN {
    fn translate(&mut self, v: f64) {
        match self.0.get_mut(0) {
            None => self.0.push(v),
            Some(x0) => *x0 += v,
        }
    }
}

impl AbsDiffEq for PolyN {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for PolyN {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-0 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
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
    fn translate(&mut self, v: f64) {
        self.0 += v
    }
}
impl HasIntegral for Poly0 {
    type IntegralOf = Poly1;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [0.0, self.0];
        Poly1(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
    }
}

impl Mul<f64> for Poly0 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly0(self.0 * rhs)
    }
}

impl MulAssign<f64> for Poly0 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0 *= rhs;
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

impl AbsDiffEq for Poly0 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly0 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-1 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
    }
}
impl HasIntegral for Poly1 {
    type IntegralOf = Poly2;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [0.0, self.0[0], self.0[1] / 2.0];
        Poly2(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
    }
}
impl Mul<f64> for Poly1 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly1([self.0[0] * rhs, self.0[1] * rhs])
    }
}

impl MulAssign<f64> for Poly1 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly1 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly1 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-2 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
    }
}
impl HasIntegral for Poly2 {
    type IntegralOf = Poly3;
    fn indefinite(&self) -> Self::IntegralOf {
        let dst = [0.0, self.0[0], self.0[1] / 2.0, self.0[2] / 3.0];
        Poly3(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
    }
}
impl Mul<f64> for Poly2 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Poly2([self.0[0] * rhs, self.0[1] * rhs, self.0[2] * rhs])
    }
}

impl MulAssign<f64> for Poly2 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly2 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly2 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-3 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
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
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
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

impl MulAssign<f64> for Poly3 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly3 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly3 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-4 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
pub struct Poly4(pub [f64; 5]);
impl Evaluate for Poly4 {
    fn evaluate(&self, x: f64) -> f64 {
        // P4(x) = (C0 + C1x) + (C2 + C3x) x2 + C4x4
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = c[4];

        let x2 = x * x;
        let r = t1.mul_add(x2, t0);

        let x4 = x2 * x2;
        t2.mul_add(x4, r)
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
    }
}
impl HasIntegral for Poly4 {
    type IntegralOf = Poly5;
    fn indefinite(&self) -> Self::IntegralOf {
        let x0 = (wide::f64x4::new([self.0[1], self.0[2], self.0[3], self.0[4]])
            * wide::f64x4::new([1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0]))
        .to_array();
        let dst = [0.0, self.0[0], x0[0], x0[1], x0[2], x0[3]];
        Poly5(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
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

impl MulAssign<f64> for Poly4 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly4 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly4 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-5 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
pub struct Poly5(pub [f64; 6]);
impl Evaluate for Poly5 {
    fn evaluate(&self, x: f64) -> f64 {
        // P5(x) = (C0 + C1x) + (C2 + C3x) x2 + (C4 + C5x) x4
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = c[5].mul_add(x, c[4]);

        let x2 = x * x;
        let r = t1.mul_add(x2, t0);

        let x4 = x2 * x2;
        t2.mul_add(x4, r)
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
    }
}
impl HasIntegral for Poly5 {
    type IntegralOf = Poly6;
    fn indefinite(&self) -> Self::IntegralOf {
        let x0 = (wide::f64x4::new([self.0[2], self.0[3], self.0[4], self.0[5]])
            * wide::f64x4::new([1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0]))
        .to_array();

        let dst = [
            0.0,
            self.0[0],
            self.0[1] * (1.0 / 2.0),
            x0[0],
            x0[1],
            x0[2],
            x0[3],
        ];
        Poly6(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
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

impl MulAssign<f64> for Poly5 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly5 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly5 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-6 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
    }
}
impl HasIntegral for Poly6 {
    type IntegralOf = Poly7;
    fn indefinite(&self) -> Self::IntegralOf {
        let x0 = (wide::f64x4::new([0.0, self.0[0], self.0[1], self.0[2]])
            * wide::f64x4::new([1.0, 1.0, 1.0 / 2.0, 1.0 / 3.0]))
        .to_array();
        let x1 = (wide::f64x4::new([self.0[3], self.0[4], self.0[5], self.0[6]])
            * wide::f64x4::new([1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0, 1.0 / 7.0]))
        .to_array();

        let dst = [x0[0], x0[1], x0[2], x0[3], x1[0], x1[1], x1[2], x1[3]];
        Poly7(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
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

impl MulAssign<f64> for Poly6 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly6 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly6 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-7 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
pub struct Poly7(pub [f64; 8]);
impl Evaluate for Poly7 {
    fn evaluate(&self, x: f64) -> f64 {
        // P7(x) = (C0 + C1x) + (C2 + C3x) x2 + ((C4 + C5x) + (C6 + C7x) x2)x4
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = c[5].mul_add(x, c[4]);
        let t3 = c[7].mul_add(x, c[6]);

        let x2 = x * x;
        let left = t1.mul_add(x2, t0);
        let right = t3.mul_add(x2, t2);

        let x4 = x2 * x2;
        right.mul_add(x4, left)
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
    }
}
impl HasIntegral for Poly7 {
    type IntegralOf = Poly8;
    fn indefinite(&self) -> Self::IntegralOf {
        let x1 = (wide::f64x4::new([self.0[0], self.0[1], self.0[2], self.0[3]])
            * wide::f64x4::new([1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0]))
        .to_array();
        let x2 = (wide::f64x4::new([self.0[4], self.0[5], self.0[6], self.0[7]])
            * wide::f64x4::new([1.0 / 5.0, 1.0 / 6.0, 1.0 / 7.0, 1.0 / 8.0]))
        .to_array();

        let dst = [0.0, x1[0], x1[1], x1[2], x1[3], x2[0], x2[1], x2[2], x2[3]];
        Poly8(dst)
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let mut indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x));
        indef
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

impl MulAssign<f64> for Poly7 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly7 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly7 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

/// Rank-8 polynomal
///
/// Evaluated using manually unrolled Estin's scheme.
#[derive(Debug, PartialEq, Clone, Copy, Default, Arbitrary, Serialize, Deserialize)]
#[cfg_attr(
    feature = "borsh",
    derive(borsh::BorshDeserialize, borsh::BorshSerialize)
)]
pub struct Poly8(pub [f64; 9]);
impl Evaluate for Poly8 {
    fn evaluate(&self, x: f64) -> f64 {
        // P8(x) = (C0 + C1x) + (C2 + C3x) x2 + ((C4 + C5x) + (C6 + C7x) x2)x4 + C8x8
        let c = self.0;

        let t0 = c[1].mul_add(x, c[0]);
        let t1 = c[3].mul_add(x, c[2]);
        let t2 = c[5].mul_add(x, c[4]);
        let t3 = c[7].mul_add(x, c[6]);

        // Left side
        let x2 = x * x;
        let left = t1.mul_add(x2, t0);
        let right2 = t3.mul_add(x2, t2);

        // Reduce right and add to left
        let x4 = x2 * x2;
        let right4 = right2.mul_add(x4, left);

        // Finally add the x^8
        let x8 = x4 * x4;
        c[8].mul_add(x8, right4)
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
    fn translate(&mut self, v: f64) {
        self.0[0] += v;
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

impl MulAssign<f64> for Poly8 {
    fn mul_assign(&mut self, rhs: f64) {
        self.0.iter_mut().for_each(|x| *x *= rhs);
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

impl AbsDiffEq for Poly8 {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        <Self::Epsilon as AbsDiffEq>::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, eps: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, eps)
    }
}

impl RelativeEq for Poly8 {
    fn default_max_relative() -> Self::Epsilon {
        <Self::Epsilon as RelativeEq>::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, eps: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
        self.0.relative_eq(&other.0, eps, max_relative)
    }
}

#[cfg(test)]
mod tests {
    // Test results taken from the reference Haskell implementation.

    use super::*;

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
        let mut poly = Poly0(7.0);
        let result = Poly0(14.0);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly1([7.0, 3.0]);
        let result = Poly1([14.0, 3.0]);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly2([7.0, 3.0, 9.0]);
        let result = Poly2([14.0, 3.0, 9.0]);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly3([7.0, 3.0, 9.0, 8.0]);
        let result = Poly3([14.0, 3.0, 9.0, 8.0]);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly4([7.0, 3.0, 9.0, 8.0, 6.0]);
        let result = Poly4([14.0, 3.0, 9.0, 8.0, 6.0]);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly5([7.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        let result = Poly5([14.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly6([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        let result = Poly6([14.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly7([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        let result = Poly7([14.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
        poly.translate(7.0);
        assert_eq!(poly, result);
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
        let mut poly = Poly8([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
        let result = Poly8([14.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
        poly.translate(7.0);
        assert_eq!(poly, result);
    }
}
