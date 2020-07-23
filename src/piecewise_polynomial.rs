use crate::polynomial::*;
use arbitrary::Arbitrary;
use std::cmp::Ordering;
use std::iter;
use std::ops::{Add, Mul, Neg};

/// A segment of a piecewise polynomial.
///
/// A segment is valid for _begin_ ≤ `x` < `end`, where _begin_ is implied
/// by the previous segment if it exists, or for `x` < `end` if this is the
/// first segment.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Segment<T> {
    pub end: f64,
    pub poly: T,
}

impl<T: Evaluate> Evaluate for Segment<T> {
    fn evaluate(&self, v: f64) -> f64 {
        self.poly.evaluate(v)
    }
}

impl<T> Mul<f64> for Segment<T>
where
    T: Mul<f64, Output = T> + Copy,
{
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Segment{
            end: self.end,
            poly: self.poly * rhs,
        }
    }
}


impl<T: HasDerivative> HasDerivative for Segment<T> {
    type DerivativeOf = Segment<T::DerivativeOf>;
    fn derivative(&self) -> Self::DerivativeOf {
        Segment {
            end: self.end,
            poly: self.poly.derivative(),
        }
    }
}

impl<T: Translate> Translate for Segment<T> {
    fn translate(&self, v: f64) -> Self {
        Segment {
            end: self.end,
            poly: self.poly.translate(v),
        }
    }
}

impl<T> HasIntegral for Segment<T>
where
    T: HasIntegral,
    T::IntegralOf: Translate,
{
    type IntegralOf = Segment<T::IntegralOf>;
    fn indefinite(&self) -> Self::IntegralOf {
        Segment {
            end: self.end,
            poly: self.poly.indefinite(),
        }
    }
    fn integral(&self, knot: Knot) -> Self::IntegralOf {
        let indef = self.indefinite();
        indef.translate(knot.y - indef.evaluate(knot.x))
    }
}

impl<T> Segment<T> {
    #[inline]
    pub fn integral_iter<'a, I>(
        segments: I,
        knot0: Knot,
    ) -> impl Iterator<Item = Segment<T::IntegralOf>> + 'a
    where
        T: HasIntegral + 'a,
        T::IntegralOf: Translate,
        I: IntoIterator<Item = &'a Segment<T>> + 'a,
    {
        let mut knot = knot0;
        segments.into_iter().map(move |seg| {
            let int = seg.integral(knot);
            knot = Knot {
                x: int.end,
                y: int.evaluate(int.end),
            };
            int
        })
    }
}

/// A piecewise-defined function.
///
/// The pieces are functions of type `T`.  This type doesn't provide any
/// guarantees that the function is smooth at the knots - the user must
/// ensure this for themself.
///
/// The following invariant must hold for all i: `segments[i].end <
/// segments[i+1].end`.
#[derive(Debug, PartialEq, Clone)]
pub struct Piecewise<T> {
    pub segments: Vec<Segment<T>>,
}

impl<T: Arbitrary> Arbitrary for Piecewise<T> {
    fn arbitrary(u: &mut arbitrary::Unstructured) -> arbitrary::Result<Self> {
        let mut ends = Vec::<f64>::arbitrary(u)?;
        if !ends.iter().all(|x| x.is_normal()) {
            return Err(arbitrary::Error::IncorrectFormat);
        }
        ends.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let segments = ends
            .into_iter()
            .map(|end| {
                Ok(Segment {
                    end,
                    poly: T::arbitrary(u)?,
                })
            })
            .collect::<arbitrary::Result<Vec<_>>>()?;
        Ok(Piecewise { segments })
    }
}

impl<T> Default for Piecewise<T> {
    fn default() -> Self {
        Piecewise {
            segments: Vec::new(),
        }
    }
}

impl<T> Mul<f64> for Piecewise<T>
where
    Segment<T>: Mul<f64, Output = Segment<T>> + Copy,
{
    type Output = Self;
    fn mul(mut self, rhs: f64) -> Self::Output {
        for seg in &mut self.segments {
            *seg = *seg * rhs;
        }
        self
    }
}

impl<T> Neg for Piecewise<T>
where
    T: Neg<Output = T> + Copy,
{
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        for seg in &mut self.segments {
            seg.poly = seg.poly.neg();
        }
        self
    }
}

impl<T: Add<Output = T> + Copy> Add for Piecewise<T> {
    type Output = Self;
    fn add(self, other: Piecewise<T>) -> Self::Output {
        let mut res = Vec::with_capacity(self.segments.len() + other.segments.len() - 1);

        let mut i = 0;
        let mut j = 0;
        let i_max = self.segments.len() - 1;
        let j_max = other.segments.len() - 1;

        loop {
            let a = &self.segments[i];
            let b = &other.segments[j];
            let a_last = i >= i_max;
            let b_last = j >= j_max;

            // Unlike the Haskell implementation, panic on NaN &c.
            let end = match a.end.partial_cmp(&b.end).unwrap() {
                Ordering::Less => {
                    if a_last {
                        j += 1;
                        b.end
                    } else {
                        i += 1;
                        a.end
                    }
                }
                Ordering::Greater => {
                    if b_last {
                        i += 1;
                        a.end
                    } else {
                        j += 1;
                        b.end
                    }
                }
                Ordering::Equal => {
                    i = i_max.min(i + 1);
                    j = j_max.min(j + 1);
                    a.end
                }
            };

            let ab = Segment {
                end,
                // I feel that copying here is silly: maybe we want to
                // implement add on the references of polynomials...
                poly: a.poly + b.poly,
            };

            res.push(ab);
            if a_last && b_last {
                break;
            }
        }

        Piecewise { segments: res }
    }
}

impl<T: Evaluate> Evaluate for Piecewise<T> {
    fn evaluate(&self, x: f64) -> f64 {
        let segment_ix = {
            assert!(!self.segments.is_empty(), "no segments to pick from");
            let close_ix = match self
                .segments
                .binary_search_by(|e| e.end.partial_cmp(&x).unwrap())
            {
                // We found exact end. The Haskell code doesn't use this:
                // it uses the next segment instead. Though if we're
                // exactly at the end and there isn't a next segment, we
                // do use the last one.
                Ok(exact_match_ix) => exact_match_ix + 1,
                // We haven't found the segment, we want to use the
                // "closest" one in some sense. The result tells us where
                // we could insert to preserve ordering: in other words,
                // the index of of the first element that's larger, or the
                // end of the vector. The first element that's larger is
                // precisely what we want. We just have to be slightly
                // careful to not run off the end.
                Err(closest) => closest,
            };
            let max_ix = self.segments.len() - 1;
            // Make sure we don't run off the end if the closest
            // segment was smaller equal, not larger.
            max_ix.min(close_ix)
        };
        self.segments[segment_ix].evaluate(x)
    }
}

impl<T: Evaluate> Piecewise<T> {
    /// Evaluate a piecewise at multiple successively increasing
    /// points. When the latter is sufficiently densely wrt the
    /// former, this is faster than doing a binary search of the
    /// segments every time, as evaluate does.
    pub fn evaluate_v(&self, xs: &[f64]) -> Vec<f64> {
        assert!(!self.segments.is_empty(), "no segments to pick from!");

        let find_seg = |start: usize, x: f64| -> usize {
            let seg_ix = self.segments[start..].iter().position(|seg| x < seg.end);
            // x is larger than all segments, use last segment. If we
            // did find a segment to fit into, remember to add the
            // offset we started looking from.
            seg_ix.map_or(self.segments.len() - 1, |i| i + start)
        };

        let mut prev_seg = 0;
        xs.iter()
            .map(|&x| {
                let seg_ix = find_seg(prev_seg, x);
                prev_seg = seg_ix;
                self.segments[seg_ix].poly.evaluate(x)
            })
            .collect()
    }
}

impl<T> HasDerivative for Piecewise<T>
where
    T: HasDerivative,
{
    type DerivativeOf = Piecewise<T::DerivativeOf>;
    fn derivative(&self) -> Self::DerivativeOf {
        Piecewise {
            segments: self.segments.iter().map(|seg| seg.derivative()).collect(),
        }
    }
}

impl<T> HasIntegral for Piecewise<T>
where
    T: HasIntegral + Clone,
    T::IntegralOf: Translate,
{
    type IntegralOf = Piecewise<T::IntegralOf>;
    fn indefinite(&self) -> Self::IntegralOf {
        if self.segments.is_empty() {
            return Default::default();
        }
        let indef_0 = self.segments[0].indefinite();
        // Everything here seems very inefficient. I think Piecewise
        // should also be able to accept slices.
        let res_vec = {
            let p_tail_int = Segment::integral_iter(
                &self.segments[1..],
                Knot {
                    x: indef_0.end,
                    y: indef_0.evaluate(indef_0.end),
                },
            );

            iter::once(indef_0).chain(p_tail_int).collect()
        };

        Piecewise { segments: res_vec }
    }

    // TODO: slow when using iter!
    fn integral(&self, knot0: Knot) -> Self::IntegralOf {
        // Slightly slower than manually inlined version... but only slightly;
        // perhaps come back to this.
        Piecewise {
            segments: Segment::integral_iter(&self.segments, knot0).collect(),
        }
    }
}

impl<T: Translate> Translate for Piecewise<T> {
    fn translate(&self, v: f64) -> Self {
        Piecewise {
            segments: self.segments.iter().map(|seg| seg.translate(v)).collect(),
        }
    }
}

mod tests {
    // We start tests from Poly1 because original implementation has
    // no Unbox instance for Poly0 for some reason so it can't do
    // piecewise evaluate that requires the instance.
    #[test]
    fn evaluate_piecewise_poly1() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly1([2.0, 3.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly1([4.0, 5.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 39.0);
    }

    #[test]
    fn evaluate_piecewise_log_poly1() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly1([2.0, 3.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly1([4.0, 5.0])),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 13.729550745276565);
    }

    #[test]
    fn integral_piecewise_poly1() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly1([2.0, 3.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly1([4.0, 5.0]),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly2([-1.5, 2.0, 1.5]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly2([-4.5, 4.0, 2.5]),
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn integral_piecewise_log_poly1() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly1([2.0, 3.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly1([4.0, 5.0])),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: IntOfLog {
                        k: 3.0,
                        poly: Poly1([-1.0, 3.0]),
                    },
                },
                Segment {
                    end: 2.0,
                    poly: IntOfLog {
                        k: 3.0,
                        poly: Poly1([-1.0, 5.0]),
                    },
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_piecewise_poly2() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly2([2.0, 3.0, 4.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly2([4.0, 5.0, 6.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 333.0);
    }

    #[test]
    fn evaluate_piecewise_log_poly2() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly2([2.0, 3.0, 4.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly2([4.0, 5.0, 6.0])),
                },
            ],
        };
        assert_approx_eq!(poly.evaluate(7.0), 36.44894859445539, 1e-12);
    }

    #[test]
    fn integral_piecewise_poly2() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly2([2.0, 3.0, 4.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly2([4.0, 5.0, 6.0]),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly3([-2.833333333333333, 2.0, 1.5, 1.3333333333333333]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly3([-6.5, 4.0, 2.5, 2.0]),
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn integral_piecewise_log_poly2() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly2([2.0, 3.0, 4.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly2([4.0, 5.0, 6.0])),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: IntOfLog {
                        k: -5.0,
                        poly: Poly2([7.0, -5.0, 4.0]),
                    },
                },
                Segment {
                    end: 2.0,
                    poly: IntOfLog {
                        k: -9.0,
                        poly: Poly2([11.0, -7.0, 6.0]),
                    },
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_piecewise_poly3() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly3([2.0, 3.0, 4.0, 5.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly3([4.0, 5.0, 6.0, 7.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 2734.0);
    }

    #[test]
    fn evaluate_piecewise_log_poly3() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly3([2.0, 3.0, 4.0, 5.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly3([4.0, 5.0, 6.0, 7.0])),
                },
            ],
        };
        assert_approx_eq!(poly.evaluate(7.0), 88.02717325878835, 1e-12);
    }

    #[test]
    fn integral_piecewise_poly3() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly3([2.0, 3.0, 4.0, 5.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly3([4.0, 5.0, 6.0, 7.0]),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly4([-4.083333333333333, 2.0, 1.5, 1.3333333333333333, 1.25]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly4([-8.25, 4.0, 2.5, 2.0, 1.75]),
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn integral_piecewise_log_poly3() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly3([2.0, 3.0, 4.0, 5.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly3([4.0, 5.0, 6.0, 7.0])),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: IntOfLog {
                        k: 25.0,
                        poly: Poly3([-23.0, 25.0, -11.0, 5.0]),
                    },
                },
                Segment {
                    end: 2.0,
                    poly: IntOfLog {
                        k: 33.0,
                        poly: Poly3([-31.0, 35.0, -15.0, 7.0]),
                    },
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_piecewise_poly4() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly4([2.0, 3.0, 4.0, 5.0, 6.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly4([4.0, 5.0, 6.0, 7.0, 8.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 21942.0);
    }

    #[test]
    fn evaluate_piecewise_log_poly4() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly4([2.0, 3.0, 4.0, 5.0, 6.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly4([4.0, 5.0, 6.0, 7.0, 8.0])),
                },
            ],
        };
        assert_approx_eq!(poly.evaluate(7.0), 202.73184850973757, 1e-12);
    }

    #[test]
    fn integral_piecewise_poly4() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly4([2.0, 3.0, 4.0, 5.0, 6.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly4([4.0, 5.0, 6.0, 7.0, 8.0]),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly5([-5.283333333333333, 2.0, 1.5, 1.3333333333333333, 1.25, 1.2]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly5([-9.85, 4.0, 2.5, 2.0, 1.75, 1.6]),
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn integral_piecewise_log_poly4() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly4([2.0, 3.0, 4.0, 5.0, 6.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly4([4.0, 5.0, 6.0, 7.0, 8.0])),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: IntOfLogPoly4 {
                        k: 2.0,
                        coeffs: [-2.0, 0.5, -1.1666666666666665, 0.9583333333333334],
                        u: -121.0,
                    },
                },
                Segment {
                    end: 2.0,
                    poly: IntOfLogPoly4 {
                        k: 2.0,
                        coeffs: [-4.0, 0.5, -1.8333333333333333, 1.2916666666666667],
                        u: -161.0,
                    },
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_piecewise_poly5() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly5([4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 173205.0);
    }

    #[test]
    fn evaluate_piecewise_log_poly5() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly5([4.0, 5.0, 6.0, 7.0, 8.0, 9.0])),
                },
            ],
        };
        assert_approx_eq!(poly.evaluate(7.0), 453.837464189018, 1e-12);
    }

    #[test]
    fn integral_piecewise_poly5() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly5([4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly6([
                        -6.449999999999999,
                        2.0,
                        1.5,
                        1.3333333333333333,
                        1.25,
                        1.2,
                        1.1666666666666667,
                    ]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly6([-11.349999999999998, 4.0, 2.5, 2.0, 1.75, 1.6, 1.5]),
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn integral_piecewise_log_poly5() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly5([4.0, 5.0, 6.0, 7.0, 8.0, 9.0])),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: IntOfLog {
                        k: 721.0,
                        poly: Poly5([-719.0, 721.0, -359.0, 121.0, -29.0, 7.0]),
                    },
                },
                Segment {
                    end: 2.0,
                    poly: IntOfLog {
                        k: 921.0,
                        poly: Poly5([-919.0, 923.0, -459.0, 155.0, -37.0, 9.0]),
                    },
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_piecewise_poly6() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly6([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 1349695.0);
    }

    #[test]
    fn evaluate_piecewise_log_poly6() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly6([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])),
                },
            ],
        };
        assert_approx_eq!(poly.evaluate(7.0), 996.7585375613455, 1e-12);
    }

    #[test]
    fn integral_piecewise_poly6() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly6([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly7([
                        -7.592857142857142,
                        2.0,
                        1.5,
                        1.3333333333333333,
                        1.25,
                        1.2,
                        1.1666666666666667,
                        1.1428571428571428,
                    ]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly7([
                        -12.778571428571428,
                        4.0,
                        2.5,
                        2.0,
                        1.75,
                        1.6,
                        1.5,
                        1.4285714285714286,
                    ]),
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn integral_piecewise_log_poly6() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly6([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: IntOfLog {
                        k: -5039.0,
                        poly: Poly6([5041.0, -5039.0, 2521.0, -839.0, 211.0, -41.0, 8.0]),
                    },
                },
                Segment {
                    end: 2.0,
                    poly: IntOfLog {
                        k: -6279.0,
                        poly: Poly6([6281.0, -6277.0, 3141.0, -1045.0, 263.0, -51.0, 10.0]),
                    },
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_piecewise_poly7() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 1.0408668e7);
    }

    #[test]
    fn evaluate_piecewise_log_poly7() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0])),
                },
            ],
        };
        assert_approx_eq!(poly.evaluate(7.0), 2158.8817270536842, 1e-12);
    }

    #[test]
    fn integral_piecewise_poly7() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly8([
                        -8.717857142857142,
                        2.0,
                        1.5,
                        1.3333333333333333,
                        1.25,
                        1.2,
                        1.1666666666666667,
                        1.1428571428571428,
                        1.125,
                    ]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly8([
                        -14.153571428571428,
                        4.0,
                        2.5,
                        2.0,
                        1.75,
                        1.6,
                        1.5,
                        1.4285714285714286,
                        1.375,
                    ]),
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn integral_piecewise_log_poly7() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0])),
                },
            ],
        };
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: IntOfLog {
                        k: 40321.0,
                        poly: Poly7([
                            -40319.0, 40321.0, -20159.0, 6721.0, -1679.0, 337.0, -55.0, 9.0,
                        ]),
                    },
                },
                Segment {
                    end: 2.0,
                    poly: IntOfLog {
                        k: 49161.0,
                        poly: Poly7([
                            -49159.0, 49163.0, -24579.0, 8195.0, -2047.0, 411.0, -67.0, 11.0,
                        ]),
                    },
                },
            ],
        };

        let knot = Knot { x: 1.0, y: 2.0 };
        assert_eq!(poly.integral(knot), result);
    }

    #[test]
    fn evaluate_piecewise_poly8() {
        use crate::piecewise_polynomial::*;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly8([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]),
                },
            ],
        };
        assert_eq!(poly.evaluate(7.0), 7.958628e7);
    }

    #[test]
    fn evaluate_piecewise_log_poly8() {
        use crate::log_polynomial::*;
        use crate::piecewise_polynomial::*;
        use assert_approx_eq::assert_approx_eq;
        let poly = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Log(Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])),
                },
                Segment {
                    end: 2.0,
                    poly: Log(Poly8([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])),
                },
            ],
        };
        assert_approx_eq!(poly.evaluate(7.0), 4625.849700383507, 1e-10);
    }
}
