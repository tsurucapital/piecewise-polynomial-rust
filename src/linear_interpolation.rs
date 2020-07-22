use crate::piecewise_polynomial::*;
use crate::polynomial::*;
use std::vec::Vec;

// Seems slow.
pub fn linear(knots: Vec<Knot>) -> Piecewise<Poly1> {
    assert!(knots.len() >= 2, "need at least 2 knots");

    let mut knots_iter = knots.into_iter();
    // We checked that there are some knots, this should never fail.
    // Note that the knot is now no longer in the stream: that's fine,
    // we'll "insert" it back in a scan soon enough.
    let first_knot = knots_iter.next().unwrap();

    let segments = knots_iter
        // Force X co-ordinates to have increasing values.
        .scan(first_knot, incr_linear)
        .collect();
    Piecewise { segments }
}

fn incr_linear(prev_knot: &mut Knot, current_knot: Knot) -> Option<Segment<Poly1>> {
    // Make sure the X coordinate is increasing.
    let knot = Knot {
        x: prev_knot.x.max(current_knot.x),
        y: current_knot.y,
    };

    let seg = segment(*prev_knot, knot);
    *prev_knot = knot;
    Some(seg)
}

fn segment(knot0: Knot, knot1: Knot) -> Segment<Poly1> {
    let end = knot1.x;
    let indef = {
        let pv = if knot1.x == knot0.x {
            0.0
        } else {
            (knot1.y - knot0.y) / (knot1.x - knot0.x)
        };
        Poly0(pv).indefinite()
    };
    let poly = indef.translate(knot0.y - indef.evaluate(knot0.x));
    Segment { end, poly }
}

mod tests {
    #[test]
    fn test_linear_simple() {
        use crate::linear_interpolation::linear;
        use crate::piecewise_polynomial::{Piecewise, Segment};
        use crate::polynomial::{Knot, Poly1};
        let knots = vec![
            Knot { x: 0.0, y: 0.0 },
            Knot { x: 1.0, y: 1.0 },
            Knot { x: 2.0, y: 2.0 },
            Knot { x: 3.0, y: 3.0 },
        ];
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 1.0,
                    poly: Poly1([0.0, 1.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly1([0.0, 1.0]),
                },
                Segment {
                    end: 3.0,
                    poly: Poly1([0.0, 1.0]),
                },
            ],
        };
        assert_eq!(linear(knots), result);
    }

    #[test]
    fn test_linear_small_fracs() {
        use crate::linear_interpolation::linear;
        use crate::piecewise_polynomial::{Piecewise, Segment};
        use crate::polynomial::{Knot, Poly1};

        let knots = vec![
            Knot {
                x: 7.807257773555076e-2,
                y: 0.9738453165629335,
            },
            Knot {
                x: 0.6947124479037923,
                y: 0.35869674342227553,
            },
            Knot {
                x: 0.6844348417809908,
                y: 0.8995724066083576,
            },
            Knot {
                x: 0.7267839023721823,
                y: 0.6997825656440388,
            },
        ];
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 0.6947124479037923,
                    poly: Poly1([1.0517290816563807, -0.9975815754061172]),
                },
                Segment {
                    end: 0.6947124479037923,
                    poly: Poly1([0.35869674342227553, 0.0]),
                },
                Segment {
                    end: 0.7267839023721823,
                    poly: Poly1([5.227299096379497, -6.229522304990376]),
                },
            ],
        };
        assert_eq!(linear(knots), result);
    }

    #[test]
    fn test_linear_big_fracs() {
        use crate::linear_interpolation::linear;
        use crate::piecewise_polynomial::{Piecewise, Segment};
        use crate::polynomial::{Knot, Poly1};
        let knots = vec![
            Knot {
                x: -7.679272597449861e18,
                y: 1.930746322207704e18,
            },
            Knot {
                x: 6.358929964150921e18,
                y: -7.235865377340728e18,
            },
            Knot {
                x: 1.3625011620979218e18,
                y: 7.384377884237804e18,
            },
            Knot {
                x: 7.408886893918922e18,
                y: -6.623845605108912e18,
            },
        ];
        let result = Piecewise {
            segments: vec![
                Segment {
                    end: 6.358929964150921e18,
                    poly: Poly1([-3.0836356634429317e18, -0.652976166950476]),
                },
                Segment {
                    end: 6.358929964150921e18,
                    poly: Poly1([-7.235865377340728e18, 0.0]),
                },
                Segment {
                    end: 7.408886893918922e18,
                    poly: Poly1([9.222339324328898e19, -13.341712495237292]),
                },
            ],
        };
        assert_eq!(linear(knots), result);
    }
}
