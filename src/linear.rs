use crate::piecewise::*;
use crate::poly::*;

/// Create a piecewise-linear function
#[inline]
pub fn linear(knots: &[Knot]) -> Piecewise<Poly1> {
    assert!(knots.len() >= 2, "need at least 2 knots");

    let mut knots_iter = knots.iter().cloned();
    // We checked that there are some knots, this should never fail.
    let mut prev_knot = knots_iter.next().unwrap();

    let segments = knots_iter
        // Force X co-ordinates to have increasing values.
        .map(|k| incr_linear(&mut prev_knot, k))
        .collect();
    Piecewise { segments }
}

fn incr_linear(prev_knot: &mut Knot, current_knot: Knot) -> Segment<Poly1> {
    // Make sure the X coordinate is increasing.
    let knot = Knot {
        x: prev_knot.x.max(current_knot.x),
        y: current_knot.y,
    };

    let seg = segment(*prev_knot, knot);
    *prev_knot = knot;
    seg
}

fn segment(knot0: Knot, knot1: Knot) -> Segment<Poly1> {
    let end = knot1.x;
    let mut poly = {
        let dx = knot1.x - knot0.x;
        let pv = if dx < f64::EPSILON {
            0.0
        } else {
            (knot1.y - knot0.y) / dx
        };
        Poly0(pv).indefinite()
    };
    poly.translate(knot0.y - poly.evaluate(knot0.x));
    Segment { end, poly }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_simple() {
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
        assert_eq!(linear(&knots), result);
    }

    #[test]
    fn test_linear_small_fracs() {
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
        assert_eq!(linear(&knots), result);
    }

    #[test]
    fn test_linear_big_fracs() {
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
        assert_eq!(linear(&knots), result);
    }
}
