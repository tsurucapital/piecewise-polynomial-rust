use crate::piecewise::*;
use crate::poly::*;
use std::iter;

/// Create a constrained cubic spline
///
/// Taken from “Constrained Cubic Spline Interpolation for Chemical
/// Engineering Applications” by CJC Kruger.
pub fn constrained_spline(ks0n: &[Knot]) -> Piecewise<Poly3> {
    assert!(ks0n.len() >= 3, "need at least 3 knots");

    // the knots are numbered 0, 1, 2, …, l, m, n
    let ks1n = &ks0n[1..];
    let ks2n = &ks1n[1..];
    let ks0m = &ks0n[..ks0n.len() - 1];
    let ks0l = &ks0m[..ks0m.len() - 1];
    let ks1m = &ks0m[1..];

    let Knot { x: x0, y: y0 } = ks0n.first().unwrap();
    let Knot { x: x1, y: y1 } = ks1n.first().unwrap();
    let Knot { x: xm, y: ym } = ks0m.last().unwrap();
    let Knot { x: xn, y: yn } = ks0n.last().unwrap();

    // derivatives: at intermediate knots, and endpoints
    // TODO: It's possible to re-write this with no intermediate allocations
    let f_mid: Vec<f64> = ks0l
        .iter()
        .zip(ks1m.iter())
        .zip(ks2n.iter())
        .map(|((&k0, &k1), &k2)| f_dx(k0, k1, k2))
        .collect(); // f' at x1 to xm              {-7a-}
    let f_x1 = f_mid[0];
    let f_x0 = (3.0 / 2.0) * (y1 - y0) / (x1 - x0) - (1.0 / 2.0) * f_x1; // {-7b-}
    let f_xm = f_mid.last().unwrap();
    let f_xn = (3.0 / 2.0) * (yn - ym) / (xn - xm) - (1.0 / 2.0) * f_xm; //{-7c-}
    let f_all = iter::once(f_x0).chain(f_mid).chain(iter::once(f_xn));

    let segments = f_all
        .clone()
        .zip(ks0m.iter())
        .zip(f_all.skip(1))
        .zip(ks1n.iter())
        .map(|(((f_x0_dx, &knot_0), f_x1_dx), &knot_1)| segment(f_x0_dx, knot_0, f_x1_dx, knot_1))
        .collect();

    Piecewise { segments }
}

/// The specified first-order derivative at some intermediate knot.
fn f_dx(knot_0: Knot, knot_1: Knot, knot_2: Knot) -> f64 {
    let (x0, y0) = (knot_0.x, knot_0.y);
    let (x1, y1) = (knot_1.x, knot_1.y);
    let (x2, y2) = (knot_2.x, knot_2.y);
    let slope01 = (y1 - y0) / (x1 - x0);
    let slope12 = (y2 - y1) / (x2 - x1);

    if slope01 * slope12 <= 0.0 {
        0.0 // slopes change sign, or at least one is flat
    } else {
        2.0 / (slope01.recip() + slope12.recip())
    }
}

/// Takes two knots and the first-order derivatives at each respective
/// /x/-coordinate, and produces the ‘constrained’ cubic segment
/// in-between.
fn segment(f_x0_dx: f64, knot_0: Knot, f_x1_dx: f64, knot_1: Knot) -> Segment<Poly3> {
    let (x0, y0) = (knot_0.x, knot_0.y);
    let (x1, y1) = (knot_1.x, knot_1.y);
    let end = x1;

    let dy = y1 - y0;
    let dx = x1 - x0;
    let slope = dy / dx;
    let x0x0 = x0 * x0;

    let f_x0_dx_dx = 2.0 * (3.0 * slope - (f_x1_dx + 2.0 * f_x0_dx)) / dx;
    let f_x1_dx_dx = 2.0 * ((2.0 * f_x1_dx + f_x0_dx) - 3.0 * slope) / dx;

    let d = (1.0 / 6.0) * (f_x1_dx_dx - f_x0_dx_dx) / dx;
    let c = (1.0 / 2.0) * (x1 * f_x0_dx_dx - x0 * f_x1_dx_dx) / dx;
    #[allow(clippy::suspicious_operation_groupings)]
    let b = slope - c * (x1 + x0) - d * (x1 * x1 + x1 * x0 + x0x0);
    let a = y0 - b * x0 - c * x0x0 - d * x0x0 * x0;

    let poly = Poly3([a, b, c, d]);

    Segment { end, poly }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constrained_spline_simple() {
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
                    poly: Poly3([0.0, 1.0, 0.0, 0.0]),
                },
                Segment {
                    end: 2.0,
                    poly: Poly3([0.0, 1.0, 0.0, 0.0]),
                },
                Segment {
                    end: 3.0,
                    poly: Poly3([0.0, 1.0, 0.0, 0.0]),
                },
            ],
        };
        assert_eq!(constrained_spline(&knots), result);
    }

    #[test]
    fn test_constrained_spline_small_fracs() {
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
                    poly: Poly3([
                        1.0148370992926183,
                        -0.5404433957934524,
                        0.29580682602721886,
                        -1.2629565745400926,
                    ]),
                },
                Segment {
                    end: 0.6844348417809908,
                    poly: Poly3([
                        -293566.6226738795,
                        1277545.8637799306,
                        -1853097.146507532,
                        895927.1016853832,
                    ]),
                },
                Segment {
                    end: 0.7267839023721823,
                    poly: Poly3([
                        424.542356105709,
                        -1744.0325950061097,
                        2395.8780877750582,
                        -1098.8493645108013,
                    ]),
                },
            ],
        };
        assert_eq!(constrained_spline(&knots), result);
    }

    #[test]
    fn test_constrained_spline_big_fracs() {
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
                    poly: Poly3([
                        -1.967768156634624e18,
                        -0.6317664571641676,
                        -2.4240440686379168e-20,
                        -1.0522021184848919e-39,
                    ]),
                },
                Segment {
                    end: 1.3625011620979218e18,
                    poly: Poly3([
                        9.741915586438009e18,
                        -0.7110587228107716,
                        -8.680586100321307e-19,
                        8.806686789053447e-38,
                    ]),
                },
                Segment {
                    end: 7.408886893918922e18,
                    poly: Poly3([
                        1.1041203706868632e19,
                        -2.78859360087304,
                        8.185200459035407e-20,
                        -3.6826046774330205e-39,
                    ]),
                },
            ],
        };
        assert_eq!(constrained_spline(&knots), result);
    }
}
