use criterion::{black_box, criterion_group, criterion_main, Criterion};
use piecewise_polynomial::*;

pub fn bench_evaluate_0(c: &mut Criterion) {
    let poly = Poly0(7.0);
    let v = 3.0;
    c.bench_function("bench_evaluate_0", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_1(c: &mut Criterion) {
    let poly = Poly1([7.0, 3.0]);
    let v = 3.0;
    c.bench_function("bench_evaluate_1", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_2(c: &mut Criterion) {
    let poly = Poly2([7.0, 3.0, 9.0]);
    let v = 3.0;
    c.bench_function("bench_evaluate_2", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_3(c: &mut Criterion) {
    let poly = Poly3([7.0, 3.0, 9.0, 8.0]);
    let v = 3.0;
    c.bench_function("bench_evaluate_3", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_4(c: &mut Criterion) {
    let poly = Poly4([7.0, 3.0, 9.0, 8.0, 6.0]);
    let v = 3.0;
    c.bench_function("bench_evaluate_4", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_5(c: &mut Criterion) {
    let poly = Poly5([7.0, 3.0, 9.0, 8.0, 6.0, 1.5]);
    let v = 3.0;
    c.bench_function("bench_evaluate_5", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_6(c: &mut Criterion) {
    let poly = Poly6([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5]);
    let v = 3.0;
    c.bench_function("bench_evaluate_6", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_7(c: &mut Criterion) {
    let poly = Poly7([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5]);
    let v = 3.0;
    c.bench_function("bench_evaluate_7", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_8(c: &mut Criterion) {
    let poly = Poly8([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
    let v = 3.0;
    c.bench_function("bench_evaluate_8", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

// We have special case for this and we use it in trader, make sure we
// benchmark it.
//
// Checks the taylor exansion version.
pub fn bench_evaluate_int_of_log_taylor_poly4(c: &mut Criterion) {
    let poly = IntOfLogPoly4 {
        k: 2.0,
        coeffs: [-2.0, 0.5, -1.1666666666666665, 0.9583333333333334],
        u: -121.0,
    };
    let v = 3.0;
    c.bench_function("bench_evaluate_int_of_log_taylor_poly4", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

// We have special case for this and we use it in trader, make sure we
// benchmark it.
//
// Checks the analytical exansion version.
pub fn bench_evaluate_int_of_log_anal_poly4(c: &mut Criterion) {
    let poly = IntOfLogPoly4 {
        k: 2.0,
        coeffs: [-2.0, 0.5, -1.1666666666666665, 0.9583333333333334],
        u: -121.0,
    };
    // We choose 10 as it falls out of the taylor threshold.
    let v = 10.0;
    c.bench_function("bench_evaluate_int_of_log_anal_poly4", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_evaluate_piecewise_of_log_poly4(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Log(Poly8([
                    14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 110.0, 111.0, 112.0,
                ])),
            },
            Segment {
                end: 4.0,
                poly: Log(Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])),
            },
            Segment {
                end: 5.0,
                poly: Log(Poly8([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])),
            },
            Segment {
                end: 6.0,
                poly: Log(Poly8([
                    14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 110.0, 111.0, 112.0,
                ])),
            },
            Segment {
                end: 7.0,
                poly: Log(Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])),
            },
            Segment {
                end: 8.0,
                poly: Log(Poly8([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])),
            },
            Segment {
                end: 9.0,
                poly: Log(Poly8([
                    14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 110.0, 111.0, 112.0,
                ])),
            },
        ],
    };

    // Pick some segment near end, don't match exactly.
    let v = poly.segments[poly.segments.len() - 2].end - 0.1;
    c.bench_function("bench_evaluate_piecewise_of_log_poly4", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

criterion_group!(
    benches,
    bench_evaluate_0,
    bench_evaluate_1,
    bench_evaluate_2,
    bench_evaluate_3,
    bench_evaluate_4,
    bench_evaluate_5,
    bench_evaluate_6,
    bench_evaluate_7,
    bench_evaluate_8,
    bench_evaluate_int_of_log_taylor_poly4,
    bench_evaluate_int_of_log_anal_poly4,
    bench_evaluate_piecewise_of_log_poly4,
);
criterion_main!(benches);
