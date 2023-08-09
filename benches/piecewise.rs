use criterion::{black_box, criterion_group, criterion_main, Criterion};
use piecewise_polynomial::*;

pub fn bench_evaluate_piecewise_poly8(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly8([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly8([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly8([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]),
            },
        ],
    };
    let v = 10.0;
    c.bench_function("bench_evaluate_piecewise_poly8", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

pub fn bench_indefinite_piecewise_poly1(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly1([2.0, 3.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly1([4.0, 5.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly1([2.0, 3.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly1([4.0, 5.0]),
            },
        ],
    };
    c.bench_function("bench_indefinite_piecewise_poly1", |b| {
        b.iter(|| black_box(&poly).indefinite())
    });
}

pub fn bench_indefinite_piecewise_poly2(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly2([2.0, 3.0, 4.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly2([4.0, 5.0, 6.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly2([2.0, 3.0, 4.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly2([4.0, 5.0, 6.0]),
            },
        ],
    };
    c.bench_function("bench_indefinite_piecewise_poly2", |b| {
        b.iter(|| black_box(&poly).indefinite())
    });
}

pub fn bench_indefinite_piecewise_poly3(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly3([2.0, 3.0, 4.0, 5.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly3([4.0, 5.0, 6.0, 7.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly3([2.0, 3.0, 4.0, 5.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly3([4.0, 5.0, 6.0, 7.0]),
            },
        ],
    };
    c.bench_function("bench_indefinite_piecewise_poly3", |b| {
        b.iter(|| black_box(&poly).indefinite())
    });
}

pub fn bench_indefinite_piecewise_poly4(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly4([2.0, 3.0, 4.0, 5.0, 6.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly4([4.0, 5.0, 6.0, 7.0, 8.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly4([2.0, 3.0, 4.0, 5.0, 6.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly4([4.0, 5.0, 6.0, 7.0, 8.0]),
            },
        ],
    };
    c.bench_function("bench_indefinite_piecewise_poly4", |b| {
        b.iter(|| black_box(&poly).indefinite())
    });
}

pub fn bench_indefinite_piecewise_poly5(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly5([4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly5([2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly5([4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
            },
        ],
    };
    c.bench_function("bench_indefinite_piecewise_poly5", |b| {
        b.iter(|| black_box(&poly).indefinite())
    });
}

pub fn bench_indefinite_piecewise_poly6(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly6([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly6([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly6([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
            },
        ],
    };
    c.bench_function("bench_indefinite_piecewise_poly6", |b| {
        b.iter(|| black_box(&poly).indefinite())
    });
}

pub fn bench_indefinite_piecewise_poly7(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]),
            },
        ],
    };
    c.bench_function("bench_indefinite_piecewise_poly7", |b| {
        b.iter(|| black_box(&poly).indefinite())
    });
}

pub fn bench_integral_piecewise_poly7(c: &mut Criterion) {
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
            Segment {
                end: 3.0,
                poly: Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
            },
            Segment {
                end: 4.0,
                poly: Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]),
            },
            Segment {
                end: 5.0,
                poly: Poly7([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
            },
            Segment {
                end: 6.0,
                poly: Poly7([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]),
            },
        ],
    };
    let v = Knot { x: 3.0, y: 5.0 };
    c.bench_function("bench_integral_piecewise_poly7", |b| {
        b.iter(|| poly.integral(black_box(v)))
    });
}

criterion_group!(
    piecewise_benches,
    bench_evaluate_piecewise_poly8,
    bench_indefinite_piecewise_poly1,
    bench_indefinite_piecewise_poly2,
    bench_indefinite_piecewise_poly3,
    bench_indefinite_piecewise_poly4,
    bench_indefinite_piecewise_poly5,
    bench_indefinite_piecewise_poly6,
    bench_indefinite_piecewise_poly7,
    bench_integral_piecewise_poly7,
);
criterion_main!(piecewise_benches);
