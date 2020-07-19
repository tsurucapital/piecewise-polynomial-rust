use criterion::{black_box, criterion_group, criterion_main, Criterion};
use piecewise_polynomial::polynomial::*;

pub fn bench_evaluate_8(c: &mut Criterion) {
    let poly = Poly8([7.0, 3.0, 9.0, 8.0, 6.0, 1.5, 3.5, 4.5, 9.0]);
    let v = 3.0;
    c.bench_function("bench_evaluate_8", |b| {
        b.iter(|| poly.evaluate(black_box(v)))
    });
}

criterion_group!(benches, bench_evaluate_8,);
criterion_main!(benches);
