use criterion::Throughput::Elements;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use piecewise_polynomial::*;

fn xs(elems: u32) -> Vec<f64> {
    (0..elems).map(f64::from).collect()
}

fn p1(elems: u32) -> Piecewise<Poly1> {
    let e: Vec<_> = xs(elems).iter().map(|&x| Knot { x: x, y: x }).collect();
    linear(&e)
}

fn p3(elems: u32) -> Piecewise<Poly3> {
    let knots: Vec<_> = xs(elems).iter().map(|&x| Knot { x: x, y: x }).collect();
    constrained_spline(&knots)
}

pub fn bench_map_eval_linear(c: &mut Criterion) {
    let mut group = c.benchmark_group("map-eval linear");
    for size in [50].iter() {
        group.throughput(Elements(u64::from(*size)));
        let input = (p1(*size), xs(*size));
        group.bench_with_input(
            BenchmarkId::from_parameter(size),
            &input,
            |b, (p1_, xs_)| {
                b.iter_with_large_drop(|| xs_.iter().map(|&x| p1_.evaluate(x)).collect::<Vec<_>>())
            },
        );
    }
}

pub fn bench_pp_eval_v_linear(c: &mut Criterion) {
    let mut group = c.benchmark_group("PP.evalV linear");
    for size in [50].iter() {
        group.throughput(Elements(u64::from(*size)));
        let input = (p1(*size), xs(*size));
        group.bench_with_input(
            BenchmarkId::from_parameter(size),
            &input,
            |b, (p1_, xs_)| b.iter_with_large_drop(|| p1_.evaluate_v(xs_.iter().copied())),
        );
    }
}

pub fn bench_map_eval_cubic(c: &mut Criterion) {
    let p3_ = p3(50);
    let xs_: Vec<f64> = black_box(xs(50));
    let run = |inp: &Vec<f64>| -> Vec<f64> { inp.iter().map(|x| p3_.evaluate(*x)).collect() };

    c.bench_function("map-eval cubic", |b| b.iter(|| run(&xs_)));
}

pub fn bench_pp_eval_v_cubic(c: &mut Criterion) {
    let p3_ = p3(50);
    let xs_: Vec<f64> = black_box(xs(50));

    c.bench_function("PP.evalV cubic", |b| {
        b.iter(|| p3_.evaluate_v(xs_.iter().copied()))
    });
}

criterion_group!(
    benches,
    bench_map_eval_linear,
    bench_pp_eval_v_linear,
    bench_map_eval_cubic,
    bench_pp_eval_v_cubic,
);
criterion_main!(benches);
