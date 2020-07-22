use criterion::{black_box, criterion_group, criterion_main, Criterion};
use piecewise_polynomial::constrained_spline_interpolation::*;
use piecewise_polynomial::linear_interpolation::*;
use piecewise_polynomial::piecewise_polynomial::*;
use piecewise_polynomial::polynomial::*;

fn xs() -> Vec<f64> {
    (0..50).map(f64::from).collect()
}

fn p1() -> Piecewise<Poly1> {
    linear(xs().iter().map(|&x| Knot { x: x, y: x }).collect())
}

fn p3() -> Piecewise<Poly3> {
    constrained_spline(xs().iter().map(|&x| Knot { x: x, y: x }).collect())
}

pub fn bench_map_eval_linear(c: &mut Criterion) {
    let p1_ = p1();
    let xs_: Vec<f64> = black_box(xs());
    let run = |inp: &Vec<f64>| -> Vec<f64> { inp.iter().map(|x| p1_.evaluate(*x)).collect() };

    c.bench_function("map-eval linear", |b| b.iter(|| run(&xs_)));
}

pub fn bench_pp_eval_v_linear(c: &mut Criterion) {
    let p1_ = p1();
    let xs_: Vec<f64> = black_box(xs());

    c.bench_function("PP.evalV linear", |b| b.iter(|| p1_.evaluate_v(&xs_)));
}

pub fn bench_map_eval_cubic(c: &mut Criterion) {
    let p3_ = p3();
    let xs_: Vec<f64> = black_box(xs());
    let run = |inp: &Vec<f64>| -> Vec<f64> { inp.iter().map(|x| p3_.evaluate(*x)).collect() };

    c.bench_function("map-eval cubic", |b| b.iter(|| run(&xs_)));
}

pub fn bench_pp_eval_v_cubic(c: &mut Criterion) {
    let p3_ = p3();
    let xs_: Vec<f64> = black_box(xs());

    c.bench_function("PP.evalV cubic", |b| b.iter(|| p3_.evaluate_v(&xs_)));
}

criterion_group!(
    benches,
    bench_map_eval_linear,
    bench_pp_eval_v_linear,
    bench_map_eval_cubic,
    bench_pp_eval_v_cubic,
);
criterion_main!(benches);
