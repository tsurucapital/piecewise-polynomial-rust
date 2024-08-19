use criterion::{black_box, criterion_group, criterion_main, Criterion};
use piecewise_polynomial::{Evaluate, IntOfLogPoly4, Piecewise, PiecewiseEvaluator, Segment};

fn evaluator(c: &mut Criterion) {
    let pp = Piecewise {
        segments: vec![
            Segment {
                end: 0.7908480990608203,
                poly: IntOfLogPoly4 {
                    k: 0.0,
                    coeffs: [0.0, 0.0, 0.0, 0.0],
                    u: 0.0,
                },
            },
            Segment {
                end: 0.7994702132819229,
                poly: IntOfLogPoly4 {
                    k: 0.0009370669680788741,
                    coeffs: [
                        -0.020752555105357416,
                        0.16650513012182652,
                        -0.6983098621266544,
                        1.4316735639121247,
                    ],
                    u: -129.92761978655466,
                },
            },
            Segment {
                end: 0.910236694826471,
                poly: IntOfLogPoly4 {
                    k: 4.177090756985528e-5,
                    coeffs: [
                        -9.632396946221638e-7,
                        -0.008561747393575922,
                        0.0719211900653713,
                        -0.226895682382893,
                    ],
                    u: 28.773672598908398,
                },
            },
            Segment {
                end: 0.9142027654012315,
                poly: IntOfLogPoly4 {
                    k: 0.00073292583408895,
                    coeffs: [
                        -0.03732198263996348,
                        0.7664146956563239,
                        -8.108156444533632,
                        42.58905496265079,
                    ],
                    u: -10391.344826496734,
                },
            },
            Segment {
                end: 0.9266792892609994,
                poly: IntOfLogPoly4 {
                    k: 7.332497880029383e-5,
                    coeffs: [
                        -5.285020782943302e-6,
                        -0.04693355843544987,
                        0.895864760994689,
                        -6.859151889173978,
                    ],
                    u: 2254.017969032283,
                },
            },
            Segment {
                end: 0.9299475163439032,
                poly: IntOfLogPoly4 {
                    k: 0.0023574348614676752,
                    coeffs: [
                        -0.1518909760249743,
                        3.8663606212051373,
                        -50.18782525984276,
                        324.3599226897312,
                    ],
                    u: -98214.38677984221,
                },
            },
            Segment {
                end: 0.9428689047374835,
                poly: IntOfLogPoly4 {
                    k: 0.00017797277291226545,
                    coeffs: [
                        -2.7008345333335962e-5,
                        -0.2397250758690307,
                        6.025512134429468,
                        -58.009084328949925,
                    ],
                    u: 23608.514218945078,
                },
            },
            Segment {
                end: 0.9454374179407244,
                poly: IntOfLogPoly4 {
                    k: 0.012766789103908128,
                    coeffs: [
                        -1.0804999197042418,
                        35.95330987760918,
                        -606.328031674692,
                        5096.048951590037,
                    ],
                    u: -2017844.3474780002,
                },
            },
            Segment {
                end: 0.9514855282456544,
                poly: IntOfLogPoly4 {
                    k: 0.0007569456648494242,
                    coeffs: [
                        -0.0002270356112813636,
                        -2.0137513470462443,
                        67.32633030227622,
                        -851.5561761106,
                    ],
                    u: 455538.13101202965,
                },
            },
            Segment {
                end: 0.9536766090786556,
                poly: IntOfLogPoly4 {
                    k: 0.05269397703660037,
                    coeffs: [
                        -5.265376146601179,
                        207.0996714811443,
                        -4120.813224418409,
                        40910.330943810295,
                    ],
                    u: -19201689.226010703,
                },
            },
            Segment {
                end: 0.9698602345262166,
                poly: IntOfLogPoly4 {
                    k: 0.0031511089755259015,
                    coeffs: [
                        -0.0013800621012309926,
                        -12.23434501367307,
                        485.8767529538376,
                        -7271.0901996702405,
                    ],
                    u: 4604597.224311776,
                },
            },
            Segment {
                end: 0.9712349811855948,
                poly: IntOfLogPoly4 {
                    k: 0.01838820278048891,
                    coeffs: [
                        -2.5035391216215768,
                        150.03671655158175,
                        -4803.3164319729585,
                        78705.73437709178,
                    ],
                    u: -61794462.67965331,
                },
            },
            Segment {
                end: 0.9772407044399701,
                poly: IntOfLogPoly4 {
                    k: 0.0038581069276785195,
                    coeffs: [
                        -0.0022704481507243393,
                        -20.10966258386586,
                        1012.3752868729623,
                        -20440.596373734945,
                    ],
                    u: 18548444.85657758,
                },
            },
            Segment {
                end: 0.9782832168474337,
                poly: IntOfLogPoly4 {
                    k: 0.01069995327916557,
                    coeffs: [
                        -1.4938931984391965,
                        108.72529227028951,
                        -4573.172245243078,
                        100403.10273958222,
                    ],
                    u: -105982701.50220917,
                },
            },
            Segment {
                end: 0.9876102189526036,
                poly: IntOfLogPoly4 {
                    k: 0.004176901369478019,
                    coeffs: [
                        -0.0029759907473936014,
                        -26.3384059404131,
                        1567.3032684266625,
                        -38922.567511457986,
                    ],
                    u: 44647303.905345015,
                },
            },
            Segment {
                end: 0.9881817944739316,
                poly: IntOfLogPoly4 {
                    k: 0.0036287745061995224,
                    coeffs: [
                        0.21730894299674777,
                        -61.56665344656205,
                        4390.07531030529,
                        -151895.99382136576,
                    ],
                    u: 260774533.05931947,
                },
            },
            Segment {
                end: 0.9910140593836662,
                poly: IntOfLogPoly4 {
                    k: 0.004151057615943841,
                    coeffs: [
                        -0.002783608841918936,
                        -24.650871240358484,
                        1287.9817344264904,
                        -21688.89628374626,
                    ],
                    u: -520533.51080991037,
                },
            },
            Segment {
                end: 0.9914299986261047,
                poly: IntOfLogPoly4 {
                    k: 0.0037216979896956544,
                    coeffs: [
                        0.235405587890172,
                        -77.30697733989999,
                        7117.088252254459,
                        -344090.3239850666,
                    ],
                    u: 852824115.7191182,
                },
            },
            Segment {
                end: 0.995009707384727,
                poly: IntOfLogPoly4 {
                    k: 0.004130632538155925,
                    coeffs: [
                        -0.0024963369579158153,
                        -22.144477980140625,
                        712.5795277345991,
                        27431.862963423453,
                    ],
                    u: -178694714.42229274,
                },
            },
            Segment {
                end: 0.9952422932381757,
                poly: IntOfLogPoly4 {
                    k: 0.004098405728249484,
                    coeffs: [
                        0.029739386937992676,
                        -35.01546973791965,
                        3284.275932968975,
                        -229380.44762049938,
                    ],
                    u: 1050238793.2899655,
                },
            },
            Segment {
                end: 1.0007933905967592,
                poly: IntOfLogPoly4 {
                    k: 0.0041290612995395345,
                    coeffs: [
                        -0.002426213925760131,
                        -21.542271082982225,
                        460.26839670649747,
                        66460.27059709422,
                    ],
                    u: -435013877.74513334,
                },
            },
            Segment {
                end: 1.0008291114115246,
                poly: IntOfLogPoly4 {
                    k: 0.004129065846414274,
                    coeffs: [
                        -0.002397551643877411,
                        -21.469975455114195,
                        551.4329333674855,
                        123943.15055456046,
                    ],
                    u: 1305218080.7041214,
                },
            },
            Segment {
                end: 1.0039321066108888,
                poly: IntOfLogPoly4 {
                    k: 0.0041290611568166255,
                    coeffs: [
                        -0.0024258403237478235,
                        -21.538256433300877,
                        469.0386975147955,
                        74227.41016686936,
                    ],
                    u: -135079363.0534059,
                },
            },
            Segment {
                end: 1.0041202266536582,
                poly: IntOfLogPoly4 {
                    k: 0.0041449406427057555,
                    coeffs: [
                        0.01779268538340848,
                        -11.224128275954737,
                        3098.108473127306,
                        409411.5290930405,
                    ],
                    u: 1918794340.278316,
                },
            },
            Segment {
                end: 1.0059277417349894,
                poly: IntOfLogPoly4 {
                    k: 0.00412834648917308,
                    coeffs: [
                        -0.002372372483551977,
                        -21.042686600654125,
                        709.3676968072226,
                        118735.77650730543,
                    ],
                    u: 218654457.0862798,
                },
            },
            Segment {
                end: 1.0062130184205806,
                poly: IntOfLogPoly4 {
                    k: 0.004391617075432078,
                    coeffs: [
                        0.22013230641350873,
                        54.363176874180375,
                        13474.1866940517,
                        1199690.3070026082,
                    ],
                    u: 4621122169.356427,
                },
            },
            Segment {
                end: 1.0084486937681607,
                poly: IntOfLogPoly4 {
                    k: 0.004116187172016964,
                    coeffs: [
                        -0.0019815363062930134,
                        -17.469260558472854,
                        1870.6931616345082,
                        262019.61295395604,
                    ],
                    u: 976530524.5763628,
                },
            },
            Segment {
                end: 1.0088569893712094,
                poly: IntOfLogPoly4 {
                    k: 0.0026081132431403833,
                    coeffs: [
                        -0.8969793596967209,
                        -230.6770577867382,
                        -23489.233438015406,
                        -1247247.1607592634,
                    ],
                    u: -3347029490.7171645,
                },
            },
            Segment {
                end: 1.013058684936052,
                poly: IntOfLogPoly4 {
                    k: 0.004186902470982606,
                    coeffs: [
                        -0.0030861403299772493,
                        -27.487177818175166,
                        -429.55313294854255,
                        62212.47240168345,
                    ],
                    u: 232672228.4441932,
                },
            },
            Segment {
                end: 1.013692769155078,
                poly: IntOfLogPoly4 {
                    k: 6.963623793947404e-5,
                    coeffs: [
                        -1.5863745500369262,
                        -272.346872629825,
                        -19323.037301111413,
                        -667485.4219936345,
                    ],
                    u: -1125921154.7424178,
                },
            },
            Segment {
                end: 1.0196161510690187,
                poly: IntOfLogPoly4 {
                    k: 0.0043824649440889025,
                    coeffs: [
                        -0.00435349350993047,
                        -38.90348728811344,
                        -2138.2673778270437,
                        -34251.04352495476,
                    ],
                    u: -822025.0445989143,
                },
            },
            Segment {
                end: 1.020573243057079,
                poly: IntOfLogPoly4 {
                    k: -0.021983483192136038,
                    coeffs: [
                        -6.768562789311658,
                        -738.6850381762822,
                        -38219.92994118669,
                        -965951.6366258537,
                    ],
                    u: -1163101377.941975,
                },
            },
            Segment {
                end: 1.0242771928058114,
                poly: IntOfLogPoly4 {
                    k: 0.005645427078081154,
                    coeffs: [
                        -0.007969993392688732,
                        -71.34520935596399,
                        -5393.666328976879,
                        -157237.78734892214,
                    ],
                    u: -200275966.63694373,
                },
            },
            Segment {
                end: 1.0254651826567727,
                poly: IntOfLogPoly4 {
                    k: 0.0763333447764266,
                    coeffs: [
                        14.667673590917119,
                        1159.6162198244301,
                        46028.2068591004,
                        918926.0767427449,
                    ],
                    u: 889428379.4691724,
                },
            },
            Segment {
                end: 1.0350440674827703,
                poly: IntOfLogPoly4 {
                    k: 0.0022472279573794644,
                    coeffs: [
                        -0.001615302278735928,
                        -14.431765310928544,
                        -760.0751785557769,
                        -15309.242775816765,
                    ],
                    u: -13475770.31744868,
                },
            },
            Segment {
                end: 1.0367695397994656,
                poly: IntOfLogPoly4 {
                    k: 0.029338694748671696,
                    coeffs: [
                        3.9085097110436413,
                        214.56572331680297,
                        5907.896086658175,
                        82044.05502622618,
                    ],
                    u: 55533311.132293984,
                },
            },
            Segment {
                end: 1.0479622538581264,
                poly: IntOfLogPoly4 {
                    k: 0.0009377708173767603,
                    coeffs: [
                        -0.0004361360571207046,
                        -3.8928692774436726,
                        -160.66189491308756,
                        -2494.321256621172,
                    ],
                    u: -1675007.5807777978,
                },
            },
            Segment {
                end: 1.0503401159254657,
                poly: IntOfLogPoly4 {
                    k: 0.01433274914871854,
                    coeffs: [
                        1.4180482169685862,
                        57.373811221789985,
                        1152.4088297295668,
                        11630.27091555109,
                    ],
                    u: 5731837.001373212,
                },
            },
            Segment {
                end: 1.0504716845944173,
                poly: IntOfLogPoly4 {
                    k: 0.00028750678699678765,
                    coeffs: [
                        -0.00012161588004725257,
                        -1.085362020550481,
                        -42.91441262313119,
                        -639.0376976537061,
                    ],
                    u: -412104.8053963409,
                },
            },
            Segment {
                end: 1.0529772202144343,
                poly: IntOfLogPoly4 {
                    k: 0.005783027899578761,
                    coeffs: [
                        0.5533464894691567,
                        21.672124382312997,
                        421.2327079180752,
                        4113.147594516959,
                    ],
                    u: 1961679.0816127567,
                },
            },
            Segment {
                end: 1.0643517510139924,
                poly: IntOfLogPoly4 {
                    k: 2.0534581806259356e-5,
                    coeffs: [
                        -4.856007686237794e-6,
                        -0.04330362533088677,
                        -1.311618006412137,
                        -15.08612429409524,
                    ],
                    u: -7591.31263054291,
                },
            },
            Segment {
                end: 1.067568962785786,
                poly: IntOfLogPoly4 {
                    k: 0.0003347268373371404,
                    coeffs: [
                        0.024923122951712715,
                        0.7685705299767573,
                        11.777065774691998,
                        90.95096330938,
                    ],
                    u: 34500.138001774576,
                },
            },
            Segment {
                end: 1.0919354438077176,
                poly: IntOfLogPoly4 {
                    k: 5.217450663303876e-6,
                    coeffs: [
                        -7.039854666126831e-7,
                        -0.0062736700376568155,
                        -0.1412779484821124,
                        -1.1945540573221736,
                    ],
                    u: -440.77122317437335,
                },
            },
            Segment {
                end: 1.0965942195299572,
                poly: IntOfLogPoly4 {
                    k: 0.00011161997847400183,
                    coeffs: [
                        0.00595972848902827,
                        0.13224522325681548,
                        1.4459520665110568,
                        7.963062112439504,
                    ],
                    u: 2169.633564141928,
                },
            },
            Segment {
                end: 1.1939095568985207,
                poly: IntOfLogPoly4 {
                    k: 7.843078602086309e-13,
                    coeffs: [
                        -4.1429237870537866e-12,
                        -1.7037894228660877e-11,
                        -1.5137051586224945e-11,
                        -3.784262896556236e-12,
                    ],
                    u: -9.081750334126259e-11,
                },
            },
        ],
    };

    // values to eval at; sampled directly from the knots above so that it's interesting
    let knots = {
        let mut knots = vec![];
        for segment in pp.segments.iter() {
            // Push in a bunch of duplicates so that we don't just go to the
            // next knot every time.
            knots.push(segment.end);
            knots.push(segment.end);
            knots.push(segment.end);
            knots.push(segment.end);
            knots.push(segment.end);
        }
        knots
    };

    // Bad case. Though we still give some sequential evals as we push
    // duplicates above.
    let shuffled_knots = {
        let mut k = knots.clone();
        k.reverse();
        k
    };

    c.bench_function("evaluate_increasing_order", |b| {
        b.iter_batched_ref(
            || PiecewiseEvaluator::new(&pp.segments),
            |pp| {
                for x in knots.iter() {
                    pp.evaluate(black_box(*x));
                }
            },
            criterion::BatchSize::SmallInput,
        );
    });

    c.bench_function("evaluate_increasing_order_vanilla", |b| {
        b.iter_batched_ref(
            || &pp,
            |pp| {
                for x in knots.iter() {
                    pp.evaluate(black_box(*x));
                }
            },
            criterion::BatchSize::SmallInput,
        );
    });

    c.bench_function("evaluate_shuffled_order", |b| {
        b.iter_batched_ref(
            || PiecewiseEvaluator::new(&pp.segments),
            |pp| {
                for x in shuffled_knots.iter() {
                    pp.evaluate(black_box(*x));
                }
            },
            criterion::BatchSize::SmallInput,
        );
    });

    c.bench_function("evaluate_shuffled_order_vanilla", |b| {
        b.iter_batched_ref(
            || &pp,
            |pp| {
                for x in shuffled_knots.iter() {
                    pp.evaluate(black_box(*x));
                }
            },
            criterion::BatchSize::SmallInput,
        );
    });
}

criterion_group!(piecewise_evaluator_benches, evaluator);
criterion_main!(piecewise_evaluator_benches);
