use crate::fasta::*;
use color_eyre::eyre::{Result, WrapErr};
use std::arch::x86_64::*;
use std::io::{BufWriter, Read, Write};

pub enum Distance {
    Levenshtein,
    Kimura,
}

fn kimura(s1: &[u8], s2: &[u8]) -> f64 {
    let dayhoff_pams = vec![
        195, /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
        196, /* 75.1% observed d; 196 PAMs estimated */
        197, 198, 199, 200, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 209, 210, 211, 212,
        213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 226, 227, 228, 229, 230, 231,
        232, 233, 234, 236, 237, 238, 239, 240, 241, 243, 244, 245, 246, 248, 249,
        250, /* 250 PAMs = 80.3% observed d */
        252, 253, 254, 255, 257, 258, 260, 261, 262, 264, 265, 267, 268, 270, 271, 273, 274, 276,
        277, 279, 281, 282, 284, 285, 287, 289, 291, 292, 294, 296, 298, 299, 301, 303, 305, 307,
        309, 311, 313, 315, 317, 319, 321, 323, 325, 328, 330, 332, 335, 337, 339, 342, 344, 347,
        349, 352, 354, 357, 360, 362, 365, 368, 371, 374, 377, 380, 383, 386, 389, 393, 396, 399,
        403, 407, 410, 414, 418, 422, 426, 430, 434, 438, 442, 447, 451, 456, 461, 466, 471, 476,
        482, 487, 493, 498, 504, 511, 517, 524, 531, 538, 545, 553, 560, 569, 577, 586, 595, 605,
        615, 626, 637, 649, 661, 675, 688, 703, 719, 736, 754, 775, 796, 819, 845, 874, 907, 945,
        /* 92.9% observed; 945 PAMs */
        988, /* 93.0% observed; 988 PAMs */
    ];
    let (matches, scored) = {
        #[cfg(target_feature = "avx2")]
        unsafe {
            kimura_avx2(s1, s2)
        }
        #[cfg(not(target_feature = "avx2"))]
        {
            kimura_x86(s1, s2)
        }
    };
    let d = 1. - matches as f64 / scored as f64;
    if d < 0.75 {
        -(1. - d - 0.2 * d * d).ln()
    } else if d > 0.930 {
        10.
    } else {
        dayhoff_pams[(d * 1000. - 750. + 0.5) as usize] as f64 / 100.0
    }
}

fn kimura_x86(s1: &[u8], s2: &[u8]) -> (usize, usize) {
    s1.iter()
        .zip(s2.iter())
        .filter(|(&n1, &n2)| n1 != b'-' && n2 != b'-') // Kimura: gaps are ignored
        .fold((0, 0), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1 } else { 0 }, scored + 1)
        })
}

unsafe fn print_mm256(a: __m256i) {
    let mut o = vec![b'_'; 32];
    _mm256_store_si256(o.as_mut_ptr() as *mut __m256i, a);
    println!("{}", std::str::from_utf8(&o).unwrap());
}

unsafe fn print_bool256(a: __m256i) {
    let mut o = vec![b'_'; 32];
    _mm256_store_si256(o.as_mut_ptr() as *mut __m256i, a);
    o.iter_mut().for_each(|b| {
        if *b == 0xff {
            *b = b'.'
        }
        if *b == 0x0 {
            *b = b' '
        }
    });
    println!("{}", std::str::from_utf8(&o).unwrap());
}


#[target_feature(enable = "avx2")]
unsafe fn kimura_avx2(s1: &[u8], s2: &[u8]) -> (usize, usize) {
    const SIMD_WIDTH: usize = 32;
    let n = s1.len();
    let lasts = n % SIMD_WIDTH;
    let (mut matches, mut scored) = s1[n - lasts..n]
        .iter()
        .zip(s2[n - lasts..n].iter())
        .filter(|(&n1, &n2)| n1 != b'-' && n2 != b'-')
        .fold((0, 0), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1 } else { 0 }, scored + 1)
        });

    let dashesx32 = _mm256_set1_epi8(b'-'.try_into().unwrap());
    for i in (0..n - lasts).step_by(SIMD_WIDTH) {
        let a = _mm256_loadu_si256(s1.as_ptr().offset(i.try_into().unwrap()) as *const __m256i);
        let b = _mm256_loadu_si256(s2.as_ptr().offset(i.try_into().unwrap()) as *const __m256i);

        let a_dashes = _mm256_cmpeq_epi8(a, dashesx32);
        let b_dashes = _mm256_cmpeq_epi8(b, dashesx32);

        let either_dashes = _mm256_or_si256(a_dashes, b_dashes); // Kimura: gaps are ignored

        let a_eq_b = _mm256_andnot_si256(either_dashes, _mm256_cmpeq_epi8(a, b));

        matches += _mm256_movemask_epi8(a_eq_b).count_ones();
        scored += SIMD_WIDTH - _mm256_movemask_epi8(either_dashes).count_ones() as usize;
    }

    (matches.try_into().unwrap(), scored)
}

fn levenshtein(s1: &[u8], s2: &[u8]) -> f64 {
        #[cfg(target_feature = "avx2")]
        unsafe {
            levenshtein_avx2(s1, s2)
        }
        #[cfg(not(target_feature = "avx2"))]
        {
            levenshtein_x86(s1, s2)
        }
}

#[target_feature(enable = "avx2")]
unsafe fn levenshtein_avx2(s1: &[u8], s2: &[u8]) -> f64 {
    const SIMD_WIDTH: usize = 32;
    let n = s1.len();
    let lasts = n % SIMD_WIDTH;
    let (mut matches, mut scored) = s1[n - lasts..n]
        .iter()
        .zip(s2[n - lasts..n].iter())
        .filter(|(&n1, &n2)| !(n1 == b'-' && n2 == b'-'))
        .fold((0, 0), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1 } else { 0 }, scored + 1)
        });

    let dashesx32 = _mm256_set1_epi8(b'-'.try_into().unwrap());
    for i in (0..n - lasts).step_by(SIMD_WIDTH) {
        let a = _mm256_loadu_si256(s1.as_ptr().offset(i.try_into().unwrap()) as *const __m256i);
        let b = _mm256_loadu_si256(s2.as_ptr().offset(i.try_into().unwrap()) as *const __m256i);

        let a_dashes = _mm256_cmpeq_epi8(a, dashesx32);
        let b_dashes = _mm256_cmpeq_epi8(b, dashesx32);

        let both_dashes = _mm256_and_si256(a_dashes, b_dashes); // Levenshtein: gaps are *not* ignored

        let a_eq_b = _mm256_andnot_si256(both_dashes, _mm256_cmpeq_epi8(a, b));

        matches += _mm256_movemask_epi8(a_eq_b).count_ones();
        scored += SIMD_WIDTH - _mm256_movemask_epi8(both_dashes).count_ones() as usize;
    }

    1.0 - matches as f64/scored as f64
}

fn levenshtein_x86(s1: &[u8], s2: &[u8]) -> f64 {
    let (matches, scored) = s1
        .iter()
        .zip(s2.iter())
        .filter(|(&n1, &n2)| !(n1 == b'-' && n2 == b'-'))
        .fold((0, 0), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1 } else { 0 }, scored + 1)
        });
    1.0 - matches as f64 / scored as f64
}

pub fn distance<T: Read>(fasta: FastaReader<T>, distance: Distance) -> (Vec<String>, Vec<f64>) {
    let fragments = fasta
        .map(|f| (f.id.clone(), f.seq.unwrap()))
        .collect::<Vec<_>>();
    let n = fragments.len();
    let mut r = vec![0f64; n * n];

    for i in 0..fragments.len() {
        for j in (i + 1)..fragments.len() {
            let d = match distance {
                Distance::Kimura => kimura(&fragments[i].1, &fragments[j].1),
                Distance::Levenshtein => levenshtein(&fragments[i].1, &fragments[j].1),
            };
            r[i * n + j] = d;
            r[j * n + i] = d;
        }
    }
    (
        fragments
            .into_iter()
            .map(|f| f.0.to_string())
            .collect::<Vec<_>>(),
        r,
    )
}
