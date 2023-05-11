use crate::fasta::*;
use aa_similarity::{AminoAcid, Blosum62, Similarity};
use once_cell::sync::OnceCell;
use rayon::prelude::*;
use smartstring::{LazyCompact, SmartString};
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use std::arch::x86_64::*;
use std::collections::HashMap;
use std::io::Read;
use std::sync::Mutex;

const ALL_AAS: [AminoAcid; 20] = [
    AminoAcid::Alanine,
    AminoAcid::Arginine,
    AminoAcid::Asparagine,
    AminoAcid::AsparticAcid,
    AminoAcid::Cysteine,
    AminoAcid::GlutamicAcid,
    AminoAcid::Glutamine,
    AminoAcid::Glycine,
    AminoAcid::Histidine,
    AminoAcid::Isoleucine,
    AminoAcid::Leucine,
    AminoAcid::Lysine,
    AminoAcid::Methionine,
    AminoAcid::Phenylalanine,
    AminoAcid::Proline,
    AminoAcid::Serine,
    AminoAcid::Threonine,
    AminoAcid::Tryptophan,
    AminoAcid::Tyrosine,
    AminoAcid::Valine,
];

static SIGMA_R_BASE: OnceCell<f32> = OnceCell::new();

pub enum Distance {
    Levenshtein,
    Kimura,
    ScoreDist,
}

#[derive(Clone, Copy)]
pub enum FastaMode {
    Msa,      // the fasta file contains a global alignment, each sequence appearing only once
    Pairwise, // the fasta file is a succession of pairwise alignment (e.g. from needleall)
}

fn kimura(s1: &[u8], s2: &[u8]) -> f32 {
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
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                unsafe { kimura_avx2(s1, s2) }
            } else {
                kimura_cpu(s1, s2)
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        kimura_cpu(s1, s2)
    };
    if scored > 0 {
        let d = 1. - matches as f32 / scored as f32;
        if d < 0.75 {
            -(1. - d - 0.2 * d * d).ln()
        } else if d > 0.930 {
            10.
        } else {
            dayhoff_pams[(d * 1000. - 750. + 0.5) as usize] as f32 / 100.0
        }
    } else {
        f32::INFINITY
    }
}

fn kimura_cpu(s1: &[u8], s2: &[u8]) -> (usize, usize) {
    s1.iter()
        .zip(s2.iter())
        .filter(|(&n1, &n2)| n1 != b'-' && n2 != b'-') // Kimura: gaps are ignored
        .fold((0, 0), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1 } else { 0 }, scored + 1)
        })
}

// from https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-108
fn scoredist(s1: &[u8], s2: &[u8]) -> f32 {
    fn sigma(s1: &[u8], s2: &[u8]) -> i16 {
        s1.iter()
            .zip(s2.iter())
            .filter(|(&n1, &n2)| n1 != b'-' && n2 != b'-' && n1 != b'X' && n2 != b'X') // gaps are ignored
            .fold(0, |score, (n1, n2)| {
                score
                    + Blosum62::similarity(
                        AminoAcid::try_from(*n1 as char).unwrap(),
                        AminoAcid::try_from(*n2 as char).unwrap(),
                    )
            })
    }

    assert!(s1.len() == s2.len());
    let l = s1.len() as f32;
    let sigma_r: f32 = SIGMA_R_BASE.get_or_init(|| {
        let mut s = 0f32;
        for i in 0..ALL_AAS.len() {
            for j in 0..=i {
                s += Blosum62::similarity(ALL_AAS[i].clone(), ALL_AAS[j].clone()) as f32;
            }
        }
        s / (ALL_AAS.len().pow(2) as f32)
    }) * l;

    let sigma_s1_s2 = sigma(s1, s2) as f32;
    let sigma_n: f32 = sigma_s1_s2 - sigma_r;
    let sigma_un: f32 = (sigma(s1, s1) as f32 + sigma(s2, s2) as f32) / 2.0 - sigma_r;

    -(sigma_n / sigma_un).ln()
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
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

fn levenshtein(s1: &[u8], s2: &[u8]) -> f32 {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { levenshtein_avx2(s1, s2) }
        } else {
            println!("Not using AVX2");
            levenshtein_cpu(s1, s2)
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    levenshtein_cpu(s1, s2)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn levenshtein_avx2(s1: &[u8], s2: &[u8]) -> f32 {
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

    if scored > 0 {
        1.0 - matches as f32 / scored as f32
    } else {
        f32::INFINITY
    }
}

fn levenshtein_cpu(s1: &[u8], s2: &[u8]) -> f32 {
    let (matches, scored) = s1
        .iter()
        .zip(s2.iter())
        .filter(|(&n1, &n2)| !(n1 == b'-' && n2 == b'-'))
        .fold((0, 0), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1 } else { 0 }, scored + 1)
        });
    if scored > 0 {
        1.0 - matches as f32 / scored as f32
    } else {
        f32::INFINITY
    }
}

pub fn distance<T: Read>(
    fasta: FastaReader<T>,
    distance: Distance,
    mode: FastaMode,
) -> (Vec<String>, Vec<f32>) {
    let fragments = fasta
        .map(|f| (f.id.clone(), f.seq.unwrap()))
        .collect::<Vec<_>>();

    match mode {
        FastaMode::Msa => {
            let n = fragments.len();
            let r = Mutex::new(vec![0f32; n * n]);

            (0..fragments.len()).into_par_iter().for_each(|i| {
                for j in (i + 1)..fragments.len() {
                    let d = match distance {
                        Distance::Kimura => kimura(&fragments[i].1, &fragments[j].1),
                        Distance::ScoreDist => scoredist(&fragments[i].1, &fragments[j].1),
                        Distance::Levenshtein => levenshtein(&fragments[i].1, &fragments[j].1),
                    };
                    r.lock().unwrap()[i * n + j] = d;
                    r.lock().unwrap()[j * n + i] = d;
                }
            });

            (
                fragments
                    .into_iter()
                    .map(|f| f.0.to_string())
                    .collect::<Vec<_>>(),
                r.into_inner().unwrap(),
            )
        }
        FastaMode::Pairwise => {
            let r = Mutex::new(HashMap::<
                SmartString<LazyCompact>,
                HashMap<SmartString<LazyCompact>, f32>,
            >::new());

            fragments.par_iter().chunks(2).for_each(|f| {
                let d = match distance {
                    Distance::Kimura => kimura(&f[0].1, &f[1].1),
                    Distance::ScoreDist => scoredist(&f[0].1, &f[1].1),
                    Distance::Levenshtein => levenshtein(&f[0].1, &f[1].1),
                };
                r.lock()
                    .unwrap()
                    .entry(f[0].0.clone())
                    .or_default()
                    .insert(f[1].0.clone(), d);
            });

            let r = r.into_inner().unwrap();
            let ids = r.keys().cloned().collect::<Vec<_>>();
            let n = ids.len();
            let mut rr = vec![0f32; n * n];
            for i in 0..n {
                for j in (i + 1)..n {
                    rr[i * n + j] = r[&ids[i]][&ids[j]];
                    rr[j * n + i] = r[&ids[j]][&ids[i]];
                }
            }

            (ids.into_iter().map(|s| s.to_string()).collect(), rr)
        }
    }
}
