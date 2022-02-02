use crate::fasta::*;
use color_eyre::eyre::{Result, WrapErr};
use log::*;
use std::fs::File;
use std::io::{BufWriter, Read, Write};

pub enum Distance {
    Levenshtein,
    Kimura,
}

fn kimura(s1: &[u8], s2: &[u8]) -> f64 {
    let (matches, scored) = s1
        .iter()
        .zip(s2.iter())
        .filter(|(&n1, &n2)| n1 != b'-' && n2 != b'-')
        .fold((0f64, 0f64), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1. } else { 0. }, scored + 1.)
        });
    let d = 1.0 - matches / scored;
    -(1.0 - d - 0.2 * d.powi(2)).ln()
}

fn levenshtein(s1: &[u8], s2: &[u8]) -> f64 {
    let (matches, scored) = s1
        .iter()
        .zip(s2.iter())
        .filter(|(&n1, &n2)| !(n1 == b'-' && n2 == b'-'))
        .fold((0f64, 0f64), |(matches, scored), (n1, n2)| {
            (matches + if n1 == n2 { 1. } else { 0. }, scored + 1.)
        });
    1.0 - matches / scored
}

pub fn distance<T: Read>(fasta: FastaReader<T>, distance: Distance) -> (Vec<String>, Vec<f64>) {
    let fragments = fasta
        .map(|f| (f.id.clone(), f.seq.unwrap().clone()))
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
