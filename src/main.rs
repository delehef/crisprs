use anyhow::Result;
use clap::Parser;
use std::fs::File;
use std::io::prelude::*;
mod fasta;
use fasta::*;
use std::collections::HashMap;

fn backtrans(
    proteins_filename: &str,
    dna_filename: &str,
    mut out: impl std::io::Write,
) -> Result<()> {
    let dna = FastaReader::new(File::open(dna_filename)?, true)
        .map(|r| (r.id, r.seq.unwrap()))
        .collect::<HashMap<_, _>>();
    for protein in FastaReader::new(File::open(proteins_filename)?, true) {
        let dna_sequence = dna.get(&protein.id).unwrap();
        let prot_sequence = protein.seq.as_ref().unwrap();
        assert!(dna_sequence.len() == 3 * prot_sequence.len());

        writeln!(
            out,
            ">{} {}",
            protein.id,
            protein.name.unwrap_or("".to_string())
        )?;
        let mut pos = 0;
        for n in prot_sequence.iter() {
            match n {
                b'-' => {
                    write!(out, "---")?;
                }
                _ => {
                    out.write(&dna_sequence[pos..pos + 3])?;
                    pos += 3;
                }
            }
        }
    }
    Ok(())
}

fn main() {
    backtrans("pipo.pp", "pipo.fa", File::open("out.test").unwrap());
}
