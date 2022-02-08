use crate::fasta::*;
use color_eyre::eyre::{Result, WrapErr};
use log::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

const PER_LINE: usize = 20;
lazy_static::lazy_static! {
    static ref PROT_DICT: HashMap<&'static str, u8> = maplit::hashmap! {
        "ATA" => b'I', "ATC" => b'I', "ATT" => b'I', "ATG" => b'M', "ACA" => b'T', "ACC" => b'T', "ACG" => b'T', "ACT" => b'T', "AAC" => b'N', "AAT" => b'N', "AAA" => b'K', "AAG" => b'K', "AGC" => b'S', "AGT" => b'S', "AGA" => b'R', "AGG" => b'R',
        "CTA" => b'L', "CTC" => b'L', "CTG" => b'L', "CTT" => b'L', "CCA" => b'P', "CCC" => b'P', "CCG" => b'P', "CCT" => b'P', "CAC" => b'H', "CAT" => b'H', "CAA" => b'Q', "CAG" => b'Q', "CGA" => b'R', "CGC" => b'R', "CGG" => b'R', "CGT" => b'R',
        "GTA" => b'V', "GTC" => b'V', "GTG" => b'V', "GTT" => b'V', "GCA" => b'A', "GCC" => b'A', "GCG" => b'A', "GCT" => b'A', "GAC" => b'D', "GAT" => b'D', "GAA" => b'E', "GAG" => b'E', "GGA" => b'G', "GGC" => b'G', "GGG" => b'G', "GGT" => b'G',
        "TCA" => b'S', "TCC" => b'S', "TCG" => b'S', "TCT" => b'S', "TTC" => b'F', "TTT" => b'F', "TTA" => b'L', "TTG" => b'L', "TAC" => b'Y', "TAT" => b'Y',                               "TGC" => b'C', "TGT" => b'C',                "TGG" => b'W',
        "NNN" => b'X',

        "ata" => b'i', "atc" => b'i', "att" => b'i', "atg" => b'm', "aca" => b't', "acc" => b't', "acg" => b't', "act" => b't', "aac" => b'n', "aat" => b'n', "aaa" => b'k', "aag" => b'k', "agc" => b's', "agt" => b's', "aga" => b'r', "agg" => b'r',
        "cta" => b'l', "ctc" => b'l', "ctg" => b'l', "ctt" => b'l', "cca" => b'p', "ccc" => b'p', "ccg" => b'p', "cct" => b'p', "cac" => b'h', "cat" => b'h', "caa" => b'q', "cag" => b'q', "cga" => b'r', "cgc" => b'r', "cgg" => b'r', "cgt" => b'r',
        "gta" => b'v', "gtc" => b'v', "gtg" => b'v', "gtt" => b'v', "gca" => b'a', "gcc" => b'a', "gcg" => b'a', "gct" => b'a', "gac" => b'd', "gat" => b'd', "gaa" => b'e', "gag" => b'e', "gga" => b'g', "ggc" => b'g', "ggg" => b'g', "ggt" => b'g',
        "tca" => b's', "tcc" => b's', "tcg" => b's', "tct" => b's', "ttc" => b'f', "ttt" => b'f', "tta" => b'l', "ttg" => b'l', "tac" => b'y', "tat" => b'y', "taa" => b'x', "tag" => b'x', "tgc" => b'c', "tgt" => b'c', "tga" => b'x', "tgg" => b'w',
        "nnn" => b'x',
    };
}

pub fn trans<'a>(
    dna_filename: &str,
    mut out: BufWriter<Box<dyn std::io::Write + 'a>>,
) -> Result<()> {
    for fragment in FastaReader::new(
        File::open(dna_filename).wrap_err(format!("while opening {}", dna_filename))?,
        true,
    ) {
        let seq = fragment.seq.unwrap();
        if seq.len() % 3 != 0 {
            return Err(eyre::eyre!(format!(
                "Length of '{}' ({}bp) is not a multiple of 3",
                fragment.id,
                seq.len()
            )));
        }

        writeln!(
            out,
            ">{} {}",
            fragment.id,
            fragment.name.unwrap_or_else(|| "".to_string())
        )?;

        let mut i = 0;
        for nnn in seq.as_slice().chunks(3) {
            let nnn = std::str::from_utf8(nnn).wrap_err(format!("invalid UTF-8: {:?}", nnn))?;
            let aa = PROT_DICT.get(&nnn).unwrap_or_else(|| {
                warn!("unknown triplet: {}", nnn);
                &b'X'
            });
            out.write_all(&[*aa])?;

            i += 1;
            if i >= PER_LINE {
                i = 0;
                writeln!(out)?;
            }
        }
        writeln!(out)?;
    }
    Ok(out.flush()?)
}

pub fn backtrans<'a>(
    proteins_filename: &str,
    dna_filename: &str,
    mut out: BufWriter<Box<dyn std::io::Write + 'a>>,
) -> Result<()> {
    let dna = FastaReader::new(
        File::open(dna_filename).wrap_err(format!("while opening {}", dna_filename))?,
        true,
    )
    .map(|r| (r.id, r.seq.unwrap()))
    .collect::<HashMap<_, _>>();
    for protein in FastaReader::new(
        File::open(proteins_filename).wrap_err(format!("while opening {}", dna_filename))?,
        true,
    ) {
        let dna_sequence = dna.get(&protein.id).unwrap();
        let prot_sequence = protein.seq.as_ref().unwrap();
        assert!(dna_sequence.len() >= 3 * prot_sequence.iter().filter(|x| **x != b'-').count());

        writeln!(
            out,
            ">{} {}",
            protein.id,
            protein.name.unwrap_or_else(|| "".to_string())
        )?;
        let mut pos = 0;
        let mut i = 0;
        for n in prot_sequence.iter() {
            match n {
                b'-' => {
                    write!(out, "---")?;
                }
                _ => {
                    out.write_all(&dna_sequence[pos..pos + 3])?;
                    pos += 3;
                }
            }

            i += 1;
            if i >= PER_LINE {
                i = 0;
                writeln!(out)?;
            }
        }
        writeln!(out)?;
    }
    Ok(out.flush()?)
}
