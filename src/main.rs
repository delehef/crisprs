use clap::app_from_crate;
use clap::*;
use color_eyre::eyre::{Result, WrapErr};
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
mod fasta;
use fasta::*;
use std::collections::HashMap;

fn backtrans<'a>(
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
        const PER_LINE: usize = 20;
        let dna_sequence = dna.get(&protein.id).unwrap();
        let prot_sequence = protein.seq.as_ref().unwrap();
        assert!(dna_sequence.len() >= 3 * prot_sequence.iter().filter(|x| **x != b'-').count());

        writeln!(
            out,
            ">{} {}",
            protein.id,
            protein.name.unwrap_or("".to_string())
        )?;
        let mut pos = 0;
        let mut i = 0;
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

fn main() -> Result<()> {
    color_eyre::install()?;

    let main_args = app_from_crate!()
        .global_setting(AppSettings::PropagateVersion)
        .global_setting(AppSettings::UseLongFormatForHelpSubcommand)
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommand(
            App::new("backtrans")
                .about("Convert a protein FASTA file to a DNA FASTA file")
                .arg(
                    Arg::new("PROTEIN_FILE")
                        .long("from")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::new("DNA_FILE")
                        .long("with")
                        .takes_value(true)
                        .required(true),
                )
                .arg(Arg::new("OUT_FILE").long("to").takes_value(true)),
        )
        .get_matches();

    match main_args.subcommand() {
        Some(("backtrans", args)) => {
            let stdout = std::io::stdout();
            let out: BufWriter<Box<dyn std::io::Write>> =
                if let Some(out_filename) = args.value_of("OUT_FILE") {
                    BufWriter::with_capacity(30_000_000, Box::new(File::create(out_filename)?))
                } else {
                    BufWriter::new(Box::new(stdout.lock()))
                };

            backtrans(
                args.value_of("PROTEIN_FILE").unwrap(),
                args.value_of("DNA_FILE").unwrap(),
                out,
            )
        }
        _ => unreachable!(),
    }
}
