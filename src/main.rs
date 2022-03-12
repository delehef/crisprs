use clap::command;
use clap::*;
use color_eyre::eyre::Result;
use log::*;
use simplelog::*;
use std::fs::File;
use std::io::BufWriter;

mod distance;
mod fasta;
mod io;
mod simd;
mod transformation;

fn main() -> Result<()> {
    color_eyre::install()?;
    TermLogger::init(
        LevelFilter::Info,
        ConfigBuilder::new()
            .set_time_level(LevelFilter::Off)
            .build(),
        TerminalMode::Stderr,
        simplelog::ColorChoice::Auto,
    )?;

    let main_args = command!()
        .propagate_version(true)
        .subcommand_required(true)
        .subcommand(
            Command::new("trans")
                .about("Convert a DNA FASTA file to a protein FASTA file")
                .arg(
                    Arg::new("DNA_FILE")
                        .long("from")
                        .takes_value(true)
                        .required(true),
                )
                .arg(Arg::new("OUT_FILE").long("to").takes_value(true)),
        )
        .subcommand(
            Command::new("backtrans")
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
        .subcommand(
            Command::new("distance")
                .about("Compute various distances on FASTA alignments")
                .arg(
                    Arg::new("METRIC")
                        .takes_value(true)
                        .required(true)
                        .possible_values(["kimura", "levenshtein"]),
                )
                .arg(Arg::new("FASTA_FILE").takes_value(true).required(true))
                .arg(Arg::new("OUT_FILE").long("to").takes_value(true)),
        )
        .get_matches();

    match main_args.subcommand() {
        Some(("trans", args)) => {
            let stdout = std::io::stdout();
            let out: BufWriter<Box<dyn std::io::Write>> =
                if let Some(out_filename) = args.value_of("OUT_FILE") {
                    BufWriter::with_capacity(30_000_000, Box::new(File::create(out_filename)?))
                } else {
                    BufWriter::new(Box::new(stdout.lock()))
                };

            transformation::trans(args.value_of("DNA_FILE").unwrap(), out)
        }
        Some(("backtrans", args)) => {
            let stdout = std::io::stdout();
            let out: BufWriter<Box<dyn std::io::Write>> =
                if let Some(out_filename) = args.value_of("OUT_FILE") {
                    BufWriter::with_capacity(30_000_000, Box::new(File::create(out_filename)?))
                } else {
                    BufWriter::new(Box::new(stdout.lock()))
                };

            transformation::backtrans(
                args.value_of("PROTEIN_FILE").unwrap(),
                args.value_of("DNA_FILE").unwrap(),
                out,
            )
        }
        Some(("distance", args)) => {
            let stdout = std::io::stdout();
            let out: BufWriter<Box<dyn std::io::Write>> =
                if let Some(out_filename) = args.value_of("OUT_FILE") {
                    BufWriter::with_capacity(30_000_000, Box::new(File::create(out_filename)?))
                } else {
                    BufWriter::new(Box::new(stdout.lock()))
                };

            let distance = match args.value_of("METRIC").unwrap() {
                "kimura" => distance::Distance::Kimura,
                "levenshtein" => distance::Distance::Levenshtein,
                _ => unreachable!(),
            };
            let (ids, m) = distance::distance(
                io::open_fasta(args.value_of("FASTA_FILE").unwrap())?,
                distance,
            );
            io::write_dist_matrix(&m, &ids, out)
        }
        _ => unreachable!(),
    }
}
