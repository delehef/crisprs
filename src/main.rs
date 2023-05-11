use clap::command;
use clap::*;
use eyre::Result;
use log::*;
use std::fs::File;
use std::io::BufWriter;

mod distance;
mod fasta;
mod io;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod simd;
mod transformation;

fn main() -> Result<()> {
    let main_args = command!()
        .propagate_version(true)
        .subcommand_required(true)
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .takes_value(false)
                .multiple_occurrences(true)
                .help("verbosity level (-v, -vv, -vvv)"),
        )
        .arg(
            Arg::new("parallel")
                .short('P')
                .long("parallel")
                .takes_value(true)
                .default_missing_value("0")
                .require_equals(true)
                .min_values(0),
        )
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
                        .possible_values(["kimura", "scoredist", "levenshtein"]),
                )
                .arg(
                    Arg::new("mode")
                        .help("whether the FASTA represents an MSA or a succession of pairwise alignments")
                        .takes_value(true)
                        .short('M')
                        .long("mode")
                        .value_parser(["msa", "pairwise"])
                        .default_value("msa"),
                )
                .arg(Arg::new("FASTA_FILE").takes_value(true).required(true))
                .arg(Arg::new("OUT_FILE").long("to").takes_value(true)),
        )
        .get_matches();

    buche::new()
        .timestamp(buche::Timestamp::Off)
        .verbosity(main_args.occurrences_of("verbose") as usize)
        .init()
        .unwrap();

    if main_args.is_present("parallel") {
        rayon::ThreadPoolBuilder::new()
            .num_threads(main_args.value_of_t("parallel").unwrap())
            .build_global()
            .unwrap();
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()
            .unwrap();
    }
    info!("Using {} threads", rayon::current_num_threads());

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if is_x86_feature_detected!("avx2") {
        info!("Using AVX2 acceleration");
    }

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
                "scoredist" => distance::Distance::ScoreDist,
                "levenshtein" => distance::Distance::Levenshtein,
                _ => unreachable!(),
            };
            let (ids, m) = distance::distance(
                io::open_fasta(args.value_of("FASTA_FILE").unwrap())?,
                distance,
                match args.value_of("mode").unwrap() {
                    "msa" => distance::FastaMode::Msa,
                    "pairwise" => distance::FastaMode::Pairwise,
                    _ => unreachable!(),
                },
            );
            io::write_dist_matrix(&m, &ids, out)
        }
        _ => unreachable!(),
    }
}
