mod design;
mod minimizer_index;
mod coverage;

use std::path::PathBuf;

use clap::ArgAction;
use coverage::compute_coverage;
use jseqio::reader::*;
use jseqio::writer::*;
use jseqio::record::*;
use clap::{Command, Arg};

use log::{info, error};

use jseqio::seq_db::SeqDB;

fn main() {
    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = Command::new("syotti2")
        .about("Bait design for targeted sequencing")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .subcommand(Command::new("design")
            .arg(Arg::new("targets")
                .help("Input FASTA or FASTQ file, possibly gzipped")
                .long("targets")
                .short('t')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))
            )
            .arg(Arg::new("output")
                .help("FASTA outputfile for the baits")
                .long("output")
                .short('o')
                .required(true)
                .value_parser(clap::value_parser!(PathBuf))
            )
            .arg(Arg::new("bait-length")
                .help("Bait length")
                .short('L')
                .long("bait-len")
                .default_value("120")
                .value_parser(clap::value_parser!(usize))
            )
            .arg(Arg::new("hamming-distance")
                .help("Maximum number of mismatches allowed in matching")
                .short('d')
                .long("hamming-distance")
                .default_value("40")
                .value_parser(clap::value_parser!(usize))
            )
            .arg(Arg::new("seed-len")
                .help("Length of seeds in matching")
                .short('g')
                .long("seed-len")
                .default_value("20")
                .value_parser(clap::value_parser!(usize))
            )
            .arg(Arg::new("minimizer-len")
                .help("Length of minimizers in indexing. Must be less or equal to seed-len")
                .short('m')
                .long("minimizer-len")
                .default_value("12")
                .value_parser(clap::value_parser!(usize))
            )
            .arg(Arg::new("cutoff")
                .help("Stop the algorithm when this coverage fraction is reached")
                .short('c')
                .long("cutoff")
                .default_value("1.0")
                .value_parser(clap::value_parser!(f64))
            )
            .arg(Arg::new("randomize")
                .help("Randomize the processing order in the greedy algorithm")
                .short('r')
                .long("randomize")
                .action(ArgAction::SetTrue)
            )
        )
    .subcommand(Command::new("coverage")
        .arg(Arg::new("baits")
            .help("FASTA file containing the baits")
            .long("baits")
            .short('b')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("targets")
            .help("FASTA or FASTQ file, possibly gzipped")
            .long("targets")
            .short('t')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("hamming-distance")
            .help("Maximum number of mismatches allowed in matching")
            .short('d')
            .long("hamming-distance")
            .default_value("40")
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("seed-len")
            .help("Length of seeds in matching")
            .short('g')
            .long("seed-len")
            .default_value("20")
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("minimizer-len")
            .help("Length of minimizers in indexing. Must be less or equal to seed-len")
            .short('m')
            .long("minimizer-len")
            .default_value("12")
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("output")
            .help("Output csv file for the coverage data.")
            .long_help("The output file will have one file per sequence in the target file, containing n comma-separated integers, where n is the length of the sequence. The i'th integer is the number of baits covering the i'th position in the sequence.")
            .long("output")
            .short('o')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        )
    );

    let cli_matches = cli.get_matches();
    match cli_matches.subcommand() {
        Some(("design", sub_matches)) => {
            let infile: &PathBuf = sub_matches.get_one("targets").unwrap();
            let outfile: &PathBuf = sub_matches.get_one("output").unwrap();
            let L: usize = *sub_matches.get_one("bait-length").unwrap();
            let d: usize = *sub_matches.get_one("hamming-distance").unwrap();
            let g: usize = *sub_matches.get_one("seed-len").unwrap();
            let m: usize = *sub_matches.get_one("minimizer-len").unwrap();
            let cutoff: f64 = *sub_matches.get_one("cutoff").unwrap();
            let randomize: bool = sub_matches.get_flag("randomize");
            if randomize { // TODO
                std::unimplemented!("Randomization not implemented yet");
            }
        
            let reader = DynamicFastXReader::from_file(infile).unwrap();
            let mut writer = std::io::BufWriter::new(std::fs::File::create(outfile).unwrap()); // Let's open this right away to crash early if there's a problem
        
            info!("Reading sequences from {}", infile.display());
            let seq_db = Box::new(reader.into_db().unwrap());
        
            info!("Indexing the sequences");
            let index = minimizer_index::MinimizerIndex::new(&seq_db, g, m);
        
            info!("Designing baits");
            design::run_algorithm(&seq_db, &index, L, d, g, cutoff, &mut writer);
        }
        Some(("coverage", sub_matches)) => {
            let targetfile: &PathBuf = sub_matches.get_one("targets").unwrap();
            let outfile: &PathBuf = sub_matches.get_one("output").unwrap();
            let baitfile: &PathBuf = sub_matches.get_one("baits").unwrap();
            let d: usize = *sub_matches.get_one("hamming-distance").unwrap();
            let g: usize = *sub_matches.get_one("seed-len").unwrap();
            let m: usize = *sub_matches.get_one("minimizer-len").unwrap();
            compute_coverage(targetfile, baitfile, outfile, d, g, m);
        }
        _ => {
            log::error!("Unknown subcommand");
        }
    }

    

}
