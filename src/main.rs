mod design;
mod minimizer_index;

use std::path::PathBuf;

use clap::ArgAction;
use jseqio::reader::*;
use jseqio::writer::*;
use jseqio::record::*;
use clap::{Command, Arg};

use log::info;

use jseqio::seq_db::SeqDB;

fn main() {
    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = Command::new("syotti2")
        .about("Bait design for targeted sequencing")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg(Arg::new("sequences")
            .help("Input FASTA or FASTQ file, possibly gzipped")
            .long("sequences")
            .short('s')
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
            .help("Length of minimizers in indexing. Must be less of equal to seed-len")
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
        );

    let cli_matches = cli.get_matches();
    let infile: &PathBuf = cli_matches.get_one("sequences").unwrap();
    let outfile: &PathBuf = cli_matches.get_one("output").unwrap();
    let L: usize = *cli_matches.get_one("bait-length").unwrap();
    let d: usize = *cli_matches.get_one("hamming-distance").unwrap();
    let g: usize = *cli_matches.get_one("seed-len").unwrap();
    let m: usize = *cli_matches.get_one("minimizer-len").unwrap();
    let cutoff: f64 = *cli_matches.get_one("cutoff").unwrap();
    let randomize: bool = cli_matches.get_flag("randomize");

    let reader = DynamicFastXReader::from_file(infile).unwrap();
    let writer = DynamicFastXWriter::new_to_file(outfile).unwrap(); // Let's open this right away to crash early if there's a problem

    info!("Reading sequences from {}", infile.display());
    let seq_db = reader.into_db().unwrap();

    info!("Indexing the sequences");
    let index = minimizer_index::MinimizerIndex::new(seq_db, m, g);

    info!("Running the bait design algorithm");
    //design::run_algorithm(&fw_db, &rc_db, L, d, g, cutoff);
    

}
