mod design;
mod minimizer_index;
mod coverage;

use std::path::PathBuf;

use clap::ArgAction;
use coverage::compute_coverage;
use coverage::into_moving_average;
use jseqio::reader::*;
use jseqio::writer::*;
use jseqio::record::*;
use clap::{Command, Arg};

use log::{info, error};

use jseqio::seq_db::SeqDB;

fn concat_into_rows<T: Copy>(vecs: Vec<Vec<T>>, n_rows: usize) -> Vec<Vec<T>> {
    let concat = vecs.into_iter().reduce(|mut acc, v| {acc.extend(v); acc}).unwrap();
    let row_len = concat.len().div_ceil(n_rows);
    let rows: Vec<Vec<T>> = concat.chunks(row_len).map(|v| v.to_vec()).collect();
    rows
}

fn pad_to_equal_length<T: Copy>(v: &mut Vec<Vec<T>>, pad_element: T) {
    let max_len = v.iter().fold(0_usize, |acc, row| acc.max(row.len()));
    for row in v.iter_mut() {
        while row.len() < max_len {
            row.push(pad_element);
        }
    }
}

fn write_picture<T: Into<f64> + Copy + Clone>(mut coverages: Vec<Vec<T>>, pad_element: T, image_outpath: &std::path::Path) {
    pad_to_equal_length(&mut coverages, pad_element);
    let mut image_out = std::fs::File::create(image_outpath).unwrap();
    coverage::write_as_png(coverages, &mut image_out);
}

fn main() {
    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = Command::new("syotti2")
        .about("Bait design for targeted sequencing")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg_required_else_help(true)
        .subcommand_required(true)
        .arg(Arg::new("threads")
            .help("Number of threads to use")
            .long("threads")
            .short('j')
            .default_value("8")
            .global(true)
            .value_parser(clap::value_parser!(usize))
        )
        .subcommand(Command::new("design")
            .arg_required_else_help(true)
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
        .arg_required_else_help(true)
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
        .arg(Arg::new("resolution")
            .help("Computes a moving average over the coverage vectors and samples those at evenly spaced intervals to get this many data points per sequence.")
            .short('r')
            .long("resolution")
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
        .arg(Arg::new("coverage-out")
            .help("Output csv file for the coverage data.")
            .long_help("The output file will have one file per sequence in the target file, containing n comma-separated integers, where n is the length of the sequence. The i'th integer is the number of baits covering the i'th position in the sequence.")
            .long("coverage-out")
            .short('o')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("mismatch-out")
            .help("Output csv file for the mismatch data.")
            .long_help("Output csv file that is like the output file, but instead of coverage, it reports the number of mismatches in the best bait covering each position, or an asterisk if there is no bait coovering the position.")
            .long("mismatch-out")
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("coverage-out-picture")
            .help("Output png file for the coverage data.")
            .long_help("Visualize the coverage depth into a png image with this path. Beware: the height is equal to the number of targets, and the width is equal to the length of the longest target. Use --resolution to reduce the width to avoid having an enormous image. To reduce the height, check out concat-into-rows.")
            .long("coverage-out-picture")
            .short('p')
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("concat-into-rows")
            .help("Let the total number of bases be n, and let this parameter be r. Concatenates the coverage vectors and splits the concatenation into ceil(n/r) rows.")
            .long("concat-into-rows")
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("variable-resolution")
            .help("Relevant if --resolution is given. In enabled, makes it so that the resolution of each coverage vector is proportional to its length, such that the longest target has resolution equal to the value passed to --resolution. This is especially useful with --coverage-out-picture for better visualization.")
            .long("variable-resolution")
            .action(clap::ArgAction::SetTrue)
        )
    );

    
    let cli_matches = cli.get_matches();

    let n_threads = *cli_matches.get_one::<usize>("threads").unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(n_threads).build_global().unwrap();

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
        
            info!("Indexing the sequences"); // TODO: move to inside design
            let index = minimizer_index::MinimizerIndex::new(&seq_db, g, m);
        
            info!("Designing baits");
            design::run_algorithm(&seq_db, &index, L, d, g, cutoff, &mut writer);
        }
        Some(("coverage", sub_matches)) => {
            let targetfile: &PathBuf = sub_matches.get_one("targets").unwrap();
            let outfile: &PathBuf = sub_matches.get_one("coverage-out").unwrap();
            let baitfile: &PathBuf = sub_matches.get_one("baits").unwrap();
            let d: usize = *sub_matches.get_one("hamming-distance").unwrap();
            let g: usize = *sub_matches.get_one("seed-len").unwrap();
            let m: usize = *sub_matches.get_one("minimizer-len").unwrap();
            let resolution = sub_matches.get_one::<usize>("resolution");
            let variable_resolution = sub_matches.get_flag("variable-resolution");
            let concat_rows = sub_matches.get_one::<usize>("concat-into-row");
            let mismatch_outfile: Option<&PathBuf> = sub_matches.get_one("mismatch-out");
            let coverage_picture_outfile: Option<&PathBuf> = sub_matches.get_one("coverage-out-picture");

            let mismatch_out = mismatch_outfile.map(|path| std::io::BufWriter::new(std::fs::File::create(path).unwrap()));
            let mut out = std::io::BufWriter::new(std::fs::File::create(outfile).unwrap());
            let bait_db = DynamicFastXReader::from_file(&baitfile).unwrap().into_db().unwrap(); // TODO: print info log
            let targets_db = DynamicFastXReader::from_file(&targetfile).unwrap().into_db().unwrap();
        
            let (mut coverages, mut mismatches) = compute_coverage(&targets_db, &bait_db, d, g, m);
            if let Some(n_rows) = concat_rows {
                coverages = concat_into_rows(coverages, *n_rows);
                mismatches = concat_into_rows(mismatches, *n_rows);
            }

            // Write coverage numbers, and possibly a picture
            if let Some(reso) = resolution {
                let cov_averages = coverage::into_resolution(coverages, *reso, variable_resolution);
                coverage::write_as_csv(&cov_averages, &mut out, |x| format!("{}", x));
                if let Some(image_outpath) = coverage_picture_outfile {
                    write_picture(cov_averages, -1.0, image_outpath);
                }
            } else { // Write coverages as-is without smoothing to resolution
                coverage::write_as_csv(&coverages, &mut out, |x| format!("{}", x));

                // Reinterpert as i32 to be able to use -1 as a padding value:
                let icoverages: Vec<Vec<i32>> = coverages.into_iter().map(|v| v.into_iter().map(|x| i32::try_from(x).unwrap()).collect()).collect();
                if let Some(image_outpath) = coverage_picture_outfile {
                    write_picture(icoverages, -1, image_outpath);
                }
            }

            // Write mismatches if needed
            if let Some(mut mismatch_out) = mismatch_out{
                if let Some(_) = resolution {
                    unimplemented!("Resolution not implemented for mismatch output");
                } else{
                    let formatter = 
                    |x: &u32| if *x == u32::MAX { 
                        "*".to_string() 
                    } else {
                        format!("{}", x)
                    };
                    coverage::write_as_csv(&mismatches, &mut mismatch_out, formatter);
                }
            }

        }
        _ => {
            log::error!("Unknown subcommand");
        }
    }

    

}
