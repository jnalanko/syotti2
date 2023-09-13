use std::path::Path;
use jseqio::reader::*;
use std::io::Write;
use crate::minimizer_index::MinimizerIndex;

pub fn compute_coverage(targetfile: &Path, baitfile: &Path, outfile: &Path, d: usize, g: usize, m: usize) {
    
    let mut out = std::io::BufWriter::new(std::fs::File::create(outfile).unwrap());
    let mut bait_reader = DynamicFastXReader::from_file(&baitfile).unwrap();
    let targets_db = DynamicFastXReader::from_file(&targetfile).unwrap().into_db().unwrap();
    let index = MinimizerIndex::new(&targets_db, g, m);

    // Initialize coverage depth vectors
    // coverages[i][j] = depth in target i at position j
    let mut coverages = Vec::<Vec::<u32>>::new(); // 32 bits ought to be enough for anyone

    while let Some(bait) = bait_reader.read_next().unwrap(){
        for (target_id, target_start) in index.get_exact_alignment_candidates(bait.seq){
            let target = &targets_db.get(target_id).seq[target_start .. target_start + bait.seq.len()];
            if syotti2::hamming_distance_not_matching_N(bait.seq, target) <= d{
                for i in 0..bait.seq.len(){
                    coverages[target_id][target_start + i] += 1;
                }
            }
        }
    }

    for v in coverages.iter() {
        // Write v as a line of comma-separated values
        for (i, x) in v.iter().enumerate(){
            if i == 0 {
                write!(out, ",").unwrap();
                out.write_all(b",").unwrap();
            }
            write!(&mut out, "{}", x).unwrap();
        }
    }

}