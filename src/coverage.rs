use std::path::Path;
use jseqio::reader::*;
use jseqio::seq_db::SeqDB;
use std::io::Write;
use crate::minimizer_index::MinimizerIndex;

pub fn compute_coverage(targets_db: &SeqDB, bait_db: &SeqDB, out: &mut impl Write, d: usize, g: usize, m: usize) {
    
    let index = MinimizerIndex::new(&targets_db, g, m);

    // Initialize coverage depth vectors
    // coverages[i][j] = depth in target i at position j
    let mut coverages = Vec::<Vec::<u32>>::new(); // 32 bits ought to be enough for anyone

    for bait in bait_db.iter() {
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
            write!(out, "{}", x).unwrap();
        }
    }

}

#[cfg(test)]
mod tests{

    #[test]
    fn coverage_basic_test(){
        
    }
}