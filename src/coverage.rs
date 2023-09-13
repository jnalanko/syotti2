use std::path::Path;
use jseqio::{reader::*, reverse_complement};
use jseqio::seq_db::SeqDB;
use std::io::Write;
use crate::minimizer_index::MinimizerIndex;

fn update_coverage(coverages: &mut Vec<Vec<u32>>, bait: &[u8], index: &MinimizerIndex, targets_db: &SeqDB, hamming_distance: usize, k: usize){
    let candidates = index.get_exact_alignment_candidates(bait);

    for (target_id, target_start) in candidates{
        let target = &targets_db.get(target_id).seq[target_start .. target_start + bait.len()];
        if syotti2::hamming_distance_not_matching_N(bait, target) <= hamming_distance{
            for i in 0..bait.len(){
                // Add 1 to the coverage using saturating add so we don't overflow
                coverages[target_id][target_start + i] = coverages[target_id][target_start + i].saturating_add(1);
            }
        }
    }
}

pub fn compute_coverage(targets_db: &SeqDB, bait_db: &SeqDB, out: &mut impl Write, d: usize, g: usize, m: usize) {
    
    let index = MinimizerIndex::new(&targets_db, g, m);

    // Initialize coverage depth vectors
    // coverages[i][j] = depth in target i at position j
    let mut coverages = Vec::<Vec::<u32>>::new(); // 32 bits ought to be enough for anyone
    for i in 0..targets_db.sequence_count(){
        coverages.push(vec![0; targets_db.get(i).seq.len()]);
    }

    for bait in bait_db.iter() {
        update_coverage(&mut coverages, bait.seq, &index, targets_db, d, g);
        update_coverage(&mut coverages, reverse_complement(bait.seq).as_slice(), &index, targets_db, d, g);
    }

    for v in coverages.iter() {
        // Write v as a line of comma-separated values
        for (i, x) in v.iter().enumerate(){
            if i > 0 {
                write!(out, ",").unwrap();
                
            }
            write!(out, "{}", x).unwrap();
        }
        write!(out, "\n").unwrap();
    }

}

#[cfg(test)]
mod tests{
    use super::compute_coverage;
    use super::*;

    #[test]
    fn coverage_basic_test(){
        let targets = vec![
        b"ACTCGTAGCACGCTATCTATCGATCGTAGCTAGCTACCACATGC".to_vec(),
                 b"ACGCTATCTATCGATCGTAGCTAGCTAC".to_vec(),
                   b"GCTATCTATCGATCGTAGCTAGCTAC".to_vec()];

        let mut target_db = SeqDB::new();
        for seq in targets.iter() {target_db.push_seq(seq);};

        let baits = vec![
            b"TCTCGTAGCA".to_vec(),
            b"TGCTATCTAT".to_vec(),
            b"TGATCGTAGC".to_vec(),
            b"AAGCTACCAC".to_vec(),
            b"GTGGTAGCTT".to_vec() // Reverse complement of the previous bait
        ]; // These tile the first sequence with at most 1 mismatch per bait, except for the last 4 bases

        let mut bait_db = SeqDB::new();
        for seq in baits.iter() {bait_db.push_seq(seq);};

        let mut t0_answer = vec![0; targets[0].len()];
        for i in 0..40 {t0_answer[i] += 1;}
        for i in 30..40 {t0_answer[i] += 1;}

        let mut t1_answer = vec![0; targets[1].len()];
        for i in 1..21 {t1_answer[i] += 1;}

        let mut t2_answer = vec![0; targets[2].len()];
        for i in 9..19 {t2_answer[i] += 1;}

        let t0_csv = t0_answer.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        let t1_csv = t1_answer.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        let t2_csv = t2_answer.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        let answer = format!("{}\n{}\n{}\n", t0_csv, t1_csv, t2_csv);

        let mut csv_out = Vec::<u8>::new();
        compute_coverage(&target_db, &bait_db, &mut csv_out, 1, 5, 3);

        eprintln!("csv_out:\n{}", String::from_utf8(csv_out.clone()).unwrap());
        eprintln!("answer:\n{}", answer);

        assert_eq!(csv_out, answer.as_bytes());
    }
}

