use std::collections::HashSet;

use jseqio::{seq_db::SeqDB, writer::DynamicFastXWriter};
use log::info;
use crate::minimizer_index::MinimizerIndex;

fn hamming(a: &[u8], b: &[u8]) -> usize{
    assert_eq!(a.len(), b.len());
    let mut dist = 0;
    for (x,y) in a.iter().zip(b.iter()){
        if x != y{
            dist += 1;
        }
    }
    dist
}

fn mark_all_that_are_covered_by(bait: &[u8], cover_marks: &mut Vec<Vec<bool>>, index: &MinimizerIndex, db: &SeqDB, hamming_distance: usize, k: usize){
    let mut align_starts = HashSet::<(usize,usize)>::new(); // Pairs (seq_id, seq_pos)
    for (bait_pos, kmer) in bait.windows(k).enumerate(){
        for (seq_id, seq_pos) in index.lookup(kmer){
            
            let align_start = seq_pos as isize - bait_pos as isize;
            if align_start >= 0 && (align_start + bait.len() as isize) <= db.get(seq_id).seq.len() as isize{
                // Is within bounds
                align_starts.insert((seq_id, align_start as usize));
            }
        }
    }

    for (seq_id, seq_pos) in align_starts{
        if hamming(bait, &db.get(seq_id).seq[seq_pos..seq_pos+bait.len()]) <= hamming_distance{
            for i in seq_pos..seq_pos+bait.len(){
                cover_marks[seq_id][i] = true; // This is within bounds because it was checked above
            }
        }
    }
}

pub fn run_algorithm(db: &SeqDB, index: &MinimizerIndex, bait_len: usize, hamming_distance: usize, k: usize, fasta_out: &mut impl std::io::Write){

    // Initialize the cover marks to falses. False means not covered.
    let mut cover_marks = Vec::<Vec::<bool>>::new();
    for rec in db.iter(){
        cover_marks.push(vec![false; rec.seq.len()]);
    }

    let mut n_baits = 0_usize;
    for (seq_id, rec) in db.iter().enumerate(){
        let mut prev_end = 0_usize;

        // First first position in cover marks that is not yet covered
        while let Some(mut bait_start) = cover_marks[seq_id][prev_end..].iter().position(|b| !*b){
            let mut bait_start = prev_end + bait_start;
            let mut bait_end = bait_start + bait_len;
            if bait_end > rec.seq.len(){
                let excess = bait_end - rec.seq.len();
                if (bait_start as i64 - excess as i64) < 0{
                    panic!("Sequence is shorter than bait length"); // TODO: handle
                }

                bait_start -= excess;
                bait_end -= excess;
            }
            let bait = &rec.seq[bait_start..bait_end];
            mark_all_that_are_covered_by(bait, &mut cover_marks, index, db, hamming_distance, k);

            n_baits += 1;
            prev_end = bait_end;

            fasta_out.write_all(format!(">{}\n", n_baits).as_bytes()).unwrap();
            fasta_out.write_all(bait).unwrap();
            fasta_out.write_all(b"\n").unwrap();
        }
    }

    info!("Selected {} baits", n_baits);
}