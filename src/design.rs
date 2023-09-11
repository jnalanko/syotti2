use std::collections::HashSet;

use jseqio::{seq_db::SeqDB, writer::DynamicFastXWriter};
use log::info;
use crate::minimizer_index::MinimizerIndex;

fn hamming_distance_not_matching_N(a: &[u8], b: &[u8]) -> usize{
    assert_eq!(a.len(), b.len());
    let mut dist = 0;
    for (x,y) in a.iter().zip(b.iter()){
        if x != y || (x == &b'N' && y == &b'N'){
            dist += 1;
        }
    }
    dist
}

// Returns the number of new bases covered
fn mark_all_that_are_covered_by(bait: &[u8], cover_marks: &mut Vec<Vec<bool>>, index: &MinimizerIndex, db: &SeqDB, hamming_distance: usize, k: usize) -> usize{
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

    let mut new_covered_bases = 0_usize;
    for (seq_id, seq_pos) in align_starts{
        if hamming_distance_not_matching_N(bait, &db.get(seq_id).seq[seq_pos..seq_pos+bait.len()]) <= hamming_distance{
            for i in seq_pos..seq_pos+bait.len(){
                new_covered_bases += !cover_marks[seq_id][i] as usize;
                cover_marks[seq_id][i] = true; // This is within bounds because it was checked above
            }
        }
    }

    new_covered_bases
}

pub fn run_algorithm(db: &SeqDB, index: &MinimizerIndex, bait_len: usize, hamming_distance: usize, k: usize, cutoff: f64, fasta_out: &mut impl std::io::Write){

    // Initialize the cover marks to falses. False means not covered.
    let mut cover_marks = Vec::<Vec::<bool>>::new();
    let mut total_seq_len = 0_usize;
    for rec in db.iter(){
        cover_marks.push(vec![false; rec.seq.len()]);
        total_seq_len += rec.seq.len();
    }

    let mut total_covered = 0_usize;
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
            total_covered += mark_all_that_are_covered_by(bait, &mut cover_marks, index, db, hamming_distance, k);
            total_covered += mark_all_that_are_covered_by(&jseqio::reverse_complement(bait), &mut cover_marks, index, db, hamming_distance, k);

            n_baits += 1;
            prev_end = bait_end;

            fasta_out.write_all(format!(">{}\n", n_baits).as_bytes()).unwrap();
            fasta_out.write_all(bait).unwrap();
            fasta_out.write_all(b"\n").unwrap();

            if (total_covered as f64) / (total_seq_len as f64) >= cutoff{
                info!("Reached coverage cutoff of {}% at {} baits", cutoff*100.0, n_baits);
                return;
            }
        }
    }

    info!("Selected {} baits", n_baits);
}

#[cfg(test)]
mod tests{
    use jseqio::reader;

    use super::*;
    
    #[test]
    fn test_hamming_distance_not_matching_N(){
        let s = b"AACCGGTTNN";
        let t = b"ATCTGTTANA";
        assert_eq!(hamming_distance_not_matching_N(s,t), 6);
    }

    #[test]
    fn basic_testcase(){
        // Ported from Syotti 1

        let d = 1;
        let g = 2;
        let bait_length = 5;

        let seqs = [
        "AAAAAACCCCCCATATATAGTTTTTTTT",
        "NNNNNNNNNNNNN",
        "AAAAAAAACTATATATGGGGGGTTTTTT", // First again
        "AAAAAAAACTATATATGGGGGGTTTTTT", // RC of the first
        "AANNAANNCCNNATNNATNNTTNNTTNN", // The first but N's such that there is no common 5-mer with 1 mismatch
        "AANNAANNCTNNATNNGGNNGGNNTTNN", // The RC of the first but N's such that there is no common 5-mer with 1 mismatch
        "NNNNNATATATANNNNNNNNN", // Island in the middle should be covered by bait TATAT
        "TACGT", // Unique
        "ACGTA", // RC of above
        ].map(|s| s.as_bytes());

        // Covering the first sequence should happen like this:
        // AAAAA covers the prefix AAAAAAC and by reverse complement the suffix GTTTTTTTT.
        // CCCCC covers CCCCCCA.
        // TATAT covers TATATA (also the last A because of reverse complement ATATA).
        // The Ns don't match to each other so they are all covered separately.
        // The third input sequence is just the reverse complement of the first, so it is automatically covered.
        let expected_baits = ["AAAAA","CCCCC","TATAT", // First sequence
                                        "NNNNN","NNNNN","NNNNN", // Second sequence
                                                                // 3. sequence: already covered
                                                                // 4. sequence: already covered
                                        "AANNA","ANNCC","NNATN","NATNN","TTNNT","NTTNN",  // 5. sequence
                                        "AANNA","ANNCT","NNATN","NGGNN","GGNNT","NTTNN", // 6. sequence
                                        "NNNNN","NNNNN","NNNNN", // 7. sequence
                                        "TACGT", // 8. sequence
                                                // 9. sequence: already covered as RC of 8.
        ].map(|s| s.as_bytes());

        let mut db = SeqDB::new();
        for s in seqs.iter(){
            db.push_record(jseqio::record::RefRecord{seq: s, head: b"", qual: None});
        }

        let index = MinimizerIndex::new(&db, g, 1);
        let mut fasta_out = Vec::<u8>::new();
        run_algorithm(&db, &index, bait_length, d, g, 1.0, &mut fasta_out);

        let reader = jseqio::reader::DynamicFastXReader::new(std::io::Cursor::new(fasta_out)).unwrap();
        let bait_db = reader.into_db().unwrap();
        let baits = bait_db.iter().map(|r| r.seq).collect::<Vec<&[u8]>>();

        for b in baits.iter(){
            println!("{}", String::from_utf8_lossy(b));
        }

        assert_eq!(baits, expected_baits);

    }

    /*
    vector<string> seqs = {"AAAAAACCCCCCATATATAGTTTTTTTT",
                           "NNNNNNNNNNNNN",
                           "AAAAAAAACTATATATGGGGGGTTTTTT", // First again
                           "AAAAAAAACTATATATGGGGGGTTTTTT", // RC of the first
                           "AANNAANNCCNNATNNATNNTTNNTTNN", // The first but N's such that there is no common 5-mer with 1 mismatch
                           "AANNAANNCTNNATNNGGNNGGNNTTNN", // The RC of the first but N's such that there is no common 5-mer with 1 mismatch
                           "NNNNNATATATANNNNNNNNN", // Island in the middle should be covered by bait TATAT
                           "TACGT", // Unique
                           "ACGTA", // RC of above
                           };
    LL d = 1;
    LL g = 2;
    LL bait_length = 5;

    FM_index fmi;
    fmi.construct(seqs);
    FM_NeighborCandidateFunction FM_NCF;
    FM_NCF.init(&fmi, g);
    NeighborFunction FM_NF;
    FM_NF.init(&FM_NCF, &seqs, d, bait_length);

    Greedy G_FM;
    G_FM.init(&FM_NF, &seqs, bait_length, d, g, false, 1);
    Greedy::Result result = G_FM.run();

    // Covering the first sequence should happen like this:
    // AAAAA covers the prefix AAAAAAC and by reverse complement the suffix GTTTTTTTT.
    // CCCCC covers CCCCCCA.
    // TATAT covers TATATA (also the last A because of reverse complement ATATA).
    // The Ns don't match to each other so they are all covered separately.
    // The third input sequence is just the reverse complement of the first, so it is automatically covered.
    vector<string> expected_baits = {"AAAAA","CCCCC","TATAT", // First sequence
                                     "NNNNN","NNNNN","NNNNN", // Second sequence
                                                              // 3. sequence: already covered
                                                              // 4. sequence: already covered
                                     "AANNA","ANNCC","NNATN","NATNN","TTNNT","NTTNN",  // 5. sequence
                                     "AANNA","ANNCT","NNATN","NGGNN","GGNNT","NTTNN", // 6. sequence
                                     "NNNNN","NNNNN","NNNNN", // 7. sequence
                                     "TACGT", // 8. sequence
                                              // 9. sequence: already covered as RC of 8.
                                      };

 */   
}