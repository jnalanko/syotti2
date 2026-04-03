use jseqio::seq_db::SeqDB;
use log::info;
use crate::minimizer_index::MinimizerIndex;


// Returns the number of new bases covered, and updates per-sequence coverage counts.
fn mark_all_that_are_covered_by(bait: &[u8], cover_marks: &mut Vec<Vec<bool>>, n_covered_by_seq: &mut Vec<usize>, index: &MinimizerIndex, db: &SeqDB, hamming_distance: usize) -> usize{
    let align_starts = index.get_exact_alignment_candidates(bait);
    let mut new_covered_bases = 0_usize;
    for (seq_id, seq_pos) in align_starts{
        if syotti2::hamming_distance_not_matching_N(bait, &db.get(seq_id).seq[seq_pos..seq_pos+bait.len()]) <= hamming_distance{
            for i in seq_pos..seq_pos+bait.len(){
                if !cover_marks[seq_id][i] {
                    new_covered_bases += 1;
                    n_covered_by_seq[seq_id] += 1;
                    cover_marks[seq_id][i] = true; // This is within bounds because it was checked above
                }
            }
        }
    }

    new_covered_bases
}

pub fn run_algorithm(db: &SeqDB, index: &MinimizerIndex, bait_len: usize, hamming_distance: usize, cutoff: f64, require_cutoff_for_every_sequence: bool, fasta_out: &mut impl std::io::Write){

    // Initialize the cover marks to falses. False means not covered.
    let mut cover_marks = Vec::<Vec::<bool>>::new();
    let mut n_covered_by_seq = Vec::<usize>::new();
    let mut total_seq_len = 0_usize;
    for rec in db.iter(){
        cover_marks.push(vec![false; rec.seq.len()]);
        n_covered_by_seq.push(0);
        total_seq_len += rec.seq.len();
    }

    let mut total_covered = 0_usize;
    let mut n_seqs_with_cutoff_reached = 0_usize;
    let mut cutoff_reached_by_seq = vec![false; db.sequence_count()];

    let mut n_baits = 0_usize;
    'outer: for (seq_id, rec) in db.iter().enumerate(){
        let mut prev_end = 0_usize;

        // Find first position in cover marks that is not yet covered
        while let Some(bait_start) = cover_marks[seq_id][prev_end..].iter().position(|b| !*b){
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
            total_covered += mark_all_that_are_covered_by(bait, &mut cover_marks, &mut n_covered_by_seq, index, db, hamming_distance);
            total_covered += mark_all_that_are_covered_by(&jseqio::reverse_complement(bait), &mut cover_marks, &mut n_covered_by_seq, index, db, hamming_distance);

            // Update per-sequence cutoff tracking
            for i in 0..db.sequence_count() {
                if !cutoff_reached_by_seq[i] && n_covered_by_seq[i] as f64 >= db.get(i).seq.len() as f64 * cutoff {
                    cutoff_reached_by_seq[i] = true;
                    n_seqs_with_cutoff_reached += 1;
                }
            }

            n_baits += 1;
            prev_end = bait_end;

            fasta_out.write_all(format!(">{}\n", n_baits).as_bytes()).unwrap();
            fasta_out.write_all(bait).unwrap();
            fasta_out.write_all(b"\n").unwrap();

            let cutoff_reached = if require_cutoff_for_every_sequence {
                n_seqs_with_cutoff_reached == db.sequence_count()
            } else {
                (total_covered as f64) / (total_seq_len as f64) >= cutoff
            };

            if cutoff_reached {
                info!("Reached coverage cutoff of {}% at {} baits", cutoff*100.0, n_baits);
                break 'outer;
            }
        }
    }

    info!("Selected {} baits", n_baits);
}

#[cfg(test)]
mod tests{

    use super::*;

    #[test]
    fn test_reverse_complement(){
        // C++ test: get_rc("ATGNAC") == "GTNCAT" (reverse complement of N is N)
        let s = b"ATGNAC";
        let rc = jseqio::reverse_complement(s);
        assert_eq!(rc, b"GTNCAT");
    }

    #[allow(non_snake_case)]
    #[test]
    fn test_hamming_distance_not_matching_N(){
        let s = b"AACCGGTTNN";
        let t = b"ATCTGTTANA";
        assert_eq!(syotti2::hamming_distance_not_matching_N(s,t), 6);
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
        run_algorithm(&db, &index, bait_length, d, 1.0, false, &mut fasta_out);

        let reader = jseqio::reader::DynamicFastXReader::new(std::io::Cursor::new(fasta_out)).unwrap();
        let bait_db = reader.into_db().unwrap();
        let baits = bait_db.iter().map(|r| r.seq).collect::<Vec<&[u8]>>();

        for b in baits.iter(){
            println!("{}", String::from_utf8_lossy(b));
        }

        assert_eq!(baits, expected_baits);

    }

    #[test]
    fn test_require_cutoff_for_every_sequence(){
        // Two sequences with no sequence or RC similarity.
        // With cutoff=0.5 and the default (flag=false), covering seq1 fully satisfies
        // the global 50% threshold (10/20 bases) and the algorithm stops — seq2 is untouched.
        // With flag=true, seq2 must also individually reach 50% before stopping.
        let d = 0;
        let g = 5;
        let bait_length = 5;
        let cutoff = 0.5;

        let mut db = SeqDB::new();
        db.push_record(jseqio::record::RefRecord{seq: b"AAAAAAAAAA", head: b"", qual: None}); // 10 A's
        db.push_record(jseqio::record::RefRecord{seq: b"CCCCCCCCCC", head: b"", qual: None}); // 10 C's; RC is GGGGGGGGGG, no overlap with seq1

        let index = MinimizerIndex::new(&db, g, 1);

        // Without the flag: one bait covers all of seq1 (10/20 = 50% global) → stop.
        let mut out = Vec::<u8>::new();
        run_algorithm(&db, &index, bait_length, d, cutoff, false, &mut out);
        let bait_db = jseqio::reader::DynamicFastXReader::new(std::io::Cursor::new(out)).unwrap().into_db().unwrap();
        assert_eq!(bait_db.sequence_count(), 1);
        assert_eq!(bait_db.get(0).seq, b"AAAAA");

        // With the flag: seq1 and seq2 must each individually reach 50%, requiring one bait each.
        let mut out = Vec::<u8>::new();
        run_algorithm(&db, &index, bait_length, d, cutoff, true, &mut out);
        let bait_db = jseqio::reader::DynamicFastXReader::new(std::io::Cursor::new(out)).unwrap().into_db().unwrap();
        assert_eq!(bait_db.sequence_count(), 2);
        assert_eq!(bait_db.get(0).seq, b"AAAAA");
        assert_eq!(bait_db.get(1).seq, b"CCCCC");
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