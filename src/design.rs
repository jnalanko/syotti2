use jseqio::seq_db::SeqDB;
use log::info;
use crate::minimizer_index::MinimizerIndex;

struct CoverageState {
    cover_marks: Vec<Vec<bool>>,
    n_covered_by_seq: Vec<usize>,
    cutoff_reached_by_seq: Vec<bool>,
    n_seqs_with_cutoff_reached: usize,
    cutoff: f64,
}

// Scans forward from `from` and returns the start position for a bait covering the next uncovered
// position, or None if everything from `from` onward is already covered.
fn find_bait_start(from: usize, seq_len: usize, cover_marks: &[bool], bait_len: usize, overhang: usize) -> Option<usize> {
    assert!(seq_len >= bait_len, "Sequence is shorter than bait length"); // TODO: handle
    let gap_start = cover_marks[from..].iter().position(|b| !*b).map(|rel| from + rel)?;
    let p = if overhang > 0 {
        // Scan at most 2*overhang + bait_len ahead for the first already-covered position.
        // If the gap fits within the full coverage width, center the bait in it so
        // both overhangs reach the covered regions on either side.
        let scan_end = (gap_start + 2 * overhang + bait_len).min(seq_len);
        let v = cover_marks[gap_start..scan_end].iter().position(|b| *b).map(|rel| gap_start + rel);
        if let Some(v) = v {
            (gap_start + v).saturating_sub(bait_len) / 2
        } else {
            gap_start + overhang
        }
    } else {
        gap_start
    };
    Some(p.min(seq_len - bait_len))
}

// Marks positions [mark_start, mark_end) in seq_id as covered, updating all coverage state.
// Returns the number of newly covered bases.
fn mark_range(seq_id: usize, mark_start: usize, mark_end: usize, seq_len: usize, cov: &mut CoverageState) -> usize {
    let mut new_covered_bases = 0_usize;
    let length_threshold = (seq_len as f64 * cov.cutoff).ceil() as usize;
    for i in mark_start..mark_end {
        if !cov.cover_marks[seq_id][i] {
            new_covered_bases += 1;
            cov.n_covered_by_seq[seq_id] += 1;
            cov.cover_marks[seq_id][i] = true;
            if !cov.cutoff_reached_by_seq[seq_id] && cov.n_covered_by_seq[seq_id] >= length_threshold {
                cov.cutoff_reached_by_seq[seq_id] = true;
                cov.n_seqs_with_cutoff_reached += 1;
            }
        }
    }
    new_covered_bases
}

// Returns whether bait and target share an exactly matching window of length seed_len.
// Runs in O(bait.len()) using a sliding mismatch count.
fn has_exact_seed_window(bait: &[u8], target: &[u8], seed_len: usize) -> bool {
    let mut mismatches = 0usize;
    for i in 0..seed_len {
        if bait[i] != target[i] { mismatches += 1; }
    }
    if mismatches == 0 { return true; }
    for i in seed_len..bait.len() {
        if bait[i - seed_len] != target[i - seed_len] { mismatches -= 1; }
        if bait[i] != target[i] { mismatches += 1; }
        if mismatches == 0 { return true; }
    }
    false
}

// Returns the number of new bases covered, and updates all coverage state including per-sequence cutoff tracking.
fn mark_all_that_are_covered_by(bait: &[u8], cov: &mut CoverageState, index: &MinimizerIndex, db: &SeqDB, hamming_distance: usize, overhang: usize, seed_len: usize) -> usize{
    let align_starts = index.get_exact_alignment_candidates(bait);
    let mut new_covered_bases = 0_usize;
    for (seq_id, seq_pos) in align_starts{
        let target = &db.get(seq_id).seq[seq_pos..seq_pos+bait.len()];
        if syotti2::hamming_distance_not_matching_N(bait, target) <= hamming_distance
            && has_exact_seed_window(bait, target, seed_len)
        {
            let seq_len = db.get(seq_id).seq.len();
            let mark_start = seq_pos.saturating_sub(overhang);
            let mark_end = (seq_pos + bait.len() + overhang).min(seq_len);
            new_covered_bases += mark_range(seq_id, mark_start, mark_end, seq_len, cov);
        }
    }

    new_covered_bases
}

pub fn run_algorithm(db: &SeqDB, index: &MinimizerIndex, bait_len: usize, hamming_distance: usize, cutoff: f64, require_cutoff_for_every_sequence: bool, overhang: usize, seed_len: usize, fasta_out: &mut impl std::io::Write){

    // Initialize the cover marks to falses. False means not covered.
    let mut total_seq_len = 0_usize;
    let mut cover_marks = Vec::<Vec::<bool>>::new();
    let mut n_covered_by_seq = Vec::<usize>::new();
    for rec in db.iter(){
        cover_marks.push(vec![false; rec.seq.len()]);
        n_covered_by_seq.push(0);
        total_seq_len += rec.seq.len();
    }

    let mut cov = CoverageState {
        cover_marks,
        n_covered_by_seq,
        cutoff_reached_by_seq: vec![false; db.sequence_count()],
        n_seqs_with_cutoff_reached: 0,
        cutoff,
    };

    let mut total_covered = 0_usize;

    let mut n_baits = 0_usize;
    for (seq_id, rec) in db.iter().enumerate(){
        let mut prev_end = 0_usize;

        while let Some(bait_start) = find_bait_start(prev_end, rec.seq.len(), &cov.cover_marks[seq_id], bait_len, overhang) {
            let bait_end = bait_start + bait_len;
            let bait = &rec.seq[bait_start..bait_end];

            // Always mark the source region (plus overhang) as covered, even if the bait
            // doesn't align back to itself (e.g. due to Ns not matching in the index).
            let mark_start = bait_start.saturating_sub(overhang);
            let mark_end = (bait_end + overhang).min(rec.seq.len());
            total_covered += mark_range(seq_id, mark_start, mark_end, rec.seq.len(), &mut cov);

            total_covered += mark_all_that_are_covered_by(bait, &mut cov, index, db, hamming_distance, overhang, seed_len);
            total_covered += mark_all_that_are_covered_by(&jseqio::reverse_complement(bait), &mut cov, index, db, hamming_distance, overhang, seed_len);

            n_baits += 1;
            prev_end = bait_end;

            fasta_out.write_all(format!(">{}\n", n_baits).as_bytes()).unwrap();
            fasta_out.write_all(bait).unwrap();
            fasta_out.write_all(b"\n").unwrap();

            let cutoff_reached = if require_cutoff_for_every_sequence {
                cov.n_seqs_with_cutoff_reached == db.sequence_count()
            } else {
                (total_covered as f64) / (total_seq_len as f64) >= cutoff
            };

            if cutoff_reached {
                info!("Reached coverage cutoff of {}% at {} baits", cutoff*100.0, n_baits);
                info!("Selected {} baits", n_baits);
                return;
            }
        }
    }
    
    // All sequences have been processed, so coverage must be 100%.
    panic!("This part of the code should never be reached");
}

#[cfg(test)]
mod tests{

    use super::*;

    // Helpers: build a cover_marks slice from a string where '.' = uncovered, 'X' = covered.
    fn marks(s: &str) -> Vec<bool> {
        s.chars().map(|c| c == 'X').collect()
    }

    #[test]
    fn find_bait_start_no_overhang_places_at_u() {
        // Without overhang the bait always starts at the first uncovered position.
        let m = marks("XX....XX");
        assert_eq!(find_bait_start(2, 8, &m, 3, 0), Some(2));
    }

    #[test]
    fn find_bait_start_no_overhang_small_gap_still_places_at_u() {
        // Gap [2,4) is shorter than bait_len=5, but with no overhang we don't center — bait starts at u=2.
        let m = marks("XX..XXXX");
        assert_eq!(find_bait_start(2, 8, &m, 5, 0), Some(2));
    }

    #[test]
    fn find_bait_start_no_overhang_clamps_to_end() {
        // u=6, bait_len=4, seq_len=8 → p=6 would overshoot, clamped to 4.
        let m = marks("XXXXXX..");
        assert_eq!(find_bait_start(6, 8, &m, 4, 0), Some(4));
    }

    #[test]
    fn find_bait_start_overhang_large_gap_shifts_right() {
        // Gap is larger than 2w+bait_len so no v is found; p = u + overhang.
        // seq: 20 uncovered bases, w=2, bait_len=5 → 2w+bait_len=9, gap=20 > 9 → p = 0+2 = 2.
        let m = marks("....................");
        assert_eq!(find_bait_start(0, 20, &m, 5, 2), Some(2));
    }

    #[test]
    fn find_bait_start_overhang_small_gap_centers() {
        // Gap [3,8), covered elsewhere. w=2, bait_len=5 → 2w+bait_len=9.
        // Gap length=5 < 9, so v=8 is found. p = (3+8-5)/2 = 3.
        let mut m = marks("XXX.....XXXXXXXXX");
        m[8] = true;
        assert_eq!(find_bait_start(3, m.len(), &m, 5, 2), Some(3));
    }

    #[test]
    fn find_bait_start_overhang_tiny_gap_centers() {
        // Gap [5,7) (length 2), w=3, bait_len=5 → p = (5+7-5)/2 = 3 (clamped fine).
        let m = marks("XXXXX..XXXXXXXXXX");
        assert_eq!(find_bait_start(5, m.len(), &m, 5, 3), Some(3));
    }

    #[test]
    fn find_bait_start_overhang_exact_span() {
        // Gap [2,9), w=2, bait_len=3 → 2w+bait_len=7 = gap length exactly.
        // v=9 is found; p = (2+9-3)/2 = 4 = u+w, so the left overhang just reaches u and the right just reaches v.
        let m = marks("XX.......XXXX");
        assert_eq!(find_bait_start(2, m.len(), &m, 3, 2), Some(4));

        let m = marks("XXX.......XXXX");
        assert_eq!(find_bait_start(0, m.len(), &m, 3, 2), Some(5));
    }

    #[test]
    fn find_bait_start_overhang_clamped_to_seq_end() {
        // Large gap near end of seq: u=8, w=3, bait_len=5, seq_len=12 → u+w=11 > 12-5=7, clamped to 7.
        let m = marks("XXXXXXXX....");
        assert_eq!(find_bait_start(8, 12, &m, 5, 3), Some(7));
    }

    #[test]
    fn find_bait_start_returns_none_when_fully_covered() {
        let m = marks("XXXXXXXX");
        assert_eq!(find_bait_start(0, 8, &m, 3, 0), None);
        assert_eq!(find_bait_start(0, 8, &m, 3, 2), None);
    }

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
        run_algorithm(&db, &index, bait_length, d, 1.0, false, 0, g, &mut fasta_out);

        let reader = jseqio::reader::DynamicFastXReader::new(std::io::Cursor::new(fasta_out)).unwrap();
        let bait_db = reader.into_db().unwrap();
        let baits = bait_db.iter().map(|r| r.seq).collect::<Vec<&[u8]>>();

        for b in baits.iter(){
            println!("{}", String::from_utf8_lossy(b));
        }

        assert_eq!(baits, expected_baits);

    }

    #[test]
    fn test_overhang() {
        let d = 0;
        let g = 8;
        let bait_length = 8;
        let overhang = 2;

        let mut db = SeqDB::new();

        //                                                                   ********    ********
        //                                                ********    ********                 ********
        //                                              01234567890123456789012345678
        db.push_record(jseqio::record::RefRecord{seq: b"ACGTATTCGTGATTCTGTAGTCAGCGTAC", head: b"", qual: None});
        db.push_record(jseqio::record::RefRecord{seq: b"ACGTATTCGTGATTCTGTAGTCAGCGTACAAATTCTGTATGCTAGCA", head: b"", qual: None});

        let index = MinimizerIndex::new(&db, g, 1);
        let mut fasta_out = Vec::<u8>::new();
        run_algorithm(&db, &index, bait_length, d, 1.0, false, overhang, g, &mut fasta_out);

        let bait_db = jseqio::reader::DynamicFastXReader::new(std::io::Cursor::new(fasta_out)).unwrap().into_db().unwrap();
        let baits = bait_db.iter().map(|r| r.seq).collect::<Vec<&[u8]>>();

        assert_eq!(baits, vec![b"GTATTCGT" as &[u8], b"CTGTAGTC", b"CAGCGTAC", b"TCTGTATG", b"TGCTAGCA"]);
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
        run_algorithm(&db, &index, bait_length, d, cutoff, false, 0, g, &mut out);
        let bait_db = jseqio::reader::DynamicFastXReader::new(std::io::Cursor::new(out)).unwrap().into_db().unwrap();
        assert_eq!(bait_db.sequence_count(), 1);
        assert_eq!(bait_db.get(0).seq, b"AAAAA");

        // With the flag: seq1 and seq2 must each individually reach 50%, requiring one bait each.
        let mut out = Vec::<u8>::new();
        run_algorithm(&db, &index, bait_length, d, cutoff, true, 0, g, &mut out);
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