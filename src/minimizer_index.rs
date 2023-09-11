use std::collections::HashSet;

use rayon::prelude::*;
use indicatif::ProgressIterator;
use jseqio::seq_db::SeqDB;
pub struct MinimizerIndex<'a>{
    seq_storage: &'a jseqio::seq_db::SeqDB,
    mphf: boomphf::Mphf<Kmer>, // Minimal perfect hash function
    locations: Vec<(u32, u32)>,
    bucket_starts: Vec<usize>,
    k: usize, // k-mer length
    m: usize, // Minimizer length
    n_mmers: usize, // Number of distinct m-mers stored in the mphf
}

// Keeps only parts of length at least k
fn split_at_non_ACGT(v: &[u8], k: usize) -> Vec<Vec<u8>>{
    let mut parts = Vec::<Vec::<u8>>::new();
    let mut start = 0 as usize;
    for i in 0 .. v.len() + 1 {
        if i == v.len() || !is_dna(v[i]) { // One-past-the-end of a part
            if i - start >= k {
                parts.push(v[start..i].to_owned());
            }
            start = i+1;
        }
    }
    parts
}

// TODO: UPPER-CASE SEQUENCES
fn is_dna(c: u8) -> bool{
    if c == b'a' || c == b'c' || c == b'g' || c == b't' {
        panic!("ERROR: lower case nucleotides found");
    }
    c == b'A' || c == b'C' || c == b'G' || c == b'T'
}

fn get_minimizer_position(kmer: &[u8], m: usize) -> usize{
    let mut minimizer = &kmer[0 .. 0+m];
    let mut min_pos = 0;
    for j in 1 .. (kmer.len() as i64) - (m as i64) + 1 {
        let j = j as usize;
        if kmer[j..j+m] < *minimizer {
            minimizer = &kmer[j..j+m];
            min_pos = j;
        }
    }
    return min_pos;
}

// The output will be stored to the positions-vector
// Will not return minimizers for k-mers with non-DNA characters
fn get_minimizer_positions(seq: &[u8], positions: &mut Vec<usize>, k: usize, m: usize){

    positions.clear();

    // Store minimizer mappings
    for i in 0 .. (seq.len() as i64) - (k as i64) + 1 {
        let i = i as usize;
        let kmer = &seq[i..i+k];
        if kmer.iter().any(|&c| !is_dna(c)){
            continue;
        }

        let min_pos = i + get_minimizer_position(kmer, m);

        if positions.len() == 0 || positions[positions.len()-1] != min_pos {
            positions.push(min_pos)
        }

    }
}

// The output will be stored to the positions-vector
fn get_minimizer_positions_with_return(seq: &[u8], k: usize, m: usize) -> Vec<usize>{

    let mut positions = Vec::<usize>::new();
    get_minimizer_positions(seq, &mut positions, k, m);
    positions
}

#[derive(Copy, Clone, PartialEq, Eq, Ord, Hash, Debug)]
struct Kmer{
    data: u64
}

impl Kmer{
    fn from_ascii(ascii: &[u8]) -> Result<Self, ()>{
        if ascii.len() > 32{
            return Err(()); // Does not fit in u64
        }
        let mut data = 0 as u64;
        for i in 0 .. ascii.len(){
            data <<= 2;
            data |= match ascii[i]{
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {return Err(());}
            };
        }
        Ok(Self{data})
    }
}

impl std::cmp::PartialOrd for Kmer{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering>{
        Some(self.cmp(other))
    }
}

impl<'a> MinimizerIndex<'a>{

    fn compress_sorted_position_list(mut L: Vec::<(Kmer, u32, u32)>, h: &boomphf::Mphf<Kmer>, n_minimizers: usize) -> (Vec<(u32, u32)>, Vec<usize>){
        
        let mut bucket_sizes: Vec::<usize> = vec![0; n_minimizers]; // Bucket sizes in left-to-right order of buckets
        for (seq, _, _) in L.iter(){
            bucket_sizes[h.hash(&seq) as usize] += 1;
        }

        // Get the starting positions of buckets
        let mut sum = 0 as usize;
        let mut bucket_starts = vec![0];
        for i in 0 .. bucket_sizes.len(){
            sum += bucket_sizes[i];
            bucket_starts.push(sum)
        }
        bucket_starts.shrink_to_fit();

        // Store the locations
        let mut locations: Vec::<(u32, u32)> = vec![(0,0); *bucket_starts.last().unwrap()]; // Will have end sentinel
        for (minmer, seq_id, pos) in L{
            let bucket = h.hash(&minmer) as usize;
            locations[bucket_starts[bucket]] = (seq_id, pos);
            bucket_starts[bucket] += 1;
        }

        // Rewind back the bucket starts
        for (i,s) in bucket_sizes.iter().enumerate(){
            bucket_starts[i] -= s;
        }

        (locations, bucket_starts)
    }

    pub fn new(db: &'a SeqDB, k: usize, m: usize) -> Self{
        if m > k {
            panic!("m > k");
        }
        
        log::info!("Finding minimizers");

        let parts = (0..db.sequence_count()).into_par_iter()
            .map(|i|{
                (i, get_minimizer_positions_with_return(db.get(i).seq, k, m))
            })
            .map(|(i, pos_list)|{
                let mut tuples = Vec::<(Kmer, u32, u32)>::new();
                for p in pos_list.iter(){
                    let minmer = Kmer::from_ascii(&db.get(i).seq[*p..*p+m]).unwrap();
                    tuples.push((minmer, i as u32, *p as u32));
                }
                tuples
            }) // Now we have lists of (minimizer, seq_id, seq_pos) tuples
            .fold(Vec::new, |mut x, y| {
                x.extend(y); x
            })
            .collect::<Vec::<Vec::<(Kmer, u32, u32)>>>();

        let mut position_list = parts.into_iter().fold(Vec::new(), |mut x, y| {
            x.extend(y); x
        });

        log::info!("Sorting tuples (minimizer, seq_id, seq_pos)");
        position_list.par_sort_unstable();

        log::info!("Collecting distinct minimizers");
        let mut minimizer_list: Vec<Kmer> = vec![];
        for (minmer, _, _) in position_list.iter(){
            match minimizer_list.last(){
                Some(last) => {
                    if minmer != last{
                        minimizer_list.push(*minmer);
                    };
                },
                None => { minimizer_list.push(*minmer); }
            }
        }

        log::info!("Removing duplicate minimizers");
        minimizer_list.dedup();
        minimizer_list.shrink_to_fit();

        log::info!("Found {} distinct minimizers", minimizer_list.len());

        log::info!("Building an MPHF for the minimizers");
        let n_mmers = minimizer_list.len();
        let mphf = boomphf::Mphf::<Kmer>::new_parallel(1.7, minimizer_list.as_slice(), None);
        drop(minimizer_list);
        
        log::info!("Compressing position lists");
        let (locations, bucket_starts) = Self::compress_sorted_position_list(position_list, &mphf, n_mmers);

        log::info!("Stored {} location pairs", locations.len());
    
        Self{seq_storage: db, mphf, locations, bucket_starts, k, m, n_mmers}
    }

    // Returns all occurrences of the query k-mer
    pub fn lookup(&self, kmer: &[u8]) -> Vec<(usize, usize)>{
        assert!(kmer.len() == self.k);
        let min_pos = get_minimizer_position(kmer, self.m);
        let minimizer = &kmer[min_pos .. min_pos + self.m];

        let mut ans: Vec<(usize,usize)> = vec![];
        match self.mphf.try_hash(&Kmer::from_ascii(&minimizer).unwrap()){
            Some(bucket) => {
                let bucket_range = self.bucket_starts[bucket as usize]..self.bucket_starts[bucket as usize + 1];
                for (seq_id, seq_pos) in self.locations[bucket_range].iter(){
                    
                    // Start of the k-mer that contains this minimizer occurrence:
                    let start = *seq_pos as i64 - min_pos as i64;

                    // Check if this occurrence is real
                    if start >= 0 && start + self.k as i64 <= self.seq_storage.get(*seq_id as usize).seq.len() as i64 {
                        // k-mer is within bounds of the sequence
                        let start = start as usize;
                        let candidate = &self.seq_storage.get(*seq_id as usize).seq[start .. start + self.k];
                        if candidate == kmer{
                            ans.push((*seq_id as usize, start));
                        }
                    }
                }
            },
            None => (), // No matches
        }

        ans
    }

    /*
    fn print_space_breakdown(&self){
        let seqs_bytes = bincode::serialize(&self.seq_storage).unwrap().len();
        let mphf_bytes = bincode::serialize(&self.mphf).unwrap().len();
        let locations_bytes = bincode::serialize(&self.locations).unwrap().len();

        eprintln!("Seqs storage:\t\t{} bytes\nMPHF:\t\t{} bytes\nLocations:\t{} bytes\n", seqs_bytes, mphf_bytes, locations_bytes);

        let mut total_kmers = 0 as usize;
        let mut total_nucleotides = 0 as usize;
        for i in 0 .. self.seq_storage.number_of_sequences() {
            let seq =  self.seq_storage.get(i);
            total_nucleotides += seq.len();
            if seq.len() >= self.k{
                total_kmers += seq.len() - self.k + 1;
            }
        }

        eprintln!("Number of k-mers including duplicates: {}", total_kmers);
        eprintln!("Total nucleotides: {}", total_nucleotides);
        eprintln!("Seq storage Bits/k-mer: {}", (seqs_bytes * 8) as f64 / total_kmers as f64);
        eprintln!("Seq storage Bits/nucleotide: {}", (seqs_bytes * 8) as f64 / total_nucleotides as f64);
        eprintln!("MPHF Bits/k-mer: {}", (mphf_bytes * 8) as f64 / total_kmers as f64);
        eprintln!("MPHF Bits/m-mer: {}", (mphf_bytes * 8) as f64 / self.n_mmers as f64);
        eprintln!("Locations Bits/k-mer: {}", (locations_bytes * 8) as f64 / total_kmers as f64);  
    }
    */

}

mod tests{

    use std::io::BufReader;

    use super::*;
    use jseqio::reader::*;

    fn number_of_kmers(seq_len: usize, k: usize) -> usize { // TODO: use everywhere
        std::cmp::max(0, (seq_len as i64) - (k as i64) + 1) as usize
    }

    fn to_ascii(S: &[u8]) -> String{
        std::str::from_utf8(S).unwrap().to_owned()
    }

    #[test]
    fn test_get_minimizer_positions(){
        let seq = "ATAGCTAGTCGATGCTGATCGTAGGTTCGTAGCTGTATGCTGACCCTGATGTCTGTAGTCGTGACTGACT";
        let k: usize = 31;
        let m: usize = 10;
        let mut positions: Vec<usize> = vec![];
        get_minimizer_positions(seq.as_bytes(), &mut positions, k, m);

        assert_eq!(positions, vec![2,22,30,42])
    }

    #[test]
    fn test_index_lookup(){

        let input = "\
>seq1
ATAGCTAGTCGATGCTGATCGTAGGTTCGTAGCTGTATGCTGACCCTGATGTCTGTAGTCGTGACTGACT
>seq2 (Substring of seq1)
GTCGATGCTGATCGTAGGTTCGTAGCTGTATGCTGACCCTGATGTCTTGACT
>seq3 (seq2 with a single change in the middle)
GTCGATGCTGATCGTAGGTTCGAAGCTGTATGCTGACCCTGATGTCTTGACT
>seq4 (contains N's)
GNCGATGCTGATCGTAGGTTCGAAGCTATTCGATGCGTATGCTGACNCCTGATGTCTTGACTATATGTCGTAGTTTCGATCGAGAGAGTATAGAANGNA"
.as_bytes();

        let k: usize = 31;
        let m: usize = 10;
        
        // Build index
        let mut reference = DynamicFastXReader::new(BufReader::new(input)).unwrap();
        let db = reference.into_db().unwrap();
        let index = MinimizerIndex::new(&db, k, m);
        
        // Read sequences
        let mut reader = DynamicFastXReader::new(BufReader::new(input)).unwrap();
        let mut seqs : Vec<Vec<u8>> = vec![];
        while let Some(rec) = reader.read_next().unwrap(){
            seqs.push(rec.seq.to_owned());
        }

        // Build true k-mer occurrences map
        let mut true_kmer_occurrences = std::collections::HashMap::<Vec<u8>, Vec<(usize, usize)>>::new();
        for (seq_id, seq) in seqs.iter().enumerate(){
            for i in 0 .. number_of_kmers(seq.len(), k){
                let key = &seq[i..i+k];
                let new_entry = (seq_id, i);
                if let Some(vec) = true_kmer_occurrences.get_mut(key){ // Existing k-mer
                    vec.push(new_entry);
                } else{ // New k-mer
                    true_kmer_occurrences.insert(key.to_owned(), vec![new_entry]);
                }
            }
        }

        // Look up all k-mers in input sequences
        for seq in seqs.iter(){
            for i in 0 .. number_of_kmers(seq.len(), k){
                let kmer = &seq[i..i+k];
                let occs = index.lookup(kmer);
                eprintln!("{} {:?} {:?}", to_ascii(&kmer), &occs, &true_kmer_occurrences[kmer]);
                assert_eq!(occs, true_kmer_occurrences[kmer]);
            }
        }

        // Look up a random k-mer (should not be found)
        let random_kmer = "ATCTTATCTGGGGCTATTGCTAGGGCTTACA".as_bytes();
        assert_eq!(index.lookup(random_kmer).len(), 0);
    }

}