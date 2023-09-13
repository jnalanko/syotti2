use std::path::Path;
use jseqio::reader::*;
use crate::minimizer_index::MinimizerIndex;

fn compute_coverage(targetfile: &Path, baitfile: &Path, outfile: &Path, d: usize, g: usize, m: usize) {
    
    let mut bait_reader = DynamicFastXReader::from_file(&baitfile).unwrap();
    let targets_db = DynamicFastXReader::from_file(&targetfile).unwrap().into_db().unwrap();
    let index = MinimizerIndex::new(&targets_db, g, m);

    while let Some(bait) = bait_reader.read_next().unwrap(){
        
    }

}