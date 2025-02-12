use jseqio::reverse_complement;
use jseqio::seq_db::SeqDB;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{io::{BufWriter, Write}, sync::atomic::{AtomicI16, AtomicU32, Ordering::Relaxed}};
use crate::minimizer_index::MinimizerIndex;

fn update_coverage(coverages: &Vec<Vec<AtomicU32>>, mismatches: &Vec<Vec<AtomicU32>>, bait: &[u8], index: &MinimizerIndex, targets_db: &SeqDB, hamming_distance: usize){
    let candidates = index.get_exact_alignment_candidates(bait);

    for (target_id, target_start) in candidates{
        let target = &targets_db.get(target_id).seq[target_start .. target_start + bait.len()];
        let distance = syotti2::hamming_distance_not_matching_N(bait, target);
        if distance <= hamming_distance{
            for i in 0..bait.len(){
                // Relaxed atomic ordering is okay for addition because addition is commutative
                coverages[target_id][target_start + i].fetch_add(1, Relaxed); // TODO: this may overflow: is there a saturating add?

                // Relaxed atomic ordering is okay for min because min is commutative
                mismatches[target_id][target_start + i].fetch_min(distance as u32, Relaxed);
            }
        }
    }
}

pub fn write_as_csv<T: std::fmt::Display, F: Fn(&T) -> String>(lines: &Vec<Vec<T>>, out: &mut impl Write, formatter: F){
    for v in lines.iter() {
        // Write v as a line of comma-separated values
        for (i, x) in v.iter().enumerate(){
            if i > 0 {
                write!(out, ",").unwrap();
                
            }
            write!(out, "{}", formatter(x)).unwrap();
        }
        write!(out, "\n").unwrap();
    }
}

// T is Generic over any type that is convertible to f64
pub fn write_as_png<T: Into<f64> + Clone>(coverages: Vec<Vec<T>>, color_scale: Option<usize>, out: &mut impl Write){
    let max_len = coverages.iter().map(|x| x.len()).max().unwrap();
    let max_coverage = match color_scale {
        Some(s) => s as f64,
        None => {
            let mut max_cov = 0.0;
            for v in coverages.iter() {
                for x in v.iter(){
                    let val: f64 = x.clone().into();
                    max_cov = f64::max(max_cov, val);
                }
            }
            max_cov
        }
    };

    let width = max_len;
    let height = coverages.len();

    let mut encoder = png::Encoder::new(BufWriter::new(out), width as u32, height as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);

    let mut pixels = vec![0_u8; width*height*3]; // Three color components per pixel
    let red_threshold = 0.0; // Coverage at or below this are drawn in red.
    //let log_max = f64::ln(max_coverage as f64); // This will be fully white
    let brigtness_scaling_factor = 255.0_f64 / max_coverage; // max_coverage * c = 255

    for (row, v) in coverages.iter().enumerate() {
        for (col, x) in v.iter().enumerate(){
            let x: f64 = x.clone().into();
            if x < 0.0 { // Padding
                pixels[row*width*3 + col*3] = 251; // Red
                pixels[row*width*3 + col*3 + 1] = 255; // Green
                pixels[row*width*3 + col*3 + 2] = 43; // Blue
            } else if x <= red_threshold {
                pixels[row*width*3 + col*3] = 255; // Red
                pixels[row*width*3 + col*3 + 1] = 0; // Green
                pixels[row*width*3 + col*3 + 2] = 0; // Blue
            } else {
                let brightness = f64::min(x, max_coverage) * brigtness_scaling_factor;
                assert!(brightness.round() <= 255.0);
                let brightness: u8 = brightness.round() as u8;
                //let mut val = f64::ln(x / coverage_floor); // Now this should be >= 0
                //val = f64::max(val, 0.0); // Clamp to zero in case of floating point inaccuracies which could make val slightly negative
                //val *= brigtness_scaling_factor; // Now the brightest pixel should have val = 255
                //let brightness = f64::min(val, 255.0) as u8; // Again clamp in case of floating point weirdness

                pixels[row*width*3 + col*3] = brightness; // Red
                pixels[row*width*3 + col*3 + 1] = brightness; // Green
                pixels[row*width*3 + col*3 + 2] = brightness; // Blue
            }
        }
    }

    let mut writer = encoder.write_header().unwrap();
    writer.write_image_data(&pixels).unwrap();


}

// Note: searches both forward and reverse complement. This means that if a bait overlaps with its own
// reverse complement, it could contribute 2 to the coverage depth at the overlapping positions.
pub fn compute_coverage(targets_db: &SeqDB, bait_db: &SeqDB, d: usize, g: usize, m: usize, n_threads: usize) -> (Vec<Vec<u32>>, Vec<Vec<u32>>){

    let index = MinimizerIndex::new(&targets_db, g, m);

    // Initialize coverage depth vectors
    // coverages[i][j] = depth in target i at position j
    let mut coverages = Vec::<Vec::<AtomicU32>>::new();  // Let's hope 32 bits are enough
    let mut mismatches = Vec::<Vec::<AtomicU32>>::new(); // Let's hope 32 bits are enough
    for i in 0..targets_db.sequence_count(){
        let len = targets_db.get(i).seq.len();
        let mut cov = Vec::<AtomicU32>::with_capacity(len);
        let mut mis = Vec::<AtomicU32>::with_capacity(len);
        for _ in 0..len {
            cov.push(AtomicU32::new(0));
            mis.push(AtomicU32::new(u32::MAX));
        }
        coverages.push(cov);
        mismatches.push(mis);
    }

    let progress = indicatif::ProgressBar::new(bait_db.sequence_count() as u64);

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(|| {
        (0..bait_db.sequence_count()).into_par_iter().for_each(|i| {
            let bait = bait_db.get(i);
            update_coverage(&coverages, &mismatches, bait.seq, &index, targets_db, d);
            update_coverage(&coverages, &mismatches, reverse_complement(bait.seq).as_slice(), &index, targets_db, d);
            progress.inc(1);
        });
    });

    // This seems to compile to basically a no-op! At least this takes no time at all.
    let non_atomic_coverages = coverages.into_iter().map(|v| v.into_iter().map(|x| x.into_inner()).collect()).collect();
    let non_atomic_mismatches = mismatches.into_iter().map(|v| v.into_iter().map(|x| x.into_inner()).collect()).collect();

    log::info!("Coverage computation finished");
    (non_atomic_coverages, non_atomic_mismatches)

}

pub fn into_moving_average(v: Vec<u32>, window_size: usize) -> Vec<f32>{
    assert!(window_size > 0);
    
    if v.len() < window_size{
        return vec![];
    }

    // Compute the average of the first window
    let mut sum = v[0..window_size].iter().fold(0_usize, |s,x| s + *x as usize);
    let mut avgs: Vec<f32> = vec![(sum as f32) / (window_size as f32)];

    // Compute averages of the rest
    for i in window_size..v.len(){
        sum = sum - v[i-window_size] as usize + v[i] as usize;
        avgs.push((sum as f32) / (window_size as f32));
    }

    avgs
}

// Samples the moving average of the coverage vectors to get exactly 'points' points for each. If variable_resolution is
// given, then every coverage vector gets a number of points proportional to its length, such that the longest
// one gets 'points' points.
pub fn into_resolution(coverages: Vec<Vec<u32>>, points: usize, variable_resolution: bool) -> Vec<Vec<f32>>{
    let mut new_covs = Vec::<Vec::<f32>>::new();
    let max_len = coverages.iter().fold(0_usize, |acc, v| acc.max(v.len()));
    for cov in coverages.into_iter(){
        let r = if variable_resolution { // Number of points in result
            let r = ((cov.len() as f64 / max_len as f64) * points as f64).round() as usize;
            r.max(1) // No zero points because that's weird
        } else {
            points
        };

        let window_len = std::cmp::max(cov.len() / r, 1);
        let sampling_step = cov.len() as f64 / r as f64;
        let avgs = into_moving_average(cov, window_len);

        let mut sampled_points: Vec<f32> = vec![0.0; r];
        
        for i in 0..r{
            sampled_points[i] = avgs[i*sampling_step as usize];
        }
        new_covs.push(sampled_points);
    }
    new_covs
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

        let mut t0_answer: Vec<u32> = vec![0; targets[0].len()];
        for i in 0..40 {t0_answer[i] += 1;}
        for i in 30..40 {t0_answer[i] += 1;}

        let mut t1_answer: Vec<u32> = vec![0; targets[1].len()];
        for i in 1..21 {t1_answer[i] += 1;}

        let mut t2_answer: Vec<u32> = vec![0; targets[2].len()];
        for i in 9..19 {t2_answer[i] += 1;}

        let answer = vec![t0_answer, t1_answer, t2_answer];

        let (coverages, _) = compute_coverage(&target_db, &bait_db,  1, 5, 3, 10);

        assert_eq!(answer, coverages);
    }

    #[test]
    fn test_moving_avg(){
        let data: Vec<Vec<u32>> = vec![vec![1,2,3,4,5], vec![1,2], vec![]];

        let window_len = 3;
        let answer: Vec<Vec<f32>> = vec![vec![2.0, 3.0, 4.0], vec![], vec![]];

        assert_eq!(into_moving_average(data[0].clone(), window_len), answer[0]);
        assert_eq!(into_moving_average(data[1].clone(), window_len), answer[1]);
        assert_eq!(into_moving_average(data[2].clone(), window_len), answer[2]);
    }
}

