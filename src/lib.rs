#[allow(non_snake_case)]
pub fn hamming_distance_not_matching_N(a: &[u8], b: &[u8]) -> usize{
    assert_eq!(a.len(), b.len());
    let mut dist = 0;
    for (x,y) in a.iter().zip(b.iter()){
        if x != y || (x == &b'N' && y == &b'N'){
            dist += 1;
        }
    }
    dist
}
