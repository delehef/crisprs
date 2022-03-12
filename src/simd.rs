use std::arch::x86_64::*;

pub unsafe fn print_mm256(a: __m256i) {
    let mut o = vec![b'_'; 32];
    _mm256_store_si256(o.as_mut_ptr() as *mut __m256i, a);
    println!("{}", std::str::from_utf8(&o).unwrap());
}

pub unsafe fn print_bool256(a: __m256i) {
    let mut o = vec![b'_'; 32];
    _mm256_store_si256(o.as_mut_ptr() as *mut __m256i, a);
    o.iter_mut().for_each(|b| {
        if *b == 0xff {
            *b = b'.'
        }
        if *b == 0x0 {
            *b = b' '
        }
    });
    println!("{}", std::str::from_utf8(&o).unwrap());
}
