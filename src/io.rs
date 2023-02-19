use crate::fasta;
use eyre::{Result, WrapErr};
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn open_fasta(filename: &str) -> Result<fasta::FastaReader<File>> {
    Ok(fasta::FastaReader::new(
        File::open(filename).wrap_err(format!("while opening {}", filename))?,
        true,
    ))
}

pub fn write_dist_matrix<T: std::fmt::Display, L: std::fmt::Display>(
    m: &[T],
    ids: &[L],
    mut out: impl std::io::Write,
) -> Result<()> {
    let n = ids.len();
    assert!(m.len() == n * n);

    // 1. Write the number of elements
    writeln!(out, "{}", ids.len())?;

    // 2. Write the matrix itself
    for (i, id) in ids.iter().enumerate() {
        write!(out, "{}", id)?;
        for j in 0..ids.len() {
            write!(out, "\t{:.5}", m[n * i + j])?;
        }
        writeln!(out)?;
    }

    // 3. Flush the output
    Ok(out.flush()?)
}
