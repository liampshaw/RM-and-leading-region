use clap::{Parser, Subcommand};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Single {
        #[arg(short, long)]
        input: String,
        #[arg(short, long)]
        output: String,
        #[arg(short, long)]
        k: usize,
        #[arg(short = 'w', long, default_value_t = 5000)]
        window_size: usize,
        #[arg(short = 's', long, default_value_t = 1000)]
        step_size: usize,
    },
    Batch {
        #[arg(short, long)]
        files: String,
        #[arg(short, long)]
        output: String,
        #[arg(short, long)]
        k: usize,
        #[arg(short = 'w', long, default_value_t = 5000)]
        window_size: usize,
        #[arg(short = 's', long, default_value_t = 1000)]
        step_size: usize,
    },
}

#[inline(always)]
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => b'N',
    }
}

fn read_fasta(path: &str) -> io::Result<Vec<u8>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            seq.extend_from_slice(line.trim().to_uppercase().as_bytes());
        }
    }
    Ok(seq)
}

/// Check palindrome in circular sequence
fn is_palindrome_circular(seq: &[u8], start: usize, k: usize) -> bool {
    let n = seq.len();
    for i in 0..k / 2 {
        let a = seq[(start + i) % n];
        let b = seq[(start + k - 1 - i) % n];
        if complement(a) != b {
            return false;
        }
    }
    true
}

/// Compute circular palindrome density
fn compute_density_circular(
    seq: &[u8],
    k: usize,
    window_size: usize,
    step_size: usize,
) -> Vec<f64> {
    let n = seq.len();

    // Compute palindrome flags at all n positions
    let mut pal_flags = vec![0u32; n];
    for i in 0..n {
        if is_palindrome_circular(seq, i, k) {
            pal_flags[i] = 1;
        }
    }

    let num_windows = n / step_size;
    let mut densities = Vec::with_capacity(num_windows);

    for w in 0..num_windows {
        let start = w * step_size;
        let mut window_sum = 0u32;

        for i in 0..window_size {
            window_sum += pal_flags[(start + i) % n];
        }

        densities.push(window_sum as f64 / window_size as f64);
    }

    densities
}

fn write_output(
    path: &str,
    densities: &[f64],
    step_size: usize,
) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);

    for (i, val) in densities.iter().enumerate() {
        let position = i * step_size;
        writeln!(writer, "{}\t{}", position, val)?;
    }

    Ok(())
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Single {
            input,
            output,
            k,
            window_size,
            step_size,
        } => {
            let seq = read_fasta(&input)?;
            let densities =
                compute_density_circular(&seq, k, window_size, step_size);
            write_output(&output, &densities, step_size)?;
        }

        Commands::Batch {
            files,
            output,
            k,
            window_size,
            step_size,
        } => {
            let file = File::open(files)?;
            let reader = BufReader::new(file);
            let filenames: Vec<String> =
                reader.lines().filter_map(Result::ok).collect();

            // Compute in parallel
            let all_results: Vec<(String, Vec<f64>)> = filenames
                .par_iter()
                .map(|filename| {
                    let seq =
                        read_fasta(filename).expect("Failed to read FASTA");
                    let densities =
                        compute_density_circular(&seq, k, window_size, step_size);

                    (filename.clone(), densities)
                })
                .collect();

            // Write output sequentially
            let mut writer = BufWriter::new(File::create(output)?);

            for (plasmid, densities) in all_results {
                for (i, val) in densities.iter().enumerate() {
                    let position = i * step_size;
                    writeln!(writer, "{}\t{}\t{}", plasmid, position, val)?;
            }
        }
    }

    }

    Ok(())
}

