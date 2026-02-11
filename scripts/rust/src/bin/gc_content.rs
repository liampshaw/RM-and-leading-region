use clap::{Parser, Subcommand};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};

#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute GC content for a single plasmid
    Single {
        #[arg(short, long)]
        input: String,
        #[arg(short, long)]
        output: String,
        #[arg(short = 'w', long, default_value_t = 5000)]
        window_size: usize,
        #[arg(short = 's', long, default_value_t = 1000)]
        step_size: usize,
    },
    /// Compute GC content for multiple plasmids listed in a file
    Batch {
        #[arg(short, long)]
        files: String, // file containing FASTA paths
        #[arg(short, long)]
        output: String,
        #[arg(short = 'w', long, default_value_t = 5000)]
        window_size: usize,
        #[arg(short = 's', long, default_value_t = 1000)]
        step_size: usize,
    },
}

/// Read a FASTA file into a sequence (uppercase)
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

/// Compute GC content in circular sliding windows
fn gc_content_circular(seq: &[u8], window_size: usize, step_size: usize) -> Vec<f64> {
    let n = seq.len();
    let num_windows = n / step_size;
    let mut densities = Vec::with_capacity(num_windows);

    for w in 0..num_windows {
        let start = w * step_size;
        let mut gc = 0u32;
        for i in 0..window_size {
            let base = seq[(start + i) % n];
            if base == b'G' || base == b'C' {
                gc += 1;
            }
        }
        densities.push(gc as f64 / window_size as f64);
    }

    densities
}

/// Write GC output for a plasmid
fn write_output(path: &str, plasmid: &str, densities: &[f64], step_size: usize) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (i, val) in densities.iter().enumerate() {
        let position = i * step_size;
        writeln!(writer, "{}\t{}\t{}", plasmid, position, val)?;
    }
    Ok(())
}

/// Write GC output to an existing writer (for batch mode)
fn write_batch_output(
    writer: &mut BufWriter<File>,
    plasmid: &str,
    densities: &[f64],
    step_size: usize,
) -> io::Result<()> {
    for (i, val) in densities.iter().enumerate() {
        let position = i * step_size;
        writeln!(writer, "{}\t{}\t{}", plasmid, position, val)?;
    }
    Ok(())
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Single {
            input,
            output,
            window_size,
            step_size,
        } => {
            let seq = read_fasta(&input)?;
            let densities = gc_content_circular(&seq, window_size, step_size);
            let plasmid = input.split('/').last().unwrap_or(&input);
            write_output(&output, plasmid, &densities, step_size)?;
        }

        Commands::Batch {
            files,
            output,
            window_size,
            step_size,
        } => {
            let file = File::open(files)?;
            let reader = BufReader::new(file);
            let filenames: Vec<String> = reader.lines().filter_map(Result::ok).collect();

            // Shared writer wrapped in Arc<Mutex<>> for thread-safe access
            let writer = Arc::new(Mutex::new(BufWriter::new(File::create(&output)?)));

            filenames.par_iter().for_each(|filename| {
                let seq = read_fasta(filename).expect("Failed to read FASTA");
                let densities = gc_content_circular(&seq, window_size, step_size);
                let plasmid = filename.split('/').last().unwrap_or(filename);

                // Lock the writer to write safely
                let mut w = writer.lock().unwrap();
                write_batch_output(&mut *w, plasmid, &densities, step_size)
                    .expect("Failed to write output");
            });
        }
    }

    Ok(())
}

