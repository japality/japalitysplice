// JapalitySplice CLI — sparse k-mer RNA-Seq splice-aware aligner
// Copyright (C) 2026 Japality Limited
// SPDX-License-Identifier: GPL-3.0-or-later

use rayon::ThreadPoolBuilder;
use std::env;
use std::process;

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn usage() -> ! {
    eprintln!(
        "JapalitySplice v{} — sparse k-mer RNA-Seq splice-aware aligner

Usage:
  japalitysplice build-index --fasta <path> --out <dir> [--kmer-size 11] [--threads 16]
  japalitysplice align --index <path> --r1 <path> [--r2 <path>] --out <dir> \
      [--max-mismatches 3] [--known-junctions <gtf/gff3>] [--threads 16]

Commands:
  build-index   Build a sparse k-mer index (.nsix) from a FASTA reference genome
  align         Align single-end or paired-end FASTQ reads against the index

Options:
  --fasta             Path to reference genome FASTA (plain or .gz)
  --out               Output directory
  --kmer-size         K-mer size for indexing (9, 11, or 13; default: 11)
  --index             Path to .nsix index file
  --r1                Path to R1 FASTQ (plain or .gz)
  --r2                Path to R2 FASTQ for paired-end alignment
  --max-mismatches    Maximum mismatches allowed (default: 3)
  --known-junctions   GTF/GFF3 annotation for annotation-guided splice detection
  --threads           Number of parallel threads (default: all available)
  --version           Print version and exit",
        VERSION
    );
    process::exit(1);
}

fn take_required(args: &[String], flag: &str) -> String {
    args.iter()
        .position(|arg| arg == flag)
        .and_then(|idx| args.get(idx + 1))
        .cloned()
        .unwrap_or_else(|| {
            eprintln!("Missing required flag: {}", flag);
            usage();
        })
}

fn take_optional(args: &[String], flag: &str) -> Option<String> {
    args.iter()
        .position(|arg| arg == flag)
        .and_then(|idx| args.get(idx + 1))
        .cloned()
}

fn parse_usize(args: &[String], flag: &str, default: usize) -> usize {
    take_optional(args, flag)
        .map(|value| {
            value.parse::<usize>().unwrap_or_else(|_| {
                eprintln!("Invalid numeric value for {}: {}", flag, value);
                process::exit(1);
            })
        })
        .unwrap_or(default)
}

fn configure_threads(args: &[String]) {
    let threads = parse_usize(args, "--threads", 0);
    if threads == 0 {
        return;
    }

    if let Err(error) = ThreadPoolBuilder::new().num_threads(threads).build_global() {
        eprintln!(
            "Failed to configure Rayon thread pool ({} threads): {}",
            threads, error
        );
        process::exit(1);
    }
}

fn load_known_junctions_if_present(annotation_path: Option<String>) {
    japalitysplice_core::clear_known_junctions();
    if let Some(path) = annotation_path {
        let count = japalitysplice_core::init_known_junctions(&path).unwrap_or_else(|error| {
            eprintln!("Failed to load known junctions {}: {}", path, error);
            process::exit(1);
        });
        eprintln!("Loaded {} known junctions from {}", count, path);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        usage();
    }

    configure_threads(&args);

    match args[1].as_str() {
        "build-index" => {
            let fasta = take_required(&args, "--fasta");
            let out_dir = take_required(&args, "--out");
            let kmer_size = parse_usize(&args, "--kmer-size", 11);

            let index_path =
                japalitysplice_core::build_index(&fasta, Some(&out_dir), kmer_size)
                    .unwrap_or_else(|error| {
                        eprintln!("Index build failed: {}", error);
                        process::exit(1);
                    });

            println!("{}", index_path);
        }
        "align" => {
            let index = take_required(&args, "--index");
            let r1 = take_required(&args, "--r1");
            let r2 = take_optional(&args, "--r2");
            let out_dir = take_required(&args, "--out");
            let max_mismatches = parse_usize(&args, "--max-mismatches", 3);
            let known_junctions = take_optional(&args, "--known-junctions");

            load_known_junctions_if_present(known_junctions);

            let sam_path = japalitysplice_core::run_alignment(
                &r1,
                r2.as_deref(),
                &index,
                &out_dir,
                max_mismatches,
            )
            .unwrap_or_else(|error| {
                eprintln!("Alignment failed: {}", error);
                process::exit(1);
            });

            println!("{}", sam_path);
        }
        "--help" | "-h" | "help" => usage(),
        "--version" | "-V" => {
            println!("japalitysplice {}", VERSION);
            process::exit(0);
        }
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            usage();
        }
    }
}
