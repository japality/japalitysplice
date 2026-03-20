# JapalitySplice

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19104744.svg)](https://doi.org/10.5281/zenodo.19104744)

**A sparse k-mer, splice-aware RNA-Seq aligner designed for resource-constrained environments.**

JapalitySplice is a novel RNA-Seq aligner written in Rust that combines sparse k-mer indexing
with empirical splice-junction detection. It is specifically designed to operate within
limited memory and compute budgets — making it suitable for edge devices, laptops, and
teaching environments — while still providing competitive mapping quality on small-to-medium
genomes.

---

## License & Commercial Use

This project is licensed under the
[GNU General Public License v3.0 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.html).

**For academic and non-commercial use**, you are free to use, modify, and distribute this
software under the terms of the GPLv3.

**For commercial use** (including but not limited to integrating this algorithm into
proprietary software, closed-source commercial apps, or internal business tools), you must
obtain a separate commercial license.
Please contact **info@japality.com** for commercial licensing inquiries.

Copyright © 2026 Japality Limited. All rights reserved.

---

## Features

| Feature | Description |
|---|---|
| **Sparse k-mer indexing** | Builds a compact in-memory index from reference genomes using configurable k-mer sizes (9 / 11 / 13) |
| **Splice-aware alignment** | Detects GT-AG / GC-AG splice junctions via empirical motif scanning |
| **Annotation-guided mode** | Optionally loads known junctions from GTF / GFF3 annotations to boost junction precision |
| **Paired-end support** | Aligns both single-end and paired-end FASTQ reads and emits paired SAM fields |
| **Low memory footprint** | Yeast index: ~0.16 GB; Arabidopsis: ~0.64 GB; Mouse: ~10.3 GB |
| **Parallel alignment** | Scales across all available cores via Rayon work-stealing |
| **Compressed I/O** | Natively reads gzipped FASTQ input; plain-text FASTA is used for index building |

## Quick Start

### Build

```bash
# Requires Rust ≥ 1.80 (uses std::sync::LazyLock)
cargo build --release
```

The binary is produced at `target/release/japalitysplice`.

### Build an Index

```bash
japalitysplice build-index \
    --fasta genome.fa \
    --out index_dir/ \
    --kmer-size 11 \
    --threads 8
```

### Align Reads

**Single-end, no annotation:**
```bash
japalitysplice align \
    --index index_dir/genome.nsix \
    --r1 reads.fastq.gz \
    --out results/ \
    --threads 8
```

**Paired-end with annotation-guided junction detection:**
```bash
japalitysplice align \
    --index index_dir/genome.nsix \
    --r1 reads_R1.fastq.gz \
    --r2 reads_R2.fastq.gz \
    --out results/ \
    --known-junctions annotation.gtf \
    --threads 8
```

Output is a standard SAM file in the specified output directory.

## Benchmark Summary

Benchmarks were conducted on an AMD Ryzen 7 9700X, 92 GB DDR5 RAM, Linux,
using 16 threads and comparing JapalitySplice (JP) against STAR 2.7.11b and HISAT2 2.2.1.

| Species | Tool | Annotation | Map% | Precision% | Recall% | Index RAM (GB) |
|---|---|---|---|---|---|---|
| Yeast (12 Mb) | STAR | off | 89.5 | 100.0 | 100.0 | 1.61 |
| Yeast | STAR | on | 89.5 | 99.7 | 99.1 | 1.61 |
| Yeast | HISAT2 | off | 85.5 | 13.1 | 14.9 | 0.11 |
| Yeast | HISAT2 | on | 85.5 | 13.1 | 15.2 | 0.74 |
| Yeast | JP | off | 91.0 | 8.2 | 1.9 | 0.16 |
| Yeast | JP | on | 91.2 | 12.9 | 3.1 | 0.16 |
| Arabidopsis (120 Mb) | STAR | off | 96.7 | 100.0 | 100.0 | 2.61 |
| Arabidopsis | STAR | on | 96.7 | 98.6 | 93.2 | 4.19 |
| Arabidopsis | HISAT2 | off | 95.5 | 94.3 | 87.4 | 0.26 |
| Arabidopsis | HISAT2 | on | 95.6 | 93.4 | 86.6 | 5.95 |
| Arabidopsis | JP | off | 88.1 | 23.0 | 11.0 | 0.64 |
| Arabidopsis | JP | on | 90.1 | 65.4 | 44.1 | 0.64 |
| Mouse (2.7 Gb) | STAR | off | 81.0 | 100.0 | 100.0 | 26.96 |
| Mouse | STAR | on | 81.3 | 97.0 | 98.3 | 31.77 |
| Mouse | HISAT2 | off | 79.2 | 62.4 | 86.4 | 4.72 |
| Mouse | JP | off | 66.0 | 55.6 | 0.2 | 10.30 |
| Mouse | JP | on | 66.0 | 90.9 | 0.8 | 10.30 |

**Strengths:** Fast index construction across all tested genomes; low RAM on yeast and
Arabidopsis; strong annotation-guided precision gains; single portable binary.

**Limitations:** Junction recall remains weak on large mammalian genomes even with annotation,
and alignment throughput trails STAR/HISAT2 on complex genomes.

Mouse HISAT2 annotation-on did not yield a complete SAM output in the benchmark workflow and is
therefore excluded from the summary table above, matching the manuscript.

For full methodology and discussion, see the accompanying paper (below).

## Citation

If you use JapalitySplice in your research, please cite:

> Chen, Y.H. & Kojima, K. (2026). JapalitySplice: A Sparse K-mer RNA-Seq Splice-Aware
> Aligner for Edge-Oriented Computing. Japality Limited, Hong Kong.
> DOI: [10.5281/zenodo.19104744](https://doi.org/10.5281/zenodo.19104744)

## Project Structure

```
japalitysplice-public/
├── Cargo.toml                 # Rust package manifest
├── LICENSE                    # GPLv3
├── README.md                  # This file
├── .gitignore
└── src/
    ├── lib.rs                 # Core library: indexing, alignment, k-mer engine
    ├── splice.rs              # Splice junction detection (motif + annotation)
    ├── dp.rs                  # Dynamic programming local alignment
    └── bin/
        └── japalitysplice.rs  # CLI entry point
```

## Mobile App

JapalitySplice is also available as a mobile application with a graphical interface,
allowing on-device RNA-Seq alignment without any server infrastructure:

- **Android**: [Google Play](https://play.google.com/store/apps/details?id=com.japality.splice)
- **iOS**: [App Store](https://apps.apple.com/app/id6760461634)

## Contributing

Bug reports and pull requests are welcome on [GitHub](https://github.com/japality/japalitysplice).
Please note that contributions are accepted under the GPLv3 license.

## Contact

- **Japality Limited** — Hong Kong
- Email: info@japality.com
