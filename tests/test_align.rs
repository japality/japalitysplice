// Integration tests for JapalitySplice
// Copyright (C) 2026 Japality Limited
// SPDX-License-Identifier: GPL-3.0-or-later

use japalitysplice_core::{build_index, run_alignment};
use serde_json::Value;
use std::fs;
use std::time::{SystemTime, UNIX_EPOCH};

fn temp_dir(prefix: &str) -> std::path::PathBuf {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("{}_{}", prefix, stamp));
    fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn test_single_end_alignment_produces_valid_sam() {
    let tmp = temp_dir("japalitysplice_se_test");
    let fasta = tmp.join("reference.fa");
    let fastq = tmp.join("reads.fastq");

    fs::write(
        &fasta,
        ">chr1\nACGTACGTTTAACCGGATCCTTAGGCCAATTCGATGGCCTAAGT\n",
    )
    .unwrap();
    fs::write(&fastq, "@read1\nACGTACGTTTAA\n+\nFFFFFFFFFFFF\n").unwrap();

    let index_path = build_index(fasta.to_str().unwrap(), Some(tmp.to_str().unwrap()), 11).unwrap();
    let sam_path =
        run_alignment(fastq.to_str().unwrap(), None, &index_path, tmp.to_str().unwrap(), 0)
            .unwrap();

    let sam = fs::read_to_string(&sam_path).unwrap();
    let headers: Vec<&str> = sam.lines().filter(|l| l.starts_with('@')).collect();
    let records: Vec<&str> = sam.lines().filter(|l| !l.starts_with('@')).collect();

    assert!(headers.iter().any(|h| h.starts_with("@HD")));
    assert!(headers.iter().any(|h| h.starts_with("@SQ")));
    assert!(headers.iter().any(|h| h.starts_with("@PG")));
    assert_eq!(records.len(), 1);

    let fields: Vec<&str> = records[0].split('\t').collect();
    assert_eq!(fields[2], "chr1");
    assert_eq!(fields[3], "1");

    fs::remove_dir_all(&tmp).unwrap();
}

#[test]
fn test_paired_end_sam_fields_are_populated() {
    let tmp = temp_dir("japalitysplice_pair_test");
    let fasta = tmp.join("reference.fa");
    let fastq_r1 = tmp.join("reads_r1.fastq");
    let fastq_r2 = tmp.join("reads_r2.fastq");

    fs::write(
        &fasta,
        ">chr1\nACGTACGTTTAACCGGATCCTTAGGCCAATTCGATGGCCTAAGT\n",
    )
    .unwrap();
    fs::write(&fastq_r1, "@pair1/1\nACGTACGTTTAA\n+\nFFFFFFFFFFFF\n").unwrap();
    fs::write(&fastq_r2, "@pair1/2\nCATCGAATTGGC\n+\nFFFFFFFFFFFF\n").unwrap();

    let index_path = build_index(fasta.to_str().unwrap(), Some(tmp.to_str().unwrap()), 11).unwrap();
    let sam_path = run_alignment(
        fastq_r1.to_str().unwrap(),
        Some(fastq_r2.to_str().unwrap()),
        &index_path,
        tmp.to_str().unwrap(),
        0,
    )
    .unwrap();

    let sam = fs::read_to_string(&sam_path).unwrap();
    let records: Vec<&str> = sam.lines().filter(|l| !l.starts_with('@')).collect();
    assert_eq!(records.len(), 2);

    let r1_fields: Vec<&str> = records[0].split('\t').collect();
    let r2_fields: Vec<&str> = records[1].split('\t').collect();

    // Both reads share the same QNAME (stripped of /1 /2)
    assert_eq!(r1_fields[0], "pair1");
    assert_eq!(r2_fields[0], "pair1");

    // Paired-end flags: 99 = paired+proper+mate_rev+read1, 147 = paired+proper+rev+read2
    assert_eq!(r1_fields[1], "99");
    assert_eq!(r2_fields[1], "147");

    // Mate reference name "=" means same contig
    assert_eq!(r1_fields[6], "=");
    assert_eq!(r2_fields[6], "=");

    // Mate positions are populated (non-zero)
    assert_eq!(r1_fields[7], "25");
    assert_eq!(r2_fields[7], "1");

    // TLEN is non-zero and opposite sign between mates
    assert_eq!(r1_fields[8], "36");
    assert_eq!(r2_fields[8], "-36");

    fs::remove_dir_all(&tmp).unwrap();
}

#[test]
fn test_built_index_uses_compact_seed_format() {
    let tmp = temp_dir("japalitysplice_index_test");
    let fasta = tmp.join("reference.fa");

    fs::write(
        &fasta,
        ">chr1\nACGTACGTTTAACCGGATCCTTAGGCCAATTCGATGGCCTAAGT\n",
    )
    .unwrap();

    let index_path = build_index(fasta.to_str().unwrap(), Some(tmp.to_str().unwrap()), 11).unwrap();

    let index_raw = fs::read_to_string(&index_path).unwrap();
    let index: Value = serde_json::from_str(&index_raw).unwrap();

    assert_eq!(index["format"], "japality_splice.nsix.v2");
    assert!(index["kmer_size"].as_u64().unwrap() == 11);
    assert!(index["seed_stride"].as_u64().unwrap() == 4);
    assert!(index["contig_count"].as_u64().unwrap() >= 1);

    let seed_index = &index["seed_index"];
    assert!(seed_index["kmer_size"].as_u64().unwrap() == 11);

    fs::remove_dir_all(&tmp).unwrap();
}

#[test]
fn test_pg_header_contains_version() {
    let tmp = temp_dir("japalitysplice_pg_test");
    let fasta = tmp.join("reference.fa");
    let fastq = tmp.join("reads.fastq");

    fs::write(
        &fasta,
        ">chr1\nACGTACGTTTAACCGGATCCTTAGGCCAATTCGATGGCCTAAGT\n",
    )
    .unwrap();
    fs::write(&fastq, "@read1\nACGTACGTTTAA\n+\nFFFFFFFFFFFF\n").unwrap();

    let index_path = build_index(fasta.to_str().unwrap(), Some(tmp.to_str().unwrap()), 11).unwrap();
    let sam_path =
        run_alignment(fastq.to_str().unwrap(), None, &index_path, tmp.to_str().unwrap(), 3)
            .unwrap();

    let sam = fs::read_to_string(&sam_path).unwrap();
    let pg_line = sam
        .lines()
        .find(|l| l.starts_with("@PG"))
        .expect("SAM output must contain @PG header");

    assert!(pg_line.contains("PN:japalitysplice"));
    assert!(pg_line.contains("VN:"));

    fs::remove_dir_all(&tmp).unwrap();
}
