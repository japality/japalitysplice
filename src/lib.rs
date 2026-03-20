// JapalitySplice — sparse k-mer RNA-Seq splice-aware aligner
// Copyright (C) 2026 Japality Limited
// SPDX-License-Identifier: GPL-3.0-or-later

use chrono::Utc;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::de::{Deserializer, MapAccess, Visitor};
use serde::ser::{SerializeMap, Serializer};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::RwLock;

mod dp;
mod splice;

static KNOWN_JUNCTIONS: std::sync::LazyLock<RwLock<Option<HashMap<String, Vec<KnownJunction>>>>> =
    std::sync::LazyLock::new(|| RwLock::new(None));

#[derive(Debug, Clone, Deserialize)]
struct Contig {
    name: String,
    sequence: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
struct SeedHit {
    contig_index: u32,
    position: u32,
}

impl SeedHit {
    fn new(contig_index: usize, position: usize) -> Result<Self, String> {
        Ok(Self {
            contig_index: u32::try_from(contig_index)
                .map_err(|_| format!("contig index {} exceeds u32 range", contig_index))?,
            position: u32::try_from(position)
                .map_err(|_| format!("position {} exceeds u32 range", position))?,
        })
    }

    fn contig_index_usize(self) -> usize {
        self.contig_index as usize
    }

    fn position_usize(self) -> usize {
        self.position as usize
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct KnownJunction {
    intron_start: usize,
    intron_end: usize,
}

#[derive(Debug, Clone)]
struct Alignment {
    mapped: bool,
    contig_name: String,
    one_based_pos: usize,
    mismatches: usize,
    cigar: String,
    is_reverse: bool,
}

struct CandidateCluster {
    contig_index: usize,
    seed_offset: usize,
    ref_pos: usize,
    start_sum: usize,
    support: usize,
}

#[derive(Debug, Clone)]
struct FastqRecord {
    header: String,
    sequence: String,
    quality: String,
}

#[derive(Debug, Clone, Default)]
struct BaseCounts {
    a: usize,
    t: usize,
    c: usize,
    g: usize,
    n: usize,
}

impl BaseCounts {
    fn push_base(&mut self, b: u8) {
        match b {
            b'A' => self.a += 1,
            b'T' => self.t += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'N' => self.n += 1,
            _ => {}
        }
    }

    fn into_hash_map(self) -> HashMap<String, usize> {
        let mut m = HashMap::new();
        m.insert("A".to_string(), self.a);
        m.insert("T".to_string(), self.t);
        m.insert("C".to_string(), self.c);
        m.insert("G".to_string(), self.g);
        m.insert("N".to_string(), self.n);
        m
    }
}

struct FastaLoadResult {
    contigs: Vec<Contig>,
    total_bases: usize,
    base_counts: HashMap<String, usize>,
}

// --- Serde helpers for index format ---

#[derive(Serialize)]
struct ContigJsonRef<'a> {
    name: &'a str,
    length: usize,
    sequence: &'a str,
}

#[derive(Serialize)]
struct NsixIndex<'a> {
    format: &'static str,
    created_at: String,
    source_fasta: String,
    kmer_size: usize,
    seed_stride: usize,
    contig_count: usize,
    total_bases: usize,
    base_counts: HashMap<String, usize>,
    contigs: Vec<ContigJsonRef<'a>>,
    seed_index: SeedIndexJson<'a>,
}

#[derive(Serialize)]
struct SeedIndexJson<'a> {
    kmer_size: usize,
    #[serde(serialize_with = "serialize_seed_index")]
    hits_by_kmer: &'a FxHashMap<u64, Vec<SeedHit>>,
}

fn serialize_seed_index<S>(
    map: &&FxHashMap<u64, Vec<SeedHit>>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut s = serializer.serialize_map(Some(map.len()))?;
    for (key, value) in *map {
        s.serialize_entry(&key.to_string(), value)?;
    }
    s.end()
}

#[derive(Deserialize)]
struct LoadedNsixIndex {
    contigs: Vec<Contig>,
    kmer_size: usize,
    seed_stride: usize,
    seed_index: LoadedSeedIndex,
}

#[derive(Deserialize)]
struct LoadedSeedIndex {
    #[allow(dead_code)]
    kmer_size: usize,
    #[serde(deserialize_with = "deserialize_seed_index")]
    hits_by_kmer: FxHashMap<u64, Vec<SeedHit>>,
}

impl std::ops::Deref for LoadedSeedIndex {
    type Target = FxHashMap<u64, Vec<SeedHit>>;
    fn deref(&self) -> &Self::Target {
        &self.hits_by_kmer
    }
}

#[derive(Deserialize)]
#[serde(untagged)]
enum SeedHitsCompat {
    Flat(Vec<u32>),
    Structured(Vec<SeedHit>),
}

impl SeedHitsCompat {
    fn into_seed_hits(self) -> Vec<SeedHit> {
        match self {
            SeedHitsCompat::Structured(v) => v,
            SeedHitsCompat::Flat(flat) => flat
                .chunks(2)
                .filter_map(|pair| {
                    if pair.len() == 2 {
                        Some(SeedHit {
                            contig_index: pair[0],
                            position: pair[1],
                        })
                    } else {
                        None
                    }
                })
                .collect(),
        }
    }
}

// --- Core algorithm ---

pub fn build_index(
    fasta_path: &str,
    output_dir: Option<&str>,
    kmer_size: usize,
) -> Result<String, String> {
    if !(kmer_size == 9 || kmer_size == 11 || kmer_size == 13) {
        return Err("kmer size must be 9/11/13".to_string());
    }

    let fasta = read_fasta_contigs(fasta_path)?;
    if fasta.contigs.is_empty() {
        return Err("no contigs".to_string());
    }

    let seed_stride = 4;
    let mut seed_index: FxHashMap<u64, Vec<SeedHit>> = FxHashMap::default();
    for (contig_index, contig) in fasta.contigs.iter().enumerate() {
        if contig.sequence.len() < kmer_size {
            continue;
        }
        let seq = contig.sequence.as_bytes();
        let max_pos = contig.sequence.len() - kmer_size;
        let mut pos = 0usize;
        while pos <= max_pos {
            if let Some(kmer) = kmer_to_u64(&seq[pos..pos + kmer_size]) {
                seed_index
                    .entry(kmer)
                    .or_default()
                    .push(SeedHit::new(contig_index, pos)?);
            }
            pos += seed_stride;
        }
    }

    let contig_json: Vec<ContigJsonRef<'_>> = fasta
        .contigs
        .iter()
        .map(|c| ContigJsonRef {
            name: &c.name,
            length: c.sequence.len(),
            sequence: &c.sequence,
        })
        .collect();

    let index = NsixIndex {
        format: "japality_splice.nsix.v2",
        created_at: Utc::now().to_rfc3339(),
        source_fasta: fasta_path.to_string(),
        kmer_size,
        seed_stride,
        contig_count: contig_json.len(),
        total_bases: fasta.total_bases,
        base_counts: fasta.base_counts,
        contigs: contig_json,
        seed_index: SeedIndexJson {
            kmer_size,
            hits_by_kmer: &seed_index,
        },
    };

    let source_name = Path::new(fasta_path)
        .file_name()
        .and_then(|x| x.to_str())
        .unwrap_or("reference.fasta")
        .replace(".fasta", ".nsix")
        .replace(".fa", ".nsix")
        .replace(".fna", ".nsix");

    let out_path = if let Some(dir) = output_dir {
        let mut p = PathBuf::from(dir);
        p.push(source_name);
        p
    } else {
        let mut p = PathBuf::from(fasta_path);
        p.set_extension("nsix");
        p
    };

    if let Some(parent) = out_path.parent() {
        let _ = fs::create_dir_all(parent);
    }

    let out_file = File::create(&out_path).map_err(|e| e.to_string())?;
    let writer = BufWriter::with_capacity(1024 * 1024, out_file);
    serde_json::to_writer(writer, &index).map_err(|e| e.to_string())?;
    Ok(out_path.to_string_lossy().to_string())
}

pub fn run_alignment(
    fastq_path_r1: &str,
    fastq_path_r2: Option<&str>,
    index_path: &str,
    output_dir: &str,
    max_mismatches: usize,
) -> Result<String, String> {
    let index_file = File::open(index_path).map_err(|e| e.to_string())?;
    let reader = BufReader::with_capacity(1024 * 1024, index_file);
    let mut index: LoadedNsixIndex = serde_json::from_reader(reader).map_err(|e| e.to_string())?;
    sanitize_contig_names(&mut index.contigs);
    let contigs = index.contigs;
    let seed_index = index.seed_index;
    let kmer_size = index.kmer_size;
    let seed_stride = index.seed_stride;

    let mut out_dir = PathBuf::from(output_dir);
    let _ = fs::create_dir_all(&out_dir);
    let stamp = Utc::now().format("%Y%m%d_%H%M%S").to_string();
    out_dir.push(format!("alignment_{}.sam", stamp));

    let out_file = fs::File::create(&out_dir).map_err(|e| e.to_string())?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);

    writeln!(writer, "@HD\tVN:1.6\tSO:unknown").map_err(|e| e.to_string())?;
    for contig in &contigs {
        writeln!(
            writer,
            "@SQ\tSN:{}\tLN:{}",
            contig.name,
            contig.sequence.len()
        )
        .map_err(|e| e.to_string())?;
    }
    writeln!(
        writer,
        "@PG\tID:japalitysplice\tPN:japalitysplice\tVN:{}\tCL:japalitysplice align",
        env!("CARGO_PKG_VERSION")
    )
    .map_err(|e| e.to_string())?;

    if let Some(r2) = fastq_path_r2 {
        let mut process_pair_batch =
            |batch: Vec<(FastqRecord, FastqRecord)>| -> Result<(), String> {
                let lines: Vec<Result<(String, String), String>> = batch
                    .into_par_iter()
                    .map(|(record_r1, record_r2)| {
                        let qname = paired_qname(&record_r1, &record_r2)?;
                        let alignment_r1 = align_read(
                            &record_r1.sequence,
                            &contigs,
                            &seed_index,
                            kmer_size,
                            seed_stride,
                            max_mismatches,
                        );
                        let alignment_r2 = align_read(
                            &record_r2.sequence,
                            &contigs,
                            &seed_index,
                            kmer_size,
                            seed_stride,
                            max_mismatches,
                        );
                        Ok((
                            to_sam_line(
                                &record_r1,
                                &alignment_r1,
                                Some(&alignment_r2),
                                true,
                                true,
                                &qname,
                            ),
                            to_sam_line(
                                &record_r2,
                                &alignment_r2,
                                Some(&alignment_r1),
                                true,
                                false,
                                &qname,
                            ),
                        ))
                    })
                    .collect();

                for pair_result in lines {
                    let (line_r1, line_r2) = pair_result?;
                    writer
                        .write_all(line_r1.as_bytes())
                        .and_then(|_| writer.write_all(b"\n"))
                        .map_err(|e| e.to_string())?;
                    writer
                        .write_all(line_r2.as_bytes())
                        .and_then(|_| writer.write_all(b"\n"))
                        .map_err(|e| e.to_string())?;
                }
                Ok(())
            };

        process_paired_fastq_batches(fastq_path_r1, r2, 4096, &mut process_pair_batch)?;
    } else {
        let mut process_batch = |batch: Vec<FastqRecord>| -> Result<(), String> {
            let lines: Vec<String> = batch
                .into_par_iter()
                .map(|record| {
                    let alignment = align_read(
                        &record.sequence,
                        &contigs,
                        &seed_index,
                        kmer_size,
                        seed_stride,
                        max_mismatches,
                    );
                    let qname = extract_qname(&record.header).to_string();
                    to_sam_line(&record, &alignment, None, false, true, &qname)
                })
                .collect();

            for line in lines {
                writer
                    .write_all(line.as_bytes())
                    .and_then(|_| writer.write_all(b"\n"))
                    .map_err(|e| e.to_string())?;
            }
            Ok(())
        };

        process_fastq_batches(fastq_path_r1, 4096, &mut process_batch)?;
    }

    writer.flush().map_err(|e| e.to_string())?;
    Ok(out_dir.to_string_lossy().to_string())
}

pub fn init_known_junctions(annotation_path: &str) -> Result<usize, String> {
    let file = fs::File::open(annotation_path).map_err(|e| e.to_string())?;
    let reader = BufReader::with_capacity(1024 * 1024, file);
    let mut transcripts: HashMap<String, (String, String, Vec<(usize, usize)>)> = HashMap::new();

    for raw_line in reader.lines() {
        let line = raw_line.map_err(|e| e.to_string())?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if !fields[2].eq_ignore_ascii_case("exon") {
            continue;
        }

        let start = match fields[3].parse::<usize>() {
            Ok(value) => value,
            Err(_) => continue,
        };
        let end = match fields[4].parse::<usize>() {
            Ok(value) => value,
            Err(_) => continue,
        };
        if start > end {
            continue;
        }

        let chrom = normalize_reference_name(fields[0]);
        let strand = fields[6].to_string();
        let attributes = parse_annotation_attributes(fields[8]);
        let transcript_id = attributes
            .get("transcript_id")
            .or_else(|| attributes.get("Parent"))
            .or_else(|| attributes.get("transcript"))
            .or_else(|| attributes.get("ID"))
            .or_else(|| attributes.get("gene_id"))
            .cloned();
        let Some(transcript_id) = transcript_id else {
            continue;
        };

        let transcript_key = format!("{}|{}|{}", chrom, strand, transcript_id);
        let entry = transcripts
            .entry(transcript_key)
            .or_insert_with(|| (chrom.clone(), strand.clone(), Vec::new()));
        entry.2.push((start, end));
    }

    let mut junctions_by_chrom: HashMap<String, Vec<KnownJunction>> = HashMap::new();
    for (_, (chrom, _strand, mut exons)) in transcripts {
        exons.sort_unstable();
        exons.dedup();
        if exons.len() < 2 {
            continue;
        }

        for pair in exons.windows(2) {
            let left_exon = pair[0];
            let right_exon = pair[1];
            let intron_start_one_based = left_exon.1.saturating_add(1);
            let intron_end_one_based = right_exon.0.saturating_sub(1);
            if intron_start_one_based > intron_end_one_based {
                continue;
            }

            let intron_start = intron_start_one_based.saturating_sub(1);
            let intron_end = intron_end_one_based;
            if intron_end <= intron_start {
                continue;
            }

            junctions_by_chrom
                .entry(chrom.clone())
                .or_default()
                .push(KnownJunction {
                    intron_start,
                    intron_end,
                });
        }
    }

    let mut total_junctions = 0usize;
    for junctions in junctions_by_chrom.values_mut() {
        junctions.sort_by(|lhs, rhs| {
            lhs.intron_start
                .cmp(&rhs.intron_start)
                .then(lhs.intron_end.cmp(&rhs.intron_end))
        });
        junctions.dedup_by(|lhs, rhs| {
            lhs.intron_start == rhs.intron_start && lhs.intron_end == rhs.intron_end
        });
        total_junctions += junctions.len();
    }

    let mut guard = KNOWN_JUNCTIONS.write().unwrap();
    *guard = Some(junctions_by_chrom);
    Ok(total_junctions)
}

pub fn clear_known_junctions() {
    let mut guard = KNOWN_JUNCTIONS.write().unwrap();
    *guard = None;
}

pub(crate) fn get_known_junction_candidates(
    contig_name: &str,
    est_intron_start: usize,
    est_intron_end: usize,
    window: usize,
    limit: usize,
) -> Vec<(usize, usize)> {
    let guard = KNOWN_JUNCTIONS.read().unwrap();
    let Some(map) = guard.as_ref() else {
        return Vec::new();
    };
    let normalized_name = normalize_reference_name(contig_name);
    let Some(junctions) = map.get(&normalized_name) else {
        return Vec::new();
    };

    let lo = est_intron_start.saturating_sub(window);
    let hi = est_intron_start.saturating_add(window);
    let start_idx = junctions.partition_point(|junction| junction.intron_start < lo);

    let mut out = Vec::new();
    for junction in junctions.iter().skip(start_idx) {
        if junction.intron_start > hi {
            break;
        }
        if junction.intron_end.abs_diff(est_intron_end) <= window {
            out.push((junction.intron_start, junction.intron_end));
            if out.len() >= limit {
                break;
            }
        }
    }
    out
}

// --- FASTQ I/O ---

fn process_fastq_batches<F>(path: &str, batch_size: usize, process: &mut F) -> Result<(), String>
where
    F: FnMut(Vec<FastqRecord>) -> Result<(), String>,
{
    let reader = open_fastq_reader(path)?;
    let mut lines = reader.lines();
    let mut batch: Vec<FastqRecord> = Vec::with_capacity(batch_size);

    loop {
        let Some(record) = read_fastq_record(&mut lines)? else {
            break;
        };
        batch.push(record);
        if batch.len() >= batch_size {
            process(std::mem::take(&mut batch))?;
        }
    }

    if !batch.is_empty() {
        process(batch)?;
    }

    Ok(())
}

fn process_paired_fastq_batches<F>(
    path_r1: &str,
    path_r2: &str,
    batch_size: usize,
    process: &mut F,
) -> Result<(), String>
where
    F: FnMut(Vec<(FastqRecord, FastqRecord)>) -> Result<(), String>,
{
    let reader_r1 = open_fastq_reader(path_r1)?;
    let reader_r2 = open_fastq_reader(path_r2)?;
    let mut lines_r1 = reader_r1.lines();
    let mut lines_r2 = reader_r2.lines();
    let mut batch: Vec<(FastqRecord, FastqRecord)> = Vec::with_capacity(batch_size);

    loop {
        let record_r1 = read_fastq_record(&mut lines_r1)?;
        let record_r2 = read_fastq_record(&mut lines_r2)?;
        match (record_r1, record_r2) {
            (Some(r1), Some(r2)) => {
                batch.push((r1, r2));
                if batch.len() >= batch_size {
                    process(std::mem::take(&mut batch))?;
                }
            }
            (None, None) => break,
            _ => return Err("paired FASTQ files have different read counts".to_string()),
        }
    }

    if !batch.is_empty() {
        process(batch)?;
    }

    Ok(())
}

fn read_fastq_record<I>(lines: &mut I) -> Result<Option<FastqRecord>, String>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    let Some(header) = lines.next() else {
        return Ok(None);
    };
    let header = header.map_err(|e| e.to_string())?;
    if !header.starts_with('@') {
        return Err("invalid FASTQ header".to_string());
    }

    let sequence_line = lines
        .next()
        .ok_or_else(|| "incomplete FASTQ record".to_string())?
        .map_err(|e| e.to_string())?;
    let separator = lines
        .next()
        .ok_or_else(|| "incomplete FASTQ record".to_string())?
        .map_err(|e| e.to_string())?;
    if !separator.starts_with('+') {
        return Err("invalid FASTQ separator".to_string());
    }
    let quality_line = lines
        .next()
        .ok_or_else(|| "incomplete FASTQ record".to_string())?
        .map_err(|e| e.to_string())?;
    let sequence = normalize_dna(&sequence_line);
    if sequence.len() != quality_line.len() {
        return Err("sequence/quality mismatch".to_string());
    }

    Ok(Some(FastqRecord {
        header,
        sequence,
        quality: quality_line,
    }))
}

fn open_fastq_reader(path: &str) -> Result<Box<dyn BufRead>, String> {
    let file = fs::File::open(path).map_err(|e| e.to_string())?;
    if path.to_ascii_lowercase().ends_with(".gz") {
        let decoder = MultiGzDecoder::new(file);
        return Ok(Box::new(BufReader::with_capacity(1024 * 1024, decoder)));
    }
    Ok(Box::new(BufReader::with_capacity(1024 * 1024, file)))
}

// --- K-mer encoding ---

pub fn kmer_to_u64(kmer: &[u8]) -> Option<u64> {
    let mut val: u64 = 0;
    for &b in kmer {
        let bits = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => return None,
        };
        val = (val << 2) | bits;
    }
    Some(val)
}

fn parse_seed_key(raw: &str) -> Option<u64> {
    raw.parse::<u64>()
        .ok()
        .or_else(|| kmer_to_u64(raw.as_bytes()))
}

fn deserialize_seed_index<'de, D>(deserializer: D) -> Result<FxHashMap<u64, Vec<SeedHit>>, D::Error>
where
    D: Deserializer<'de>,
{
    struct SeedIndexVisitor;

    impl<'de> Visitor<'de> for SeedIndexVisitor {
        type Value = FxHashMap<u64, Vec<SeedHit>>;

        fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            formatter.write_str("a JSON object mapping k-mers to hit arrays")
        }

        fn visit_map<M>(self, mut access: M) -> Result<Self::Value, M::Error>
        where
            M: MapAccess<'de>,
        {
            let mut out = FxHashMap::default();
            while let Some(kmer) = access.next_key::<String>()? {
                let hits = access.next_value::<SeedHitsCompat>()?.into_seed_hits();
                let Some(key) = parse_seed_key(&kmer) else {
                    continue;
                };
                if !hits.is_empty() {
                    out.insert(key, hits);
                }
            }
            Ok(out)
        }
    }

    deserializer.deserialize_map(SeedIndexVisitor)
}

// --- FASTA I/O ---

fn sanitize_contig_names(contigs: &mut [Contig]) {
    for contig in contigs {
        let sanitized = contig
            .name
            .split_whitespace()
            .next()
            .unwrap_or(&contig.name)
            .to_string();
        contig.name = sanitized;
    }
}

fn append_normalized_dna(
    raw: &str,
    sequence: &mut String,
    base_counts: &mut BaseCounts,
    total_bases: &mut usize,
) {
    for byte in raw.bytes() {
        let normalized = match byte {
            b'A' | b'a' => b'A',
            b'T' | b't' => b'T',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'N' | b'n' => b'N',
            _ => continue,
        };
        sequence.push(normalized as char);
        base_counts.push_base(normalized);
        *total_bases += 1;
    }
}

fn flush_contig(
    contigs: &mut Vec<Contig>,
    current_name: &mut Option<String>,
    current_seq: &mut String,
) {
    if let Some(name) = current_name.take() {
        if !current_seq.is_empty() {
            contigs.push(Contig {
                name,
                sequence: std::mem::take(current_seq),
            });
        }
    }
}

fn read_fasta_contigs(path: &str) -> Result<FastaLoadResult, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let reader = BufReader::new(file);
    let mut out = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();
    let mut total_bases = 0usize;
    let mut base_counts = BaseCounts::default();

    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            flush_contig(&mut out, &mut current_name, &mut current_seq);
            current_name = Some(
                rest.trim()
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
            );
        } else {
            append_normalized_dna(line, &mut current_seq, &mut base_counts, &mut total_bases);
        }
    }

    flush_contig(&mut out, &mut current_name, &mut current_seq);

    Ok(FastaLoadResult {
        contigs: out,
        total_bases,
        base_counts: base_counts.into_hash_map(),
    })
}

// --- Alignment core ---

fn extend_seed_local(
    contig: &[u8],
    ref_pos: usize,
    query: &[u8],
    seed_offset: usize,
    kmer_size: usize,
    max_miss: usize,
) -> Option<(usize, usize, String)> {
    let allowed_mm = max_miss.min(3);

    let mut r_right = ref_pos + kmer_size;
    let mut q_right = seed_offset + kmer_size;
    let mut score: i32 = 0;
    let mut best_score: i32 = 0;
    let mut best_r_right = r_right;
    let mut best_q_right = q_right;

    while r_right < contig.len() && q_right < query.len() {
        if contig[r_right] == query[q_right] {
            score += 1;
        } else {
            score -= 2;
        }
        r_right += 1;
        q_right += 1;
        if score >= best_score {
            best_score = score;
            best_r_right = r_right;
            best_q_right = q_right;
        }
        if score < best_score - 10 {
            break;
        }
    }

    let mut r_left = ref_pos;
    let mut q_left = seed_offset;
    let mut score_l: i32 = 0;
    let mut best_score_l: i32 = 0;
    let mut best_r_left = r_left;
    let mut best_q_left = q_left;

    while r_left > 0 && q_left > 0 {
        if contig[r_left - 1] == query[q_left - 1] {
            score_l += 1;
        } else {
            score_l -= 2;
        }
        r_left -= 1;
        q_left -= 1;
        if score_l >= best_score_l {
            best_score_l = score_l;
            best_r_left = r_left;
            best_q_left = q_left;
        }
        if score_l < best_score_l - 10 {
            break;
        }
    }

    let total_aligned = best_r_right - best_r_left;

    let mut actual_mm = 0;
    for i in 0..total_aligned {
        if contig[best_r_left + i] != query[best_q_left + i] {
            actual_mm += 1;
        }
    }

    let min_aligned = (query.len() * 70) / 100;
    if total_aligned < min_aligned {
        return None;
    }
    if actual_mm > allowed_mm {
        return None;
    }

    let mut cigar = String::new();
    if best_q_left > 0 {
        cigar.push_str(&format!("{}S", best_q_left));
    }
    cigar.push_str(&format!("{}M", total_aligned));
    let soft_clip_end = query.len() - best_q_right;
    if soft_clip_end > 0 {
        cigar.push_str(&format!("{}S", soft_clip_end));
    }
    Some((best_r_left, actual_mm, cigar))
}

fn cigar_soft_clip_bases(cigar: &str) -> usize {
    let mut total = 0usize;
    let mut value = 0usize;

    for ch in cigar.chars() {
        if let Some(digit) = ch.to_digit(10) {
            value = value.saturating_mul(10).saturating_add(digit as usize);
        } else {
            if ch == 'S' {
                total = total.saturating_add(value);
            }
            value = 0;
        }
    }

    total
}

fn should_replace_alignment(
    best_alignment: &Option<Alignment>,
    best_support: usize,
    candidate: &Alignment,
    candidate_support: usize,
) -> bool {
    let Some(best) = best_alignment.as_ref() else {
        return true;
    };

    if candidate.mismatches != best.mismatches {
        return candidate.mismatches < best.mismatches;
    }

    let candidate_soft_clip = cigar_soft_clip_bases(&candidate.cigar);
    let best_soft_clip = cigar_soft_clip_bases(&best.cigar);
    if candidate_soft_clip != best_soft_clip {
        return candidate_soft_clip < best_soft_clip;
    }

    if candidate_support != best_support {
        return candidate_support > best_support;
    }

    candidate.one_based_pos < best.one_based_pos
}

fn try_splice_mapping(
    sequence: &str,
    contigs: &[Contig],
    seed_index: &FxHashMap<u64, Vec<SeedHit>>,
    kmer_size: usize,
    is_reverse: bool,
) -> Option<Alignment> {
    if let Some((contig_index, pos, cigar, mismatches)) =
        splice::find_splice(sequence, contigs, seed_index, kmer_size)
    {
        let contig = &contigs[contig_index];
        return Some(Alignment {
            mapped: true,
            contig_name: contig.name.clone(),
            one_based_pos: pos + 1,
            mismatches,
            cigar,
            is_reverse,
        });
    }
    None
}

fn align_read_directional(
    sequence: &str,
    contigs: &[Contig],
    seed_index: &FxHashMap<u64, Vec<SeedHit>>,
    kmer_size: usize,
    seed_stride: usize,
    max_mismatches: usize,
    is_reverse: bool,
) -> Alignment {
    let mut best_alignment: Option<Alignment> = None;
    let mut best_support = 0usize;
    let seq_bytes = sequence.as_bytes();

    if sequence.len() >= kmer_size {
        let max_seed_offset = sequence.len() - kmer_size;
        let center_offset = sequence.len() / 2;
        let candidate_bin = seed_stride.max(1) * 2;
        let max_hits_per_seed = 64usize;
        let repetitive_seed_cutoff = 256usize;
        let mut candidates: FxHashMap<(usize, usize), CandidateCluster> = FxHashMap::default();

        for seed_offset in 0..=max_seed_offset {
            let seed_key = match kmer_to_u64(&seq_bytes[seed_offset..seed_offset + kmer_size]) {
                Some(k) => k,
                None => continue,
            };
            if let Some(hits) = seed_index.get(&seed_key) {
                if hits.len() > repetitive_seed_cutoff {
                    continue;
                }

                for hit in hits.iter().take(max_hits_per_seed) {
                    let contig_index = hit.contig_index_usize();
                    let ref_pos = hit.position_usize();
                    if contigs.get(contig_index).is_none() {
                        continue;
                    }

                    let approx_start = ref_pos.saturating_sub(seed_offset);
                    let bucket = approx_start / candidate_bin;
                    let key = (contig_index, bucket);
                    let entry = candidates.entry(key).or_insert_with(|| CandidateCluster {
                        contig_index,
                        seed_offset,
                        ref_pos,
                        start_sum: 0,
                        support: 0,
                    });
                    entry.support += 1;
                    entry.start_sum = entry.start_sum.saturating_add(approx_start);

                    if seed_offset.abs_diff(center_offset)
                        < entry.seed_offset.abs_diff(center_offset)
                    {
                        entry.seed_offset = seed_offset;
                        entry.ref_pos = ref_pos;
                    }
                }
            }
        }

        let mut ranked_candidates: Vec<CandidateCluster> = candidates.into_values().collect();
        ranked_candidates.sort_by(|lhs, rhs| {
            rhs.support
                .cmp(&lhs.support)
                .then_with(|| {
                    lhs.seed_offset
                        .abs_diff(center_offset)
                        .cmp(&rhs.seed_offset.abs_diff(center_offset))
                })
                .then_with(|| lhs.contig_index.cmp(&rhs.contig_index))
                .then_with(|| lhs.ref_pos.cmp(&rhs.ref_pos))
        });

        let min_support = if sequence.len() >= 75 { 2 } else { 1 };
        let mut examined = 0usize;
        for candidate in ranked_candidates {
            if candidate.support < min_support && examined >= 6 {
                break;
            }
            if examined >= 16 {
                break;
            }
            examined += 1;

            let Some(contig) = contigs.get(candidate.contig_index) else {
                continue;
            };

            if let Some((ref_start, miss, cigar)) = extend_seed_local(
                contig.sequence.as_bytes(),
                candidate.ref_pos,
                seq_bytes,
                candidate.seed_offset,
                kmer_size,
                max_mismatches,
            ) {
                return Alignment {
                    mapped: true,
                    contig_name: contig.name.clone(),
                    one_based_pos: ref_start + 1,
                    mismatches: miss,
                    cigar,
                    is_reverse,
                };
            }

            let approx_start = candidate.start_sum / candidate.support.max(1);
            let pad = 20 + (max_mismatches * 4);
            let expected_start = approx_start.saturating_sub(pad);
            let expected_end = (approx_start + sequence.len() + pad).min(contig.sequence.len());
            if expected_start >= expected_end {
                continue;
            }

            let window = &contig.sequence.as_bytes()[expected_start..expected_end];
            if let Some((start_ref, _end_ref, cigar, mismatches)) =
                dp::semi_global_align(window, seq_bytes, max_mismatches + 3)
            {
                if mismatches <= (max_mismatches + 3) {
                    let actual_pos = expected_start + start_ref;
                    let aln = Alignment {
                        mapped: true,
                        contig_name: contig.name.clone(),
                        one_based_pos: actual_pos + 1,
                        mismatches,
                        cigar,
                        is_reverse,
                    };
                    if mismatches <= max_mismatches
                        && cigar_soft_clip_bases(&aln.cigar) <= sequence.len() / 5
                    {
                        return aln;
                    }
                    if should_replace_alignment(
                        &best_alignment,
                        best_support,
                        &aln,
                        candidate.support,
                    ) {
                        best_support = candidate.support;
                        best_alignment = Some(aln);
                    }
                }
            }
        }
    }

    if best_alignment.is_none() && sequence.len() >= kmer_size * 2 {
        if let Some(aln) = try_splice_mapping(sequence, contigs, seed_index, kmer_size, is_reverse)
        {
            return aln;
        }
    }

    best_alignment.unwrap_or_else(|| Alignment {
        mapped: false,
        contig_name: "*".to_string(),
        one_based_pos: 0,
        mismatches: 0,
        cigar: "*".to_string(),
        is_reverse: false,
    })
}

fn align_read(
    sequence: &str,
    contigs: &[Contig],
    seed_index: &FxHashMap<u64, Vec<SeedHit>>,
    kmer_size: usize,
    seed_stride: usize,
    max_mismatches: usize,
) -> Alignment {
    let fwd = align_read_directional(
        sequence,
        contigs,
        seed_index,
        kmer_size,
        seed_stride,
        max_mismatches,
        false,
    );
    if fwd.mapped {
        return fwd;
    }

    let rev_seq = reverse_complement(sequence);
    let rev = align_read_directional(
        &rev_seq,
        contigs,
        seed_index,
        kmer_size,
        seed_stride,
        max_mismatches,
        true,
    );

    if rev.mapped {
        return rev;
    }

    fwd
}

// --- SAM output ---

fn extract_qname(header: &str) -> &str {
    header
        .trim_start_matches('@')
        .split_whitespace()
        .next()
        .unwrap_or("read")
}

fn normalize_paired_qname(qname: &str) -> &str {
    qname
        .strip_suffix("/1")
        .or_else(|| qname.strip_suffix("/2"))
        .unwrap_or(qname)
}

fn paired_qname(record_r1: &FastqRecord, record_r2: &FastqRecord) -> Result<String, String> {
    let qname_r1 = normalize_paired_qname(extract_qname(&record_r1.header));
    let qname_r2 = normalize_paired_qname(extract_qname(&record_r2.header));
    if qname_r1 != qname_r2 {
        return Err(format!(
            "paired FASTQ names do not match: {} vs {}",
            qname_r1, qname_r2
        ));
    }
    Ok(qname_r1.to_string())
}

fn cigar_reference_bases(cigar: &str) -> usize {
    let mut value = 0usize;
    let mut count = 0usize;
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + ch.to_digit(10).unwrap_or(0) as usize;
            continue;
        }
        if matches!(ch, 'M' | 'D' | 'N' | '=' | 'X') {
            value += count;
        }
        count = 0;
    }
    value
}

fn alignment_end_one_based(alignment: &Alignment) -> Option<usize> {
    if !alignment.mapped {
        return None;
    }
    let span = cigar_reference_bases(&alignment.cigar).max(1);
    Some(alignment.one_based_pos + span - 1)
}

fn is_proper_fr_pair(alignment: &Alignment, mate: &Alignment, is_read1: bool) -> bool {
    if !alignment.mapped
        || !mate.mapped
        || alignment.contig_name != mate.contig_name
        || alignment.is_reverse == mate.is_reverse
    {
        return false;
    }

    if is_read1 {
        !alignment.is_reverse && mate.is_reverse && alignment.one_based_pos <= mate.one_based_pos
    } else {
        alignment.is_reverse && !mate.is_reverse && mate.one_based_pos <= alignment.one_based_pos
    }
}

fn template_length(alignment: &Alignment, mate: &Alignment, is_read1: bool) -> i64 {
    if !alignment.mapped || !mate.mapped || alignment.contig_name != mate.contig_name {
        return 0;
    }
    let this_end = alignment_end_one_based(alignment).unwrap_or(alignment.one_based_pos);
    let mate_end = alignment_end_one_based(mate).unwrap_or(mate.one_based_pos);
    let left_start = alignment.one_based_pos.min(mate.one_based_pos);
    let right_end = this_end.max(mate_end);
    let span = right_end.saturating_sub(left_start) as i64 + 1;
    let this_is_leftmost = alignment.one_based_pos < mate.one_based_pos
        || (alignment.one_based_pos == mate.one_based_pos && is_read1);
    if this_is_leftmost {
        span
    } else {
        -span
    }
}

fn to_sam_line(
    record: &FastqRecord,
    alignment: &Alignment,
    mate: Option<&Alignment>,
    is_paired: bool,
    is_read1: bool,
    qname: &str,
) -> String {
    let mut flag = 0usize;
    if is_paired {
        flag |= 0x1;
        if is_read1 {
            flag |= 0x40;
        } else {
            flag |= 0x80;
        }
        if let Some(mate_alignment) = mate {
            if !mate_alignment.mapped {
                flag |= 0x8;
            } else if mate_alignment.is_reverse {
                flag |= 0x20;
            }
            if is_proper_fr_pair(alignment, mate_alignment, is_read1) {
                flag |= 0x2;
            }
        }
    }
    if !alignment.mapped {
        flag |= 0x4;
    } else if alignment.is_reverse {
        flag |= 0x10;
    }

    let mapq = if alignment.mapped {
        let soft_clip_bases = cigar_soft_clip_bases(&alignment.cigar);
        let read_len = record.sequence.len().max(1);
        let clip_frac = soft_clip_bases as f64 / read_len as f64;

        if alignment.mismatches == 0 && clip_frac < 0.05 {
            60
        } else if alignment.mismatches <= 1 && clip_frac < 0.10 {
            40
        } else if alignment.mismatches <= 2 && clip_frac < 0.15 {
            20
        } else {
            3
        }
    } else {
        0
    };

    let mate_rname = if let Some(mate_alignment) = mate {
        if mate_alignment.mapped {
            if alignment.mapped && alignment.contig_name == mate_alignment.contig_name {
                "=".to_string()
            } else {
                mate_alignment.contig_name.clone()
            }
        } else {
            "*".to_string()
        }
    } else {
        "*".to_string()
    };
    let mate_pos = mate
        .filter(|mate_alignment| mate_alignment.mapped)
        .map(|mate_alignment| mate_alignment.one_based_pos.to_string())
        .unwrap_or_else(|| "0".to_string());
    let tlen = mate
        .map(|mate_alignment| template_length(alignment, mate_alignment, is_read1).to_string())
        .unwrap_or_else(|| "0".to_string());

    let mut fields = vec![
        qname.to_string(),
        flag.to_string(),
        if alignment.mapped {
            alignment.contig_name.clone()
        } else {
            "*".to_string()
        },
        if alignment.mapped {
            alignment.one_based_pos.to_string()
        } else {
            "0".to_string()
        },
        mapq.to_string(),
        if alignment.mapped {
            alignment.cigar.clone()
        } else {
            "*".to_string()
        },
        mate_rname,
        mate_pos,
        tlen,
        if alignment.is_reverse {
            reverse_complement(&record.sequence)
        } else {
            record.sequence.clone()
        },
        if alignment.is_reverse {
            record.quality.chars().rev().collect()
        } else {
            record.quality.clone()
        },
    ];

    if alignment.mapped {
        fields.push(format!("NM:i:{}", alignment.mismatches));
    }

    fields.join("\t")
}

// --- Utilities ---

fn normalize_dna(sequence: &str) -> String {
    sequence
        .chars()
        .filter_map(|c| {
            let uc = c.to_ascii_uppercase();
            match uc {
                'A' | 'C' | 'G' | 'T' | 'N' => Some(uc),
                _ => None,
            }
        })
        .collect()
}

fn normalize_reference_name(name: &str) -> String {
    let trimmed = name.trim();
    if trimmed.len() >= 3 && trimmed[..3].eq_ignore_ascii_case("chr") {
        trimmed[3..].to_string()
    } else {
        trimmed.to_string()
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            _ => 'N',
        })
        .collect()
}

fn parse_annotation_attributes(attributes_raw: &str) -> HashMap<String, String> {
    let mut parsed = HashMap::new();
    for token in attributes_raw.split(';') {
        let item = token.trim();
        if item.is_empty() {
            continue;
        }
        if let Some((key, value)) = item.split_once('=') {
            let cleaned = value
                .trim()
                .trim_matches('"')
                .split(',')
                .next()
                .unwrap_or("")
                .to_string();
            parsed.insert(key.trim().to_string(), cleaned);
            continue;
        }
        if let Some((key, value)) = item.split_once(' ') {
            let cleaned = value
                .trim()
                .trim_matches('"')
                .split(',')
                .next()
                .unwrap_or("")
                .to_string();
            parsed.insert(key.trim().to_string(), cleaned);
        }
    }
    parsed
}
