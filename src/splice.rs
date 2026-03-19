// Splice junction detection via empirical motif scanning
// Copyright (C) 2024-2026 Japality Limited
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::{Contig, SeedHit};
use rustc_hash::FxHashMap;

#[derive(Clone, Copy)]
struct JunctionCandidate {
    intron_start: usize,
    intron_end: usize,
    motif_bonus: f64,
}

/// Search for nearest canonical splice motif (GT-AG, GC-AG) within a window
fn find_canonical_junction(
    contig_seq: &[u8],
    est_intron_start: usize,
    est_intron_end: usize,
    window: usize,
) -> Option<(usize, usize)> {
    let seq_len = contig_seq.len();
    let lo_start = est_intron_start.saturating_sub(window);
    let hi_start = (est_intron_start + window).min(seq_len.saturating_sub(2));
    let lo_end = est_intron_end.saturating_sub(window);
    let hi_end = (est_intron_end + window).min(seq_len);

    let mut best: Option<(usize, usize, usize)> = None;

    for s in lo_start..=hi_start {
        if s + 2 > seq_len {
            break;
        }
        let d1 = contig_seq[s].to_ascii_uppercase();
        let d2 = contig_seq[s + 1].to_ascii_uppercase();
        let is_donor = (d1 == b'G' && d2 == b'T') || (d1 == b'G' && d2 == b'C');
        if !is_donor {
            continue;
        }

        for e in lo_end..=hi_end {
            if e < 2 || e > seq_len {
                continue;
            }
            let a1 = contig_seq[e - 2].to_ascii_uppercase();
            let a2 = contig_seq[e - 1].to_ascii_uppercase();
            let is_acceptor = a1 == b'A' && a2 == b'G';
            if !is_acceptor {
                continue;
            }

            let intron_len = e.saturating_sub(s);
            if intron_len < 30 || intron_len > 10_000 {
                continue;
            }

            let dist = (s as isize - est_intron_start as isize).unsigned_abs()
                + (e as isize - est_intron_end as isize).unsigned_abs();

            if best.is_none() || dist < best.unwrap().2 {
                best = Some((s, e, dist));
            }
        }
    }
    best.map(|(s, e, _)| (s, e))
}

fn upsert_candidate(candidates: &mut Vec<JunctionCandidate>, candidate: JunctionCandidate) {
    if let Some(existing) = candidates.iter_mut().find(|existing| {
        existing.intron_start == candidate.intron_start
            && existing.intron_end == candidate.intron_end
    }) {
        existing.motif_bonus = existing.motif_bonus.max(candidate.motif_bonus);
    } else {
        candidates.push(candidate);
    }
}

pub fn find_splice(
    sequence: &str,
    contigs: &[Contig],
    seed_index: &FxHashMap<u64, Vec<SeedHit>>,
    kmer_size: usize,
) -> Option<(usize, usize, String, usize)> {
    let seq_b = sequence.as_bytes();
    let q_len = sequence.len();

    if q_len < kmer_size * 3 {
        return None;
    }

    let min_exon = kmer_size + 6;
    let step = ((q_len - 2 * min_exon) / 10).max(3);
    let mut split_points: Vec<usize> = Vec::new();
    let mut s = min_exon;
    while s <= q_len - min_exon {
        split_points.push(s);
        s += step;
    }

    let mut best: Option<(usize, usize, String, usize, f64)> = None;

    for &split in &split_points {
        let left_positions: Vec<usize> = vec![2, split.saturating_sub(kmer_size + 2), split / 2];
        let right_positions: Vec<usize> = vec![
            split + 2,
            split + (q_len - split) / 2,
            q_len.saturating_sub(kmer_size + 2),
        ];

        for &left_pos in &left_positions {
            if left_pos + kmer_size > split {
                continue;
            }
            let left_key = match crate::kmer_to_u64(&seq_b[left_pos..left_pos + kmer_size]) {
                Some(k) => k,
                None => continue,
            };
            let left_hits = match seed_index.get(&left_key) {
                Some(h) if h.len() <= 20 => h,
                _ => continue,
            };

            for &right_pos in &right_positions {
                if right_pos < split || right_pos + kmer_size > q_len {
                    continue;
                }
                let right_key = match crate::kmer_to_u64(&seq_b[right_pos..right_pos + kmer_size]) {
                    Some(k) => k,
                    None => continue,
                };
                let right_hits = match seed_index.get(&right_key) {
                    Some(h) if h.len() <= 20 => h,
                    _ => continue,
                };

                for l_hit in left_hits.iter().take(8) {
                    for r_hit in right_hits.iter().take(8) {
                        if l_hit.contig_index != r_hit.contig_index {
                            continue;
                        }

                        let contig_index = l_hit.contig_index as usize;
                        let left_ref_pos = l_hit.position as usize;
                        let right_ref_pos = r_hit.position as usize;
                        let contig = &contigs[contig_index];
                        let contig_seq_bytes = contig.sequence.as_bytes();

                        let ref_start = left_ref_pos.saturating_sub(left_pos);
                        let est_intron_start = ref_start + split;
                        let est_intron_end_raw =
                            right_ref_pos.saturating_sub(right_pos.saturating_sub(split));

                        if est_intron_end_raw <= est_intron_start {
                            continue;
                        }
                        let est_intron_size = est_intron_end_raw - est_intron_start;
                        if est_intron_size < 30 || est_intron_size > 10_000 {
                            continue;
                        }

                        let mut candidates: Vec<JunctionCandidate> = Vec::new();

                        // Canonical motif-based junction
                        if let Some((junction_start, junction_end)) = find_canonical_junction(
                            contig_seq_bytes,
                            est_intron_start,
                            est_intron_end_raw,
                            10,
                        ) {
                            upsert_candidate(
                                &mut candidates,
                                JunctionCandidate {
                                    intron_start: junction_start,
                                    intron_end: junction_end,
                                    motif_bonus: 0.05,
                                },
                            );
                        }

                        // Known annotation-guided junctions
                        for (junction_start, junction_end) in crate::get_known_junction_candidates(
                            &contig.name,
                            est_intron_start,
                            est_intron_end_raw,
                            12,
                            16,
                        ) {
                            let boundary_dist = junction_start.abs_diff(est_intron_start)
                                + junction_end.abs_diff(est_intron_end_raw);
                            let known_bonus = (0.12 - (boundary_dist as f64 / 250.0)).max(0.03);
                            upsert_candidate(
                                &mut candidates,
                                JunctionCandidate {
                                    intron_start: junction_start,
                                    intron_end: junction_end,
                                    motif_bonus: known_bonus,
                                },
                            );
                        }

                        if candidates.is_empty() {
                            continue;
                        }

                        for candidate in candidates {
                            let intron_start = candidate.intron_start;
                            let intron_end = candidate.intron_end;
                            let intron_size = intron_end - intron_start;

                            let delta = intron_start as isize - est_intron_start as isize;
                            let exon1_len = (split as isize + delta) as usize;
                            let exon2_len = q_len.saturating_sub(exon1_len);
                            if exon1_len < kmer_size || exon2_len < kmer_size || exon1_len >= q_len
                            {
                                continue;
                            }

                            // Verify exon alignment quality
                            let e1_check = exon1_len.min(25);
                            let mut matches_e1 = 0usize;
                            for i in 0..e1_check {
                                let rp = intron_start.checked_sub(e1_check).unwrap_or(0) + i;
                                let qp = exon1_len - e1_check + i;
                                if rp < contig_seq_bytes.len()
                                    && qp < q_len
                                    && contig_seq_bytes[rp].to_ascii_uppercase()
                                        == seq_b[qp].to_ascii_uppercase()
                                {
                                    matches_e1 += 1;
                                }
                            }
                            let e2_check = exon2_len.min(25);
                            let mut matches_e2 = 0usize;
                            for i in 0..e2_check {
                                let rp = intron_end + i;
                                let qp = exon1_len + i;
                                if rp < contig_seq_bytes.len()
                                    && qp < q_len
                                    && contig_seq_bytes[rp].to_ascii_uppercase()
                                        == seq_b[qp].to_ascii_uppercase()
                                {
                                    matches_e2 += 1;
                                }
                            }
                            let match_rate =
                                (matches_e1 + matches_e2) as f64 / (e1_check + e2_check) as f64;
                            if match_rate < 0.70 {
                                continue;
                            }

                            let combined_score = match_rate + candidate.motif_bonus;

                            if best.is_none() || combined_score > best.as_ref().unwrap().4 {
                                let nm = (e1_check - matches_e1) + (e2_check - matches_e2);
                                let cigar = format!("{}M{}N{}M", exon1_len, intron_size, exon2_len);
                                best = Some((contig_index, ref_start, cigar, nm, combined_score));
                            }
                        }
                    }
                }
            }
        }
    }

    best.map(|(ci, pos, cigar, nm, _score)| (ci, pos, cigar, nm))
}
