// Dynamic programming local alignment
// Copyright (C) 2024-2026 Japality Limited
// SPDX-License-Identifier: GPL-3.0-or-later

use std::cell::RefCell;

thread_local! {
    static DP_BUFFER: RefCell<Vec<i32>> = RefCell::new(Vec::new());
}

pub fn semi_global_align(
    ref_seq: &[u8],
    query_seq: &[u8],
    max_mismatches: usize,
) -> Option<(usize, usize, String, usize)> {
    DP_BUFFER.with(|buf| {
        let mut dp = buf.borrow_mut();
        semi_global_align_impl(&mut dp, ref_seq, query_seq, max_mismatches)
    })
}

fn semi_global_align_impl(
    dp: &mut Vec<i32>,
    ref_seq: &[u8],
    query_seq: &[u8],
    max_mismatches: usize,
) -> Option<(usize, usize, String, usize)> {
    let m = query_seq.len();
    let n = ref_seq.len();
    if m == 0 || n == 0 {
        return None;
    }

    let cols = m + 1;
    let required = (n + 1) * cols;
    dp.resize(required, 0);
    dp[..required].fill(0);
    let match_score = 1i32;
    let mismatch_penalty = -1i32;
    let gap_penalty = -1i32;

    for i in 0..=n {
        dp[i * cols] = 0; // Free start
    }
    for j in 1..=m {
        dp[j] = dp[j - 1] + gap_penalty;
    }

    let mut max_score = -99999;
    let mut end_ref = 0;
    let mut end_query = m; // how many query bases are consumed
    let min_query_aligned = (m * 70) / 100; // allow right soft-clip up to 30%

    for i in 1..=n {
        let r_char = ref_seq[i - 1];
        let row_idx = i * cols;
        let prev_row_idx = (i - 1) * cols;

        for j in 1..=m {
            let q_char = query_seq[j - 1];
            let cost = if r_char == q_char {
                match_score
            } else {
                mismatch_penalty
            };

            let score_diag = dp[prev_row_idx + j - 1] + cost;
            let score_up = dp[prev_row_idx + j] + gap_penalty;
            let score_left = dp[row_idx + j - 1] + gap_penalty;

            let current_score = score_diag.max(score_up).max(score_left);
            dp[row_idx + j] = current_score;

            // Allow ending at j >= min_query_aligned (right soft-clip)
            if j >= min_query_aligned && current_score > max_score {
                max_score = current_score;
                end_ref = i;
                end_query = j;
            }
        }
    }

    let min_score_required = (end_query as i32 * match_score) - (max_mismatches as i32 * 2);
    if max_score < min_score_required.min((end_query as i32) / 2) {
        return None;
    }

    let mut i = end_ref;
    let mut j = end_query;
    let mut align_ops = Vec::new();
    let mut actual_mismatches = 0;

    // Add right soft-clip if we didn't consume the full query
    let right_clip = m - end_query;

    while j > 0 && i > 0 {
        let current = dp[i * cols + j];
        let ref_c = ref_seq[i - 1];
        let q_c = query_seq[j - 1];
        let cost = if ref_c == q_c {
            match_score
        } else {
            mismatch_penalty
        };

        let diag_score = dp[(i - 1) * cols + j - 1] + cost;
        let _up_score = dp[(i - 1) * cols + j] + gap_penalty;
        let left_score = dp[i * cols + j - 1] + gap_penalty;

        if current == diag_score {
            if ref_c == q_c {
                align_ops.push("M");
            } else {
                align_ops.push("M"); // M covers mismatches in SAM
                actual_mismatches += 1;
            }
            i -= 1;
            j -= 1;
        } else if current == left_score {
            align_ops.push("I");
            actual_mismatches += 1;
            j -= 1;
        } else {
            align_ops.push("D");
            actual_mismatches += 1;
            i -= 1;
        }
    }
    while j > 0 {
        align_ops.push("S"); // Soft clip unaligned query part at start
        j -= 1;
    }

    let start_ref = i;
    align_ops.reverse();

    // Compress CIGAR
    let mut cigar = String::new();
    if align_ops.is_empty() {
        return None;
    }

    let mut current_op = align_ops[0];
    let mut count = 1;
    for k in 1..align_ops.len() {
        if align_ops[k] == current_op {
            count += 1;
        } else {
            cigar.push_str(&format!("{}{}", count, current_op));
            current_op = align_ops[k];
            count = 1;
        }
    }
    cigar.push_str(&format!("{}{}", count, current_op));

    // Add right soft-clip
    if right_clip > 0 {
        cigar.push_str(&format!("{}S", right_clip));
    }

    Some((start_ref, end_ref, cigar, actual_mismatches))
}
