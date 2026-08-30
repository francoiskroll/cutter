#!/usr/bin/env python3
"""
trim_softclips.py
=================
For each alignment in a BAM file, trims soft-clipped bases while retaining
the innermost INDEX_LEN bases of each soft clip — i.e. the in-line index.

Soft-clip structure (your amplicon sequencing setup):
  Left soft-clip:  [adapter ...] [left_index INDEX_LEN nt]  <- keep last INDEX_LEN
  Right soft-clip: [right_index INDEX_LEN nt] [adapter ...]  <- keep first INDEX_LEN

After trimming:
  - SEQ and QUAL are shortened accordingly
  - CIGAR is updated: trimmed soft-clip bases become hard-clips (H),
    retained index bases remain soft-clips (S)
  - POS (reference start) is adjusted if left-side bases are removed
  - Unmapped reads are passed through unchanged

Usage
-----
python trim_softclips.py input.bam output.bam [--index-len 8] [--threads 4]

Requirements
------------
pip install pysam
"""

import argparse
import sys
import pysam


def trim_softclips(read: pysam.AlignedSegment, index_len: int) -> pysam.AlignedSegment:
    """
    Trim a single read's soft clips, keeping index_len bases on each side.
    Returns the (modified) read. Unmapped reads are returned unchanged.
    """
    if read.is_unmapped:
        return read

    cigar = read.cigartuples
    if cigar is None:
        return read

    # pysam CIGAR operation codes
    CMATCH  = 0   # M
    CINS    = 1   # I
    CDEL    = 2   # D
    CREF    = 3   # N
    CSOFT   = 4   # S
    CHARD   = 5   # H
    CPAD    = 6   # P
    CEQUAL  = 7   # =
    CDIFF   = 8   # X

    seq  = read.query_sequence       # full SEQ including existing soft clips
    qual = read.query_qualities      # full QUAL including existing soft clips

    # ── Measure existing soft clips ───────────────────────────────────────────

    left_soft  = 0
    right_soft = 0

    # Left: first op may be H then S (hard clips always outside soft clips)
    idx = 0
    if cigar[idx][0] == CHARD:
        idx += 1
    if idx < len(cigar) and cigar[idx][0] == CSOFT:
        left_soft = cigar[idx][1]

    # Right: last op may be H then S
    idx = len(cigar) - 1
    if cigar[idx][0] == CHARD:
        idx -= 1
    if idx >= 0 and cigar[idx][0] == CSOFT:
        right_soft = cigar[idx][1]

    # Nothing to trim if soft clips are already <= index_len
    left_trim  = max(0, left_soft  - index_len)
    right_trim = max(0, right_soft - index_len)

    if left_trim == 0 and right_trim == 0:
        return read

    # ── Trim SEQ and QUAL ─────────────────────────────────────────────────────

    seq_len   = len(seq)
    new_start = left_trim
    new_end   = seq_len - right_trim

    new_seq  = seq[new_start:new_end]
    new_qual = qual[new_start:new_end] if qual is not None else None

    # ── Rebuild CIGAR ─────────────────────────────────────────────────────────
    # Strategy: rebuild op by op, adjusting the soft/hard clip ops at each end.
    #
    # Original left end (example):  H5 S12 M150 ...
    # After trimming 4 of the 12 soft:  H9 S8 M150 ...
    # (trimmed soft bases become additional hard clip)

    new_cigar = list(cigar)

    # ── Left end ──────────────────────────────────────────────────────────────
    if left_trim > 0:
        i = 0
        # Accumulate existing hard clip on left
        existing_hard_left = 0
        if new_cigar[i][0] == CHARD:
            existing_hard_left = new_cigar[i][1]
            new_cigar.pop(i)

        # Now new_cigar[0] must be CSOFT
        assert new_cigar[0][0] == CSOFT, "Expected SOFT after HARD on left"
        remaining_soft = new_cigar[0][1] - left_trim
        new_cigar.pop(0)

        # Re-insert updated ops
        ops_to_insert = []
        new_hard = existing_hard_left + left_trim
        if new_hard > 0:
            ops_to_insert.append((CHARD, new_hard))
        if remaining_soft > 0:
            ops_to_insert.append((CSOFT, remaining_soft))
        new_cigar = ops_to_insert + new_cigar

    # ── Right end ─────────────────────────────────────────────────────────────
    if right_trim > 0:
        i = len(new_cigar) - 1
        existing_hard_right = 0
        if new_cigar[i][0] == CHARD:
            existing_hard_right = new_cigar[i][1]
            new_cigar.pop(i)
            i -= 1

        assert new_cigar[i][0] == CSOFT, "Expected SOFT before HARD on right"
        remaining_soft = new_cigar[i][1] - right_trim
        new_cigar.pop(i)

        ops_to_append = []
        if remaining_soft > 0:
            ops_to_append.append((CSOFT, remaining_soft))
        new_hard = existing_hard_right + right_trim
        if new_hard > 0:
            ops_to_append.append((CHARD, new_hard))
        new_cigar = new_cigar + ops_to_append

    # ── Adjust POS if left bases were removed ─────────────────────────────────
    # Soft-clipped bases don't consume reference, so trimming them from the
    # left doesn't change reference_start. Only hard clips are relevant here
    # and they also don't consume reference. So POS is unchanged.

    # ── Write back ────────────────────────────────────────────────────────────
    read.query_sequence   = new_seq
    read.query_qualities  = new_qual
    read.cigartuples      = new_cigar

    return read


def main():
    parser = argparse.ArgumentParser(
        description="Trim soft clips in a BAM file, keeping the innermost "
                    "index_len bases (the in-line index) as soft clip.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("input_bam",  help="Input BAM file (must be sorted+indexed, or unsorted)")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument("--index-len", type=int, default=8,
                        help="Number of index bases to retain at each end (default: 8)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of compression threads for output (default: 1)")
    args = parser.parse_args()

    print(f"\nSettings:")
    print(f"  Input BAM    : {args.input_bam}")
    print(f"  Output BAM   : {args.output_bam}")
    print(f"  Index length : {args.index_len} nt retained per soft-clip end")
    print(f"  Threads      : {args.threads}\n")

    n_total      = 0
    n_trimmed    = 0
    n_unmapped   = 0
    n_no_clip    = 0

    with pysam.AlignmentFile(args.input_bam, "rb") as bam_in, \
         pysam.AlignmentFile(args.output_bam, "wb",
                             header=bam_in.header,
                             threads=args.threads) as bam_out:

        for read in bam_in:
            n_total += 1

            if read.is_unmapped:
                n_unmapped += 1
                bam_out.write(read)
                continue

            cigar = read.cigartuples
            if cigar is None:
                n_no_clip += 1
                bam_out.write(read)
                continue

            # Check if any trimming is needed
            left_soft  = cigar[0][1]  if cigar[0][0]  == 4 else 0
            right_soft = cigar[-1][1] if cigar[-1][0] == 4 else 0
            # Also account for existing hard clips on outside
            if cigar[0][0] == 5 and len(cigar) > 1:
                left_soft  = cigar[1][1] if cigar[1][0]  == 4 else 0
            if cigar[-1][0] == 5 and len(cigar) > 1:
                right_soft = cigar[-2][1] if cigar[-2][0] == 4 else 0

            needs_trim = (left_soft > args.index_len or right_soft > args.index_len)

            if not needs_trim:
                n_no_clip += 1
                bam_out.write(read)
                continue

            read = trim_softclips(read, args.index_len)
            n_trimmed += 1
            bam_out.write(read)

            if n_total % 100_000 == 0:
                print(f"  Processed {n_total:,} reads...", flush=True)

    print(f"\n{'='*50}")
    print(f"  Done — {n_total:,} reads processed")
    print(f"  Trimmed          : {n_trimmed:,}")
    print(f"  No trim needed   : {n_no_clip:,}")
    print(f"  Unmapped (passed): {n_unmapped:,}")
    print(f"{'='*50}\n")
    print(f"  Output: {args.output_bam}")
    print(f"  Don't forget to index: samtools index {args.output_bam}\n")


if __name__ == "__main__":
    main()
