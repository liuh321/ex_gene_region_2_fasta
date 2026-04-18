#!/usr/bin/env python3
"""
Extract region sequences from a FASTA file.

Usage:
    python extract_region.py <fa_file> <region_file> [-o output.fa]

region.txt format:
    Chr1    100    300
    Chr2    500    900
    ...
    (1-based, inclusive coordinates)

Supports gzipped FASTA files (.gz).
"""

import argparse
import gzip
import sys
from pathlib import Path


def parse_fasta(fa_path: str) -> dict[str, str]:
    """Read a FASTA file (plain or gzipped) into a dict {id: seq}."""
    sequences: dict[str, str] = {}
    open_fn = gzip.open if fa_path.endswith(".gz") else open
    mode = "rt" if fa_path.endswith(".gz") else "r"

    with open_fn(fa_path, mode) as fh:
        current_id: str | None = None
        buf: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(buf)
                current_id = line[1:].split()[0]   # first token after '>'
                buf = []
            else:
                buf.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(buf)

    return sequences


def parse_regions(region_path: str) -> list[tuple[str, int, int]]:
    """Read region file and return list of (chrom, start, end) — 1-based inclusive."""
    regions: list[tuple[str, int, int]] = []
    with open(region_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                print(f"[WARN] line {lineno} skipped (fewer than 3 columns): {line!r}",
                      file=sys.stderr)
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            regions.append((chrom, start, end))
    return regions


def extract_regions(sequences: dict[str, str],
                    regions: list[tuple[str, int, int]],
                    out_fh) -> None:
    """Extract each region and write to out_fh in FASTA format."""
    for chrom, start, end in regions:
        if chrom not in sequences:
            print(f"[WARN] '{chrom}' not found in FASTA — skipping.", file=sys.stderr)
            continue

        seq = sequences[chrom]
        seq_len = len(seq)

        # Convert 1-based inclusive → 0-based half-open
        pos = start - 1
        size = end - start + 1

        if pos < 0 or pos >= seq_len:
            print(f"[WARN] start={start} out of range for '{chrom}' (len={seq_len}) — skipping.",
                  file=sys.stderr)
            continue
        if end > seq_len:
            print(f"[WARN] end={end} exceeds '{chrom}' length ({seq_len}); clipping.",
                  file=sys.stderr)
            size = seq_len - pos

        subseq = seq[pos: pos + size]
        region_id = f"{chrom}_{start}_{end}"
        out_fh.write(f">{region_id}\n{subseq}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract region sequences from a FASTA file (equivalent to deal_fa.pl -format 5)."
    )
    parser.add_argument("fa_file",    help="Input FASTA file (plain or .gz)")
    parser.add_argument("region_file", help="Region file: chrom  start  end  (1-based, inclusive)")
    parser.add_argument("-o", "--output", default="-",
                        help="Output FASTA file (default: stdout)")
    args = parser.parse_args()

    # Validate inputs
    if not Path(args.fa_file).exists():
        sys.exit(f"[ERROR] FASTA file not found: {args.fa_file}")
    if not Path(args.region_file).exists():
        sys.exit(f"[ERROR] Region file not found: {args.region_file}")

    print("[INFO] Loading FASTA...", file=sys.stderr)
    sequences = parse_fasta(args.fa_file)
    print(f"[INFO] Loaded {len(sequences)} sequences.", file=sys.stderr)

    regions = parse_regions(args.region_file)
    print(f"[INFO] Loaded {len(regions)} regions.", file=sys.stderr)

    if args.output == "-":
        extract_regions(sequences, regions, sys.stdout)
    else:
        with open(args.output, "w") as out_fh:
            extract_regions(sequences, regions, out_fh)
        print(f"[INFO] Output written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
