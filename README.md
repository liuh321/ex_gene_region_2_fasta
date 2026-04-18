Extract region sequences from a FASTA file.

Usage:
    python extract_region.py <fa_file> <region_file> [-o output.fa]

region.txt format:
    Chr1    100    300
    Chr2    500    900
    ...
    (1-based, inclusive coordinates)

Supports gzipped FASTA files (.gz).
