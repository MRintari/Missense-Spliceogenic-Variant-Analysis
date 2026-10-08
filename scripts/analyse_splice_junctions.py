#!/usr/bin/env python3
import argparse
import csv
import os
import sys
import pysam


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Quantify allele-specific intron retention at a specific junction boundary."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to single BAM file (.bam) OR text file listing BAM paths (.txt/.list)."
    )
    parser.add_argument(
        "-v", "--variant", required=True,
        help="Variant in 'chr:pos:ref:alt' format (e.g., 'chr1:235458763:G:A')."
    )
    parser.add_argument(
        "-b", "--boundary", type=int, required=True,
        help="1-based coordinate of the exact Exon-Intron boundary (Donor or Acceptor site)."
    )
    parser.add_argument(
        "-g", "--group", choices=["Carrier", "Control"], required=True,
        help="Experimental group for these BAMs: 'Carrier' or 'Control'."
    )
    parser.add_argument(
        "-o", "--output-csv", required=True,
        help="Output CSV file path for read-level metadata."
    )
    parser.add_argument(
        "-w", "--min-overlap", type=int, default=2,
        help="Minimum required bp overlap on BOTH sides of boundary to count as Retained (default: 2)."
    )
    parser.add_argument(
        "-q", "--min-mapq", type=int, default=20,
        help="Minimum mapping quality (MAPQ) to include a read (default: 20)."
    )
    return parser.parse_args()


def parse_variant(var_str):
    parts = var_str.strip().split(":")
    if len(parts) != 4:
        sys.exit(f"Error: Variant '{var_str}' invalid. Must be 'chr:pos:ref:alt'.")
    return parts[0], int(parts[1]), parts[2].upper(), parts[3].upper()


def extract_path_metadata(bam_path):
    abs_path = os.path.abspath(bam_path)
    parts = abs_path.split(os.sep)

    seq_date = "Unknown_Date"
    rna_id = "Unknown_RNAID"
    sample_id = os.path.basename(bam_path).replace(".bam", "")

    if len(parts) >= 4:
        sample_id = parts[-2] if parts[-2] else parts[-1]
        rna_id = parts[-3]
        possible_date = parts[-4]
        if len(possible_date) == 10 and possible_date.count("-") == 2:
            seq_date = possible_date

    return seq_date, rna_id, sample_id


def get_bam_files(input_path):
    if not os.path.exists(input_path):
        sys.exit(f"Error: Input file '{input_path}' does not exist.")
    if input_path.endswith(".bam"):
        return [input_path]

    bam_list = []
    with open(input_path, "r") as f:
        for line in f:
            path = line.strip()
            if path and not path.startswith("#") and os.path.exists(path):
                bam_list.append(path)
    return bam_list


def process_bam(bam_file, chrom, vus_pos, ref_allele, alt_allele, boundary, min_overlap, min_mapq, group_label, csv_writer):
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error reading '{bam_file}': {e}", file=sys.stderr)
        return

    seq_date, rna_id, sample_id = extract_path_metadata(bam_file)

    # Convert 1-based CLI inputs to 0-based coordinates for internal PySAM logic
    vus_pos_0 = vus_pos - 1
    boundary_0 = boundary - 1

    # Broaden window (1000 bp) so reads starting/ending far away are fetched properly
    fetch_start = max(0, min(vus_pos_0, boundary_0) - 1000)
    fetch_end = max(vus_pos_0, boundary_0) + 1000

    # Dictionary to aggregate information at the FRAGMENT level (by query_name)
    fragments = {}

    for read in bam.fetch(chrom, fetch_start, fetch_end):
        # Filter out low quality, secondary, supplementary, and duplicate alignments
        if (read.is_secondary or
            read.is_supplementary or
            read.is_duplicate or
            read.is_qcfail or
            read.mapping_quality < min_mapq):
            continue

        read_id = read.query_name

        if read_id not in fragments:
            fragments[read_id] = {
                'allele': None,
                'has_spliced': False,
                'has_retained': False,
                'mapq': read.mapping_quality
            }

        # # 1. Check VUS position allele on this mate
        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        for r_idx, g_idx in aligned_pairs:
            if g_idx == vus_pos_0:
                base = read.query_sequence[r_idx].upper()
                if base == alt_allele:
                    fragments[read_id]['allele'] = 'ALT'
                elif base == ref_allele:
                    fragments[read_id]['allele'] = 'REF'
                break

        # # 2. Check Boundary Splicing / Retention on this mate
        left_check = boundary_0 - min_overlap
        right_check = boundary_0 + min_overlap

        # Continuous read check (Retained)
        for start, end in read.get_blocks():
            if start <= left_check and end >= right_check:
                fragments[read_id]['has_retained'] = True
                break

        # CIGAR gap check (Spliced)
        cigar_tuples = read.cigartuples or []
        curr_pos = read.reference_start
        for op, length in cigar_tuples:
            if op in (0, 7, 8):  # Match/Mismatch
                curr_pos += length
            elif op == 3:  # 'N' (Splice junction gap)
                splice_start = curr_pos
                splice_end = curr_pos + length
                if abs(splice_start - boundary_0) <= 2 or abs(splice_end - boundary_0) <= 2:
                    fragments[read_id]['has_spliced'] = True
                curr_pos += length
            elif op == 2:  # Deletion
                curr_pos += length

    bam.close()

    # Final Classification Phase (Option 1: Spliced Priority Rule)
    for read_id, info in fragments.items():
        allele_key = info['allele']

        # Skip fragments that do not cover the VUS position
        if not allele_key:
            continue

        if info['has_spliced']:
            splicing_status = "Spliced"
            retained_binary = 0
        elif info['has_retained']:
            splicing_status = "Retained"
            retained_binary = 1
        else:
            continue

        csv_writer.writerow([
            read_id,
            sample_id,
            rna_id,
            seq_date,
            group_label,
            allele_key,
            splicing_status,
            retained_binary,
            info['mapq']
        ])


def main():
    args = parse_arguments()
    chrom, vus_pos, ref_allele, alt_allele = parse_variant(args.variant)
    bam_files = get_bam_files(args.input)

    file_exists = os.path.exists(args.output_csv)

    with open(args.output_csv, "a", newline="") as out_f:
        writer = csv.writer(out_f)

        if not file_exists or os.path.getsize(args.output_csv) == 0:
            writer.writerow([
                "Fragment_ID", "Sample_ID", "RNA_ID", "Seq_Date",
                "Group", "Allele", "Splicing_Status", "Retained_Binary", "MAPQ"
            ])

        print(f"Analyzing {len(bam_files)} BAM file(s)...")
        print(f"Variant: {chrom}:{vus_pos} ({ref_allele} > {alt_allele})")
        print(f"Boundary: {chrom}:{args.boundary}")

        for b_file in bam_files:
            process_bam(
                b_file, chrom, vus_pos, ref_allele, alt_allele,
                args.boundary, args.min_overlap, args.min_mapq,
                args.group, writer
            )

    print(f"\nDone! Extracted reads written to '{args.output_csv}'")


if __name__ == "__main__":
    main()
