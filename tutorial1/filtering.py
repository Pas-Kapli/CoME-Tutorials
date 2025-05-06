#!/usr/bin/env python3

import argparse
import os

def parse_fasta(file_path):
    """Parse a FASTA file and return taxon names and sequences."""
    headers, sequences = [], []
    with open(file_path, "r") as f:
        temp_seq = ""
        for line in f:
            if line.startswith(">"):
                headers.append(line.strip())
                if temp_seq:
                    sequences.append(temp_seq)
                    temp_seq = ""
            else:
                temp_seq += line.strip()
        if temp_seq:
            sequences.append(temp_seq)
    return headers, sequences

def write_fasta(headers, sequences, output_file):
    """Write sequences to a FASTA file."""
    with open(output_file, "w") as f:
        for header, seq in zip(headers, sequences):
            f.write(f"{header}\n{seq}\n")

def remove_columns_with_missing_data(sequences, threshold):
    """
    Remove columns where the proportion of missing data ('-', 'N', 'n')
    exceeds the given threshold or contains only lowercase letters.
    """
    if not sequences:
        return []

    seq_length = len(sequences[0])
    seq_count = len(sequences)
    
    valid_columns = [
        i for i in range(seq_length)
        if sum(seq[i] in "-Nn" for seq in sequences) / seq_count <= threshold
        and not all(seq[i].islower() for seq in sequences)
    ]
    
    filtered_seqs = ["".join(seq[i] for i in valid_columns) for seq in sequences]
    
    return filtered_seqs

def check_alignment_length(sequences, min_length):
    """Check if the alignment length is greater than or equal to the given threshold."""
    if not sequences or len(sequences[0]) < min_length:
        return False
    return True

def filter_by_average_differences(headers, sequences, threshold):
    """Remove sequences exceeding average percentage differences."""
    seq_count = len(sequences)
    seq_length = len(sequences[0])
    avg_diffs = []

    for i in range(seq_count):
        diffs = []
        for j in range(seq_count):
            if i != j:
                pair_diff = sum(
                    sequences[i][k] != sequences[j][k]
                    and sequences[i][k] not in "-Nn"
                    and sequences[j][k] not in "-Nn"
                    for k in range(seq_length)
                )
                valid_bases = sum(
                    sequences[i][k] not in "-Nn" and sequences[j][k] not in "-Nn"
                    for k in range(seq_length)
                )
                diffs.append(pair_diff / valid_bases if valid_bases > 0 else 0)
        avg_diffs.append(sum(diffs) / len(diffs) if diffs else 0)
    
    filtered_headers = [
        header for header, diff in zip(headers, avg_diffs) if diff < threshold
    ]
    filtered_seqs = [
        seq for seq, diff in zip(sequences, avg_diffs) if diff < threshold
    ]
    return filtered_headers, filtered_seqs

def main():
    parser = argparse.ArgumentParser(description="Filter alignment sequences.")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output file")
    parser.add_argument("--format", choices=["fasta", "phylip"], default="fasta", help="Output format")
    parser.add_argument("--diff_threshold", type=float, required=True, help="Threshold for average differences")
    parser.add_argument("--missing_threshold", type=float, required=True, help="Threshold for missing data in sequences")
    parser.add_argument("--alignment_missing_threshold", type=float, required=True, help="Threshold for alignment removal")
    parser.add_argument("--column_missing_threshold", type=float, required=True, help="Threshold for missing data in columns")
    parser.add_argument("--min_alignment_length", type=int, required=True, help="Minimum length for output alignment")
    args = parser.parse_args()
    
    headers, sequences = parse_fasta(args.input_file)
    
    if not sequences:
        print("No sequences found in input file. Exiting.")
        return
    
    sequences = remove_columns_with_missing_data(sequences, args.column_missing_threshold)
    
    if not check_alignment_length(sequences, args.min_alignment_length):
        print("Filtered alignment is too short. No output generated.")
        return
    
    headers, sequences = filter_by_average_differences(headers, sequences, args.diff_threshold)
    
    if args.format == "fasta":
        write_fasta(headers, sequences, args.output_file)
    elif args.format == "phylip":
        max_name_length = max(len(header[1:]) for header in headers)
        with open(args.output_file, "w") as f:
            f.write(f"{len(headers)} {len(sequences[0])}\n")
            for header, seq in zip(headers, sequences):
                f.write(f"{header[1:]:<{max_name_length + 2}}{seq}\n")
    
    print(f"Filtered alignment saved to {args.output_file}.")

if __name__ == "__main__":
    main()
