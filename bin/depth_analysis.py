import pysam
import csv
import os
import glob
import sys

# Get meta.id passed as the first argument
meta_id = sys.argv[1]

# Use glob to find all BAM files in the current directory
bam_files = glob.glob("*.bam")

output_file = f"{meta_id}_merged_depth_analysis.csv"

def analyze_bam(bam, meta_id, depth_data):
    # Open BAM file
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        # Get the single reference (assuming only one reference in BAM)
        reference = bam_file.references[0]
        print(f"Analyzing BAM file '{bam}' for reference: {reference}")

        # Iterate through positions in the reference
        for pileupcolumn in bam_file.pileup(contig=reference):
            position = pileupcolumn.pos + 1  # 1-based position
            base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}

            # Count base occurrences at each position
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    if pileupread.query_position is not None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counts[base] = base_counts.get(base, 0) + 1
                else:
                    # Count reads that don't align to a specific base as 'N'
                    base_counts['N'] += 1

            # Use the sum of base-specific depths as the total depth
            total_depth = sum(base_counts.values())

            # Append data for CSV, adding the meta_id and reference name as columns
            depth_data.append([meta_id, os.path.basename(bam), reference, position, total_depth] + list(base_counts.values()))

def main(meta_id, bam_files, output_file):
    depth_data = []

    # Loop through all BAM files found with glob
    for bam in bam_files:
        if bam.endswith(".bam"):
            analyze_bam(bam, meta_id, depth_data)
        else:
            print(f"Skipping non-BAM file: {bam}")

    # Write all results to a single CSV
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Meta ID', 'BAM File', 'Reference', 'Position', 'Total Depth', 'A Depth', 'T Depth', 'C Depth', 'G Depth', 'N Depth'])
        writer.writerows(depth_data)

    print(f"All BAM analyses completed and merged into {output_file}")

# Run the function
if __name__ == "__main__":
    main(meta_id, bam_files, output_file)
