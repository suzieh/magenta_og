# usr/bin/python3
# usage: parse_bt2_output.py -i [alignment_output_path] -m [mapping_file] -o [output_table]

import os
import argparse


# Helper function: Read the donor / pre read ids
def get_unique_ids(filename):
    unique_values = set()
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split()  # Assumes whitespace delimiter
            if columns:  # Ensure the line is not empty
                unique_values.add(columns[0])
    return unique_values

# Helper function: count FASTQ reads
def count_fastq_reads(fastq_file):
    with open(fastq_file, 'r') as f:
        return sum(1 for _ in f) // 4

# Helper function: get ratio of alignments to reads (self-alignments of donors and pre samples)
def read_match_ratio(fastq_file, align_file)
    with open(fastq_file, 'r') as f1:
        total = sum(1 for _ in f1) // 4
    with open(align_file, 'r') as f2:
        match = sum(1 for _ in f2)
    return match / total

# Helper function: get the donor and pre-treatment cross-over
def crossover(fastq_file, db_id)
    with open({fastq_file.replace(".fastq","")}_{db_id}_aligned.txt, 'r') as f2:
        match = sum(1 for _ in f2)
    return match

# Main function
def main():
    # Argument parser set-up
    parser = argparse.ArgumentParser(
        description="Parse the Bowtie2 Alignment output for engraftment estimates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input alignment files.")
    parser.add_argument("-m", "--mapping", type=str, required=True, help="Mapping file (table in MAGEnTa format).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file path and name.")
    args = parser.parse_args()

    # Parse arguments
    mapping_file = args.mapping
    output_file = args.output

    # Loop over elements of mapping file table
    with open(mapping_file, 'r') as table:
        header = table.readline()  # skip header
        if not os.path.exists(output_file):
            with open(output_file, 'w') as f:
                f.write('Simulation,Total_Donor,Total_Pre,Percent_Donor,Percent_Pre,Alignment_Adjustment,Uniquely_Donor,Uniquely_Pre,Ambiguous,Amb_Donor,Amb_Pre,Unmapped\n')
        for line in table:
            cols = line.strip().split('\t')
            if len(cols) != 8:
                continue  # skip malformed lines
            simulation_pre = cols[0]
            post_fastq = cols[1]
            donor_id = cols[2]
            donor_fwd = cols[3]
            donor_rev = cols[4]
            pre_id = cols[5]
            pre_fwd = cols[6]
            pre_rev = cols[7]

            # Paths to donor and pre alignment files (assume .txt files already created from prior alignments)
            donor_file = f'{post_fastq.replace(".fastq","")}_{donor_id}_aligned.txt'  # generated alignment file
            pre_file = f'{post_fastq.replace(".fastq","")}__{pre_id}_aligned.txt'     # generated alignment file

            if not os.path.exists(donor_file) or not os.path.exists(pre_file) or not os.path.exists(post_fastq):
                print(f"Missing alignment files for {simulation_pre}, skipping.")
                continue

            # Unique read IDs
            d_reads = get_unique_ids(donor_file)
            p_reads = get_unique_ids(pre_file)

            # Total reads in post-treatment FASTQ
            total_n = count_fastq_reads(post_fastq)

            # Set operations
            intersect = d_reads.intersection(p_reads)
            d_uniq = len(d_reads) - len(intersect)
            p_uniq = len(p_reads) - len(intersect)
            ambig = len(intersect)
            unknown = total_n - len(intersect) - d_uniq - p_uniq

            # Assume base cases give an estimate of unmatched cases (alignments of donor and pre samples with themselves)
            base_donor = read_match_ratio(donor_fwd, f'{donor_fwd.replace(".fastq","")}_{donor_id}_aligned.txt')
            base_pre = read_match_ratio(pre_fwd, f'{pre_fwd.replace(".fastq","")}_{pre_id}_aligned.txt')
            mean_ratio = (base_donor + base_pre) / 2

            # Calculate probabilities
            base_p_d = crossover(pre_fwd, donor_id)
            base_d_p = crossover(donor_fwd, pre_id)
            prob_d_p = base_p_d / (base_pre * total_n)
            prob_p_d = base_d_p / (base_donor * total_n)
            prob_d = d_uniq / (d_uniq + p_uniq) if (d_uniq + p_uniq) > 0 else 0.5
            prob_p = p_uniq / (d_uniq + p_uniq) if (d_uniq + p_uniq) > 0 else 0.5

            amb_d = (prob_p_d * prob_d) / ((prob_p_d * prob_d) + (prob_d_p * prob_p)) if ((prob_p_d * prob_d) + (prob_d_p * prob_p)) != 0 else 0.5
            amb_p = 1 - amb_d

            d_tot = d_uniq + round(amb_d * ambig)
            p_tot = p_uniq + round(amb_p * ambig)

            corr = mean_ratio if (total_n * mean_ratio) > (d_tot + p_tot) else 1 # expected alignment rate
            corr_d = d_tot / (total_n * corr)
            corr_p = p_tot / (total_n * corr)

            # Write to output file
            with open(output_file, 'a') as f:
                f.write(f'{simulation_pre},{d_tot},{p_tot},{corr_d},{corr_p},{d_uniq},{p_uniq},{ambig},{round(amb_d * ambig)},{round(amb_p * ambig)},{unknown}\n')



# Run the main function
if __name__ == "__main__":
    main()

