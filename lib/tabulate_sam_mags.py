# usage : python3 tabulate_mags.py [input_dir] [output_file]

import os
import csv
from collections import defaultdict

# Directory containing the SAM files
directory = sys.argv[1]

# Dictionary to store counts of MAGIDs in each file
file_magid_counts = defaultdict(lambda: defaultdict(int))

# Set to store all unique MAGIDs
all_magids = set()

# Iterate over each SAM file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.sam'):  # Check for .sam files
        filepath = os.path.join(directory, filename)
        
        # Open and read the SAM file
        with open(filepath, 'r') as sam_file:
            reader = csv.reader(sam_file, delimiter='\t')
            
            # Iterate over each line in the SAM file
            for row in reader:
                if not row[0].startswith('@'):  # Skip header lines (lines that start with '@')
                    read_id = row[0]  # First column: ReadID
                    magid = row[2]    # Third column: MAGID
                    
                    # Count occurrences of each MAGID for each file
                    file_magid_counts[filename][magid] += 1
                    all_magids.add(magid)  # Store unique MAGIDs

# Write the results to a tab-delimited output file
with open('output_file.tsv', 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t')
    
    # Write the header row with file names and sorted MAGIDs as columns
    header = ['FileName'] + sorted(all_magids)
    writer.writerow(header)
    
    # Write each file's counts of MAGIDs
    for filename, magid_count_dict in file_magid_counts.items():
        row = [filename]
        for magid in sorted(all_magids):
            row.append(magid_count_dict.get(magid, 0))  # Append count, or 0 if MAGID not present
        writer.writerow(row)

print("Output file 'output_file.tsv' created successfully.")