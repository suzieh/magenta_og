# usage : python3 tabulate_mags.py [input_dir] [output_file] [keep_n_mags]

import os
import sys
import csv
from collections import defaultdict

# Directory containing the files
directory = sys.argv[1]

# Specify the number of top MAGIDs to keep
top_n_magids = int(sys.argv[3])  # Adjust this to the desired number of top MAGIDs

# Dictionary to store counts of MAGIDs in each file
## lambda anonymous function allows for a new defaultdict(int)
##     each time a new file is called.
file_magid_counts = defaultdict(lambda: defaultdict(int))

# Dictionary to store total counts of each MAGID across all files
total_magid_counts = defaultdict(int)

# Set to store all unique MAGIDs
all_magids = set()

# Iterate over each SAM file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.txt'):  # Check for .txt files
        filepath = os.path.join(directory, filename)
        
        # Open and read the SAM file
        with open(filepath, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            
            # Iterate over each line in the SAM file
            for row in reader:
                read_id = row[0]  # Assuming Read ID in first column
                magid = row[1]    # Assuming MAGID in second column
                    
                # Count occurrences of each MAGID for each file
                file_magid_counts[filename][magid] += 1
                total_magid_counts[magid] += 1  # Track total counts across all files

# Sort MAGIDs by their total counts and keep only the top N MAGIDs
top_magids = sorted(total_magid_counts, key=total_magid_counts.get, reverse=True)[:top_n_magids]

# Write the results to a tab-delimited output file
with open(sys.argv[2], 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t')
    
    # Write the header row with file names and the top N MAGIDs as columns
    header = ['FileName'] + top_magids
    writer.writerow(header)
    
    # Write each file's counts of the top N MAGIDs
    for filename, magid_count_dict in file_magid_counts.items():
        row = [filename]
        for magid in top_magids:
            row.append(magid_count_dict.get(magid, 0))  # Append count, or 0 if MAGID not present
        writer.writerow(row)

print(f"Output file {sys.argv[2]} created successfully, keeping the top {top_n_magids} MAGIDs.")
