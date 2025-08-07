import sys
import os

# Define a global output file where all results will be stored
final_output_file = "seq_matching_results.txt"

# Ensure output file exists and write headers only once
if not os.path.exists(final_output_file):
    with open(final_output_file, "w") as outfile:
        outfile.write("Result_File\tTotal_Percentage_Matched\tTotal_Abs_Percentage_Matched\n")  # Use tab separation

# Process each input file given in command-line arguments
for result_file in sys.argv[1:]:
    total_residue = 0
    total_residue_matched = 0
    total_residue_matched_abs = 0

    with open(result_file, 'r') as infile:
        for i, lines in enumerate(infile):
            # Skip header line (assuming the first line is the header)
            if i == 0 and "pdbname" in lines.lower():
                print(f"[INFO] Skipping header in {result_file}")
                continue
            
            # Skip empty lines
            if not lines.strip():
                continue

            line = lines.strip().split('|')

            # Ensure the first part is formatted correctly
            split_data = line[0].strip().split()

            if len(split_data) < 3:
                print(f"[WARNING] Skipping malformed line in {result_file}: {line[0]}")
                continue  # Skip this line
            
            emdb, resolution, seid = split_data[:3]  # Take only the first 3 elements

            try:
                se_length = int(seid.split('_')[-1])
            except ValueError:
                print(f"[WARNING] Could not parse SE length from: {seid} in {result_file}")
                continue  # Skip this line
            
            # Prevent IndexError
            if len(line) < 2:
                print(f"[WARNING] Skipping malformed line (missing expected values) in {result_file}: {line}")
                continue

            if 'nan' in line[1]:
                print(f"[INFO] Skipping line with 'nan' values in {result_file}: {line}")
                continue
            else:
                total_residue += se_length
                rank = line[1].strip().split()[0]  # Handle multiple spaces

                try:
                    if int(rank) <= se_length / 3:
                        total_residue_matched += se_length

                    if int(rank) == 0:
                        total_residue_matched_abs += se_length
                except ValueError:
                    print(f"[WARNING] Invalid rank value in {result_file}: {rank}")
                    continue

    # Prevent division by zero error
    if total_residue == 0:
        total_percentage_matched = 0.0
        total_abs_percentage_matched = 0.0
    else:
        total_percentage_matched = round(100 * (total_residue_matched / total_residue), 3)
        total_abs_percentage_matched = round(100 * (total_residue_matched_abs / total_residue), 3)

    # Append new result to the common output file
    with open(final_output_file, "a") as outfile:
        outfile.write(f"{result_file}\t{total_percentage_matched}\t{total_abs_percentage_matched}\n")  # Use tab separation

#print(f"✅ Results saved to {final_output_file}")

