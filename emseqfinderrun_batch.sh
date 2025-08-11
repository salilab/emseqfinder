#!/bin/bash

module load imp

# Resolution used for database generation
resolution=4

# Output file for overall results
final_output_file="batch_matching_results.txt"

# Add header if file doesn’t exist
if [[ ! -f "$final_output_file" ]]; then
    echo "Result_File Total_Percentage_Matched Total_Abs_Percentage_Matched" > "$final_output_file"
fi

# Loop through all PDB files
for pdbfile in pdb_files/*.pdb; do

    filename_with_ext=$(basename "$pdbfile")
    base="${filename_with_ext%.*}"

    mapfile="cryoem_maps/${base}.map"
    fastafile="fasta_files/${base}.fasta"

    # Check file existence
    if [[ ! -f "$mapfile" ]]; then
        echo "[WARNING] Missing map file for '$base'. Skipping..."
        continue
    fi
    if [[ ! -f "$fastafile" ]]; then
        echo "[WARNING] Missing FASTA file for '$base'. Skipping..."
        continue
    fi

    echo "==============================================="
    echo "[INFO] Processing $base"

    # Compute dynamic threshold
    threshold=$(python3 -m IMP.emseqfinder.compute_dynamic_threshold "$mapfile")
    if [[ -z "$threshold" ]]; then
        echo "[WARNING] Failed to compute threshold for $base. Skipping..."
        continue
    fi
    echo "[INFO] Using threshold: $threshold"

    # Clean and recreate project directory
    rm -rf "$base"
    mkdir -p "$base/0system"
    mkdir -p "$base/1structure_elements"

    # STRIDE assignment
    stride_file="$base/${base}.stride"
    rm -f "$stride_file"
    stride "$pdbfile" -rA -f"$stride_file"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] STRIDE failed for $base. Skipping..."
        continue
    fi

    # Create unique fragment directory
    frag_dir="fragments_${base}"
    rm -rf "$frag_dir"
    mkdir -p "$frag_dir"
    python3 -m IMP.emseqfinder.fragdb.get_fraglib_from_native \
        "$pdbfile" "$stride_file" "$frag_dir"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Fragment generation failed. Skipping $base."
        continue
    fi

    # Copy system files
    cp "$pdbfile" "$base/0system/native.pdb"
    cp "$mapfile" "$base/0system/emdb.map"
    cp "$fastafile" "$base/0system/"
    cp "$frag_dir"/*.pdb "$base/1structure_elements/"

    # Normalize map
    python3 -m IMP.emseqfinder.mldb.normalize_map_for_parts_fitting \
        "$base" --thresh "$threshold"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Normalization failed. Skipping $base."
        continue
    fi

    # Remove old database files if exist
    rm -f "${base}_ML_side.dat" "${base}_ML_side.pkl" "${base}_ML_side_ML_prob.dat" "${base}_ML_output.txt"

    # Create database
    python3 -m IMP.emseqfinder.mldb.get_database_for_one_emdb_using_parts \
        "$base" "${base}_ML_side.dat" "$resolution"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] ML DB generation failed. Skipping $base."
        continue
    fi

    # Convert to pickle
    python3 -m IMP.emseqfinder.convert_MLDB_topkl \
        "${base}_ML_side.dat" "${base}_ML_side"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] PKL conversion failed. Skipping $base."
        continue
    fi

    # Prediction
    python3 -m IMP.emseqfinder.final_ML_predict \
        "${base}_ML_side.pkl" 10000
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Prediction failed. Skipping $base."
        continue
    fi

    # Evaluate prediction
    python3 -m IMP.emseqfinder.evaluate_output_database \
        "${base}_ML_side_ML_prob.dat" "${base}_ML_output.txt"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Evaluation failed. Skipping $base."
        continue
    fi

    # Calculate and append sequence match
    python3 -m IMP.emseqfinder.calculate_seq_match_batch \
        "${base}_ML_output.txt" >> "$final_output_file"
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Match calculation failed. Skipping $base."
        continue
    fi

    # Optional cleanup
    rm -rf "$frag_dir"

    echo "[INFO] Finished processing $base ✅"
    echo "==============================================="
    echo

done

echo "[DONE] Batch processing complete. Results saved to: $final_output_file"
