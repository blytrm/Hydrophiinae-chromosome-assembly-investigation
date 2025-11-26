#!/usr/bin/env bash

# Define the target chromosome file
TARGET_CHR="Hmajchr2_0-10mb.fasta"

# Check if target chromosome file exists
if [ ! -f "$TARGET_CHR" ]; then
    echo "Error: Target chromosome file '$TARGET_CHR' not found!"
    exit 1
fi

# Loop through all available reference genome files
for ref_file in *_z.fasta; do
    # Check if the file exists (in case no files match the pattern)
    if [ ! -f "$ref_file" ]; then
        echo "No reference files found matching pattern '*_z.fasta'"
        exit 1
    fi
    
    # Extract the species name from the filename (remove _z.fasta suffix)
    species_name=$(basename "$ref_file" _z.fasta)
    
    # Create output filename
    output_file="maj2_${species_name}Z.paf"
    
    echo "Processing: $ref_file -> $output_file"
    
    # Run minimap2 alignment
    minimap2 -x asm20 -t 3 -c --MD "$ref_file" "$TARGET_CHR" > "$output_file"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully created: $output_file"
    else
        echo "Error processing: $ref_file"
    fi
    
    echo "----------------------------------------"
done

echo "All alignments completed!"
