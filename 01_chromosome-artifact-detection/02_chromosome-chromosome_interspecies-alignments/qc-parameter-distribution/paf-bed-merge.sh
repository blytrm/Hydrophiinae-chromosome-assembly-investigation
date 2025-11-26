#!/bin/bash

# Define the input/output directory
DATA_DIR="/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/2_Hmajch2-chz/Chromosome-Chromosome_Alignments/qc_paf/Interm-raw-data"

# Change to the data directory
cd "$DATA_DIR"

# Loop through all final PAF files and create merged BED files
for file in maj2_*_final_qc.paf; do
    # Check if the file exists (in case no files match the pattern)
    if [ ! -f "$file" ]; then
        echo "No PAF files found matching pattern 'maj2_*_final_qc.paf'"
        exit 1
    fi
    
    # Extract the species name from the filename
    species=$(basename "$file" _final_qc.paf)
    echo "Processing $species..."
    
    # Create merged BED file using your approach
    tail -n+2 "$file" | cut -f 6-8 | awk -v OFS="\t" '{print $0}' | sort -k2,2n -k3,3n | bedtools merge -d 1000 > "${species}_merged.bed"
    
    echo "Created ${species}_merged.bed"
done

# Verify the files were created
echo "=== Created merged BED files ==="
ls -la *_merged.bed

# Show first few lines of each file
for file in *_merged.bed; do
    echo "=== $file ==="
    head -5 "$file"
    echo ""
done
