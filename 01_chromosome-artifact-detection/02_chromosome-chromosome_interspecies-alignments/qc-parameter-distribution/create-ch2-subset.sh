#!/bin/bash

# Extract first 10Mb of chromosome 2 from H.major assembly
# Script to subset chromosome 2 (0-10Mb region)

# Set file paths
ASSEMBLY_FILE= "/Users/billytrim/Desktop/Sanders_Lab/H.major_C2_Assessment/2_Hmajch2-chz/interm_raw_files/H.maj_seq.fasta"
OUTPUT_FILE="chr2_0-10mb.fasta"

# Extract chromosome 2 and subset to first 10Mb
seqkit grep -p "chr2" "$ASSEMBLY_FILE" | seqkit subseq -r 1:10000000 > "$OUTPUT_FILE"

# Check if extraction was successful
if [ $? -eq 0 ]; then
    echo "Successfully extracted chromosome 2 (0-10Mb) to $OUTPUT_FILE"
    
    # Get some basic stats
    echo ""
    echo "File statistics:"
    seqkit stats "$OUTPUT_FILE"
else
    echo "Error: Failed to extract chromosome 2"
    exit 1
fi 
