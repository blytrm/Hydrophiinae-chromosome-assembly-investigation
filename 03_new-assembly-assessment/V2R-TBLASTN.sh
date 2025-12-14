#!/bin/bash
set -e

QUERY="HMA.7tm.cdhit.fa"
NEW_REF="new-hmaj-assembly.fa"
OLD_REF="old-filtered.fa"

echo "Making BLAST databases"
makeblastdb -input_type fasta -in "${NEW_REF}" -dbtype nucl -out "${NEW_REF}" -parse_seqids &
makeblastdb -input_type fasta -in "${OLD_REF}" -dbtype nucl -out "${OLD_REF}" -parse_seqids &
wait

# Run TBLASTN on both assemblies in parallel
echo "Running both assemblies"
tblastn -query "${QUERY}" -db "${NEW_REF}" -out "new.7tm.tblastn.outfmt6" -outfmt 6 -num_threads 16 &
tblastn -query "${QUERY}" -db "${OLD_REF}" -out "old.7tm.tblastn.outfmt6" -outfmt 6 -num_threads 16 &
wait

echo "OLD assembly complete!"
