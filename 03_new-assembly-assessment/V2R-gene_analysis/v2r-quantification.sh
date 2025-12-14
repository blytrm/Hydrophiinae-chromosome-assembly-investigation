# Old Assembly V2R Quantification
out_csv="old_v2r_summary.csv"
echo "chromosome,count" > "$out_csv"

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chrZ; do
  count=$(
    awk -v c="$chr" '$7==1 && $8==$4 && $4>=200 && $11<=1e-10 && $3>=60 && $2==c' "old.7tm.tblastn.outfmt6" | \
      awk '{start = ($9 < $10) ? $9 : $10; end = ($9 < $10) ? $10 : $9; print $2 "\t" start-1 "\t" end}' | \
      sort -k1,1 -k2,2n | \
      ./bedtools2/bin/bedtools merge -i - | wc -l
  )
  echo "$chr,$count" >> "$out_csv"
done

wg=$(awk '$7==1 && $8==$4 && $4>=200 && $11<=1e-10 && $3>=60' "old.7tm.tblastn.outfmt6" | \
    awk '{start = ($9 < $10) ? $9 : $10; end = ($9 < $10) ? $10 : $9; print $2 "\t" start-1 "\t" end}' | \
    sort -k1,1 -k2,2n | \
    ./bedtools2/bin/bedtools merge -i - | wc -l)

echo "whole_genome,$wg" >> "$out_csv"

# New Assembly V2R Quantification
out_csv="new_v2r_summary.csv"
echo "chromosome,count" > "$out_csv"

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chrZ; do
  count=$(
    awk -v c="$chr" '$7==1 && $8==$4 && $4>=200 && $11<=1e-10 && $3>=60 && $2==c' "new.7tm.tblastn.outfmt6" | \
      awk '{start = ($9 < $10) ? $9 : $10; end = ($9 < $10) ? $10 : $9; print $2 "\t" start-1 "\t" end}' | \
      sort -k1,1 -k2,2n | \
      ./bedtools2/bin/bedtools merge -i - | wc -l
  )
  echo "$chr,$count" >> "$out_csv"
done

out_csv="old_v2r_summary.csv"
echo "chromosome,count" > "$out_csv"
