# Google Cloud Shell: Syntenic Visualisation of Whole-Genome Alignments to Compare an Old and New Assembly

**Hydrophis major**

## Whole Chromosome-level Genome Assembly Synteny Plot

### Comparing New and Old Hydrophis Major Chromosome-Level Assemblies

### Start Google Cloud Shell VM + prime

```bash
gcloud compute instances create fast \
    --project=eternal-outlook-478022-s4 \
    --zone=us-central1-a \
    --machine-type=n2-standard-32 \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --boot-disk-size=500GB \
    --tags=http-server,https-server \
    --create-disk=name=fast-ssd-data,size=500GB,mode=rw,type=pd-ssd \
    --scopes=cloud-platform
```

```bash
gcloud compute ssh fast --zone=us-central1-a --project=eternal-outlook-478022-s4
```

```bash
curl -O Anaconda3-2025.12-1-Linux-x86_64.sh
bash Anaconda3-2025.12-1-Linux-x86_64.sh
```

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y git wget curl build-essential
```

### Seqkit installation

```bash
cd ~
wget https://github.com/shenwei356/seqkit/releases/download/v2.12.0/seqkit_linux_amd64.tar.gz
tar -xzf seqkit_linux_amd64.tar.gz
chmod +x seqkit
mkdir -p ~/bin
mv seqkit ~/bin
export PATH=$HOME/bin:$PATH
seqkit --version
```

### MUMmer installation

```bash
wget https://github.com/mummer4/mummer/releases/download/v4.0.0/mummer-4.0.0.tar.gz
tar -xzf mummer-4.0.0.tar.gz
cd mummer-4.0.0
./configure --prefix=$HOME/mummer
make
make install
export PATH=$HOME/mummer/bin:$PATH
nucmer --version
```

### Loading files to VM from Google Cloud Storage

```bash
mkdir syn
cd syn
gcloud auth login
```

```bash
gsutil -m cp gs://summer25-26_sanders-lab/old-filtered.fa .
gsutil -m cp gs://summer25-26_sanders-lab/fixed-new-filtered.fa .
```

### Check file format

```bash
REF=fixed-new-filtered.fa QUERY=old-filtered.fa
```

```bash
seqkit stats *.fa --all
```

```bash
seqkit seq --name --only-id $REF
seqkit seq --name --only-id $QUERY
```

### Fix formatting

```bash
cat <(seqkit grep -v -r -p '^chrZ$' $REF) <(seqkit grep -r -p '^chrZ$' $REF) > new-hmaj-assembly.fa # order chrZ last, like in QUERY
```

```bash
grep "^>" new-hmaj-assembly.fa && grep "^>" $QUERY # check
```

### Reference file

```bash
seqkit faidx new-hmaj-assembly.fa && seqkit faidx old-filtered.fa
head old-filtered.fa.fai && head new-hmaj-assembly.fa.fai
```

```bash
cat new-hmaj-assembly.fa.fai old-filtered.fai > alignment.ref.tsv # reference file
cat alignment.ref.tsv # the first two columns being the names and sizes
```

### Alignment

```bash
REF=new-hmaj-assembly.fa QUERY=old-filtered.fa NAME=alignment
nucmer --prefix $NAME $REF $QUERY
```

```bash
head $NAME.delta
```

```bash
show-coords -lTH $NAME.delta > $NAME.coords
```

```bash
NEW_HEADER="reference_start,reference_end,query_start,query_end,reference_alignment_length,query_alignment_length,percent_identity,reference_length,query_length,reference_chromosome,query_chromosome"
```

```bash
(echo "$NEW_HEADER" && awk '{$1=$1}1' OFS="," $NAME.coords) > $NAME.csv
head $NAME.csv
```

```bash
awk '{print $1 "," $2}' old-filtered.fa.fai > old_ref.csv
sed -i '1i chromosome,size' old_ref.csv
```

```bash
awk '{print $1 "," $2}' new-hmaj-assembly.fa.fai > new_ref.csv
sed -i '1i chromosome,size' new_ref.csv
```

```bash
echo "chromosome,size,assembly" > combined_ref.csv && awk '{print $1 "," $2 ",old"}' old-filtered.fa.fai >> combined_ref.csv && awk '{print $1 "," $2 ",new"}' new-hmaj-assembly.fa.fai >> combined_ref.csv
```

```bash
awk -F',' 'NR==1 || ($5>=1000 && $6>=1000)' alignment.csv > filtered-alignment.csv
```
### Download files

```bash
gsutil -m cp filtered-alignment.csv gs://summer25-26_sanders-lab/ 
gsutil -m cp combined_ref.csv gs://summer25-26_sanders-lab/
wait echo "Files Downloaded"
```
