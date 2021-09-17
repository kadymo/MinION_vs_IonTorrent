# MinION_methods

# Project description
The Oxford Nanopore MinION is one of the newer third generation sequencing technologies available, is significantly more affordable, more user friendly and results in much longer read lengths than its competitors. For this study, 16S amplicons were sequenced using both the Nanopore MinION and the Ion torrent S5 platforms. Two databases (SILVA & RDP) for taxonomic assignment and two common pipelines (QIIME & MOTHUR) were used for downstream analysis. The combinations of each were compared and assessed for their accuracy and coverage. 

# Samples used
Samples come from three different marine waters:
- Oudekraal (beach, designation ODK)
- Muizenberg (surfzone, designation AB)
- IEP (shelf - cruise, designation N)

A DNA standard was used to benchmark (mock community)

Sample numbers: ODK3b; ODK5b, ODK7b; ODK9a; (ODK11b - removed); ODK13b; ODK15a; ODK21b; AB1.2; N42; DNAs1; DNAs2

# Platforms used
**ONT - MinION sequencing (flow cells chemistry R9.4.1)**

Run 1 (FAN-53765): ODK3b; ODK5b, ODK7b; ODK9a; (ODK11b - removed); ODK13b; ODK15a; ODK21b; DNAs1 

Run 2 (FAN-41723): AB1.2; N42; DNAs2

**IT - Ion Torrent (Ion Chef Instrument, templating kits Ion 510 amp, Ion 52 amp, and Ion 530 Kit – Chef)**

Run 1: ODK3b; ODK5b, ODK7b; ODK9a; (ODK11b - removed); ODK13b; ODK15a; ODK21b; AB1.2; N42; DNAs1; DNAs2

# Pipleines used
QIIME

MOTHUR

# Databases used
SILVA

RDP

# Data processing
**Raw sequences produced (fastq)**

**Quality control - ONT**
1. qcat on all barcodes (each run done separately) to trim adapters and barcodes: qcat -f input.fastq --trim -b output_folder
2. Select relevant barcodes/samples
3. Get quality reads (each run done separately): NanoFilt -q 10 -l 1300 --maxlength 1600

**File organization = ONT V4**
1. Fastq files converted into .fasta using seqtk
2. Sequence names relabeled in each sample with samplename+number (script by Kat – bash script with sed)

**Select V4 region - for ONT only**
1. ran Batch script 1
2. ran Batch script 2
3. ran degap.seqs() on good.align file

**Dereplicate, remove chimeras & cluster**
1. Using Qiime (separately for each barcode V4 = V4_sample_script.sh | Full length = ONT_FL_script.sh)
2. Using mothur

**Assign taxonomy (as part of script above)**
1. SILVA
2. RDP

**Output files**
- biom (.biom)
- taxonomy (.tsv)
- treefile (tree.nwk)

**Convert biom to txt**

in qiime: biom convert -i feature-table.biom -o table-from-biom.txt --to-tsv

# Importing into phyloseq
