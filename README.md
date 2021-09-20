# MinION_methods

# Project description
The Oxford Nanopore MinION is one of the newer third generation sequencing technologies available, is significantly more affordable, more user friendly and results in much longer read lengths than its competitors. For this study, 16S amplicons were sequenced using both the Nanopore MinION (ONT) and the Ion torrent S5 (IT) platforms. Two databases (SILVA & RDP) for taxonomic assignment and two common pipelines (QIIME & MOTHUR) were used for downstream analysis. The combinations of each were compared and assessed for their accuracy and coverage. We compared the 16S amplicons: IT-V4; ONT-V4; & ONT-FL (full-length).

# Samples used
Samples come from three different marine waters:
- Oudekraal (beach, label ODK)
- Muizenberg (surfzone, label AB)
- IEP (shelf - cruise, label N)

A DNA standard was used to benchmark (mock community)

Sample numbers: ODK3b; ODK5b, ODK7b; ODK9a; (ODK11b - removed); ODK13b; ODK15a; ODK21b; ODK36a; AB1.2; N42; DNAs1; DNAs2

# Platforms used
**ONT - MinION sequencing (flow cells chemistry R9.4.1)**

Run 1 (FAN-53765): ODK3b; ODK5b, ODK7b; ODK9a; (ODK11b - removed); ODK13b; ODK15a; ODK21b; ODK36a; DNAs1 

Run 2 (FAN-41723): AB1.2; N42; DNAs2

**IT - Ion Torrent (Ion Chef Instrument, templating kits Ion 510 amp, Ion 52 amp, and Ion 530 Kit – Chef)**

Run 1: ODK3b; ODK5b, ODK7b; ODK9a; (ODK11b - removed); ODK13b; ODK15a; ODK21b; ODK36a; AB1.2; N42; DNAs1; DNAs2

# Pipelines used
QIIME (v2020.8)

MOTHUR (vxx)

# Databases used
SILVA (vxx)

RDP (vxx)

# Data processing - ONT
**Raw sequences produced (fastq)**

**Quality control**
1. qcat on all barcodes (each run done separately) to trim adapters and barcodes: qcat -f input.fastq --trim -b output_folder
2. Select relevant barcodes/samples
3. Get quality reads (each run done separately): NanoFilt -q 10 -l 1300 --maxlength 1600

**File organization = for ONT V4**
1. Fastq files converted into .fasta using seqtk
2. Sequence names relabeled in each sample with samplename+number (script by Kat – relabel_seqs.sh)

**Select V4 region in Mothur**
1. ran MinION_Batch1.batch
2. ran MinION_Batch2.batch
3. ran degap.seqs() on good.align file

**Dereplicate, remove chimeras & cluster**
1. Using Qiime (separately for each barcode V4 = ONT_V4_sample_{database}.sh | Full length = ONT_FL_{database}.sh)
2. Using Mothur

**Assign taxonomy (as part of script above)**
1. SILVA (separately for each barcode V4 = ONT_V4_sample_S.sh | Full length = ONT_FL_SILVA.sh)
2. RDP (separately for each barcode V4 = ONT_V4_sample_R.sh | Full length = ONT_FL_RDP.sh)

**Output files**
1.Mothur
- taxonomy(.taxonomy)
- count (.count_table)

2.Qiime
- biom (.biom)
- taxonomy (.tsv)
- treefile (tree.nwk)

**Convert biom to txt**

in qiime: biom convert -i feature-table.biom -o table-from-biom.txt --to-tsv

# Data processing - IT

# Importing into phyloseq
Impoting table-from-biom.txt and taxonomy.tsv for each sample

**Before importing we need to fix file headers to be readable into R**

To do this follow the following commands in terminal/unix environment:

sed -i -e "1d" *.txt

sed -i -e "s/ //" *.txt

sed -i -e "s/#OTUID/OTUID/" *.txt

sed -i -e "s/#OTUID/OTUID/" *.tsv

**Data import into R**

Use the data_importer.R function called in the script importing_data.R

**Making phyloseq objects**

Create phyloseq object using the script data_to_phyloseq.R
