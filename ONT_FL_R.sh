#!/bin/sh
#SBATCH --account=biosci
#SBATCH --partition=grace
#SBATCH --time=130:00:00
#SBATCH --nodes=1 --ntasks=24
#SBATCH --job-name="ONT_FL_R"
#SBATCH --mail-user=kathryn.morrissey@uct.ac.za
#SBATCH --mail-type=ALL

module load python/miniconda3-qiime2

qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path /scratch/kathryn/ONT_FL/4_QC \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path /scratch/kathryn/ONT_FL_R/4_single-end-demux.qza

qiime demux summarize \
--i-data /scratch/kathryn/ONT_FL_R/4_single-end-demux.qza \
--o-visualization /scratch/kathryn/ONT_FL_R/4_single-end-demux.qzv

qiime vsearch dereplicate-sequences \
--i-sequences /scratch/kathryn/ONT_FL_R/4_single-end-demux.qza \
--o-dereplicated-table /scratch/kathryn/ONT_FL_R/5_derep-table.qza \
--o-dereplicated-sequences /scratch/kathryn/ONT_FL_R/5_derep-seqs.qza

qiime feature-table tabulate-seqs \
--i-data /scratch/kathryn/ONT_FL_R/5_derep-seqs.qza \
--o-visualization /scratch/kathryn/ONT_FL_R/5_derep-seqs.qzv

qiime feature-table summarize \
--i-table /scratch/kathryn/ONT_FL_R/5_derep-table.qza \
--o-visualization /scratch/kathryn/ONT_FL_R/5_derep-table.qzv

qiime vsearch uchime-ref \
--i-table /scratch/kathryn/ONT_FL_R/5_derep-table.qza \
--i-sequences /scratch/kathryn/ONT_FL_R/5_derep-seqs.qza \
--i-reference-sequences /home/kathryn/RDP/RDP_ext_ref-seqs.qza \
--output-dir /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out \
--p-threads 24

qiime metadata tabulate \
--m-input-file /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/stats.qza \
--o-visualization /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/stats.qzv

#filtering chimeric seqeuences
qiime feature-table filter-features \
--i-table /scratch/kathryn/ONT_FL_R/5_derep-table.qza \
--m-metadata-file /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/nonchimeras.qza \
--o-filtered-table /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
--i-data /scratch/kathryn/ONT_FL_R/5_derep-seqs.qza \
--m-metadata-file /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/nonchimeras.qza \
--o-filtered-data /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/rep-seqs-nonchimeric-wo-borderline.qza

#visualizing non-chimeric data...
qiime feature-table summarize \
--i-table /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qza \
--o-visualization /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qzv

#OTU clustering
qiime vsearch cluster-features-open-reference \
--i-table /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qza \
--i-sequences /scratch/kathryn/ONT_FL_R/6.1_uchime-ref-out/rep-seqs-nonchimeric-wo-borderline.qza \
--i-reference-sequences /home/kathryn/RDP/RDP_ext_ref-seqs.qza \
--p-perc-identity 0.80 \
--o-clustered-table /scratch/kathryn/ONT_FL_R/6.2_table-op_ref-85.qza \
--o-clustered-sequences /scratch/kathryn/ONT_FL_R/6.2_rep-seqs-op_ref-85.qza \
--o-new-reference-sequences /scratch/kathryn/ONT_FL_R/6.2_new-ref-seqs-op_ref-85.qza \
--p-threads 24

#Finally assigning taxonomy
qiime feature-classifier classify-sklearn \
--i-classifier /home/kathryn/RDP/RDP_uniform-classifier.qza \
--i-reads /scratch/kathryn/ONT_FL_R/6.2_rep-seqs-op_ref-85.qza \
--o-classification /scratch/kathryn/ONT_FL_R/11_taxonomy-sklearn.qza \
--p-reads-per-batch 5000

#  make a directory for files to be exported 
mkdir /scratch/kathryn/ONT_FL_R/exported

#  converting all necessary files for Phyloseq input
qiime tools export --input-path /scratch/kathryn/ONT_FL_R/6.2_table-op_ref-85.qza --output-path /scratch/kathryn/ONT_FL_R/exported/

#export 11_taxonomy-sklearn.qza to text file (taxonomy.tsv)"
qiime tools export --input-path /scratch/kathryn/ONT_FL_R/11_taxonomy-sklearn.qza --output-path /scratch/kathryn/ONT_FL_R/exported/

#export 10_rooted-tree.qza to newick file format (tree.nwk)"

qiime tools export --input-path /scratch/kathryn/ONT_FL_R/6.2_rep-seqs-op_ref-85.qza --output-path /scratch/kathryn/ONT_FL_R/exported/


#edit the taxonomy.tsv file for compatible header line
sed -i -e 's/Feature ID/#OTUID/g; s/Taxon/taxonomy/g; s/Confidence/confidence/g' /scratch/kathryn/ONT_FL_R/exported/taxonomy.tsv

#Finally, add taxonomy to biom file
biom add-metadata -i /scratch/kathryn/ONT_FL_R/exported/feature-table.biom -o /scratch/kathryn/ONT_FL_R/exported/table-with-taxonomy.biom --observation-metadata-fp /scratch/kathryn/ONT_FL_R/exported/taxonomy.tsv --sc-separated taxonomy

