#!/bin/sh
#SBATCH --account=biosci
#SBATCH --partition=ada
#SBATCH --time=130:00:00
#SBATCH --nodes=1 --ntasks=40
#SBATCH --job-name="ONT_Q_V4_s_sampleX"
#SBATCH --mail-user=kathryn.morrissey@uct.ac.za
#SBATCH --mail-type=ALL

module load python/miniconda3-qiime2

qiime vsearch dereplicate-sequences \
--i-sequences /scratch/kathryn/ONT_V4_qual/sampleX.qza \
--o-dereplicated-table /scratch/kathryn/ONT_V4/sampleX/5_derep-table.qza \
--o-dereplicated-sequences /scratch/kathryn/ONT_V4/sampleX/5_derep-seqs.qza

qiime feature-table tabulate-seqs \
--i-data /scratch/kathryn/ONT_V4/sampleX/5_derep-seqs.qza \
--o-visualization /scratch/kathryn/ONT_V4/sampleX/5_derep-seqs.qzv

qiime feature-table summarize \
--i-table /scratch/kathryn/ONT_V4/sampleX/5_derep-table.qza \
--o-visualization /scratch/kathryn/ONT_V4/sampleX/5_derep-table.qzv

qiime vsearch uchime-ref \
--i-table /scratch/kathryn/ONT_V4/sampleX/5_derep-table.qza \
--i-sequences /scratch/kathryn/ONT_V4/sampleX/5_derep-seqs.qza \
--i-reference-sequences /home/kathryn/SILVA/99_otus.qza \
--output-dir /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out \
--p-threads 40

qiime metadata tabulate \
--m-input-file /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/stats.qza \
--o-visualization /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/stats.qzv

#filtering chimeric seqeuences
qiime feature-table filter-features \
--i-table /scratch/kathryn/ONT_V4/sampleX/5_derep-table.qza \
--m-metadata-file /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/nonchimeras.qza \
--o-filtered-table /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
--i-data /scratch/kathryn/ONT_V4/sampleX/5_derep-seqs.qza \
--m-metadata-file /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/nonchimeras.qza \
--o-filtered-data /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/rep-seqs-nonchimeric-wo-borderline.qza

#visualizing non-chimeric data...
qiime feature-table summarize \
--i-table /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qza \
--o-visualization /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qzv

#OTU clustering
qiime vsearch cluster-features-open-reference \
--i-table /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/table-nonchimeric-wo-borderline.qza \
--i-sequences /scratch/kathryn/ONT_V4/sampleX/6.1_uchime-ref-out/rep-seqs-nonchimeric-wo-borderline.qza \
--i-reference-sequences /home/kathryn/SILVA/99_otus.qza \
--p-perc-identity 0.80 \
--o-clustered-table /scratch/kathryn/ONT_V4/sampleX/6.2_table-op_ref-85.qza \
--o-clustered-sequences /scratch/kathryn/ONT_V4/sampleX/6.2_rep-seqs-op_ref-85.qza \
--o-new-reference-sequences /scratch/kathryn/ONT_V4/sampleX/6.2_new-ref-seqs-op_ref-85.qza \
--p-threads 40

#Finally assigning taxonomy
qiime feature-classifier classify-sklearn \
--i-classifier /home/kathryn/SILVA/classifier.qza \
--i-reads /scratch/kathryn/ONT_V4/sampleX/6.2_rep-seqs-op_ref-85.qza \
--o-classification /scratch/kathryn/ONT_V4/sampleX/11_taxonomy-sklearn.qza \
--p-reads-per-batch 5000
--p-n-jobs -1

#  make a directory for files to be sampleX_exported 
mkdir /scratch/kathryn/ONT_V4/sampleX/sampleX_exported

#  converting all necessary files for Phyloseq input
qiime tools export --input-path /scratch/kathryn/ONT_V4/sampleX/6.2_table-op_ref-85.qza --output-path /scratch/kathryn/ONT_V4/sampleX/sampleX_exported/

#export 11_taxonomy-sklearn.qza to text file (taxonomy.tsv)"
qiime tools export --input-path /scratch/kathryn/ONT_V4/sampleX/11_taxonomy-sklearn.qza --output-path /scratch/kathryn/ONT_V4/sampleX/sampleX_exported/

#export 10_rooted-tree.qza to newick file format (tree.nwk)"

qiime tools export --input-path /scratch/kathryn/ONT_V4/sampleX/6.2_rep-seqs-op_ref-85.qza --output-path /scratch/kathryn/ONT_V4/sampleX/sampleX_exported/


#edONT the taxonomy.tsv file for compatible header line
sed -i -e 's/Feature ID/#OTUID/g; s/Taxon/taxonomy/g; s/Confidence/confidence/g' /scratch/kathryn/ONT_V4/sampleX/sampleX_exported/taxonomy.tsv

#Finally, add taxonomy to biom file
biom add-metadata -i /scratch/kathryn/ONT_V4/sampleX/sampleX_exported/feature-table.biom -o /scratch/kathryn/ONT_V4/sampleX/sampleX_exported/table-wONTh-taxonomy.biom --observation-metadata-fp /scratch/kathryn/ONT_V4/sampleX/sampleX_exported/taxonomy.tsv --sc-separated taxonomy

