library("DESeq2")
library("phyloseq")
library("ggplot2")
library("plyr")
library("ape")
library("microbiome")
library("microbiomeSeq") 
library("adespatial")
library("igraph")
library("ggpubr")

setwd("~/Desktop/UCT projects/Oudekraal/Methods paper/R")

#create a tree:
raandom_tree =rtree(ntaxa(data.1), rooted=TRUE, tip.label=taxa_names(data.1))
plot(raandom_tree)

#merge tree data with data.1 phyloseq object:

data.1.1 = merge_phyloseq(data.1, raandom_tree)
data.1.1

#remove DNAs samples:
Data.1.2 = subset_samples(data.1.1, Sample_type != "None")

#Remove ONT V4 region:

Data.1.3 = subset_samples(Data.1.2, PR !="ONT_V4")
Data.1.3
sample_data(Data.1.3)

#divide into panels: 
sample_data(Data.1.3)

#new sample column: Platform_database_Pipeline:
variable3 = sample_data(Data.1.3)$Platform
variable4 = sample_data(Data.1.3)$Database
sample_data(Data.1.3)$P_DB <-  mapply(paste0, variable3, sep = "_", variable4)
variable6 = sample_data(Data.1.3)$P_DB
variable5 = sample_data(Data.1.3)$Pipeline
sample_data(Data.1.3)$PDBP <-  mapply(paste0, variable6, sep = "_", variable5)
sample_data(Data.1.3)

#Mothurobject:
Data.1.3.M <- subset_samples(Data.1.3, Pipeline =="MOTHUR")
Data.1.3.MS <- subset_samples(Data.1.3.M, Database =="SILVA")
Data.1.3.MR <- subset_samples(Data.1.3.M, Database =="RDP")
sample_sums(Data.1.3.MR)
# Discard samples that have less than 10000 reads...?
#Data.1.3.MR.filt <- prune_samples(sample_sums(Data.1.3.MR) >= 10000, Data.1.3.MR)
#sample_sums(Data.1.3.MR.filt)

#QIIME2 object:
Data.1.3.Q <- subset_samples(Data.1.3, Pipeline =="QIIME")
Data.1.3.QS <- subset_samples(Data.1.3.Q, Database == "SILVA")
Data.1.3.QR <- subset_samples(Data.1.3.Q, Database == "RDP")

#import Phyloseq object into DESeq: (works!)
diagddsQS = phyloseq_to_deseq2(Data.1.3.QS, ~ Platform)
diagddsQS = DESeq(diagddsQS, test="Wald", sfType = "poscounts")

#investigate test results table:

res = results(diagddsQS, cooksCutoff = FALSE)
res
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Data.1.2)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#Log2fold change (Genus level):

#QS:
data3QS <- subset_taxa(Data.1.3.QS, Genus != "Unclassified")
data3QS <- subset_taxa(data3QS, Genus != "Unclassified_XXXXX")
#data3QS <- subset_taxa(data3QS, Genus != "Chloroplast_uncultured.marine.eukaryote_X")
#data3QS <- subset_taxa(data3QS, Genus != "Chloroplast_XX")
data4QS = tax_glom(data3QS, "Genus")
#prune to top 15 taxa:
Genus15QS = prune_taxa(names(sort(taxa_sums(data4QS), TRUE))[1:15], data4QS)
otu_table(Genus15QS) <- t(otu_table(Genus15QS))
physeqQS<- taxa_level(Genus15QS, "Genus")
deseq_sigQS <- differential_abundance(physeqQS, grouping_column = "Platform", output_norm = NULL, pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)


#QR:
data3QR <- subset_taxa(Data.1.3.QR, Genus != "Unclassified_XXXXX")
data3QR <- subset_taxa(data3QR, Genus != "Mitochondria_X")
data3QR <- subset_taxa(Data.1.3.QR, Genus != "Unclassified")
data4QR = tax_glom(data3QR, "Genus")
#prune to top 15 taxa:
Genus15QR = prune_taxa(names(sort(taxa_sums(data4QR), TRUE))[1:15], data4QR)
otu_table(Genus15QR) <- t(otu_table(Genus15QR))
physeqQR<- taxa_level(Genus15QR, "Genus")
deseq_sigQR <- differential_abundance(physeqQR, grouping_column = "Platform", output_norm = NULL, pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

#MS:
data3MS <- subset_taxa(Data.1.3.MS, Genus != "Unclassified_XXXXX")
data3MS <- subset_taxa(data3MS, Genus != "Mitochondria_X")
data3MR <- subset_taxa(data3MR, Genus != "Bacteria_XXXXX")
data4MS = tax_glom(data3MS, "Genus")
#prune to top 15 taxa:
Genus15MS = prune_taxa(names(sort(taxa_sums(data4MS), TRUE))[1:15], data4MS)
otu_table(Genus15MS) <- t(otu_table(Genus15MS))
physeqMS<- taxa_level(Genus15MS, "Genus")
deseq_sigMS <- differential_abundance(physeqMS, grouping_column = "Platform", output_norm = NULL, pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

#MR:
data3MR <- subset_taxa(Data.1.3.MR, Genus != "Unclassified_XXXXX")
data3MR <- subset_taxa(data3MR, Genus != "Chloroplast_XX")
data3MR <- subset_taxa(data3MR, Genus != "Bacteria_XXXXX")
data3MR <- subset_taxa(data3MR, Genus != "Mitochondria_X")
data4MR = tax_glom(data3MR, "Genus")
#prune to top 15 taxa:
Genus15MR = prune_taxa(names(sort(taxa_sums(data4MR), TRUE))[1:15], data4MR)
otu_table(Genus15MR) <- t(otu_table(Genus15MR))
physeqMR<- taxa_level(Genus15MR, "Genus")
deseq_sigMR <- differential_abundance(physeqMR, grouping_column = "Platform", output_norm = NULL, pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

#create colour palette:
Database_colours <- list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))
platform_region_colour <- list(PR = c(IT_V4 = "#80A4D8", ONT_V4 = "#AAD794", ONT_V4="#FAD37B"))

#To generate a plot showing differentially abundant taxa between among compared groups , 
#corresponding adjusted p-values and rank of importance as detected by random forest classifier:
#QR:
p_QR<-plot_signif(deseq_sigQR$plotdata, top.taxa = 20)
print(p_QR)
plot_MDA(deseq_sigQR$importance, top.taxa=20)
#logfold change:
p_QR <- plot_MA(deseq_sigQR$SignFeaturesTable, label=F)
QR <- print(p_QR$lfcplot)
QR <- QR + theme(legend.position = "none") + scale_fill_manual(values=c("#80A4D8","#AAD794")) + coord_flip()
QR <- QR + coord_flip()
QR

#QS:
p_QS<-plot_signif(deseq_sigQS$plotdata, top.taxa = 20)
print(p_QS)
plot_MDA(deseq_sigQS$importance, top.taxa=20)
#logfold change:
p_QS <- plot_MA(deseq_sigQS$SignFeaturesTable, label=T)
QS <- print(p_QS$lfcplot)
QS <- QS + theme(legend.position = "none") + scale_fill_manual(values=c("#80A4D8","#AAD794")) + coord_flip()
QS <- QS + coord_flip()
QS

#MR:
p_MR<-plot_signif(deseq_sigMR$plotdata, top.taxa = 20)
print(p_MR)
plot_MDA(deseq_sigMR$importance, top.taxa=20)
#logfold change:
p_MR <- plot_MA(deseq_sigMR$SignFeaturesTable, label=T)
MR <- print(p_MR$lfcplot)
MR <- MR + theme(legend.position = "none") + scale_fill_manual(values=c("#80A4D8","#AAD794")) + coord_flip()
MR <- MR + coord_flip()
MR

#MS:
p_MS<-plot_signif(deseq_sigMS$plotdata, top.taxa = 20)
print(p_MS)
plot_MDA(deseq_sigMS$importance, top.taxa=20)
#logfold change:
p_MS <- plot_MA(deseq_sigMS$SignFeaturesTable, label=T, main="MOTHUR | SILVA")
MS <- print(p_MS$lfcplot)
MS <- MS + theme(legend.position = "none") + scale_fill_manual(values=c("#80A4D8","#AAD794")) + coord_flip()
MS <- MS + coord_flip()
MS

# 4 figures arranged in 2 rows and 2 columns

library(cowplot)
plot_grid(QS, QR, MS, MR, labels="AUTO", ncol=2)

# The information provided here is strictly limited to log fold change and basemean values annotated for each taxa/feature, it help with identification of features with extreme log folds more easily.
#How do we only plot the top 20??
#p1_QR <- plotMA(deseq_sigQS$SignFeaturesTable, label=F)

#using ggpubr:(must input results object)
ggmaplot(sigtab, main = expression("ONT" %->% "IT"),
         fdr = 0.05, fc = 2, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(sigtab$Genus),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

#  print(p$maplot)
print(p1_QR$lfcplot)


# Ordination --------------------------------------------------------------

#relative abundance:

Data_rel <- normalise_data(data.1.1, norm.method = "relative")

#sort by abundance:

myTaxa = names(sort(taxa_sums(data.1.1), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, data.1.1)
plot(phy_tree(ex1), show.node.label = TRUE)
ex1
plot_tree(ex1, color = "Platform", label.tips = "Phylum", ladderize = "left", justify = "left")

#log transformation:
xt <- transform(data.1.1, 'log10')
xt

#transform rel abundance:
GP1 = transform_sample_counts(data.1.1, function(x) 100 * x/sum(x))


#keep the most abundant 5 phyla:
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

GP1

#PCoA (unweighted unifrac) by platform:
ordu = ordinate(GP1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(GP1, ordu, color="Platform", shape="Sample_type")

#plot a tree by sample:

?plot_tree
sample_names(data.1.1)
sample_variables(ex1)

#Unifrac:


theme_set(theme_bw())

UniFrac(data.1.1, weighted =TRUE, normalized = TRUE, parallel=TRUE, fast=TRUE)

#list distance functions:
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
