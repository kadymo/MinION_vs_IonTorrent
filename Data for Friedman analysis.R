library(ggplot2)
library(phyloseq)
library("tidyverse") 
library(plyr)
library(vegan)
library(phylosmith)
library("data.table")
library(cowplot)
require("pheatmap")



setwd("~/Dropbox/MinION protocols/Methods paper/R_analysis/Data_files/Phyloseq_objects")


tax<-read.csv2("Edited_tax_table.csv", header = T)
tax<-column_to_rownames(tax, "OTU")  
tax<- as.matrix(tax)
otu<-read.csv2("Edited_otu_table.csv", header = T)
otu<-column_to_rownames(otu, "OTU")

mapfile = "Sample_data_all.csv"
map <- read.csv2(mapfile)
rownames(map) <- map$Sample_name
map<-map[,-1]              


tax_final<-tax_table(tax)
otu_final<- otu_table(otu, taxa_are_rows = TRUE)
MAP<-sample_data(map)

data<-merge_phyloseq(tax_final,otu_final,MAP)
data


##PPD column added to sample data####
variable1 = sample_data(data)$Platform
variable2 = sample_data(data)$Pipeline
variable3 = sample_data(data)$Database
sample_data(data)$PPD <-  mapply(paste0, variable1, sep = "_", variable2, sep = "_", variable3)

## SP column added to sample data####
variable1 = sample_data(data)$Sample_original
variable2 = sample_data(data)$Platform
sample_data(data)$SP <-  mapply(paste0, variable1, sep = "_", variable2)

## SD column added to sample data####
variable1 = sample_data(data)$Sample_original
variable2 = sample_data(data)$Database
sample_data(data)$SD <-  mapply(paste0, variable1, sep = "_", variable2)

## PRP column added to sample data####
variable1 = sample_data(data)$Platform
variable2 = sample_data(data)$Region
variable3 = sample_data(data)$Pipeline
sample_data(data)$PRP <-  mapply(paste0, variable1, sep = "_", variable2, sep = "_", variable3)

## PRPD column added to sample data####
variable1 = sample_data(data)$Platform
variable2 = sample_data(data)$Region
variable3 = sample_data(data)$Pipeline
variable4 = sample_data(data)$Database
sample_data(data)$PRPD <- mapply(paste0, variable1, sep = "_", variable2, sep = "_", variable3, sep = "_", variable4)

## PR column added to sample data####
variable1 = sample_data(data)$Platform
variable2 = sample_data(data)$Region
sample_data(data)$PR <- mapply(paste0, variable1, sep = "_", variable2)



###FILTER#####

data.1 <- subset_taxa(data, Kingdom == "Bacteria" & Class != "Chloroplast" & Order != "Chloroplast" & Family != "mitochondria" &  Phylum != "Cyanobacteria/Chloroplast")
data.2<- subset_taxa(data.1, Phylum != "Bacteria_X" & Phylum != "Unclassified_X")

test<-as.data.frame(tax_table(data.1))

data.3<-subset_samples(data.2, PR !="ONT_V4")

###Microbiome#####
microbiome<-subset_samples(data.3, Type == "Sample")
microbiome_phylum<-tax_glom(microbiome, taxrank="Phylum")
microbiome_phylum_relab<-transform_sample_counts(microbiome_phylum, function(x) 100 * x/sum(x))

microbiome_tax<-as.data.frame(tax_table(microbiome_phylum_relab))
#write.csv2(microbiome_tax, "microbiome_tax.csv")

microbiome_otu<-as.data.frame(otu_table(microbiome_phylum_relab))
#write.csv2(microbiome_otu, "microbiome_otu.csv")

microbiome_sam<-as.data.frame(sample_data(microbiome_phylum_relab))
#write.csv2(microbiome_sam, "microbiome_sam.csv")


###DNA Standard#####
standard<-subset_samples(data.3, Type == "Standard")
standard_phylum<-tax_glom(standard, taxrank="Phylum")
standard_phylum_relab<-transform_sample_counts(standard_phylum, function(x) 100 * x/sum(x))

standard_tax<-as.data.frame(tax_table(standard_phylum_relab))
#write.csv2(standard_tax, "standard_tax.csv")

standard_otu<-as.data.frame(otu_table(standard_phylum_relab))
#write.csv2(standard_otu, "standard_otu.csv")

standard_sam<-as.data.frame(sample_data(standard_phylum_relab))
#write.csv2(standard_sam, "standard_sam.csv")

