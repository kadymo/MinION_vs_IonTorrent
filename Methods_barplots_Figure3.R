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

###Phyloseq object after edit####

tax<-read.csv2("All_Final_TAX.csv", header = T)
tax<-column_to_rownames(tax, "OTU")  
tax<- as.matrix(tax)
otu<-read.csv2("All_Final_OTU.csv", header = T)
otu<-column_to_rownames(otu, "X")

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

###Data set up####

##Mothur##

EMU_data<-subset_samples(data.2, Pipeline == "EMU")
EMU_data_phylum<-tax_glom(EMU_data, taxrank="Phylum")
EMU_data_class<-tax_glom(EMU_data, taxrank="Class")
#EMU_data_order<-tax_glom(mothur_data, taxrank="Order")
#EMU_data_family<-tax_glom(mothur_data, taxrank="Family")
#EMU_genus<-tax_glom(mothur_data, taxrank="Genus")
#EMU_data_species<-tax_glom(mothur_data, taxrank="Species")

##Qiime##

#Qiime_data<-subset_samples(data.2, Pipeline == "QIIME")
#Qiime_data_phylum<-tax_glom(Qiime_data, taxrank="Phylum")
#Qiime_data_class<-tax_glom(Qiime_data, taxrank="Class")
#Qiime_data_order<-tax_glom(Qiime_data, taxrank="Order")
#Qiime_data_family<-tax_glom(Qiime_data, taxrank="Family")
#Qiime_data_genus<-tax_glom(Qiime_data, taxrank="Genus")
#Qiime_data_species<-tax_glom(Qiime_data, taxrank="Species")


####Barplots####

###Barplot final##### 


colourway<-c("#E91E53", "#EF9A9A", "#E67E22", "#F1C40F", "#AA8A0A",
             "#388E3C", "#56bc5b", "#c5e1a5", "#00695c", "#00668f",
             "#F2F2F2", "#444444", "#838383","#B3B3B3", "#000000")

DATA_phylum_emu_top15 = names(sort(taxa_sums(EMU_data_phylum), TRUE)[1:15])
DATA_phylum_emu_top15 = prune_taxa(DATA_phylum_emu_top15, EMU_data_phylum)
DATA_phylum_emu_top15<-transform_sample_counts(DATA_phylum_emu_top15, function(x) 100 * x/sum(x))


emu_phylum = plot_bar(DATA_phylum_emu_top15, "Sample_original", fill="Phylum")+
  facet_grid( Database~Platform,scales = "free_x",space="free_x")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size= 15,color="black"), 
        axis.text.y =element_text(size= 15,color="black"), 
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 9), 
        legend.position = "right",
        axis.title.y = element_text(size= 15),
        axis.title.x = element_text(size= 15),
        strip.text.x = element_text(size = 15, color = "black"),
        strip.text.y = element_text(size = 15, color = "black"),
        strip.background =element_rect(color="black", fill="white", 
                                       size=1, linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
scale_fill_manual(values=colourway)+
  coord_cartesian(expand = F)

emu_phylum$data$Sample_original <- factor(emu_phylum$data$Sample_original, 
                                                           levels = c("AB1", "N42","ODK3b","ODK5b",
                                                                    "ODK7b", "ODK9a","ODK13b",
                                                                      "ODK15a", "ODK21b","ODK36c"))
emu_phylum 
#ggsave("Emu_phylum.pdf",scale = 1,dpi = 400)

DATA_class_emu_top15 = names(sort(taxa_sums(EMU_data_class), TRUE)[1:15])
DATA_class_emu_top15 = prune_taxa(DATA_class_emu_top15, EMU_data_class)
DATA_class_emu_top15<-transform_sample_counts(DATA_class_emu_top15, function(x) 100 * x/sum(x))


emu_class = plot_bar(DATA_class_emu_top15, "Sample_original", fill="Class")+
  facet_grid( Database~ Platform)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size= 15,color="black"), 
        axis.text.y =element_text(size= 15,color="black"), 
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 9), 
        legend.position = "right",
        axis.title.y = element_text(size= 15),
        axis.title.x = element_text(size= 15),
        strip.text.x = element_text(size = 15, color = "black"),
        strip.text.y = element_text(size = 15, color = "black"),
        strip.background =element_rect(color="black", fill="white", 
                                       size=1, linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_fill_manual(values=colourway)+
  coord_cartesian(expand = F)

emu_class$data$Sample_original <- factor(emu_class$data$Sample_original, 
                                          levels = c("AB1", "N42","ODK3b","ODK5b",
                                                     "ODK7b", "ODK9a","ODK13b",
                                                     "ODK15a", "ODK21b","ODK36c"))
emu_class

#ggsave("Emu_class.pdf",scale = 1,dpi = 400)

