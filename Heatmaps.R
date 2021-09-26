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

###Data set up####

##Mothur##

mothur_data<-subset_samples(data.2, Pipeline == "MOTHUR")
mothur_data_phylum<-tax_glom(mothur_data, taxrank="Phylum")
mothur_data_class<-tax_glom(mothur_data, taxrank="Class")
mothur_data_order<-tax_glom(mothur_data, taxrank="Order")
mothur_data_family<-tax_glom(mothur_data, taxrank="Family")
mothur_data_genus<-tax_glom(mothur_data, taxrank="Genus")
mothur_data_species<-tax_glom(mothur_data, taxrank="Species")

##Qiime##

Qiime_data<-subset_samples(data.2, Pipeline == "QIIME")
Qiime_data_phylum<-tax_glom(Qiime_data, taxrank="Phylum")
Qiime_data_class<-tax_glom(Qiime_data, taxrank="Class")
Qiime_data_order<-tax_glom(Qiime_data, taxrank="Order")
Qiime_data_family<-tax_glom(Qiime_data, taxrank="Family")
Qiime_data_genus<-tax_glom(Qiime_data, taxrank="Genus")
Qiime_data_species<-tax_glom(Qiime_data, taxrank="Species")


###Heatmaps_phylum#####

#Exclude V4

##Mothur##

#Phylum

mothur_data_phylum_noV4<-subset_samples(mothur_data_phylum, PR !="ONT_V4")
sample_data(mothur_data_phylum_noV4)
tax_table(mothur_data_phylum_noV4)
otu_table(mothur_data_phylum_noV4)

mothur_data_phylum_top20 = names(sort(taxa_sums(mothur_data_phylum_noV4), TRUE)[1:15])
mothur_data_phylum_top20 = prune_taxa(mothur_data_phylum_top20, mothur_data_phylum_noV4)
mothur_data_phylum_top20<-transform_sample_counts(mothur_data_phylum_top20, function(x) 100 * x/sum(x))

tax_m_phyla<-as.data.frame(tax_table(mothur_data_phylum_top20))
tax_m_phyla.1<-rownames_to_column(tax_m_phyla)
tax_m_phyla.2<-as.data.frame(tax_m_phyla.1[,c(1,3)])

otu_m_phyla<-as.data.frame(otu_table(mothur_data_phylum_top20))
otu_m_phyla.1<- rownames_to_column(otu_m_phyla)

otu_m_phyla.2<-left_join(tax_m_phyla.2,otu_m_phyla.1 , by ="rowname")
m_phyla<- otu_m_phyla.2[,-1]
m_phyla.1<-column_to_rownames(m_phyla,"Phylum")
m_phyla.1[m_phyla.1 == 0] <- NA   


m_phyla.2<-select(m_phyla.1,AB1_FL_ONT_M_RDP,  AB1_V4_IT_M_RDP,AB1_FL_ONT_M_S, AB1_V4_IT_M_S,
                  N42_FL_ONT_M_RDP,  N42_V4_IT_M_RDP,N42_FL_ONT_M_S, N42_V4_IT_M_S,
                  ODK3_FL_ONT_M_RDP,  ODK3_V4_IT_M_RDP, ODK3_FL_ONT_M_S,ODK3_V4_IT_M_S,
                  ODK5_FL_ONT_M_RDP, ODK5_V4_IT_M_RDP, ODK5_FL_ONT_M_S, ODK5_V4_IT_M_S,
                  ODK7_FL_ONT_M_RDP,  ODK7_V4_IT_M_RDP, ODK7_FL_ONT_M_S,ODK7_V4_IT_M_S,
                  ODK9_FL_ONT_M_RDP,  ODK9_V4_IT_M_RDP, ODK9_FL_ONT_M_S,ODK9_V4_IT_M_S,
                  ODK13_FL_ONT_M_RDP,  ODK13_V4_IT_M_RDP,ODK13_FL_ONT_M_S, ODK13_V4_IT_M_S,
                  ODK15_FL_ONT_M_RDP,  ODK15_V4_IT_M_RDP,ODK15_FL_ONT_M_S, ODK15_V4_IT_M_S,
                  ODK21_FL_ONT_M_RDP,  ODK21_V4_IT_M_RDP,ODK21_FL_ONT_M_S, ODK21_V4_IT_M_S,
                  ODK36_FL_ONT_M_RDP,  ODK36_V4_IT_M_RDP, ODK36_FL_ONT_M_S,ODK36_V4_IT_M_S,
                  DNAs1_FL_ONT_M_RDP,  DNAs1_V4_IT_M_RDP, DNAs1_FL_ONT_M_S,DNAs1_V4_IT_M_S,
                  DNAs2_FL_ONT_M_RDP,  DNAs2_V4_IT_M_RDP, DNAs2_FL_ONT_M_S,DNAs2_V4_IT_M_S)

map_m_phyla<-data.frame(sample_data(mothur_data_phylum_top20))

col_annotation_m<-dplyr::select(map_m_phyla,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_m_p<-pheatmap(m_phyla.2, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = col_annotation_m,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                   cellwidth = 12,
                   color=color(100),
                   na_col = "grey94",
                   annotation_colors = my_colour)


##Qiime##

#Phylum

Qiime_data_phylum_noV4<-subset_samples(Qiime_data_phylum, PR !="ONT_V4")
Qiime_data_phylum_top20 = names(sort(taxa_sums(Qiime_data_phylum_noV4), TRUE)[1:15])
Qiime_data_phylum_top20 = prune_taxa(Qiime_data_phylum_top20, Qiime_data_phylum_noV4)
Qiime_data_phylum_top20<-transform_sample_counts(Qiime_data_phylum_top20, function(x) 100 * x/sum(x))

tax_Q_phyla<-as.data.frame(tax_table(Qiime_data_phylum_top20))
tax_Q_phyla.1<-rownames_to_column(tax_Q_phyla)
tax_Q_phyla.2<-as.data.frame(tax_Q_phyla.1[,c(1,3)])

otu_Q_phyla<-as.data.frame(otu_table(Qiime_data_phylum_top20))
otu_Q_phyla.1<- rownames_to_column(otu_Q_phyla)

otu_Q_phyla.2<-left_join(tax_Q_phyla.2,otu_Q_phyla.1 , by ="rowname")
q_phyla<- otu_Q_phyla.2[,-1]
q_phyla.1<-column_to_rownames(q_phyla,"Phylum")
q_phyla.1[q_phyla.1 == 0] <- NA   


q_phyla.2<-select(q_phyla.1,AB1_FL_ONT_Q_RDP,  AB1_V4_IT_Q_RDP,AB1_FL_ONT_Q_S, AB1_V4_IT_Q_S,
                  N42_FL_ONT_Q_RDP,  N42_V4_IT_Q_RDP,N42_FL_ONT_Q_S, N42_V4_IT_Q_S,
                  ODK3_FL_ONT_Q_RDP,  ODK3_V4_IT_Q_RDP, ODK3_FL_ONT_Q_S,ODK3_V4_IT_Q_S,
                  ODK5_FL_ONT_Q_RDP, ODK5_V4_IT_Q_RDP, ODK5_FL_ONT_Q_S, ODK5_V4_IT_Q_S,
                  ODK7_FL_ONT_Q_RDP,  ODK7_V4_IT_Q_RDP, ODK7_FL_ONT_Q_S,ODK7_V4_IT_Q_S,
                  ODK9_FL_ONT_Q_RDP,  ODK9_V4_IT_Q_RDP, ODK9_FL_ONT_Q_S,ODK9_V4_IT_Q_S,
                  ODK13_FL_ONT_Q_RDP,  ODK13_V4_IT_Q_RDP,ODK13_FL_ONT_Q_S, ODK13_V4_IT_Q_S,
                  ODK15_FL_ONT_Q_RDP,  ODK15_V4_IT_Q_RDP,ODK15_FL_ONT_Q_S, ODK15_V4_IT_Q_S,
                  ODK21_FL_ONT_Q_RDP,  ODK21_V4_IT_Q_RDP,ODK21_FL_ONT_Q_S, ODK21_V4_IT_Q_S,
                  ODK36_FL_ONT_Q_RDP,  ODK36_V4_IT_Q_RDP, ODK36_FL_ONT_Q_S,ODK36_V4_IT_Q_S,
                  DNAs1_FL_ONT_Q_RDP,  DNAs1_V4_IT_Q_RDP, DNAs1_FL_ONT_Q_S,DNAs1_V4_IT_Q_S,
                  DNAs2_FL_ONT_Q_RDP,  DNAs2_V4_IT_Q_RDP, DNAs2_FL_ONT_Q_S,DNAs2_V4_IT_Q_S)

map_q_phyla<-data.frame(sample_data(Qiime_data_phylum_top20))

col_annotation_q<-dplyr::select(map_q_phyla,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_q_p<-pheatmap(q_phyla.2, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = col_annotation_q,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                   cellwidth = 12,
                   color=color(100),
                   na_col = "grey94",
                   annotation_colors = my_colour)

save_pheatmap_pdf <- function(x, filename, width=15, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save_pheatmap_pdf(heat_m_p, "heat_m_p.pdf")
#save_pheatmap_pdf(heat_q_p, "heat_q_p.pdf")


###Heatmaps_class#####

##Mothur##

#Class
mothur_data_class_noV4<-subset_samples(mothur_data_class, PR !="ONT_V4")
mothur_data_class_top20 = names(sort(taxa_sums(mothur_data_class_noV4), TRUE)[1:15])
mothur_data_class_top20 = prune_taxa(mothur_data_class_top20, mothur_data_class_noV4)
mothur_data_class_top20<-transform_sample_counts(mothur_data_class_top20, function(x) 100 * x/sum(x))

tax_m_class<-as.data.frame(tax_table(mothur_data_class_top20))
tax_m_class.1<-rownames_to_column(tax_m_class)
tax_m_class.2<-as.data.frame(tax_m_class.1[,c(1,4)])

otu_m_class<-as.data.frame(otu_table(mothur_data_class_top20))
otu_m_class.1<- rownames_to_column(otu_m_class)

otu_m_class.2<-left_join(tax_m_class.2,otu_m_class.1 , by ="rowname")
m_class<- otu_m_class.2[,-1]
m_class.1<-column_to_rownames(m_class,"Class")
m_class.1[m_class.1 == 0] <- NA   


m_class.2<-select(m_class.1,AB1_FL_ONT_M_RDP,  AB1_V4_IT_M_RDP,AB1_FL_ONT_M_S, AB1_V4_IT_M_S,
                  N42_FL_ONT_M_RDP,  N42_V4_IT_M_RDP,N42_FL_ONT_M_S, N42_V4_IT_M_S,
                  ODK3_FL_ONT_M_RDP,  ODK3_V4_IT_M_RDP, ODK3_FL_ONT_M_S,ODK3_V4_IT_M_S,
                  ODK5_FL_ONT_M_RDP, ODK5_V4_IT_M_RDP, ODK5_FL_ONT_M_S, ODK5_V4_IT_M_S,
                  ODK7_FL_ONT_M_RDP,  ODK7_V4_IT_M_RDP, ODK7_FL_ONT_M_S,ODK7_V4_IT_M_S,
                  ODK9_FL_ONT_M_RDP,  ODK9_V4_IT_M_RDP, ODK9_FL_ONT_M_S,ODK9_V4_IT_M_S,
                  ODK13_FL_ONT_M_RDP,  ODK13_V4_IT_M_RDP,ODK13_FL_ONT_M_S, ODK13_V4_IT_M_S,
                  ODK15_FL_ONT_M_RDP,  ODK15_V4_IT_M_RDP,ODK15_FL_ONT_M_S, ODK15_V4_IT_M_S,
                  ODK21_FL_ONT_M_RDP,  ODK21_V4_IT_M_RDP,ODK21_FL_ONT_M_S, ODK21_V4_IT_M_S,
                  ODK36_FL_ONT_M_RDP,  ODK36_V4_IT_M_RDP, ODK36_FL_ONT_M_S,ODK36_V4_IT_M_S,
                  DNAs1_FL_ONT_M_RDP,  DNAs1_V4_IT_M_RDP, DNAs1_FL_ONT_M_S,DNAs1_V4_IT_M_S,
                  DNAs2_FL_ONT_M_RDP,  DNAs2_V4_IT_M_RDP, DNAs2_FL_ONT_M_S,DNAs2_V4_IT_M_S)

map_m_class<-data.frame(sample_data(mothur_data_class_top20))

col_annotation_m<-dplyr::select(map_m_class,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_m_c<-pheatmap(m_class.2, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = col_annotation_m,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                   cellwidth = 12,
                   color=color(100),
                   na_col = "grey94",
                   annotation_colors = my_colour)


##Qiime##

#Class
Qiime_data_class_noV4<-subset_samples(Qiime_data_class, PR !="ONT_V4")
Qiime_data_class_top20 = names(sort(taxa_sums(Qiime_data_class_noV4), TRUE)[1:15])
Qiime_data_class_top20 = prune_taxa(Qiime_data_class_top20, Qiime_data_class_noV4)
Qiime_data_class_top20<-transform_sample_counts(Qiime_data_class_top20, function(x) 100 * x/sum(x))

tax_Q_class<-as.data.frame(tax_table(Qiime_data_class_top20))
tax_Q_class.1<-rownames_to_column(tax_Q_class)
tax_Q_class.2<-as.data.frame(tax_Q_class.1[,c(1,4)])

otu_Q_class<-as.data.frame(otu_table(Qiime_data_class_top20))
otu_Q_class.1<- rownames_to_column(otu_Q_class)

otu_Q_class.2<-left_join(tax_Q_class.2,otu_Q_class.1 , by ="rowname")
q_class<- otu_Q_class.2[,-1]
q_class.1<-column_to_rownames(q_class,"Class")
q_class.1[q_class.1 == 0] <- NA   


q_class.2<-select(q_class.1, AB1_FL_ONT_Q_RDP,  AB1_V4_IT_Q_RDP,AB1_FL_ONT_Q_S, AB1_V4_IT_Q_S,
                  N42_FL_ONT_Q_RDP,  N42_V4_IT_Q_RDP,N42_FL_ONT_Q_S, N42_V4_IT_Q_S,
                  ODK3_FL_ONT_Q_RDP,  ODK3_V4_IT_Q_RDP, ODK3_FL_ONT_Q_S,ODK3_V4_IT_Q_S,
                  ODK5_FL_ONT_Q_RDP, ODK5_V4_IT_Q_RDP, ODK5_FL_ONT_Q_S, ODK5_V4_IT_Q_S,
                  ODK7_FL_ONT_Q_RDP,  ODK7_V4_IT_Q_RDP, ODK7_FL_ONT_Q_S,ODK7_V4_IT_Q_S,
                  ODK9_FL_ONT_Q_RDP,  ODK9_V4_IT_Q_RDP, ODK9_FL_ONT_Q_S,ODK9_V4_IT_Q_S,
                  ODK13_FL_ONT_Q_RDP,  ODK13_V4_IT_Q_RDP,ODK13_FL_ONT_Q_S, ODK13_V4_IT_Q_S,
                  ODK15_FL_ONT_Q_RDP,  ODK15_V4_IT_Q_RDP,ODK15_FL_ONT_Q_S, ODK15_V4_IT_Q_S,
                  ODK21_FL_ONT_Q_RDP,  ODK21_V4_IT_Q_RDP,ODK21_FL_ONT_Q_S, ODK21_V4_IT_Q_S,
                  ODK36_FL_ONT_Q_RDP,  ODK36_V4_IT_Q_RDP, ODK36_FL_ONT_Q_S,ODK36_V4_IT_Q_S,
                  DNAs1_FL_ONT_Q_RDP,  DNAs1_V4_IT_Q_RDP, DNAs1_FL_ONT_Q_S,DNAs1_V4_IT_Q_S,
                  DNAs2_FL_ONT_Q_RDP,  DNAs2_V4_IT_Q_RDP, DNAs2_FL_ONT_Q_S,DNAs2_V4_IT_Q_S)

map_q_class<-data.frame(sample_data(Qiime_data_class_top20))

col_annotation_q<-dplyr::select(map_q_class,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_q_c<-pheatmap(q_class.2, cluster_row = FALSE, cluster_cols = FALSE,
                   annotation_col = col_annotation_q,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44, 46),
                   color=color(100),
                   cellwidth = 12,
                   na_col = "grey94",
                   #border_color = "White",
                   annotation_colors = my_colour)

save_pheatmap_pdf <- function(x, filename, width=15, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save_pheatmap_pdf(heat_m_c, "heat_m_c.pdf")
#save_pheatmap_pdf(heat_q_c, "heat_q_c.pdf")

               

###Heatmaps_order#####

##Mothur##

#order


mothur_data_order_noV4<-subset_samples(mothur_data_order, PR !="ONT_V4")
mothur_data_order_top20 = names(sort(taxa_sums(mothur_data_order_noV4), TRUE)[1:15])
mothur_data_order_top20 = prune_taxa(mothur_data_order_top20, mothur_data_order_noV4)
mothur_data_order_top20<-transform_sample_counts(mothur_data_order_top20, function(x) 100 * x/sum(x))

tax_m_order<-as.data.frame(tax_table(mothur_data_order_top20))
tax_m_order.1<-rownames_to_column(tax_m_order)
tax_m_order.2<-as.data.frame(tax_m_order.1[,c(1,5)])

otu_m_order<-as.data.frame(otu_table(mothur_data_order_top20))
otu_m_order.1<- rownames_to_column(otu_m_order)

otu_m_order.2<-left_join(tax_m_order.2,otu_m_order.1 , by ="rowname")
m_order<- otu_m_order.2[,-1]
m_order.1<-column_to_rownames(m_order,"Order")
m_order.1[m_order.1 == 0] <- NA   


m_order.2<-select(m_order.1,AB1_FL_ONT_M_RDP,  AB1_V4_IT_M_RDP,AB1_FL_ONT_M_S, AB1_V4_IT_M_S,
                  N42_FL_ONT_M_RDP,  N42_V4_IT_M_RDP,N42_FL_ONT_M_S, N42_V4_IT_M_S,
                  ODK3_FL_ONT_M_RDP,  ODK3_V4_IT_M_RDP, ODK3_FL_ONT_M_S,ODK3_V4_IT_M_S,
                  ODK5_FL_ONT_M_RDP, ODK5_V4_IT_M_RDP, ODK5_FL_ONT_M_S, ODK5_V4_IT_M_S,
                  ODK7_FL_ONT_M_RDP,  ODK7_V4_IT_M_RDP, ODK7_FL_ONT_M_S,ODK7_V4_IT_M_S,
                  ODK9_FL_ONT_M_RDP,  ODK9_V4_IT_M_RDP, ODK9_FL_ONT_M_S,ODK9_V4_IT_M_S,
                  ODK13_FL_ONT_M_RDP,  ODK13_V4_IT_M_RDP,ODK13_FL_ONT_M_S, ODK13_V4_IT_M_S,
                  ODK15_FL_ONT_M_RDP,  ODK15_V4_IT_M_RDP,ODK15_FL_ONT_M_S, ODK15_V4_IT_M_S,
                  ODK21_FL_ONT_M_RDP,  ODK21_V4_IT_M_RDP,ODK21_FL_ONT_M_S, ODK21_V4_IT_M_S,
                  ODK36_FL_ONT_M_RDP,  ODK36_V4_IT_M_RDP, ODK36_FL_ONT_M_S,ODK36_V4_IT_M_S,
                  DNAs1_FL_ONT_M_RDP,  DNAs1_V4_IT_M_RDP, DNAs1_FL_ONT_M_S,DNAs1_V4_IT_M_S,
                  DNAs2_FL_ONT_M_RDP,  DNAs2_V4_IT_M_RDP, DNAs2_FL_ONT_M_S,DNAs2_V4_IT_M_S)

map_m_order<-data.frame(sample_data(mothur_data_order_top20))

col_annotation_m<-dplyr::select(map_m_order,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_m_o<-pheatmap(m_order.2, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = col_annotation_m,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                   cellwidth = 12,
                   color=color(100),
                   na_col = "grey94",
                   annotation_colors = my_colour)


##Qiime##

#order
Qiime_data_order_noV4<-subset_samples(Qiime_data_order, PR !="ONT_V4")
Qiime_data_order_top20 = names(sort(taxa_sums(Qiime_data_order_noV4), TRUE)[1:15])
Qiime_data_order_top20 = prune_taxa(Qiime_data_order_top20, Qiime_data_order_noV4)
Qiime_data_order_top20<-transform_sample_counts(Qiime_data_order_top20, function(x) 100 * x/sum(x))

tax_Q_order<-as.data.frame(tax_table(Qiime_data_order_top20))
tax_Q_order.1<-rownames_to_column(tax_Q_order)
tax_Q_order.2<-as.data.frame(tax_Q_order.1[,c(1,5)])

otu_Q_order<-as.data.frame(otu_table(Qiime_data_order_top20))
otu_Q_order.1<- rownames_to_column(otu_Q_order)

otu_Q_order.2<-left_join(tax_Q_order.2,otu_Q_order.1 , by ="rowname")
q_order<- otu_Q_order.2[,-1]
q_order.1<-column_to_rownames(q_order,"Order")
q_order.1[q_order.1 == 0] <- NA   


q_order.2<-select(q_order.1,AB1_FL_ONT_Q_RDP,  AB1_V4_IT_Q_RDP,AB1_FL_ONT_Q_S, AB1_V4_IT_Q_S,
                  N42_FL_ONT_Q_RDP,  N42_V4_IT_Q_RDP,N42_FL_ONT_Q_S, N42_V4_IT_Q_S,
                  ODK3_FL_ONT_Q_RDP,  ODK3_V4_IT_Q_RDP, ODK3_FL_ONT_Q_S,ODK3_V4_IT_Q_S,
                  ODK5_FL_ONT_Q_RDP, ODK5_V4_IT_Q_RDP, ODK5_FL_ONT_Q_S, ODK5_V4_IT_Q_S,
                  ODK7_FL_ONT_Q_RDP,  ODK7_V4_IT_Q_RDP, ODK7_FL_ONT_Q_S,ODK7_V4_IT_Q_S,
                  ODK9_FL_ONT_Q_RDP,  ODK9_V4_IT_Q_RDP, ODK9_FL_ONT_Q_S,ODK9_V4_IT_Q_S,
                  ODK13_FL_ONT_Q_RDP,  ODK13_V4_IT_Q_RDP,ODK13_FL_ONT_Q_S, ODK13_V4_IT_Q_S,
                  ODK15_FL_ONT_Q_RDP,  ODK15_V4_IT_Q_RDP,ODK15_FL_ONT_Q_S, ODK15_V4_IT_Q_S,
                  ODK21_FL_ONT_Q_RDP,  ODK21_V4_IT_Q_RDP,ODK21_FL_ONT_Q_S, ODK21_V4_IT_Q_S,
                  ODK36_FL_ONT_Q_RDP,  ODK36_V4_IT_Q_RDP, ODK36_FL_ONT_Q_S,ODK36_V4_IT_Q_S,
                  DNAs1_FL_ONT_Q_RDP,  DNAs1_V4_IT_Q_RDP, DNAs1_FL_ONT_Q_S,DNAs1_V4_IT_Q_S,
                  DNAs2_FL_ONT_Q_RDP,  DNAs2_V4_IT_Q_RDP, DNAs2_FL_ONT_Q_S,DNAs2_V4_IT_Q_S)


map_q_order<-data.frame(sample_data(Qiime_data_order_top20))

col_annotation_q<-dplyr::select(map_q_order,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_q_o<-pheatmap(q_order.2, cluster_row = FALSE, cluster_cols = FALSE,
                   annotation_col = col_annotation_q,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44, 46),
                   color=color(100),
                   cellwidth = 12,
                   na_col = "grey94",
                   #border_color = "White",
                   annotation_colors = my_colour)

save_pheatmap_pdf <- function(x, filename, width=15, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save_pheatmap_pdf(heat_m_o, "heat_m_o.pdf")
#save_pheatmap_pdf(heat_q_o, "heat_q_o.pdf")

###Heatmaps_family#####

##Mothur##

#family
mothur_data_family_noV4<-subset_samples(mothur_data_family, PR !="ONT_V4")
mothur_data_family_top20 = names(sort(taxa_sums(mothur_data_family_noV4), TRUE)[1:15])
mothur_data_family_top20 = prune_taxa(mothur_data_family_top20, mothur_data_family_noV4)
mothur_data_family_top20<-transform_sample_counts(mothur_data_family_top20, function(x) 100 * x/sum(x))

tax_m_family<-as.data.frame(tax_table(mothur_data_family_top20))
tax_m_family.1<-rownames_to_column(tax_m_family)
tax_m_family.2<-as.data.frame(tax_m_family.1[,c(1,6)])

otu_m_family<-as.data.frame(otu_table(mothur_data_family_top20))
otu_m_family.1<- rownames_to_column(otu_m_family)

otu_m_family.2<-left_join(tax_m_family.2,otu_m_family.1 , by ="rowname")
m_family<- otu_m_family.2[,-1]
m_family.1<-column_to_rownames(m_family,"Family")
m_family.1[m_family.1 == 0] <- NA   


m_family.2<-select(m_family.1,AB1_FL_ONT_M_RDP,  AB1_V4_IT_M_RDP,AB1_FL_ONT_M_S, AB1_V4_IT_M_S,
                   N42_FL_ONT_M_RDP,  N42_V4_IT_M_RDP,N42_FL_ONT_M_S, N42_V4_IT_M_S,
                   ODK3_FL_ONT_M_RDP,  ODK3_V4_IT_M_RDP, ODK3_FL_ONT_M_S,ODK3_V4_IT_M_S,
                   ODK5_FL_ONT_M_RDP, ODK5_V4_IT_M_RDP, ODK5_FL_ONT_M_S, ODK5_V4_IT_M_S,
                   ODK7_FL_ONT_M_RDP,  ODK7_V4_IT_M_RDP, ODK7_FL_ONT_M_S,ODK7_V4_IT_M_S,
                   ODK9_FL_ONT_M_RDP,  ODK9_V4_IT_M_RDP, ODK9_FL_ONT_M_S,ODK9_V4_IT_M_S,
                   ODK13_FL_ONT_M_RDP,  ODK13_V4_IT_M_RDP,ODK13_FL_ONT_M_S, ODK13_V4_IT_M_S,
                   ODK15_FL_ONT_M_RDP,  ODK15_V4_IT_M_RDP,ODK15_FL_ONT_M_S, ODK15_V4_IT_M_S,
                   ODK21_FL_ONT_M_RDP,  ODK21_V4_IT_M_RDP,ODK21_FL_ONT_M_S, ODK21_V4_IT_M_S,
                   ODK36_FL_ONT_M_RDP,  ODK36_V4_IT_M_RDP, ODK36_FL_ONT_M_S,ODK36_V4_IT_M_S,
                   DNAs1_FL_ONT_M_RDP,  DNAs1_V4_IT_M_RDP, DNAs1_FL_ONT_M_S,DNAs1_V4_IT_M_S,
                   DNAs2_FL_ONT_M_RDP,  DNAs2_V4_IT_M_RDP, DNAs2_FL_ONT_M_S,DNAs2_V4_IT_M_S)

map_m_family<-data.frame(sample_data(mothur_data_family_top20))

col_annotation_m<-dplyr::select(map_m_family,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_m_f<-pheatmap(m_family.2, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = col_annotation_m,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                   cellwidth = 12,
                   color=color(100),
                   na_col = "grey94",
                   annotation_colors = my_colour)


##Qiime##

#family
Qiime_data_family_noV4<-subset_samples(Qiime_data_family, PR !="ONT_V4")
Qiime_data_family_top20 = names(sort(taxa_sums(Qiime_data_family_noV4), TRUE)[1:15])
Qiime_data_family_top20 = prune_taxa(Qiime_data_family_top20, Qiime_data_family_noV4)
Qiime_data_family_top20<-transform_sample_counts(Qiime_data_family_top20, function(x) 100 * x/sum(x))

tax_Q_family<-as.data.frame(tax_table(Qiime_data_family_top20))
tax_Q_family.1<-rownames_to_column(tax_Q_family)
tax_Q_family.2<-as.data.frame(tax_Q_family.1[,c(1,6)])

otu_Q_family<-as.data.frame(otu_table(Qiime_data_family_top20))
otu_Q_family.1<- rownames_to_column(otu_Q_family)

otu_Q_family.2<-left_join(tax_Q_family.2,otu_Q_family.1 , by ="rowname")
q_family<- otu_Q_family.2[,-1]
q_family.1<-column_to_rownames(q_family,"Family")
q_family.1[q_family.1 == 0] <- NA   


q_family.2<-select(q_family.1,AB1_FL_ONT_Q_RDP,  AB1_V4_IT_Q_RDP,AB1_FL_ONT_Q_S, AB1_V4_IT_Q_S,
                   N42_FL_ONT_Q_RDP,  N42_V4_IT_Q_RDP,N42_FL_ONT_Q_S, N42_V4_IT_Q_S,
                   ODK3_FL_ONT_Q_RDP,  ODK3_V4_IT_Q_RDP, ODK3_FL_ONT_Q_S,ODK3_V4_IT_Q_S,
                   ODK5_FL_ONT_Q_RDP, ODK5_V4_IT_Q_RDP, ODK5_FL_ONT_Q_S, ODK5_V4_IT_Q_S,
                   ODK7_FL_ONT_Q_RDP,  ODK7_V4_IT_Q_RDP, ODK7_FL_ONT_Q_S,ODK7_V4_IT_Q_S,
                   ODK9_FL_ONT_Q_RDP,  ODK9_V4_IT_Q_RDP, ODK9_FL_ONT_Q_S,ODK9_V4_IT_Q_S,
                   ODK13_FL_ONT_Q_RDP,  ODK13_V4_IT_Q_RDP,ODK13_FL_ONT_Q_S, ODK13_V4_IT_Q_S,
                   ODK15_FL_ONT_Q_RDP,  ODK15_V4_IT_Q_RDP,ODK15_FL_ONT_Q_S, ODK15_V4_IT_Q_S,
                   ODK21_FL_ONT_Q_RDP,  ODK21_V4_IT_Q_RDP,ODK21_FL_ONT_Q_S, ODK21_V4_IT_Q_S,
                   ODK36_FL_ONT_Q_RDP,  ODK36_V4_IT_Q_RDP, ODK36_FL_ONT_Q_S,ODK36_V4_IT_Q_S,
                   DNAs1_FL_ONT_Q_RDP,  DNAs1_V4_IT_Q_RDP, DNAs1_FL_ONT_Q_S,DNAs1_V4_IT_Q_S,
                   DNAs2_FL_ONT_Q_RDP,  DNAs2_V4_IT_Q_RDP, DNAs2_FL_ONT_Q_S,DNAs2_V4_IT_Q_S)

map_q_family<-data.frame(sample_data(Qiime_data_family_top20))

col_annotation_q<-dplyr::select(map_q_family,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_q_f<-pheatmap(q_family.2, cluster_row = FALSE, cluster_cols = FALSE,
                   annotation_col = col_annotation_q,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44, 46),
                   color=color(100),
                   cellwidth = 12,
                   na_col = "grey94",
                   #bfamily_color = "White",
                   annotation_colors = my_colour)

save_pheatmap_pdf <- function(x, filename, width=15, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save_pheatmap_pdf(heat_m_f, "heat_m_f.pdf")
#save_pheatmap_pdf(heat_q_f, "heat_q_f.pdf")

###Heatmaps_genus#####

##Mothur##

#genus
mothur_data_genus_noV4<-subset_samples(mothur_data_genus, PR !="ONT_V4")
mothur_data_genus_top20 = names(sort(taxa_sums(mothur_data_genus_noV4), TRUE)[1:15])
mothur_data_genus_top20 = prune_taxa(mothur_data_genus_top20, mothur_data_genus_noV4)
mothur_data_genus_top20<-transform_sample_counts(mothur_data_genus_top20, function(x) 100 * x/sum(x))

tax_m_genus<-as.data.frame(tax_table(mothur_data_genus_top20))
tax_m_genus.1<-rownames_to_column(tax_m_genus)
tax_m_genus.2<-as.data.frame(tax_m_genus.1[,c(1,7)])

otu_m_genus<-as.data.frame(otu_table(mothur_data_genus_top20))
otu_m_genus.1<- rownames_to_column(otu_m_genus)

otu_m_genus.2<-left_join(tax_m_genus.2,otu_m_genus.1 , by ="rowname")
m_genus<- otu_m_genus.2[,-1]
m_genus.1<-column_to_rownames(m_genus,"Genus")
m_genus.1[m_genus.1 == 0] <- NA   


m_genus.2<-select(m_genus.1,AB1_FL_ONT_M_RDP,  AB1_V4_IT_M_RDP,AB1_FL_ONT_M_S, AB1_V4_IT_M_S,
                  N42_FL_ONT_M_RDP,  N42_V4_IT_M_RDP,N42_FL_ONT_M_S, N42_V4_IT_M_S,
                  ODK3_FL_ONT_M_RDP,  ODK3_V4_IT_M_RDP, ODK3_FL_ONT_M_S,ODK3_V4_IT_M_S,
                  ODK5_FL_ONT_M_RDP, ODK5_V4_IT_M_RDP, ODK5_FL_ONT_M_S, ODK5_V4_IT_M_S,
                  ODK7_FL_ONT_M_RDP,  ODK7_V4_IT_M_RDP, ODK7_FL_ONT_M_S,ODK7_V4_IT_M_S,
                  ODK9_FL_ONT_M_RDP,  ODK9_V4_IT_M_RDP, ODK9_FL_ONT_M_S,ODK9_V4_IT_M_S,
                  ODK13_FL_ONT_M_RDP,  ODK13_V4_IT_M_RDP,ODK13_FL_ONT_M_S, ODK13_V4_IT_M_S,
                  ODK15_FL_ONT_M_RDP,  ODK15_V4_IT_M_RDP,ODK15_FL_ONT_M_S, ODK15_V4_IT_M_S,
                  ODK21_FL_ONT_M_RDP,  ODK21_V4_IT_M_RDP,ODK21_FL_ONT_M_S, ODK21_V4_IT_M_S,
                  ODK36_FL_ONT_M_RDP,  ODK36_V4_IT_M_RDP, ODK36_FL_ONT_M_S,ODK36_V4_IT_M_S,
                  DNAs1_FL_ONT_M_RDP,  DNAs1_V4_IT_M_RDP, DNAs1_FL_ONT_M_S,DNAs1_V4_IT_M_S,
                  DNAs2_FL_ONT_M_RDP,  DNAs2_V4_IT_M_RDP, DNAs2_FL_ONT_M_S,DNAs2_V4_IT_M_S)

map_m_genus<-data.frame(sample_data(mothur_data_genus_top20))

col_annotation_m<-dplyr::select(map_m_genus,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_m_g<-pheatmap(m_genus.2, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = col_annotation_m,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                   cellwidth = 12,
                   color=color(100),
                   na_col = "grey94",
                   annotation_colors = my_colour)


##Qiime##

#genus
Qiime_data_genus_noV4<-subset_samples(Qiime_data_genus, PR !="ONT_V4")
Qiime_data_genus_top20 = names(sort(taxa_sums(Qiime_data_genus_noV4), TRUE)[1:15])
Qiime_data_genus_top20 = prune_taxa(Qiime_data_genus_top20, Qiime_data_genus_noV4)
Qiime_data_genus_top20<-transform_sample_counts(Qiime_data_genus_top20, function(x) 100 * x/sum(x))

tax_Q_genus<-as.data.frame(tax_table(Qiime_data_genus_top20))
tax_Q_genus.1<-rownames_to_column(tax_Q_genus)
tax_Q_genus.2<-as.data.frame(tax_Q_genus.1[,c(1,7)])

otu_Q_genus<-as.data.frame(otu_table(Qiime_data_genus_top20))
otu_Q_genus.1<- rownames_to_column(otu_Q_genus)

otu_Q_genus.2<-left_join(tax_Q_genus.2,otu_Q_genus.1 , by ="rowname")
q_genus<- otu_Q_genus.2[,-1]
q_genus.1<-column_to_rownames(q_genus,"Genus")
q_genus.1[q_genus.1 == 0] <- NA   


q_genus.2<-select(q_genus.1, AB1_FL_ONT_Q_RDP,  AB1_V4_IT_Q_RDP,AB1_FL_ONT_Q_S, AB1_V4_IT_Q_S,
                  N42_FL_ONT_Q_RDP,  N42_V4_IT_Q_RDP,N42_FL_ONT_Q_S, N42_V4_IT_Q_S,
                  ODK3_FL_ONT_Q_RDP,  ODK3_V4_IT_Q_RDP, ODK3_FL_ONT_Q_S,ODK3_V4_IT_Q_S,
                  ODK5_FL_ONT_Q_RDP, ODK5_V4_IT_Q_RDP, ODK5_FL_ONT_Q_S, ODK5_V4_IT_Q_S,
                  ODK7_FL_ONT_Q_RDP,  ODK7_V4_IT_Q_RDP, ODK7_FL_ONT_Q_S,ODK7_V4_IT_Q_S,
                  ODK9_FL_ONT_Q_RDP,  ODK9_V4_IT_Q_RDP, ODK9_FL_ONT_Q_S,ODK9_V4_IT_Q_S,
                  ODK13_FL_ONT_Q_RDP,  ODK13_V4_IT_Q_RDP,ODK13_FL_ONT_Q_S, ODK13_V4_IT_Q_S,
                  ODK15_FL_ONT_Q_RDP,  ODK15_V4_IT_Q_RDP,ODK15_FL_ONT_Q_S, ODK15_V4_IT_Q_S,
                  ODK21_FL_ONT_Q_RDP,  ODK21_V4_IT_Q_RDP,ODK21_FL_ONT_Q_S, ODK21_V4_IT_Q_S,
                  ODK36_FL_ONT_Q_RDP,  ODK36_V4_IT_Q_RDP, ODK36_FL_ONT_Q_S,ODK36_V4_IT_Q_S,
                  DNAs1_FL_ONT_Q_RDP,  DNAs1_V4_IT_Q_RDP, DNAs1_FL_ONT_Q_S,DNAs1_V4_IT_Q_S,
                  DNAs2_FL_ONT_Q_RDP,  DNAs2_V4_IT_Q_RDP, DNAs2_FL_ONT_Q_S,DNAs2_V4_IT_Q_S)

map_q_genus<-data.frame(sample_data(Qiime_data_genus_top20))

col_annotation_q<-dplyr::select(map_q_genus,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_q_g<-pheatmap(q_genus.2, cluster_row = FALSE, cluster_cols = FALSE,
                   annotation_col = col_annotation_q,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44, 46),
                   color=color(100),
                   cellwidth = 12,
                   na_col = "grey94",
                   #bgenus_color = "White",
                   annotation_colors = my_colour)

save_pheatmap_pdf <- function(x, filename, width=15, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save_pheatmap_pdf(heat_m_g, "heat_m_g.pdf")
#save_pheatmap_pdf(heat_q_g, "heat_q_g.pdf")

###Heatmaps_species#####

##Mothur##

#species
mothur_data_species_noV4<-subset_samples(mothur_data_species, PR !="ONT_V4")
mothur_data_species_top20 = names(sort(taxa_sums(mothur_data_species_noV4), TRUE)[1:15])
mothur_data_species_top20 = prune_taxa(mothur_data_species_top20, mothur_data_species_noV4)
mothur_data_species_top20<-transform_sample_counts(mothur_data_species_top20, function(x) 100 * x/sum(x))

tax_m_species<-as.data.frame(tax_table(mothur_data_species_top20))
tax_m_species.1<-rownames_to_column(tax_m_species)
tax_m_species.2<-as.data.frame(tax_m_species.1[,c(1,8)])

otu_m_species<-as.data.frame(otu_table(mothur_data_species_top20))
otu_m_species.1<- rownames_to_column(otu_m_species)

otu_m_species.2<-left_join(tax_m_species.2,otu_m_species.1 , by ="rowname")
m_species<- otu_m_species.2[,-1]
m_species.1<-column_to_rownames(m_species,"Species")
m_species.1[m_species.1 == 0] <- NA   


m_species.2<-select(m_species.1,AB1_FL_ONT_M_RDP,  AB1_V4_IT_M_RDP,AB1_FL_ONT_M_S, AB1_V4_IT_M_S,
                    N42_FL_ONT_M_RDP,  N42_V4_IT_M_RDP,N42_FL_ONT_M_S, N42_V4_IT_M_S,
                    ODK3_FL_ONT_M_RDP,  ODK3_V4_IT_M_RDP, ODK3_FL_ONT_M_S,ODK3_V4_IT_M_S,
                    ODK5_FL_ONT_M_RDP, ODK5_V4_IT_M_RDP, ODK5_FL_ONT_M_S, ODK5_V4_IT_M_S,
                    ODK7_FL_ONT_M_RDP,  ODK7_V4_IT_M_RDP, ODK7_FL_ONT_M_S,ODK7_V4_IT_M_S,
                    ODK9_FL_ONT_M_RDP,  ODK9_V4_IT_M_RDP, ODK9_FL_ONT_M_S,ODK9_V4_IT_M_S,
                    ODK13_FL_ONT_M_RDP,  ODK13_V4_IT_M_RDP,ODK13_FL_ONT_M_S, ODK13_V4_IT_M_S,
                    ODK15_FL_ONT_M_RDP,  ODK15_V4_IT_M_RDP,ODK15_FL_ONT_M_S, ODK15_V4_IT_M_S,
                    ODK21_FL_ONT_M_RDP,  ODK21_V4_IT_M_RDP,ODK21_FL_ONT_M_S, ODK21_V4_IT_M_S,
                    ODK36_FL_ONT_M_RDP,  ODK36_V4_IT_M_RDP, ODK36_FL_ONT_M_S,ODK36_V4_IT_M_S,
                    DNAs1_FL_ONT_M_RDP,  DNAs1_V4_IT_M_RDP, DNAs1_FL_ONT_M_S,DNAs1_V4_IT_M_S,
                    DNAs2_FL_ONT_M_RDP,  DNAs2_V4_IT_M_RDP, DNAs2_FL_ONT_M_S,DNAs2_V4_IT_M_S)

map_m_species<-data.frame(sample_data(mothur_data_species_top20))

col_annotation_m<-dplyr::select(map_m_species,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_m_s<-pheatmap(m_species.2, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE,
                   annotation_col = col_annotation_m,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                   cellwidth = 12,
                   color=color(100),
                   na_col = "grey94",
                   annotation_colors = my_colour)


##Qiime##

#species
Qiime_data_species_noV4<-subset_samples(Qiime_data_species, PR !="ONT_V4")
Qiime_data_species_top20 = names(sort(taxa_sums(Qiime_data_species_noV4), TRUE)[1:15])
Qiime_data_species_top20 = prune_taxa(Qiime_data_species_top20, Qiime_data_species_noV4)
Qiime_data_species_top20<-transform_sample_counts(Qiime_data_species_top20, function(x) 100 * x/sum(x))

tax_Q_species<-as.data.frame(tax_table(Qiime_data_species_top20))
tax_Q_species.1<-rownames_to_column(tax_Q_species)
tax_Q_species.2<-as.data.frame(tax_Q_species.1[,c(1,8)])

otu_Q_species<-as.data.frame(otu_table(Qiime_data_species_top20))
otu_Q_species.1<- rownames_to_column(otu_Q_species)

otu_Q_species.2<-left_join(tax_Q_species.2,otu_Q_species.1 , by ="rowname")
q_species<- otu_Q_species.2[,-1]
q_species.1<-column_to_rownames(q_species,"Species")
q_species.1[q_species.1 == 0] <- NA   


q_species.2<-select(q_species.1,AB1_FL_ONT_Q_RDP,  AB1_V4_IT_Q_RDP,AB1_FL_ONT_Q_S, AB1_V4_IT_Q_S,
                    N42_FL_ONT_Q_RDP,  N42_V4_IT_Q_RDP,N42_FL_ONT_Q_S, N42_V4_IT_Q_S,
                    ODK3_FL_ONT_Q_RDP,  ODK3_V4_IT_Q_RDP, ODK3_FL_ONT_Q_S,ODK3_V4_IT_Q_S,
                    ODK5_FL_ONT_Q_RDP, ODK5_V4_IT_Q_RDP, ODK5_FL_ONT_Q_S, ODK5_V4_IT_Q_S,
                    ODK7_FL_ONT_Q_RDP,  ODK7_V4_IT_Q_RDP, ODK7_FL_ONT_Q_S,ODK7_V4_IT_Q_S,
                    ODK9_FL_ONT_Q_RDP,  ODK9_V4_IT_Q_RDP, ODK9_FL_ONT_Q_S,ODK9_V4_IT_Q_S,
                    ODK13_FL_ONT_Q_RDP,  ODK13_V4_IT_Q_RDP,ODK13_FL_ONT_Q_S, ODK13_V4_IT_Q_S,
                    ODK15_FL_ONT_Q_RDP,  ODK15_V4_IT_Q_RDP,ODK15_FL_ONT_Q_S, ODK15_V4_IT_Q_S,
                    ODK21_FL_ONT_Q_RDP,  ODK21_V4_IT_Q_RDP,ODK21_FL_ONT_Q_S, ODK21_V4_IT_Q_S,
                    ODK36_FL_ONT_Q_RDP,  ODK36_V4_IT_Q_RDP, ODK36_FL_ONT_Q_S,ODK36_V4_IT_Q_S,
                    DNAs1_FL_ONT_Q_RDP,  DNAs1_V4_IT_Q_RDP, DNAs1_FL_ONT_Q_S,DNAs1_V4_IT_Q_S,
                    DNAs2_FL_ONT_Q_RDP,  DNAs2_V4_IT_Q_RDP, DNAs2_FL_ONT_Q_S,DNAs2_V4_IT_Q_S)

map_q_species<-data.frame(sample_data(Qiime_data_species_top20))

col_annotation_q<-dplyr::select(map_q_species,  Database,Sample_original)

my_colour = list(Database = c(RDP = "#22499C", SILVA= "#BF8AB2"))

color= colorRampPalette(c("#003c0d","#9DD15A", "#fdff9b", "#f2bc7c", "#6e081d"))

heat_q_s<-pheatmap(q_species.2, cluster_row = FALSE, cluster_cols = FALSE,
                   annotation_col = col_annotation_q,
                   gaps_col = c(2,4,6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44, 46),
                   color=color(100),
                   cellwidth = 12,
                   na_col = "grey94",
                   #bspecies_color = "White",
                   annotation_colors = my_colour)

save_pheatmap_pdf <- function(x, filename, width=15, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save_pheatmap_pdf(heat_m_s, "heat_m_s.pdf")
#save_pheatmap_pdf(heat_q_s, "heat_q_s.pdf")
