library(ggplot2)
library(phyloseq)
library("tidyverse") 
library(plyr)
library(vegan)
#install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
#devtools::install_github('schuyler-smith/phylosmith')
library("phylosmith")


setwd("~/Dropbox/MinION/Methods paper/R_analysis/Data_files/Phyloseq_objects")

###MinION upload####

#Mothur_Silva
TAXfile="ONT_S_m_tax.csv"
OTUfile="ONT_S_m_otu.csv"

otu_ont_m_s= read.table(OTUfile, header=TRUE, sep=";")
otu_ont_m_s=column_to_rownames(otu_ont_m_s, "OTU")  
tax_ont_m_s= read.table(TAXfile, header=TRUE, sep=";")
tax_ont_m_s=column_to_rownames(tax_ont_m_s, "OTU")  
tax_ont_m_s<- as.matrix(tax_ont_m_s)

OTU_ont_m_s = otu_table(otu_ont_m_s, taxa_are_rows = TRUE)
TAX_ont_m_s= tax_table(tax_ont_m_s)

ONT_M_S<- merge_phyloseq(OTU_ont_m_s, TAX_ont_m_s)
ONT_M_S


#Mothur_RDP
TAXfile="ONT_R_m_tax.csv"
OTUfile="ONT_R_m_otu.csv"

otu_ont_m_r= read.table(OTUfile, header=TRUE, sep=";")
otu_ont_m_r=column_to_rownames(otu_ont_m_r, "OTU")  
tax_ont_m_r= read.table(TAXfile, header=TRUE, sep=";")
tax_ont_m_r=column_to_rownames(tax_ont_m_r, "OTU")  
tax_ont_m_r<- as.matrix(tax_ont_m_r)

OTU_ont_m_r = otu_table(otu_ont_m_r, taxa_are_rows = TRUE)
TAX_ont_m_r= tax_table(tax_ont_m_r)

ONT_m_r<- merge_phyloseq(OTU_ont_m_r, TAX_ont_m_r)
ONT_m_r

#QIMME_Silva
TAXfile="ONT_S_q_tax.csv"
OTUfile="ONT_S_q_otu.csv"

otu_ont_q_s= read.table(OTUfile, header=TRUE, sep=";")
otu_ont_q_s=column_to_rownames(otu_ont_q_s, "OTU")  
tax_ont_q_s= read.table(TAXfile, header=TRUE, sep=";")
tax_ont_q_s=column_to_rownames(tax_ont_q_s, "OTU")  
tax_ont_q_s<- as.matrix(tax_ont_q_s)

OTU_ont_q_s = otu_table(otu_ont_q_s, taxa_are_rows = TRUE)
TAX_ont_q_s= tax_table(tax_ont_q_s)

ONT_q_s<- merge_phyloseq(OTU_ont_q_s, TAX_ont_q_s)
ONT_q_s

#QIMME_RDP
TAXfile="ONT_R_q_tax.csv"
OTUfile="ONT_R_q_otu.csv"

otu_ont_q_r= read.table(OTUfile, header=TRUE, sep=";")
otu_ont_q_r=column_to_rownames(otu_ont_q_r, "OTU")  
tax_ont_q_r= read.table(TAXfile, header=TRUE, sep=";")
tax_ont_q_r=column_to_rownames(tax_ont_q_r, "OTU")  
tax_ont_q_r<- as.matrix(tax_ont_q_r)

OTU_ont_q_r = otu_table(otu_ont_q_r, taxa_are_rows = TRUE)
TAX_ont_q_r= tax_table(tax_ont_q_r)

ONT_q_r<- merge_phyloseq(OTU_ont_q_r, TAX_ont_q_r)
ONT_q_r

###Ion-Torrent####

#Mothur_Silva
TAXfile="IT_S_m_tax.csv"
OTUfile="IT_S_m_otu.csv"

otu_IT_m_s= read.table(OTUfile, header=TRUE, sep=";")
otu_IT_m_s=column_to_rownames(otu_IT_m_s, "OTU")  
tax_IT_m_s= read.table(TAXfile, header=TRUE, sep=";")
tax_IT_m_s=column_to_rownames(tax_IT_m_s, "OTU")  
tax_IT_m_s<- as.matrix(tax_IT_m_s)

OTU_IT_m_s = otu_table(otu_IT_m_s, taxa_are_rows = TRUE)
TAX_IT_m_s= tax_table(tax_IT_m_s)

IT_M_S<- merge_phyloseq(OTU_IT_m_s, TAX_IT_m_s)
IT_M_S

#Mothur_RDP
TAXfile="IT_R_m_tax.csv"
OTUfile="IT_R_m_otu.csv"

otu_IT_m_r= read.table(OTUfile, header=TRUE, sep=";")
otu_IT_m_r=column_to_rownames(otu_IT_m_r, "OTU")  
tax_IT_m_r= read.table(TAXfile, header=TRUE, sep=";")
tax_IT_m_r=column_to_rownames(tax_IT_m_r, "OTU")  
tax_IT_m_r<- as.matrix(tax_IT_m_r)

OTU_IT_m_r = otu_table(otu_IT_m_r, taxa_are_rows = TRUE)
TAX_IT_m_r= tax_table(tax_IT_m_r)

IT_m_r<- merge_phyloseq(OTU_IT_m_r, TAX_IT_m_r)
IT_m_r

#QIMME_Silva
TAXfile="IT_S_q_tax.csv"
OTUfile="IT_S_q_otus.csv"

otu_IT_q_s= read.table(OTUfile, header=TRUE, sep=";")
otu_IT_q_s=column_to_rownames(otu_IT_q_s, "OTU")  
tax_IT_q_s= read.table(TAXfile, header=TRUE, sep=";")
tax_IT_q_s=column_to_rownames(tax_IT_q_s, "OTU")  
tax_IT_q_s<- as.matrix(tax_IT_q_s)

OTU_IT_q_s = otu_table(otu_IT_q_s, taxa_are_rows = TRUE)
TAX_IT_q_s= tax_table(tax_IT_q_s)

IT_q_s<- merge_phyloseq(OTU_IT_q_s, TAX_IT_q_s)
IT_q_s

#QIMME_RDP
TAXfile="IT_R_q_tax.csv"
OTUfile="IT_R_q_otus.csv"

otu_IT_q_r= read.table(OTUfile, header=TRUE, sep=";")
otu_IT_q_r=column_to_rownames(otu_IT_q_r, "OTU")  
tax_IT_q_r= read.table(TAXfile, header=TRUE, sep=";")
tax_IT_q_r=column_to_rownames(tax_IT_q_r, "OTU")  
tax_IT_q_r<- as.matrix(tax_IT_q_r)

OTU_IT_q_r = otu_table(otu_IT_q_r, taxa_are_rows = TRUE)
TAX_IT_q_r= tax_table(tax_IT_q_r)

IT_q_r<- merge_phyloseq(OTU_IT_q_r, TAX_IT_q_r)
IT_q_r

###Metadata####
mapfile = "Sample_data_all.csv"

map <- read.csv2(mapfile)
rownames(map) <- map$Sample_name
map<-map[,-1]
head(map)
MAP <- sample_data(map)

###Merge Phyloseq####

g<-merge_phyloseq(ONT_m_r,ONT_M_S,ONT_q_r, ONT_q_s, IT_m_r, IT_M_S, IT_q_r,IT_q_s)
g

tax_g<-as.data.frame(tax_table(g))
otu_g<- as.data.frame(otu_table(g))

#write_csv2(otu_g, "otu_g.csv")
#write_csv2(tax_g, "tax_g.csv")

g1<-merge_phyloseq(g, MAP)

tax_g1<-as.data.frame(tax_table(g1))
tax_g1<-rownames_to_column(tax_g1)
otu_g1<- as.data.frame(otu_table(g1))

#write_csv2(otu_g1, "Phyloseq_final_otu.csv")
#write_csv2(tax_g1, "Phyloseq_final_tax.csv")


#rarecurve<-rarecurve(t(otu_table(g1)))
                     

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


##FILTER
data.1 <- subset_taxa(data, Kingdom == "Bacteria" & Class != "Chloroplast" & Order != "Chloroplast" & Family != "mitochondria" &  Phylum != "Cyanobacteria/Chloroplast")


##add standard separation column
variable1 = sample_data(data.1)$Location
sample_data(data.1)$Type <-  mapply(paste0, variable1)
sample_data(data.1)$Type <- as.character(sample_data(data.1)$Type)
sample_data(data.1)$Type[sample_data(data.1)$Type == "Oudekraal"] <- "Sample"
sample_data(data.1)$Type[sample_data(data.1)$Type == "IEP"] <- "Sample"
sample_data(data.1)$Type[sample_data(data.1)$Type == "Muizenberg"] <- "Sample"
sample_data(data.1)$Type[sample_data(data.1)$Type == "None"] <- "Standard"

