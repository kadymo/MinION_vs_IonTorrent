setwd("~/Documents/MinION_Seq_data/Methods/MinION_qiime_analysis/qiime/V4_out/All")

library(dplyr)
library(plyr)
library(ggplot2)
library(tidyr)

source("~/Documents/MinION_Seq_data/Methods/MinION_qiime_analysis/qiime/V4_out/All/data_importer.R")
setwd("~/Documents/MinION_Seq_data/Methods/MinION_qiime_analysis/qiime/V4_out/All/SILVA_files")
#SILVA
AB1_S <- data_importer("AB1.txt", "AB1.taxonomy.tsv")
DNAs1_S <- data_importer("DNAs1.txt", "DNAs1.taxonomy.tsv")
DNAs2_S <- data_importer("DNAs2.txt", "DNAs2.taxonomy.tsv")
N42_S <- data_importer("N42.txt", "N42.taxonomy.tsv")
ODK3b_S <- data_importer("ODK3b.txt", "ODK3b.taxonomy.tsv")
ODK5b_S <- data_importer("ODK5b.txt", "ODK5b.taxonomy.tsv")
ODK7b_S <- data_importer("ODK7b.txt", "ODK7b.taxonomy.tsv")
ODK9a_S <- data_importer("ODK9a.txt", "ODK9a.taxonomy.tsv")
ODK11b_S <- data_importer("ODK11b.txt", "ODK11b.taxonomy.tsv")
ODK13b_S <- data_importer("ODK13b.txt", "ODK13b.taxonomy.tsv")
ODK15a_S <- data_importer("ODK15a.txt", "ODK15a.taxonomy.tsv")
ODK21b_S <- data_importer("ODK21b.txt", "ODK21b.taxonomy.tsv")
ODK36a_S <- data_importer("ODK36a.txt", "ODK36a.taxonomy.tsv")

setwd("~/Documents/MinION_Seq_data/Methods/MinION_qiime_analysis/qiime/V4_out/All/RDP_files/")
#RDP
AB1_R <- data_importer("AB1_RDP.txt", "AB1_RDP.taxonomy.tsv")
DNAs1_R <- data_importer("DNAs1_RDP.txt", "DNAs1_RDP.taxonomy.tsv")
DNAs2_R <- data_importer("DNAs2_RDP.txt", "DNAs2_RDP.taxonomy.tsv")
N42_R <- data_importer("N42_RDP.txt", "N42_RDP.taxonomy.tsv")
ODK3b_R <- data_importer("ODK3b_RDP.txt", "ODK3b_RDP.taxonomy.tsv")
ODK5b_R <- data_importer("ODK5b_RDP.txt", "ODK5b_RDP.taxonomy.tsv")
ODK7b_R <- data_importer("ODK7b_RDP.txt", "ODK7b_RDP.taxonomy.tsv")
ODK9a_R <- data_importer("ODK9a_RDP.txt", "ODK9a_RDP.taxonomy.tsv")
ODK11b_R <- data_importer("ODK11b_RDP.txt", "ODK11b_RDP.taxonomy.tsv")
ODK13b_R <- data_importer("ODK13b_RDP.txt", "ODK13b_RDP.taxonomy.tsv")
ODK15a_R <- data_importer("ODK15a_RDP.txt", "ODK15a_RDP.taxonomy.tsv")
ODK21b_R <- data_importer("ODK21b_RDP.txt", "ODK21b_RDP.taxonomy.tsv")
ODK36a_R <- data_importer("ODK36a_RDP.txt", "ODK36a_RDP.taxonomy.tsv")


all_r <- AB1_R %>%
  left_join(DNAs1_R, by='taxonomy') %>%
  left_join(DNAs2_R, by='taxonomy') %>%
  left_join(N42_R, by='taxonomy') %>%
  left_join(ODK3b_R, by='taxonomy') %>%
  left_join(ODK5b_R, by='taxonomy') %>%
  left_join(ODK7b_R, by='taxonomy') %>%
  left_join(ODK9a_R, by='taxonomy') %>%
  left_join(ODK11b_R, by='taxonomy') %>%
  left_join(ODK13b_R, by='taxonomy') %>%
  left_join(ODK15a_R, by='taxonomy') %>%
  left_join(ODK21b_R, by='taxonomy') %>%
  left_join(ODK36a_R, by='taxonomy') 

RDP_all <- ddply(all_r,"taxonomy",numcolwise(sum))
RDP_all[is.na(RDP_all)] <- 0
RDP_all  <- separate(RDP_all, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";", fill="right")


all_s <- AB1_S %>%
  left_join(DNAs1_S, by='taxonomy') %>%
  left_join(DNAs2_S, by='taxonomy') %>%
  left_join(N42_S, by='taxonomy') %>%
  left_join(ODK3b_S, by='taxonomy') %>%
  left_join(ODK5b_S, by='taxonomy') %>%
  left_join(ODK7b_S, by='taxonomy') %>%
  left_join(ODK9a_S, by='taxonomy') %>%
  left_join(ODK11b_S, by='taxonomy') %>%
  left_join(ODK13b_S, by='taxonomy') %>%
  left_join(ODK15a_S, by='taxonomy') %>%
  left_join(ODK21b_S, by='taxonomy') %>%
  left_join(ODK36a_S, by='taxonomy')

SILVA_all <- ddply(all_s,"taxonomy",numcolwise(sum))
SILVA_all[is.na(SILVA_all)] <- 0
SILVA_all  <- separate(SILVA_all, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";", fill="right")

