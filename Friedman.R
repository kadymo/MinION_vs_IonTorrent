library(tidyverse)
library(ggpubr)
library(rstatix)
library(datarium)
library(phyloseq)

setwd("~/Desktop/UCT projects/Oudekraal/Methods paper/R")
load("~/Desktop/UCT projects/Oudekraal/Methods paper/R/DESeq.RData")

#Data prep:

#Grab Phyla:
gp.Bact = subset_taxa(Data.1.3, Phylum == "Bacteroidetes")

sample_data(gp.Bact)
library(DAtest)
#try DA.fri with phyloseq object:
?DA.fri
Fri_Bacteroidetes <- DA.fri(gp.Bact, "P_DB_P", "Sample_original", relative=TRUE)
#posthoc:
Wil.Bact <- DA.wil(gp.Bact, "P_DB_P", "Sample_original", relative=TRUE, p.adj="bonferroni")
?p.adjust
