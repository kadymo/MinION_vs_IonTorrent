library(ggplot2)
library(phyloseq)
library("tidyverse") 
library(plyr)
library(vegan)
#install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
#devtools::install_github('schuyler-smith/phylosmith')
library("phylosmith")


setwd("~/Dropbox/MinION/Methods paper/R_analysis/Data_files/Phyloseq_objects")

###Mothur sum taxa####
ont_r_m_combine<-read.csv2("ONT_R_m_combine.csv")
ont_r_m_combine.1<- ddply(ont_r_m_combine,"Taxonomy",numcolwise(sum))
ont_r_m_combine.2  <- separate(ont_r_m_combine.1, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=",", fill="right")
#write.csv2(ont_r_m_combine.2, "ont_r_m_combine.2.csv")

IT_S_m_combine<-read.csv2("IT_S_m_combine.csv")
IT_S_m_combine.1<- ddply(IT_S_m_combine,"Taxonomy",numcolwise(sum))
IT_S_m_combine.2  <- separate(IT_S_m_combine.1, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";", fill="right")
#write.csv2(IT_S_m_combine.2, "IT_S_m_combine.2.csv")


ONT_S_m_combine<-read.csv2("ONT_S_m_combine.csv")
ONT_S_m_combine.1<- ddply(ONT_S_m_combine,"Taxonomy",numcolwise(sum))
ONT_S_m_combine.2  <- separate(ONT_S_m_combine.1, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=",", fill="right")
#write.csv2(ONT_S_m_combine.2, "ONT_S_m_combine.2.csv")

IT_R_m_combine<-read.csv2("IT_R_m_combine.csv")
IT_R_m_combine.1<- ddply(IT_R_m_combine,"Taxonomy",numcolwise(sum))
IT_R_m_combine.2  <- separate(IT_R_m_combine.1, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";", fill="right")
#write.csv2(IT_R_m_combine.2, "IT_R_m_combine.2.csv")


###MinION####

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


####Can be used to add PPD column to sample data####
#variable1 = sample_data(data)$Platform
#variable2 = sample_data(data)$Pipeline
#variable3 = sample_data(data)$Database
#sample_data(data)$PPD <-  mapply(paste0, variable1, sep = "_", variable2, sep = "_", variable3)

## SP column added to sample data####
variable1 = sample_data(data)$Sample_original
variable2 = sample_data(data)$Platform
sample_data(data)$SP <-  mapply(paste0, variable1, sep = "_", variable2)


##FILTER
data.1 <- subset_taxa(data, Kingdom == "Bacteria" & Class != "Chloroplast" & Family != "mitochondria" &  Phylum != "Cyanobacteria/Chloroplast")

##add standard separation column
variable1 = sample_data(data.1)$Location
sample_data(data.1)$Type <-  mapply(paste0, variable1)
sample_data(data.1)$Type <- as.character(sample_data(data.1)$Type)
sample_data(data.1)$Type[sample_data(data.1)$Type == "Oudekraal"] <- "Sample"
sample_data(data.1)$Type[sample_data(data.1)$Type == "IEP"] <- "Sample"
sample_data(data.1)$Type[sample_data(data.1)$Type == "Muizenberg"] <- "Sample"
sample_data(data.1)$Type[sample_data(data.1)$Type == "None"] <- "Standard"

#Add new variable for plotting - alpha diversity
variable1 = sample_data(data.1)$Platform
variable2 = sample_data(data.1)$Type
sample_data(data.1)$PT <-  mapply(paste0, variable1, sep = "_", variable2)

##Alpha diversity
std <- subset_samples(data.1, Type == "Standard")
samples <- subset_samples(data.1, Type == "Sample")

p <- plot_richness(samples, x="Pipeline" , color = "Platform", measures = c("Shannon"))
p3 <- plot_richness(std, x="Pipeline" , measures = c("Shannon"))+ facet_grid(~Database, scales = "free_x") + geom_jitter()
p2 = p + theme_bw() + geom_boxplot(data=p$data, aes(x=Pipeline, y=value, color=Platform), alpha=0.1) + 
  facet_grid(~Database, scales = "free_x") + geom_point(data=p$data, aes(x=Pipeline, y=value, color=Platform), position = position_dodge(width = 0.75)) +
  geom_point(data=p3$data, aes(x=Pipeline, y=value, shape=PT), colour = "black", position = position_dodge(width = 0.75))
p2$layers <- p2$layers[-1]
p2

#Alpha diversity testing
library(emmeans)
library(multcomp)
rich <- estimate_richness(samples, measures = c("Shannon"))
variable1 = sample_data(samples)$Platform
variable2 = sample_data(samples)$Pipeline
variable3 = sample_data(samples)$Database
rich$tab <-  mapply(paste0, variable1, sep = "_", variable2, sep = "_", variable3)
res.aov <- aov(Shannon ~ tab, data = rich)
summary(res.aov)
lsmeans <- lsmeans(res.aov, 'tab')


contrasts <- list(
  "ONT_r_pipe" =c(0,0,0,0,-1,0,1,0),
  "IT_r_pipe" =c(-1,0,1,0,0,0,0,0),
  "ONT_s_pipe" =c(0,0,0,0,0,-1,0,1),
  "IT_s_pipe" =c(0,-1,0,1,0,0,0,0),
  "Q_R_plat" =c(0,0,-1,0,0,0,1,0),
  "Q_S_plat" =c(0,0,0,-1,0,0,0,1),
  "M_R_plat" =c(-1,0,0,0,1,0,0,0),
  "M_S_plat" =c(0,-1,0,0,0,1,0,0),
  "ONT_Q_db" =c(0,0,0,0,0,0,-1,1),
  "IT_Q_db" =c(0,0,-1,1,0,0,0,0),
  "ONT_M_db" =c(0,0,0,0,-1,1,0,0),
  "IT_M_db" =c(-1,1,0,0,0,0,0,0)
)

sum_test <- summary(as.glht(contrast(lsmeans, contrasts)), test = adjusted("bonferroni"))
capture.output(sum_test, file = "sum_test.txt")

