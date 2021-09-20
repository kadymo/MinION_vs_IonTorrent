setwd("~/Documents/MinION_Seq_data/Methods/MinION_qiime_analysis/qiime/V4_out/All")

#Split files
#Taxonomy
TAX_S = SILVA_all %>% select((Kingdom:Species))
TAX_R = RDP_all %>% select((Kingdom:Species))

#OTU table
OTU_S = SILVA_all %>% select(!(Kingdom:Species))
OTU_R = RDP_all %>% select(!(Kingdom:Species))

#sample table:
SAM_S <- read.csv("sample_s.csv", header=TRUE, row.names = 1, ",")
SAM_R <- read.csv("sample_r.csv", header=TRUE, row.names = 1, sep = ",")

#transform into matrices:
OTU_S <- as.matrix(OTU_S)
OTU_R <- as.matrix(OTU_R)
TAX_S <- as.matrix(TAX_S)
TAX_R <- as.matrix(TAX_R)

#transform into phyloseq objects:
OTU_S = otu_table(OTU_S, taxa_are_rows = TRUE)
TAX_S = tax_table(TAX_S)
samples_s = sample_data(SAM_S)
ONT_S <- phyloseq(OTU_S, TAX_S, samples_s)

OTU_R = otu_table(OTU_R, taxa_are_rows = TRUE)
TAX_R = tax_table(TAX_R)
samples_r = sample_data(SAM_R)
ONT_R <- phyloseq(OTU_R, TAX_R, samples_r)

##Write files

write.csv(tax_table(ONT_R), "ONT_R_q_tax.csv")
write.csv(otu_table(ONT_R), "ONT_R_q_otus.csv")
write.csv(sample_data(ONT_R), "ONT_R_q_sam.csv")

write.csv(tax_table(ONT_S), "ONT_S_q_tax.csv")
write.csv(otu_table(ONT_S), "ONT_S_q_otus.csv")
write.csv(sample_data(ONT_S), "ONT_S_q_sam.csv")



