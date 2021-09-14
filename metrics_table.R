# Metrics -----------------------------------------------------------------
##Use unfiltered data!!! - phyloseq object = data

#separate out taxa groups of interest
archaea <- subset_taxa(data, Kingdom == "Archaea")
bacteria <- subset_taxa(data, Kingdom == "Bacteria")
unknown <- subset_taxa(data, Kingdom != "Bacteria"& Kingdom != "Archaea")
unclass <- subset_taxa(data, Kingdom == "Unclassified")
euks <- subset_taxa(data, Kingdom == "Eukaryota")
chloro <- subset_taxa(data, Phylum == "Cyanobacteria/Chloroplast" & Class == "Chloroplast")
mito <- subset_taxa(data, Family == "mitochondria") #none?

#find unclassified taxa; list them; convert to tax_tables
df <- as.data.frame(tax_table(bacteria))
unlcass_s <- df[grepl("X", df$Species), ]
unlcass_s <- as.matrix(unlcass_s)
unlcass_s = tax_table(unlcass_s)
unlcass_g <- df[grepl("X", df$Genus), ]
unlcass_g <- as.matrix(unlcass_g)
unlcass_g = tax_table(unlcass_g)
unlcass_f <- df[grepl("X", df$Family), ]
unlcass_f <- as.matrix(unlcass_f)
unlcass_f = tax_table(unlcass_f)
unlcass_o <- df[grepl("X", df$Order), ]
unlcass_o <- as.matrix(unlcass_o)
unlcass_o = tax_table(unlcass_o)
unlcass_c <- df[grepl("X", df$Class), ]
unlcass_c <- as.matrix(unlcass_c)
unlcass_c = tax_table(unlcass_c)
unlcass_p <- df[grepl("X", df$Phylum), ]
unlcass_p <- as.matrix(unlcass_p)
unlcass_p = tax_table(unlcass_p)

#crop unclassified taxa out for each level
#phylum
unclass_taxa_p <- subset(otu_table(bacteria), rownames(otu_table(bacteria)) %in% row.names(unlcass_p))
unclass_phyla <- merge_phyloseq(unclass_taxa_p, tax_table(bacteria), sample_data(bacteria))

#class
unclass_taxa_c <- subset(otu_table(bacteria), rownames(otu_table(bacteria)) %in% row.names(unlcass_c))
unclass_class <- merge_phyloseq(unclass_taxa_c, tax_table(bacteria), sample_data(bacteria))

#order
unclass_taxa_o <- subset(otu_table(bacteria), rownames(otu_table(bacteria)) %in% row.names(unlcass_o))
unclass_order <- merge_phyloseq(unclass_taxa_o, tax_table(bacteria), sample_data(bacteria))

#family
unclass_taxa_f <- subset(otu_table(bacteria), rownames(otu_table(bacteria)) %in% row.names(unlcass_f))
unclass_fam <- merge_phyloseq(unclass_taxa_f, tax_table(bacteria), sample_data(bacteria))

#genus
unclass_taxa_g <- subset(otu_table(bacteria), rownames(otu_table(bacteria)) %in% row.names(unlcass_g))
unclass_genus <- merge_phyloseq(unclass_taxa_g, tax_table(bacteria), sample_data(bacteria))

#species
unclass_taxa_s <- subset(otu_table(bacteria), rownames(otu_table(bacteria)) %in% row.names(unlcass_s))
unclass_spp<- merge_phyloseq(unclass_taxa_s, tax_table(bacteria), sample_data(bacteria))


##get numbers
bact <- sample_sums(bacteria)
phyla <- sample_sums(unclass_phyla)
class <- sample_sums(unclass_class)
ord <- sample_sums(unclass_order)
fam <- sample_sums(unclass_fam)
gen <- sample_sums(unclass_genus)
spp <- sample_sums(unclass_spp)

all <- rbind(bact, phyla, class, ord, fam, gen, spp)
write.csv(all, "all_uclass.csv")
