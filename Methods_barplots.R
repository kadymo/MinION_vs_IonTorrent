###DNA Standard Barplot##### 

dna_standard<-subset_samples(data.2, Type == "Standard")
mothur_dna_standard<-subset_samples(dna_standard, Pipeline == "MOTHUR")
qiime_dna_standard<-subset_samples(dna_standard, Pipeline == "QIIME")

tax_table(qiime_dna_standard)
mothur_dna_standard<-transform_sample_counts(mothur_dna_standard, function(x) 100 * x/sum(x))
qiime_dna_standard<-transform_sample_counts(qiime_dna_standard, function(x) 100 * x/sum(x))

dna_m = plot_bar(mothur_dna_standard, "SP", fill="Phylum")+
  facet_grid( Database~.,scales = "free_x",space="free_x")+
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
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dna_m

dna_q = plot_bar(qiime_dna_standard, "SP", fill="Phylum")+
  facet_grid( Database~.,scales = "free_x",space="free_x")+
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
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dna_q


###Barplots##### 
colourway<-c("#E91E53", "#EF9A9A", "#E67E22", "#F1C40F", "#AA8A0A",
             "#388E3C", "#56bc5b", "#c5e1a5", "#00695c", "#00668f",
             "#F2F2F2", "#444444", "#838383","#B3B3B3", "#000000")
##Mothur##

#Phylum
mothur_data_phylum_top20 = names(sort(taxa_sums(mothur_data_phylum), TRUE)[1:15])
mothur_data_phylum_top20 = prune_taxa(mothur_data_phylum_top20, mothur_data_phylum)
mothur_data_phylum_top20<-transform_sample_counts(mothur_data_phylum_top20, function(x) 100 * x/sum(x))

p_m = plot_bar(mothur_data_phylum_top20, "SP", fill="Phylum")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
p_m


#Class
mothur_data_class_top20 = names(sort(taxa_sums(mothur_data_class), TRUE)[1:15])
mothur_data_class_top20 = prune_taxa(mothur_data_class_top20, mothur_data_class)
mothur_data_class_top20<-transform_sample_counts(mothur_data_class_top20, function(x) 100 * x/sum(x))

c_m = plot_bar(mothur_data_class_top20, "SP", fill="Class")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
c_m

#Order
mothur_data_order_top20 = names(sort(taxa_sums(mothur_data_order), TRUE)[1:15])
mothur_data_order_top20 = prune_taxa(mothur_data_order_top20, mothur_data_order)
mothur_data_order_top20<-transform_sample_counts(mothur_data_order_top20, function(x) 100 * x/sum(x))

o_m = plot_bar(mothur_data_order_top20, "SP", fill="Order")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
o_m

#Family
mothur_data_family_top20 = names(sort(taxa_sums(mothur_data_family), TRUE)[1:15])
mothur_data_family_top20 = prune_taxa(mothur_data_family_top20, mothur_data_family)
mothur_data_family_top20<-transform_sample_counts(mothur_data_family_top20, function(x) 100 * x/sum(x))

f_m = plot_bar(mothur_data_family_top20, "SP", fill="Family")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
f_m   

#Genus
mothur_data_genus_top20 = names(sort(taxa_sums(mothur_data_genus), TRUE)[1:15])
mothur_data_genus_top20 = prune_taxa(mothur_data_genus_top20, mothur_data_genus)
mothur_data_genus_top20<-transform_sample_counts(mothur_data_genus_top20, function(x) 100 * x/sum(x))

g_m = plot_bar(mothur_data_genus_top20, "SP", fill="Genus")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
g_m    

#Species
mothur_data_species_top20 = names(sort(taxa_sums(mothur_data_species), TRUE)[1:15])
mothur_data_species_top20 = prune_taxa(mothur_data_species_top20, mothur_data_species)
mothur_data_species_top20<-transform_sample_counts(mothur_data_species_top20, function(x) 100 * x/sum(x))

s_m = plot_bar(mothur_data_species_top20, "SP", fill="Species")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
s_m    


##Qiime##

#Phylum
Qiime_data_phylum_top20 = names(sort(taxa_sums(Qiime_data_phylum), TRUE)[1:15])
Qiime_data_phylum_top20 = prune_taxa(Qiime_data_phylum_top20, Qiime_data_phylum)
Qiime_data_phylum_top20<-transform_sample_counts(Qiime_data_phylum_top20, function(x) 100 * x/sum(x))

p_q = plot_bar(Qiime_data_phylum_top20, "SP", fill="Phylum")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
p_q

#Class
Qiime_data_class_top20 = names(sort(taxa_sums(Qiime_data_class), TRUE)[1:15])
Qiime_data_class_top20 = prune_taxa(Qiime_data_class_top20, Qiime_data_class)
Qiime_data_class_top20<-transform_sample_counts(Qiime_data_class_top20, function(x) 100 * x/sum(x))

c_q = plot_bar(Qiime_data_class_top20, "SP", fill="Class")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
c_q

#Order
Qiime_data_order_top20 = names(sort(taxa_sums(Qiime_data_order), TRUE)[1:15])
Qiime_data_order_top20 = prune_taxa(Qiime_data_order_top20, Qiime_data_order)
Qiime_data_order_top20<-transform_sample_counts(Qiime_data_order_top20, function(x) 100 * x/sum(x))

o_q = plot_bar(Qiime_data_order_top20, "SP", fill="Order")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
o_q

#Family
Qiime_data_family_top20 = names(sort(taxa_sums(Qiime_data_family), TRUE)[1:15])
Qiime_data_family_top20 = prune_taxa(Qiime_data_family_top20, Qiime_data_family)
Qiime_data_family_top20<-transform_sample_counts(Qiime_data_family_top20, function(x) 100 * x/sum(x))

f_q = plot_bar(Qiime_data_family_top20, "SP", fill="Family")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
f_q   

#Genus
Qiime_data_genus_top20 = names(sort(taxa_sums(Qiime_data_genus), TRUE)[1:15])
Qiime_data_genus_top20 = prune_taxa(Qiime_data_genus_top20, Qiime_data_genus)
Qiime_data_genus_top20<-transform_sample_counts(Qiime_data_genus_top20, function(x) 100 * x/sum(x))

g_q = plot_bar(Qiime_data_genus_top20, "SP", fill="Genus")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
g_q    

#Species
Qiime_data_species_top20 = names(sort(taxa_sums(Qiime_data_species), TRUE)[1:15])
Qiime_data_species_top20 = prune_taxa(Qiime_data_species_top20, Qiime_data_species)
Qiime_data_species_top20<-transform_sample_counts(Qiime_data_species_top20, function(x) 100 * x/sum(x))

s_q = plot_bar(Qiime_data_species_top20, "SP", fill="Species")+
  facet_grid( Database~Type,scales = "free_x",space="free_x")+
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
  scale_fill_manual(values=colourway)
s_q

###Combined mothur and qiime plots##

phyla_m<- ggplotGrob(p_m) 
phyla_q <- ggplotGrob(p_q) 

phylum<-plot_grid (phyla_m, phyla_q,
                   ncol = 1, nrow = 2, align = "lrb", axis="rb")
phylum

