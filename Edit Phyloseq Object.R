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

tax<-read.csv2("Phyloseq_final_tax.csv", header=TRUE)
tax<- column_to_rownames(tax, var="rowname")

Kingdom = tax$Kingdom
Phylum = tax$Phylum
Class = tax$Class
Order= tax$Order
Family=tax$Family
Genus= tax$Genus
Species= tax$Species
variable1= "X"
variable2= "_X"


##Go column by column 

#Change unknown Kingdom to Unclassifed
tax$Kingdom <- replace(tax$Kingdom,tax$Kingdom =="unknown", "Unclassified")

#Replace NA Kingdom to Unclassifed
tax$Kingdom <-tax$Kingdom %>% replace_na("Unclassified")

#Change all "_unclassified" to "_X" in Phylum
tax$Phylum<-gsub("_unclassified", "_X", tax$Phylum)

#Replace NA Phylum to Unclassifed
tax$Phylum <-tax$Phylum %>% replace_na("Unclassified")

#Replace Unknown_x Phylum to Unclassifed_x
tax$Phylum  <- replace(tax$Phylum ,tax$Phylum =="unknown_X", "Unclassified_X")

#Replace repeated unclassifed kingdom with  Unclassifed_x
tax$Phylum  <- replace(tax$Phylum ,tax$Phylum =="Unclassified", "Unclassified_X")

#making blanks NA
tax.1<-tax %>% mutate_if(is.character, list(~na_if(.,""))) 

#Change all "_unclassified" to "_X" in Class
tax1.1<-tax.1 %>% 
  mutate(Class = if_else(str_ends(Class, "_unclassified"),mapply(paste0, Phylum, variable2), Class))

#Replace any class with phylum that end in _X with _XX, with variable1 = X
tax.2<-tax1.1 %>% 
  mutate(Class = if_else(str_ends(Phylum, "_X"),mapply(paste0, Phylum, variable1), Class))

#Replace any class that are NA with phylum + _X
tax.3<-tax.2 %>% 
  mutate(Class = if_else(is.na(Class),mapply(paste0, Phylum, variable2), Class))

#Replace any class that have _cl with _X
tax.4<-tax.3 %>% 
  mutate(Class = if_else(str_ends(Class, "_cl"),mapply(paste0, Phylum, variable2), Class))

#Replace any class that have uncultured with Phylum_uncultured
tax.4.1<-tax.4 %>% 
  mutate(Class = if_else(str_detect(Class, 'uncultured'), mapply(paste0, Phylum, sep = "_", Class), Class))

#Replace any class that have marine metagenome with Phylum_marine metagenom
tax.4.2<-tax.4.1 %>% 
  mutate(Class = if_else(Class == "marine metagenome", mapply(paste0, Phylum, sep = "_", Class), Class))

tax.4.3<-tax.4.2 %>% 
  mutate(Class = if_else(Class == "metagenome", mapply(paste0, Phylum, sep = "_", Class), Class))

tax.4.4<-tax.4.3 %>% 
  mutate(Class = if_else(Class == "Unknown_Class", mapply(paste0, Phylum, variable2), Class))

#Remove _or from orders
tax.5<-tax.4.4 %>% 
  mutate_at("Order", str_replace, "_or", "")

#Replace any order with class that end in _X with _XX, with variable1 = X
tax.6<-tax.5 %>% 
  mutate(Order = if_else(str_ends(Class, "_X"),mapply(paste0, Class, variable1), Order))

#Replace any order with class that end in _XX with _XXX, with variable1 = X
tax.7<-tax.6 %>% 
  mutate(Order = if_else(str_ends(Class, "_XX"),mapply(paste0, Class, variable1), Order))

###chnage name of one order
tax.8<-tax.7 %>% 
  mutate_at("Order", str_replace, "44524", "11-24")

tax.8.1<-tax.8 %>% 
  mutate_at("Order", str_replace, "24-Nov", "11-24")

##replace orders that are teh same as class with class+ _X
tax.9<-tax.8.1 %>%
  mutate(Order = ifelse(Order==Class,mapply(paste0, Class, variable2), Order))

#Change all "_unclassified" to "_X" in Order
tax.9$Order<-gsub("_unclassified", "_X", tax.9$Order)

#Replace any order that are NA with class + _X
tax.10<-tax.9 %>% 
  mutate(Order = if_else(is.na(Order),mapply(paste0, Class, variable2),Order))

#Replace any order with class that end in _uncultured with , with variable1 = X
tax.11<-tax.10 %>% 
  mutate(Order = if_else(str_detect(Class, 'uncultured'),mapply(paste0, Class, variable2), Order))

#Replace any order with class that end in _uncultured bacterium with , with variable1 = X
tax.12<-tax.11 %>% 
  mutate(Order = if_else(str_ends(Class, "_marine metagenome"),mapply(paste0, Class, variable2), Order))

tax.12.1<-tax.12 %>% 
  mutate(Order = if_else(str_ends(Class, "_metagenome"),mapply(paste0, Class, variable2), Order))

tax.12.2<-tax.12.1 %>% 
  mutate(Order = if_else(Order == "Unknown_Order", mapply(paste0, Class, variable2), Order))

tax.12.3<-tax.12.2 %>% 
  mutate(Order  = if_else(Order  == "metagenome", mapply(paste0, Class, variable2), Order ))

tax.12.4<-tax.12.3%>% 
  mutate(Order = if_else(str_detect(Order, 'unidentified'),mapply(paste0, Class, variable2), Order))

tax.12.5<-tax.12.4 %>% 
  mutate(Order  = if_else(str_starts(Order, "Family"), mapply(paste0, Class, sep = "_", Order), Order ))

#Replace any Order that have uncultured with Class_uncultured
tax.13<-tax.12.5 %>% 
  mutate(Order = if_else(str_detect(Order, 'uncultured'), mapply(paste0, Class, sep = "_", Order), Order))

tax.13.1<-tax.13 %>% 
  mutate(Order = if_else(Order == "marine metagenome", mapply(paste0, Class, sep = "_", Order), Order))

#Remove _fa from families
tax.14<-tax.13.1 %>% 
  mutate_at("Family", str_replace, "_fa", "")

#Replace any family with order that end in _X with _XX, with variable1 = X
tax.15<-tax.14 %>% 
  mutate(Family = if_else(str_ends(Order, "_X"),mapply(paste0, Order, variable1), Family))

#Replace any family with order that end in _XX with _XXX, with variable1 = X
tax.16<-tax.15 %>% 
  mutate(Family = if_else(str_ends(Order, "_XX"),mapply(paste0, Order, variable1), Family))

#Replace any family with order that end in _XXX with _XXXX, with variable1 = X
tax.17<-tax.16 %>% 
  mutate(Family = if_else(str_ends(Order, "_XXX"),mapply(paste0, Order, variable1), Family))

##replace families that are the same as order with order + _X
tax.18<-tax.17 %>%
  mutate(Family = ifelse(Family ==Order,mapply(paste0, Order, variable2), Family))

#Change all "_unclassified" to "_X" in Family
tax.18$Family<-gsub("_unclassified", "_X", tax.18$Family)

#Replace any order that are NA with class + _X
tax.19<-tax.18 %>% 
  mutate(Family = if_else(is.na(Family),mapply(paste0, Order, variable2),Family))

#Replace any Family with class that end in _uncultured bacterium with , with variable1 = X
tax.19.1<-tax.19 %>% 
  mutate(Family = if_else(str_detect(Order, 'uncultured'),mapply(paste0, Order, variable2), Family))

#Replace any amily with class that end in _uncultured bacterium with , with variable1 = X
tax.19.2<-tax.19.1%>% 
  mutate(Family = if_else(str_ends(Order, "_marine metagenome"),mapply(paste0, Order, variable2), Family))

tax.19.3<-tax.19.2 %>% 
  mutate(Family  = if_else(Family  == "Unknown_Family", mapply(paste0, Order, variable2), Family ))

tax.19.4<-tax.19.3 %>% 
  mutate(Family  = if_else(Family  == "Unknown Family", mapply(paste0, Order, variable2), Family ))

tax.19.5<-tax.19.4 %>% 
  mutate(Family  = if_else(Family  == "metagenome", mapply(paste0, Order, variable2), Family ))

tax.19.6<-tax.19.5%>% 
  mutate(Family = if_else(str_detect(Family, 'unidentified'),mapply(paste0, Order, variable2), Family))

tax.19.7<-tax.19.6 %>% 
  mutate(Family  = if_else(Family  == "hydrothermal vent metagenome", mapply(paste0, Order, sep = "_", Family), Family ))

tax.19.8<-tax.19.7 %>% 
  mutate(Family  = if_else(str_detect(Order, "Family"), mapply(paste0, Order, variable2), Family ))

tax.19.9<-tax.19.8 %>% 
  mutate(Family  = if_else(str_starts(Family, "Family"), mapply(paste0, Order, sep = "_", Family), Family ))

#Replace any Order that have uncultured with Class_uncultured
tax.20<-tax.19.9 %>% 
  mutate(Family  = if_else(str_detect(Family , "^uncultured"), mapply(paste0, Order, sep = "_", Family ), Family ))

tax.20.1<-tax.20 %>% 
  mutate(Family  = if_else(Family  == "marine metagenome", mapply(paste0, Order, sep = "_", Family ), Family ))

#Remove _ge from genera
tax.21<-tax.20.1 %>% 
  mutate_at("Genus", str_replace, "_ge", "")

#Replace any Genus with order that end in _X with _XX, with variable1 = X
tax.22<-tax.21 %>% 
  mutate(Genus = if_else(str_ends(Family, "_X"),mapply(paste0, Family, variable1), Genus))

#Replace any Genus with Family that end in _XX with _XXX, with variable1 = X
tax.23<-tax.22 %>% 
  mutate(Genus = if_else(str_ends(Family, "_XX"),mapply(paste0, Family, variable1), Genus))

#Replace any Genus with Family that end in _XXX with _XXXX, with variable1 = X
tax.24<-tax.23 %>% 
  mutate(Genus = if_else(str_ends(Family, "_XXX"),mapply(paste0, Family, variable1), Genus))

#Replace any Genus with Family that end in _XXXX with _XXXXX, with variable1 = X
tax.25<-tax.24 %>% 
  mutate(Genus = if_else(str_ends(Family, "_XXXX"),mapply(paste0, Family, variable1), Genus))

##replace genera that are the same as Family with Family + _X
tax.26<-tax.25 %>%
  mutate(Genus = ifelse(Genus ==Family,mapply(paste0, Family, variable2), Genus))

#Change all "_unclassified" to "_X" in Genus
tax.26$Genus<-gsub("_unclassified", "_X", tax.26$Genus)

#Replace any Family that are NA with class + _X
tax.27<-tax.26 %>% 
  mutate(Genus = if_else(is.na(Genus),mapply(paste0, Family, variable2),Genus))

#Replace any genus with class that end in _uncultured bacterium with , with variable1 = X
tax.27.1<-tax.27 %>% 
  mutate(Genus = if_else(str_detect(Family, 'uncultured'),mapply(paste0, Family, variable2), Genus))

#Replace any genus with class that end in _uncultured bacterium with , with variable1 = X
tax.27.2<-tax.27.1%>% 
  mutate(Genus = if_else(str_ends(Family, "_marine metagenome"),mapply(paste0, Family, variable2), Genus))

#Replace any genus that have uncultured with Class_uncultured
tax.28<-tax.27.2 %>% 
  mutate(Genus  = if_else(str_detect(Genus , "^uncultured"), mapply(paste0, Family, sep = "_", Genus ), Genus ))

tax.28.1<-tax.28 %>% 
  mutate(Genus  = if_else(Genus  == "marine metagenome", mapply(paste0, Family, sep = "_", Genus ), Genus ))

tax.28.2<-tax.28.1 %>% 
  mutate(Genus = if_else(str_detect(Family, "Family"), mapply(paste0, Family, variable2), Genus ))

tax.28.3<-tax.28.2 %>% 
  mutate(Genus  = if_else(str_starts(Genus, "Family"), mapply(paste0, Family, sep = "_", Genus), Genus ))

tax.28.4<-tax.28.3 %>% 
  mutate(Genus  = if_else(Genus  == "metagenome", mapply(paste0, Family, variable2), Genus ))

tax.28.5<-tax.28.4 %>% 
  mutate(Genus  = if_else(str_detect(Family, "hydrothermal vent metagenome"), mapply(paste0, Family, variable2), Genus ))

tax.28.6<-tax.28.5 %>% 
  mutate(Genus  = if_else(str_starts(Genus, "hydrothermal vent metagenome"), mapply(paste0, Family, sep = "_", Genus), Genus ))

tax.28.7<-tax.28.6 %>% 
  mutate(Genus  = if_else(Genus  == "Ambiguous_taxa", mapply(paste0, Family, sep = "_", Genus), Genus ))


#Remove _sp from species
tax.29<-tax.28.7 %>% 
  mutate_at("Species", str_replace, "_sp", "")

#Replace any Species with order that end in _X with _XX, with variable1 = X
tax.30<-tax.29 %>% 
  mutate(Species = if_else(str_ends(Genus, "_X"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XX with _XXX, with variable1 = X
tax.31<-tax.30 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XX"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XXX with _XXXX, with variable1 = X
tax.32<-tax.31 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XXX"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XXXX with _XXXXX, with variable1 = X
tax.33<-tax.32 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XXXX"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XXXXX with _XXXXXX, with variable1 = X
tax.34<-tax.33 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XXXXX"),mapply(paste0, Genus, variable1), Species))

##replace species that are the same as Genus with Genus + _X
tax.35<-tax.34 %>%
  mutate(Species = ifelse(Species ==Genus,mapply(paste0, Genus, variable2), Species))

#Change all "_unclassified" to "_X" in Species
tax.35$Species<-gsub("_unclassified", "_X", tax.35$Species)

#Replace any Species that are NA with class + _X
tax.36<-tax.35 %>% 
  mutate(Species = if_else(is.na(Species),mapply(paste0, Genus, variable2),Species))

#Replace any Species with class that end in _uncultured bacterium with , with variable1 = X
tax.36.1<-tax.36 %>% 
  mutate(Species = if_else(str_detect(Genus, 'uncultured'),mapply(paste0, Genus, variable2), Species))

#Replace any Species with class that end in _uncultured bacterium with , with variable1 = X
tax.36.2<-tax.36.1%>% 
  mutate(Species = if_else(str_ends(Genus, "_marine metagenome"),mapply(paste0, Genus, variable2), Species))

#Replace any Species that have uncultured with Class_uncultured
tax.37<-tax.36.2 %>% 
  mutate(Species  = if_else(str_detect(Species , "^uncultured"), mapply(paste0, Genus, sep = "_", Species ), Species ))

tax.37.1<-tax.37 %>% 
  mutate(Species  = if_else(str_ends(Species, 'metagenome'), mapply(paste0, Genus, sep = "_", Species ), Species ))

tax.37.2<-tax.37.1 %>% 
  mutate(Species = if_else(str_detect(Genus, "Family"), mapply(paste0, Genus, variable2), Species  ))

tax.37.3<-tax.37.2 %>% 
  mutate(Species  = if_else(str_starts(Species, "Family"), mapply(paste0, Genus, sep = "_", Species), Species ))

tax.37.4<-tax.37.3 %>% 
  mutate(Species = if_else(str_detect(Genus, "hydrothermal vent metagenome"), mapply(paste0, Genus, variable2), Species))

tax.37.5<-tax.37.4 %>% 
  mutate(Species  = if_else(str_starts(Species, "hydrothermal vent metagenome"), mapply(paste0, Genus, sep = "_", Species), Species ))

tax.37.6<-tax.37.5 %>% 
  mutate(Species  = if_else(str_detect(Genus, "Ambiguous_taxa"), mapply(paste0, Genus, variable2), Species ))

tax.37.7<-tax.37.6 %>% 
  mutate(Species  = if_else(Species == "Ambiguous_taxa", mapply(paste0, Genus, sep = "_", Species), Species ))

#Concatenate Genus+Species if species names starts lowercase except "uncultured"
tax.38<-tax.37.7 %>% 
  mutate(Species = if_else(str_detect(Species, '^[:lower:]'),mapply(paste0, Genus, sep = " ", Species),Species))


###replacing dupes
tax.39<-tax.38 %>% 
  mutate(Phylum = ifelse(Class== "Campylobacteria" & Phylum == "Epsilonbacteraeota","Campilobacterota",Phylum))

tax.39.1<-tax.39 %>% 
  mutate(Phylum = ifelse(Class== "Chlorobia" & Phylum == "Bacteroidetes","Chlorobi",Phylum))

tax.39.2<-tax.39.1 %>% 
  mutate(Phylum = ifelse(Class== "Ignavibacteria" & Phylum == "Ignavibacteriae","Bacteroidetes",Phylum))

tax.39.3<-tax.39.2 %>% 
  mutate(Phylum = ifelse(Class== "Rhodothermia" & Phylum == "Bacteroidetes","Rhodothermaeota",Phylum))

tax.39.4<-tax.39.3 %>% 
  mutate(Class = ifelse(Order== "Caldilineales" & Class  == "Anaerolineae","Caldilineae",Class ))

tax.39.5<-tax.39.4 %>% 
  mutate(Class = ifelse(Order== "Opitutales" & Class  == "Verrucomicrobiae","Opitutae",Class ))   

tax.39.6<-tax.39.5 %>% 
  mutate(Class = ifelse(Order== "Sphingobacteriales" & Class  == "Bacteroidia","Sphingobacteriia",Class ))   

tax.39.7<-tax.39.6 %>% 
  mutate(Class = ifelse(Order== "Chitinophagales" & Class  == "Bacteroidia","Chitinophagia",Class ))   

tax.39.8<-tax.39.7 %>% 
  mutate(Class = ifelse(Order== "Cytophagales" & Class  == "Bacteroidia","Cytophagia",Class ))

tax.39.9<-tax.39.8 %>% 
  mutate(Class = ifelse(Order== "Flavobacteriales" & Class  == "Flavobacteriia","Bacteroidia",Class ))

tax.39.10<-tax.39.9 %>% 
  mutate(Class = ifelse(Order== "Chlamydiales" & Class  == "Chlamydiia","Chlamydiae",Class ))

tax.39.11<-tax.39.10 %>% 
  mutate(Class = ifelse(Order== "Oligoflexales" & Class  == "Deltaproteobacteria","Oligoflexia",Class ))

tax.39.12<-tax.39.11 %>% 
  mutate(Class = ifelse(Order== "Balneolales" & Class  == "Rhodothermia","Balneolia",Class ))

tax.39.13<-tax.39.12 %>% 
  mutate(Class = ifelse(Order== "Blastocatellales" & Class  == "Blastocatellia","Blastocatellia_(Subgroup_4)",Class ))

tax.39.14<-tax.39.13 %>% 
  mutate(Class = ifelse(Class  == "Blastocatellia (Subgroup 4)","Blastocatellia_(Subgroup_4)",Class ))

tax.39.15<-tax.39.14 %>% 
  mutate(Order = ifelse(Family== "Campylobacteria_X" & Order  == "Blastocatellia","Blastocatellia_(Subgroup_4)",Order ))

tax.39.16<-tax.39.15 %>% 
  mutate(Phylum = ifelse(Class== "Balneolia" & Phylum == "Rhodothermaeota","Balneolaeota",Phylum))

tax.39.17<-tax.39.16 %>% 
  mutate(Order = ifelse(Family== "Nocardiaceae" & Order == "Corynebacteriales","Mycobacteriales",Order))

tax.39.18<-tax.39.17 %>% 
  mutate(Order = ifelse(Family== "Trueperaceae" & Order == "Deinococcales","Trueperales",Order))

tax.39.19<-tax.39.18 %>% 
  mutate(Order= ifelse(Family== "Acidaminococcaceae" & Order == "Selenomonadales","Acidaminococcales",Order))

tax.39.20<-tax.39.19 %>% 
  mutate(Order = ifelse(Family== "Ilumatobacteraceae" & Order == "Microtrichales","Acidimicrobiales",Order))

tax.39.21<-tax.39.20 %>% 
  mutate(Order = ifelse(Family== "Iamiaceae" & Order == "Microtrichales","Acidimicrobiales",Order))

tax.39.22<-tax.39.21 %>% 
  mutate(Order = ifelse(Family== "Bacteriovoracaceae" & Order == "Bdellovibrionales","Bacteriovoracales",Order))

tax.39.22<-tax.39.21 %>% 
  mutate(Order = ifelse(Family== "Lentimicrobiaceae" & Order == "Sphingobacteriales","Bacteroidales",Order))

tax.39.23<-tax.39.22 %>% 
  mutate(Order = ifelse(Family== "Leptospiraceae" & Order == "Spirochaetales","Leptospirales",Order))

tax.39.24<-tax.39.23 %>% 
  mutate(Order = ifelse(Family== "Rhodocyclaceae" & Order == "Betaproteobacteriales","Rhodocyclales",Order))

tax.39.25<-tax.39.24 %>% 
  mutate(Order = ifelse(Family== "Burkholderiaceae" & Order == "Betaproteobacteriales","Burkholderiales",Order))

tax.39.26<-tax.39.25 %>% 
  mutate(Order = ifelse(Family== "Methylophilaceae" & Order == "Betaproteobacteriales","Nitrosomonadales",Order))

tax.39.27<-tax.39.26 %>% 
  mutate(Order = ifelse(Family== "Gallionellaceae" & Order == "Betaproteobacteriales","Nitrosomonadales",Order))

tax.39.28<-tax.39.27 %>% 
  mutate(Order = ifelse(Family== "Nitrosomonadaceae" & Order == "Betaproteobacteriales","Nitrosomonadales",Order))

tax.39.29<-tax.39.28 %>% 
  mutate(Order = ifelse(Family== "Calditrichaceae" & Order == "Deferribacterales","Calditrichales",Order))

tax.39.30<-tax.39.29 %>% 
  mutate(Order = ifelse(Family== "Parvularculaceae" & Order == "Caulobacterales","Parvularculales",Order))

tax.39.31<-tax.39.30 %>% 
  mutate(Order = ifelse(Family== "Saprospiraceae" & Order == "Chitinophagales","Saprospirales",Order))

tax.39.32<-tax.39.31 %>% 
  mutate(Order = ifelse(Family== "Woeseiaceae" & Order == "Steroidobacterales","Chromatiales",Order))

tax.39.33<-tax.39.32 %>% 
  mutate(Order = ifelse(Family== "Thalassobaculaceae" & Order == "Thalassobaculales","Rhodospirillales",Order))

tax.39.34<-tax.39.33 %>% 
  mutate(Order = ifelse(Family== "Ectothiorhodospiraceae" & Order == "Ectothiorhodospirales","Chromatiales",Order))

tax.39.35<-tax.39.34 %>% 
  mutate(Order = ifelse(Family== "Enterobacteriaceae" & Order == "Enterobacteriales","Enterobacterales",Order))

tax.39.36<-tax.39.35 %>% 
  mutate(Order = ifelse(Family== "Francisellaceae" & Order == "Francisellales","Thiotrichales",Order))

tax.39.37<-tax.39.36 %>% 
  mutate(Order = ifelse(Family== "Corynebacteriaceae" & Order == "Corynebacteriales","Mycobacteriales",Order))

tax.39.38<-tax.39.37 %>% 
  mutate(Order = ifelse(Family== "Mycobacteriaceae" & Order == "Corynebacteriales","Mycobacteriales",Order))

tax.39.39<-tax.39.38 %>% 
  mutate(Order = ifelse(Family== "Gemmataceae" & Order == "Planctomycetales","Gemmatales",Order))

tax.39.40<-tax.39.39 %>% 
  mutate(Order = ifelse(Family== "Isosphaeraceae" & Order == "Planctomycetales","Isosphaerales",Order))

tax.39.41<-tax.39.40 %>% 
  mutate(Order = ifelse(Family== "Kiloniellaceae" & Order == "Rhodovibrionales","Kiloniellales",Order))

tax.39.42<-tax.39.41 %>% 
  mutate(Order = ifelse(Family== "Coxiellaceae" & Order == "Coxiellales","Legionellales",Order))

tax.39.43<-tax.39.42 %>% 
  mutate(Order = ifelse(Family== "Leptospiraceae" & Order == "Spirochaetales","Leptospirales",Order))

tax.39.44<-tax.39.43%>% 
  mutate(Order = ifelse(Family== "Steroidobacteraceae" & Order == "Steroidobacterales","Nevskiales",Order))

tax.39.45<-tax.39.44%>% 
  mutate(Order = ifelse(Family== "Veillonellaceae" & Order == "Selenomonadales","Veillonellales",Order))

tax.39.46<-tax.39.45%>% 
  mutate(Order = ifelse(Family== "Reyranellaceae" & Order == "Reyranellales","Rhodospirillales",Order))

tax.39.47<-tax.39.46%>% 
  mutate(Order = ifelse(Family== "Acetobacteraceae" & Order == "Acetobacterales","Rhodospirillales",Order))

tax.39.48<-tax.39.47%>% 
  mutate(Order = ifelse(Family== "Puniceicoccaceae" & Order == "Opitutales","Puniceicoccales",Order))

tax.39.49<-tax.39.48%>% 
  mutate(Order = ifelse(Family== "Acetobacteraceae" & Order == "Acetobacterales","Rhodospirillales",Order))

tax.39.50<-tax.39.49%>% 
  mutate(Order = ifelse(Family== "Geminicoccaceae" & Order == "Tistrellales","Rhodospirillales",Order))

tax.39.51<-tax.39.50%>% 
  mutate(Order = ifelse(Family== "Bacteriovoracaceae" & Order == "Bdellovibrionales","Bacteriovoracales",Order))

tax.39.52<-tax.39.51%>% 
  mutate(Order = ifelse(Family== "Parvibaculaceae" & Order == "Parvibaculales","Rhizobiales",Order))

tax.39.53<-tax.39.52%>% 
  mutate(Order = ifelse(Family== "Piscirickettsiaceae" & Order == "Piscirickettsiales","Thiotrichales",Order))

tax.39.54<-tax.39.53%>% 
  mutate(Family = ifelse(Genus== "Bryobacter" & Family == "Solibacteraceae (Subgroup 3)","Solibacteraceae_(Subgroup_3)",Family))

tax.39.55<-tax.39.54%>% 
  mutate(Family = ifelse(Genus== "Parabacteroides" & Family == "Porphyromonadaceae","Tannerellaceae",Family))

tax.39.56<-tax.39.55%>% 
  mutate(Family = ifelse(Genus== "Fulvivirga" & Family == "Cyclobacteriaceae","Fulvivirgaceae",Family))

tax.39.57<-tax.39.56%>% 
  mutate(Family = ifelse(Genus== "Marivirga" & Family == "Cyclobacteriaceae","Marivirgaceae",Family))

tax.39.58<-tax.39.57%>% 
  mutate(Family = ifelse(Genus== "Brumimicrobium" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

tax.39.59<-tax.39.58%>% 
  mutate(Family = ifelse(Genus== "Crocinitomix" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

tax.39.60<-tax.39.59%>% 
  mutate(Family = ifelse(Genus== "Leuconostoc" & Family == "Lactobacillaceae","Leuconostocaceae",Family))

tax.39.61<-tax.39.60%>% 
  mutate(Family = ifelse(Genus== "Mitsuokella" & Family == "Veillonellaceae","Selenomonadaceae",Family))

tax.39.62<-tax.39.61%>% 
  mutate(Family = ifelse(Genus== "Planctomicrobium" & Family == "Rubinisphaeraceae","Planctomycetaceae",Family))

tax.39.63<-tax.39.62%>% 
  mutate(Family = ifelse(Genus== "Methylobacterium" & Family == "Beijerinckiaceae","Methylobacteriaceae",Family))

tax.39.64<-tax.39.63%>% 
  mutate(Family = ifelse(Genus== "Nordella" & Family == "Rhizobiales Incertae Sedis","Rhizobiales_Incertae_Sedis",Family))

tax.39.65<-tax.39.64%>% 
  mutate(Family = ifelse(Genus== "Phreatobacter" & Family == "Rhizobiales_Incertae_Sedis","Phreatobacteraceae",Family))

tax.39.66<-tax.39.65%>% 
  mutate(Family = ifelse(Genus== "Phreatobacter" & Family == "Rhizobiales Incertae Sedis","Phreatobacteraceae",Family))

tax.39.67<-tax.39.66%>% 
  mutate(Family = ifelse(Genus== "Bradyrhizobium" & Family == "Xanthobacteraceae","Bradyrhizobiaceae",Family))

tax.39.68<-tax.39.67%>% 
  mutate(Family = ifelse(Genus== "Magnetospira" & Family == "Magnetospiraceae","Thalassospiraceae",Family))

tax.39.69<-tax.39.68%>% 
  mutate(Family = ifelse(Genus== "Marinobacter" & Family == "Marinobacteraceae","Alteromonadaceae",Family))

tax.39.70<-tax.39.69%>% 
  mutate(Family = ifelse(Genus== "Psychrobium" & Family == "Alteromonadales_incertae_sedis","Shewanellaceae",Family))

tax.39.71<-tax.39.70%>% 
  mutate(Family = ifelse(Genus== "Alcaligenes" & Family == "Burkholderiaceae","Alcaligenaceae",Family))

tax.39.72<-tax.39.71%>% 
  mutate(Family = ifelse(Genus== "Rhodoferax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.39.73<-tax.39.72%>% 
  mutate(Family = ifelse(Genus== "Xylophilus" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.39.74<-tax.39.73%>% 
  mutate(Family = ifelse(Genus== "Methyloversatilis" & Family == "Rhodocyclaceae","Sterolibacteriaceae",Family))

tax.39.75<-tax.39.74%>% 
  mutate(Family = ifelse(Genus== "Zoogloea" & Family == "Rhodocyclaceae","Zoogloeaceae",Family))

tax.39.76<-tax.39.75%>% 
  mutate(Family = ifelse(Genus== "Providencia" & Family == "Enterobacteriaceae","Morganellaceae",Family))

tax.39.77<-tax.39.76%>% 
  mutate(Family = ifelse(Genus== "Amphritea" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.39.78<-tax.39.77%>% 
  mutate(Family = ifelse(Genus== "Profundimonas" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.39.79<-tax.39.78%>% 
  mutate(Family = ifelse(Genus== "Oleibacter" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

tax.39.80<-tax.39.79%>% 
  mutate(Family = ifelse(Genus== "Nevskia" & Family == "Solimonadaceae","Nevskiaceae",Family))

tax.39.81<-tax.39.80%>% 
  mutate(Family = ifelse(Genus== "Luteolibacter" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

tax.39.82<-tax.39.81%>% 
  mutate(Family = ifelse(Genus== "Mitsuokella" & Family == "Veillonellaceae","Selenomonadaceae",Family))

tax.39.83<-tax.39.82%>% 
  mutate(Family = ifelse(Genus== "Selenomonas" & Family == "Veillonellaceae","Selenomonadaceae",Family))

tax.39.84<-tax.39.83%>% 
  mutate(Family = ifelse(Genus== "Acidovorax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.39.85<-tax.39.84%>% 
  mutate(Family = ifelse(Genus== "Pseudacidovorax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.39.86<-tax.39.85%>% 
  mutate(Family = ifelse(Genus== "Curvibacter" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.39.87<-tax.39.86%>% 
  mutate(Family = ifelse(Genus== "Hydrogenophaga" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.39.88<-tax.39.87%>% 
  mutate(Family = ifelse(Genus== "Thauera" & Family == "Rhodocyclaceae","Zoogloeaceae",Family))

tax.39.89<-tax.39.88%>% 
  mutate(Family = ifelse(Genus== "Neiella" & Family == "Alteromonadaceae","Alteromonadales_incertae_sedis",Family))

tax.39.90<-tax.39.89%>% 
  mutate(Family = ifelse(Genus== "Macellibacteroides" & Family == "Tannerellaceae","Porphyromonadaceae",Family))

tax.39.91<-tax.39.90%>% 
  mutate(Family = ifelse(Genus== "Salinirepens" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

tax.39.92<-tax.39.91%>% 
  mutate(Family = ifelse(Genus== "Dongia" & Family == "Dongiaceae","Rhodospirillaceae",Family))

tax.39.93<-tax.39.92%>% 
  mutate(Family = ifelse(Genus== "Novosphingobium" & Family == "Sphingomonadaceae","Erythrobacteraceae",Family))

tax.39.94<-tax.39.93%>% 
  mutate(Family = ifelse(Genus== "Sphaerochaeta" & Family == "Spirochaetaceae","Sphaerochaetaceae",Family))

tax.39.95<-tax.39.94%>% 
  mutate(Family = ifelse(Genus== "Thiomicrorhabdus" & Family == "Thiomicrospiraceae","Piscirickettsiaceae",Family))

tax.39.96<-tax.39.95%>% 
  mutate(Family = ifelse(Genus== "Methylophaga" & Family == "Methylophagaceae","Piscirickettsiaceae",Family))

tax.39.97<-tax.39.96%>% 
  mutate(Family = ifelse(Genus== "Pseudohongiella" & Family == "Pseudohongiellaceae","Oceanospirillales_incertae_sedis",Family))

tax.39.98<-tax.39.97%>% 
  mutate(Family = ifelse(Genus== "Reinekea" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

tax.39.99<-tax.39.98%>% 
  mutate(Family = ifelse(Genus== "Nitrincola" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.40<-tax.39.99%>% 
  mutate(Family = ifelse(Genus== "Neptunomonas" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.40.1<-tax.40%>% 
  mutate(Family = ifelse(Genus== "Neptuniibacter" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.40.2<-tax.40.1%>% 
  mutate(Family = ifelse(Genus== "Marinobacterium" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.40.3<-tax.40.2%>% 
  mutate(Family = ifelse(Genus== "Pantoea" & Family == "Enterobacteriaceae","Erwiniaceae",Family))

tax.40.4<-tax.40.3%>% 
  mutate(Family = ifelse(Genus== "Rheinheimera" & Family == "Alteromonadaceae","Chromatiaceae",Family))

tax.40.4.1<-tax.40.4%>% 
  mutate(Family = ifelse(Genus== "Desulfatiglans" & Family == "Desulfarculaceae","Desulfobacteraceae",Family))

tax.40.5<-tax.40.4.1%>% 
  mutate(Family = ifelse(Genus== "Ferribacterium" & Family == "Rhodocyclaceae","Azonexaceae",Family))

tax.40.6<-tax.40.5%>% 
  mutate(Family = ifelse(Genus== "Noviherbaspirillum" & Family == "Burkholderiaceae","Oxalobacteraceae",Family))

tax.40.7<-tax.40.6%>% 
  mutate(Family = ifelse(Genus== "Massilia" & Family == "Burkholderiaceae","Oxalobacteraceae",Family))

tax.40.8<-tax.40.7%>% 
  mutate(Family = ifelse(Genus== "Variovorax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.8.1<-tax.40.8%>% 
  mutate(Family = ifelse(Genus== "Advenella" & Family == "Burkholderiaceae","Alcaligenaceae",Family))

tax.40.9<-tax.40.8.1%>% 
  mutate(Family = ifelse(Genus== "Altererythrobacter" & Family == "Sphingomonadaceae","Erythrobacteraceae",Family))

tax.40.10<-tax.40.9%>% 
  mutate(Family = ifelse(Genus== "Mesorhizobium" & Family == "Rhizobiaceae","Phyllobacteriaceae",Family))

tax.40.11<-tax.40.10%>% 
  mutate(Family = ifelse(Genus== "Ahrensia" & Family == "Rhizobiaceae","Ahrensiaceae",Family))

tax.40.12<-tax.40.11%>% 
  mutate(Family = ifelse(Genus== "Pseudahrensia" & Family == "Rhizobiaceae","Ahrensiaceae",Family))

tax.40.13<-tax.40.12%>% 
  mutate(Family = ifelse(Genus== "Rubinisphaera" & Family == "Rubinisphaeraceae","Planctomycetaceae",Family))

tax.40.14<-tax.40.13%>% 
  mutate(Family = ifelse(Genus== "Exiguobacterium" & Family == "Bacillales_Incertae Sedis XII","Bacillales_Incertae_Sedis_XII",Family))

tax.40.15<-tax.40.14%>% 
  mutate(Family = ifelse(Genus== "Oceanirhabdus" & Family == "Clostridiaceae 1","Clostridiaceae_1",Family))

tax.40.16<-tax.40.15%>% 
  mutate(Family = ifelse(Genus== "Bacillus" & Family == "Bacillaceae 1","Bacillaceae",Family))

tax.40.17<-tax.40.16%>% 
  mutate(Family = ifelse(Genus== "Bacillus" & Family == "Bacillaceae_1","Bacillaceae",Family))

tax.40.18<-tax.40.17%>% 
  mutate(Family = ifelse(Genus== "Roseibacillus" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

tax.40.19<-tax.40.18%>% 
  mutate(Family = ifelse(Genus== "Litoribacillus" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

tax.40.20<-tax.40.19%>% 
  mutate(Family = ifelse(Genus== "Wandonia" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

tax.40.21<-tax.40.20%>% 
  mutate(Family = ifelse(Genus== "Marinoscillum" & Family == "Cyclobacteriaceae","Reichenbachiellaceae",Family))

tax.40.22<-tax.40.21%>% 
  mutate(Family = ifelse(Genus== "Egicoccus" & Family == "Nitriliruptoraceae","Egicoccaceae",Family))

tax.40.23<-tax.40.22%>% 
  mutate(Family = ifelse(Genus== "Persicirhabdus" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

tax.40.24<-tax.40.23%>% 
  mutate(Family = ifelse(Genus== "Aestuariicella" & Family == "Cellvibrionaceae","Alteromonadales_incertae_sedis",Family))

tax.40.25<-tax.40.24%>% 
  mutate(Family = ifelse(Genus== "Afipia" & Family == "Xanthobacteraceae","Bradyrhizobiaceae",Family))

tax.40.26<-tax.40.25%>% 
  mutate(Family = ifelse(Genus== "Agarivorans" & Family == "Psychromonadaceae","Alteromonadaceae",Family))

tax.40.27<-tax.40.26%>% 
  mutate(Family = ifelse(Genus== "Alkanibacter" & Family == "Solimonadaceae","Nevskiaceae",Family))

tax.40.28<-tax.40.27%>% 
  mutate(Family = ifelse(Genus== "Anderseniella" & Family == "Rhizobiales Incertae Sedis","Parvibaculaceae",Family))

tax.40.29<-tax.40.28%>% 
  mutate(Family = ifelse(Genus== "Aquabacterium" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.30<-tax.40.29%>% 
  mutate(Family = ifelse(Genus== "Aquincola" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.31<-tax.40.30%>% 
  mutate(Family = ifelse(Genus== "Bythopirellula" & Family == "Pirellulaceae","Lacipirellulaceae",Family))

tax.40.32<-tax.40.31%>% 
  mutate(Family = ifelse(Genus== "Catellicoccus" & Family == "Carnobacteriaceae","Enterococcaceae",Family))

tax.40.33<-tax.40.32%>% 
  mutate(Family = ifelse(Genus== "Comamonas" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.34<-tax.40.33%>% 
  mutate(Family = ifelse(Genus== "Corallomonas" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.40.35<-tax.40.34%>% 
  mutate(Family = ifelse(Genus== "Dechloromonas" & Family == "Rhodocyclaceae","Azonexaceae",Family))

tax.40.36<-tax.40.35%>% 
  mutate(Family = ifelse(Genus== "Defluviicoccus" & Family == "Rhodopirillaceae","Geminicoccaceae",Family))

tax.40.37<-tax.40.36%>% 
  mutate(Family = ifelse(Genus== "Dyadobacter" & Family == "Cytophagaceae","Spirosomaceae",Family))

tax.40.38<-tax.40.37%>% 
  mutate(Family = ifelse(Genus== "Eilatimonas" & Family == "Kordiimonadales_Incertae_Sedis","Temperatibacteraceae",Family))

tax.40.39<-tax.40.38%>% 
  mutate(Family = ifelse(Genus== "Ekhidna" & Family == "Cyclobacteriaceae","Reichenbachiellaceae",Family))

tax.40.40<-tax.40.39%>% 
  mutate(Family = ifelse(Genus== "Erythrobacter" & Family == "Sphingomonadaceae","Erythrobacteraceae",Family))

tax.40.41<-tax.40.40%>% 
  mutate(Family = ifelse(Genus== "Fabibacter" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

tax.40.42<-tax.40.41%>% 
  mutate(Family = ifelse(Genus== "Ferrimonas" & Family == "Shewanellaceae","Ferrimonadaceae",Family))

tax.40.43<-tax.40.42%>% 
  mutate(Family = ifelse(Genus== "Fulvitalea" & Family == "Cyclobacteriaceae","Ferrimonadaceae",Family))

tax.40.44<-tax.40.43%>% 
  mutate(Family = ifelse(Genus== "Gimesia" & Family == "Gimesiaceae","Planctomycetaceae",Family))

tax.40.45<-tax.40.44%>% 
  mutate(Family = ifelse(Genus== "Granulosicoccus" & Family == "Thiohalorhabdaceae","Granulosicoccaceae",Family))

tax.40.46<-tax.40.45%>% 
  mutate(Family = ifelse(Genus== "Halobacteriovorax" & Family == "Bacteriovoracaceae","Halobacteriovoraceae",Family))

tax.40.47<-tax.40.46%>% 
  mutate(Family = ifelse(Genus== "Haloferula" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

tax.40.48<-tax.40.47%>% 
  mutate(Family = ifelse(Genus== "Laribacter" & Family == "Aquaspirillaceae","Chromobacteriaceae",Family))

tax.40.49<-tax.40.48%>% 
  mutate(Family = ifelse(Genus== "Lewinella" & Family == "Saprospiraceae","Lewinellaceae",Family))

tax.40.50<-tax.40.49%>% 
  mutate(Family = ifelse(Genus== "Limnohabitans" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.51<-tax.40.50%>% 
  mutate(Family = ifelse(Genus== "Macromonas" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.52<-tax.40.51%>% 
  mutate(Family = ifelse(Genus== "Malikia" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.53<-tax.40.52%>% 
  mutate(Family = ifelse(Genus== "Marinomonas" & Family == "Marinomonadaceae","Oceanospirillaceae",Family))

tax.40.54<-tax.40.53%>% 
  mutate(Family = ifelse(Genus== "Methylibium" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.55<-tax.40.54%>% 
  mutate(Family = ifelse(Genus== "Motiliproteus" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

tax.40.56<-tax.40.55%>% 
  mutate(Family = ifelse(Genus== "Nisaea" & Family == "Nisaeaceae","Thalassobaculaceae",Family))

tax.40.57<-tax.40.56%>% 
  mutate(Family = ifelse(Genus== "Oceaniserpentilla" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

tax.40.58<-tax.40.57%>% 
  mutate(Family = ifelse(Genus== "Oceanococcus" & Family == "Solimonadaceae","Ectothiorhodospiraceae",Family))

tax.40.59<-tax.40.58%>% 
  mutate(Family = ifelse(Genus== "Oceanospirillum" & Family == "Halomonadaceae","Oceanospirillaceae",Family))

tax.40.60<-tax.40.59%>% 
  mutate(Family = ifelse(Genus== "Odoribacter" & Family == "Marinifilaceae","Odoribacteraceae",Family))

tax.40.61<-tax.40.60%>% 
  mutate(Family = ifelse(Genus== "Oleispira" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

tax.40.62<-tax.40.61%>% 
  mutate(Family = ifelse(Genus== "Ottowia" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.63<-tax.40.62%>% 
  mutate(Family = ifelse(Genus== "Pelagibius" & Family == "Kiloniellaceae","Rhodovibrionaceae",Family))

tax.40.64<-tax.40.63%>% 
  mutate(Family = ifelse(Genus== "Pelomonas" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.65<-tax.40.64%>% 
  mutate(Family = ifelse(Genus== "Pelosinus" & Family == "Veillonellaceae","Sporomusaceae",Family))

tax.40.66<-tax.40.65%>% 
  mutate(Family = ifelse(Genus== "Persicobacter" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

tax.40.67<-tax.40.66%>% 
  mutate(Family = ifelse(Genus== "Portibacter" & Family == "Saprospiraceae","Lewinellaceae",Family))

tax.40.68<-tax.40.67%>% 
  mutate(Family = ifelse(Genus== "Pseudarcicella" & Family == "Cytophagaceae","Spirosomaceae",Family))

tax.40.69<-tax.40.68%>% 
  mutate(Family = ifelse(Genus== "Pseudobacteriovorax" & Family == "Oligoflexaceae","Pseudobacteriovoracaceae",Family))

tax.40.70<-tax.40.69%>% 
  mutate(Family = ifelse(Genus== "Reichenbachiella" & Family == "Cyclobacteriaceae","Reichenbachiellaceae",Family))

tax.40.71<-tax.40.70%>% 
  mutate(Family = ifelse(Genus== "Rubricoccus" & Family == "Rhodothermaceae","Rubricoccaceae",Family))

tax.40.72<-tax.40.71%>% 
  mutate(Family = ifelse(Genus== "Rubrivirga" & Family == "Rhodothermaceae","Rubricoccaceae",Family))

tax.40.73<-tax.40.72%>% 
  mutate(Family = ifelse(Genus== "Rubrivivax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.74<-tax.40.73%>% 
  mutate(Family = ifelse(Genus== "Simplicispira" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.75<-tax.40.74%>% 
  mutate(Family = ifelse(Genus== "Spirosoma" & Family == "Cytophagaceae","Spirosomaceae",Family))

tax.40.76<-tax.40.75%>% 
  mutate(Family = ifelse(Genus== "Spongiimonas" & Family == "Flavobacteriaceae","Weeksellaceae",Family))

tax.40.77<-tax.40.76%>% 
  mutate(Family = ifelse(Genus== "Spongiispira" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

tax.40.78<-tax.40.77%>% 
  mutate(Family = ifelse(Genus== "Tunicatimonas" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

tax.40.79<-tax.40.78%>% 
  mutate(Family = ifelse(Genus== "Tunicatimonas" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

tax.40.80<-tax.40.79%>% 
  mutate(Family = ifelse(Genus== "Ideonella" & Family == "Burkholderiaceae","Comamonadaceae",Family))

tax.40.81<-tax.40.80%>% 
  mutate(Class= ifelse(Order== "Bacteriovoracales" & Class == "Deltaproteobacteria","Oligoflexia",Class))

tax.40.82<-tax.40.81%>% 
  mutate(Class= ifelse(Order== "Bacteroidales" & Class == "Sphingobacteriia","Bacteroidia",Class))

tax.40.83<-tax.40.82%>% 
  mutate(Class= ifelse(Order== "Calditrichales" & Class == "Deferribacteres","Calditrichia",Class))

tax.40.84<-tax.40.83%>% 
  mutate(Class= ifelse(Order== "Leptospirales" & Class == "Leptospirae","Spirochaetia",Class))

tax.40.85<-tax.40.84%>% 
  mutate(Class= ifelse(Order== "Nitrosomonadales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

tax.40.86<-tax.40.85%>% 
  mutate(Class= ifelse(Order== "Rhodocyclales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

tax.40.87<-tax.40.86%>% 
  mutate(Class= ifelse(Order== "Saprospirales" & Class == "Chitinophagia","Saprospiria",Class))

tax.40.88<-tax.40.87%>% 
  mutate(Class= ifelse(Order== "Burkholderiales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

tax.40.89<-tax.40.88%>% 
  mutate(Order= ifelse(Family== "Alteromonadales_incertae_sedis" & Order == "Cellvibrionales","Alteromonadales",Order))

tax.40.90<-tax.40.89%>% 
  mutate(Order= ifelse(Family== "Chromatiaceae" & Order == "Alteromonadales","Chromatiales",Order))

tax.40.91<-tax.40.90%>% 
  mutate(Order= ifelse(Family== "Chromobacteriaceae" & Order == "Betaproteobacteriales","Neisseriales",Order))

tax.40.92<-tax.40.91%>% 
  mutate(Order= ifelse(Family== "Desulfobacteraceae" & Order == "Desulfarculales","Desulfobacterales",Order))

tax.40.93<-tax.40.92%>% 
  mutate(Order= ifelse(Family== "Ectothiorhodospiraceae" & Order == "Salinisphaerales","Chromatiales",Order))

tax.40.94<-tax.40.93%>% 
  mutate(Order= ifelse(Family== "Egicoccaceae" & Order == "Nitriliruptorales","Egicoccales",Order))

tax.40.95<-tax.40.94%>% 
  mutate(Order= ifelse(Family== "Ferrimonadaceae" & Order == "Cytophagales","Alteromonadales",Order))

tax.40.96<-tax.40.95%>% 
  mutate(Order= ifelse(Family== "Granulosicoccaceae" & Order == "Thiohalorhabdales","Chromatiales",Order))

tax.40.97<-tax.40.96%>% 
  mutate(Order= ifelse(Family== "Nevskiaceae" & Order == "Salinisphaerales","Nevskiales",Order))

tax.40.98<-tax.40.97%>% 
  mutate(Order= ifelse(Family== "Piscirickettsiaceae" & Order == "Nitrosococcales","Thiotrichales",Order))

tax.40.99<-tax.40.98%>% 
  mutate(Order= ifelse(Family== "Piscirickettsiaceae" & Order == "Thiomicrospirales","Thiotrichales",Order))

tax.41<-tax.40.99%>% 
  mutate(Order= ifelse(Family== "Rhodospirillaceae" & Order == "Dongiales","Rhodospirillales",Order))

tax.41.1<-tax.41%>% 
  mutate(Order= ifelse(Family== "Rhodovibrionaceae" & Order == "Kiloniellales","Rhodospirillales",Order))

tax.41.2<-tax.41.1%>% 
  mutate(Order= ifelse(Family== "Selenomonadaceae" & Order == "Veillonellales","Selenomonadales",Order))

tax.41.3<-tax.41.2%>% 
  mutate(Order= ifelse(Family== "Sporomusaceae" & Order == "Veillonellales","Selenomonadales",Order))

tax.41.4<-tax.41.3%>% 
  mutate(Order= ifelse(Family== "Sterolibacteriaceae" & Order == "Rhodocyclales","Nitrosomonadales",Order))

tax.41.5<-tax.41.4%>% 
  mutate(Order= ifelse(Family== "Thalassobaculaceae" & Order == "Thalassobaculales","Rhodospirillales",Order))

tax.41.6<-tax.41.5%>% 
  mutate(Family= ifelse(Family== "Rubrobacteriaceae" ,"Rubrobacteraceae",Family))

tax.41.7<-tax.41.6%>% 
  mutate(Genus= ifelse(Genus== "Prevotella 9" ,"Prevotella",Genus))

tax.41.8<-tax.41.7%>% 
  mutate(Family= ifelse(Genus== "Fulvitalea" & Family == "Ferrimonadaceae","Flammeovirgaceae",Family))

tax.41.9<-tax.41.8%>% 
  mutate(Order= ifelse(Family== "Flammeovirgaceae" & Order == "Alteromonadales","Cytophagales",Order))

tax.41.10<-tax.41.9%>% 
  mutate(Class= ifelse(Order== "Alteromonadales" & Class == "Cytophagia","Gammaproteobacteria",Class))

tax.41.11<-tax.41.10%>% 
  mutate(Class= ifelse(Order== "Neisseriales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

tax.41.12<-tax.41.11%>% 
  mutate(Phylum= ifelse(Class== "Calditrichia" & Phylum == "Deferribacteres","Calditrichaeota",Phylum))

####Find duplicates####

TAX<-tax.41.12

#write.csv2(TAX,"Edited_tax_table.csv" )

TAX<-as.matrix(TAX)
otu<-read.csv2("Phyloseq_final_otu.csv", header=TRUE)
otu<-column_to_rownames(otu, "rowname")
mapfile = "Sample_data_all.csv"
map <- read.csv2(mapfile)
rownames(map) <- map$Sample_name
map<-map[,-1]              

tax_final<-tax_table(TAX)
otu_final<- otu_table(otu, taxa_are_rows = TRUE)
MAP<-sample_data(map)

data<-merge_phyloseq(tax_final,otu_final,MAP)
data

data.1 <- subset_taxa(data, Kingdom == "Bacteria" & Class != "Chloroplast" & Family != "mitochondria" &  Phylum != "Cyanobacteria/Chloroplast")
data.2<- subset_taxa(data.1, Phylum != "Bacteria_X")

data_phylum<-tax_glom(data.2, taxrank="Phylum")
data_class<-tax_glom(data.2, taxrank="Class")
data_order<-tax_glom(data.2, taxrank="Order")
data_family<-tax_glom(data.2, taxrank="Family")
data_genus<-tax_glom(data.2, taxrank="Genus")
data_species<-tax_glom(data.2, taxrank="Species")


#Phyla duplicates
data_phylum_t<-as.data.frame(tax_table(data_phylum))
dup_phyla<- data_phylum_t %>% filter(duplicated(data_phylum_t$Phylum)| duplicated(Phylum, fromLast = TRUE))

#Class duplicates
data_class_t<-as.data.frame(tax_table(data_class))
dup_class<- data_class_t %>% filter(duplicated(data_class_t$Class)| duplicated(Class, fromLast = TRUE))

#Order duplicates
data_order_t<-as.data.frame(tax_table(data_order))
dup_order<- data_order_t %>% filter(duplicated(data_order_t$Order)| duplicated(Order, fromLast = TRUE))

#Family duplicates
data_family_t<-as.data.frame(tax_table(data_family))
dup_Family<- data_family_t %>% filter(duplicated(data_family_t$Family)| duplicated(Family, fromLast = TRUE))

#Genus duplicates
data_genus_t<-as.data.frame(tax_table(data_genus))
dup_Genus<- data_genus_t %>% filter(duplicated(data_genus_t$Genus)| duplicated(Genus, fromLast = TRUE))

#Species duplicates
data_species_t<-as.data.frame(tax_table(data_species))
dup_species<- data_species_t %>% filter(duplicated(data_species_t$Species)| duplicated(Species, fromLast = TRUE))
