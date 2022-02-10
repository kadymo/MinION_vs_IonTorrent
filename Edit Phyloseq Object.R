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

ALL_TAX<-read.csv2("ALL_TAX.csv", header=TRUE, na.strings=c("","NA"))
#ALL_TAX <- read.delim("ALL_TAX.txt", header=TRUE)
ALL_TAX<- column_to_rownames(ALL_TAX, var="OTU")

Kingdom = ALL_TAX$Kingdom
Phylum = ALL_TAX$Phylum
Class = ALL_TAX$Class
Order= ALL_TAX$Order
Family=ALL_TAX$Family
Genus= ALL_TAX$Genus
Species= ALL_TAX$Species
variable1= "X"
variable2= "_X"


##Go column by column 

#Change unknown Kingdom to Unclassifed
ALL_TAX$Kingdom <- replace(ALL_TAX$Kingdom,ALL_TAX$Kingdom =="unknown", "Unclassified")

#Replace NA Kingdom to Unclassifed
ALL_TAX$Kingdom <-ALL_TAX$Kingdom %>% replace_na("Unclassified")

#Change all "_unclassified" to "_X" in Phylum
ALL_TAX$Phylum<-gsub("_unclassified", "_X", ALL_TAX$Phylum)

#Replace NA Phylum to Unclassifed
ALL_TAX$Phylum <-ALL_TAX$Phylum %>% replace_na("Unclassified")

#Replace Unknown_x Phylum to Unclassified_x
ALL_TAX$Phylum  <- replace(ALL_TAX$Phylum ,ALL_TAX$Phylum =="unknown_X", "Unclassified_X")

#Replace repeated unclassifed kingdom with  Unclassified_x
ALL_TAX$Phylum  <- replace(ALL_TAX$Phylum ,ALL_TAX$Phylum =="Unclassified", "Unclassified_X")

#making blanks NA
ALL_TAX.1<-ALL_TAX %>% mutate_if(is.character, list(~na_if(.,""))) 
#ALL_TAX.1 <- read.csv("data2.csv", header=T, na.strings=c("","NA"))

#Change all "_unclassified" to "_X" in Class
ALL_TAX1.1<-ALL_TAX.1 %>% 
  mutate(Class = if_else(str_ends(Class, "_Unclassified"),mapply(paste0, Phylum, variable2), Class))

#Replace any class with phylum that end in _X with _XX, with variable1 = X
ALL_TAX.2<-ALL_TAX1.1 %>% 
  mutate(Class = if_else(str_ends(Phylum, "_X"),mapply(paste0, Phylum, variable1), Class))

#Replace any class that are NA with phylum + _X
ALL_TAX.3<-ALL_TAX.2 %>% 
  mutate(Class = if_else(is.na(Class),mapply(paste0, Phylum, variable2), Class))

#Replace any class that have _cl with _X
ALL_TAX.4<-ALL_TAX.3 %>% 
  mutate(Class = if_else(str_ends(Class, "_cl"),mapply(paste0, Phylum, variable2), Class))

#Replace any class that have uncultured with Phylum_uncultured
ALL_TAX.4.1<-ALL_TAX.4 %>% 
  mutate(Class = if_else(str_detect(Class, 'uncultured'), mapply(paste0, Phylum, sep = "_", Class), Class))

#Replace any class that have marine metagenome with Phylum_marine metagenome
ALL_TAX.4.2<-ALL_TAX.4.1 %>% 
  mutate(Class = if_else(Class == "marine metagenome", mapply(paste0, Phylum, sep = "_", Class), Class))

ALL_TAX.4.3<-ALL_TAX.4.2 %>% 
  mutate(Class = if_else(Class == "metagenome", mapply(paste0, Phylum, sep = "_", Class), Class))

ALL_TAX.4.4<-ALL_TAX.4.3 %>% 
  mutate(Class = if_else(Class == "Unknown_Class", mapply(paste0, Phylum, variable2), Class))

ALL_TAX.4.5<-ALL_TAX.4.4 %>% 
  mutate(Class = if_else(Class == "uncultured bacterium", mapply(paste0, Phylum, sep = "_", Class), Class))

ALL_TAX.4.6<-ALL_TAX.4.5 %>% 
  mutate(Class = if_else(Class == "uncultured organism", mapply(paste0, Phylum, sep = "_", Class), Class))


#Remove _or from orders
ALL_TAX.5<-ALL_TAX.4.6 %>% 
  mutate_at("Order", str_replace, "_or", "")

#Replace any order with class that end in _X with _XX, with variable1 = X
ALL_TAX.6<-ALL_TAX.5 %>% 
  mutate(Order = if_else(str_ends(Class, "_X"),mapply(paste0, Class, variable1), Order))

#Replace any order with class that end in _XX with _XXX, with variable1 = X
ALL_TAX.7<-ALL_TAX.6 %>% 
  mutate(Order = if_else(str_ends(Class, "_XX"),mapply(paste0, Class, variable1), Order))

###chnage name of one order
ALL_TAX.8<-ALL_TAX.7 %>% 
  mutate_at("Order", str_replace, "44524", "11-24")

ALL_TAX.8.1<-ALL_TAX.8 %>% 
  mutate_at("Order", str_replace, "24-Nov", "11-24")

##replace orders that are teh same as class with class+ _X
ALL_TAX.9<-ALL_TAX.8.1 %>%
  mutate(Order = ifelse(Order==Class,mapply(paste0, Class, variable2), Order))

#Change all "_unclassified" to "_X" in Order
ALL_TAX.9$Order<-gsub("_unclassified", "_X", ALL_TAX.9$Order)

#Replace any order that are NA with class + _X
ALL_TAX.10<-ALL_TAX.9 %>% 
  mutate(Order = if_else(is.na(Order),mapply(paste0, Class, variable2),Order))

#Replace any order with class that end in _uncultured with , with variable1 = X
ALL_TAX.11<-ALL_TAX.10 %>% 
  mutate(Order = if_else(str_detect(Class, 'uncultured'),mapply(paste0, Class, variable2), Order))

#Replace any order with class that end in _uncultured bacterium with , with variable1 = X
ALL_TAX.12<-ALL_TAX.11 %>% 
  mutate(Order = if_else(str_ends(Class, "_marine metagenome"),mapply(paste0, Class, variable2), Order))

ALL_TAX.12.1<-ALL_TAX.12 %>% 
  mutate(Order = if_else(str_ends(Class, "_metagenome"),mapply(paste0, Class, variable2), Order))

ALL_TAX.12.2<-ALL_TAX.12.1 %>% 
  mutate(Order = if_else(Order == "Unknown_Order", mapply(paste0, Class, variable2), Order))

ALL_TAX.12.3<-ALL_TAX.12.2 %>% 
  mutate(Order  = if_else(Order  == "metagenome", mapply(paste0, Class, variable2), Order ))

ALL_TAX.12.4<-ALL_TAX.12.3%>% 
  mutate(Order = if_else(str_detect(Order, 'unidentified'),mapply(paste0, Class, variable2), Order))

ALL_TAX.12.5<-ALL_TAX.12.4 %>% 
  mutate(Order  = if_else(str_starts(Order, "Family"), mapply(paste0, Class, sep = "_", Order), Order ))

#Replace any Order that have uncultured with Class_uncultured
ALL_TAX.13<-ALL_TAX.12.5 %>% 
  mutate(Order = if_else(str_detect(Order, 'uncultured'), mapply(paste0, Class, sep = "_", Order), Order))

ALL_TAX.13.1<-ALL_TAX.13 %>% 
  mutate(Order = if_else(Order == "marine metagenome", mapply(paste0, Class, sep = "_", Order), Order))

#Remove _fa from families
ALL_TAX.14<-ALL_TAX.13.1 %>% 
  mutate_at("Family", str_replace, "_fa", "")

#Replace any family with order that end in _X with _XX, with variable1 = X
ALL_TAX.15<-ALL_TAX.14 %>% 
  mutate(Family = if_else(str_ends(Order, "_X"),mapply(paste0, Order, variable1), Family))

#Replace any family with order that end in _XX with _XXX, with variable1 = X
ALL_TAX.16<-ALL_TAX.15 %>% 
  mutate(Family = if_else(str_ends(Order, "_XX"),mapply(paste0, Order, variable1), Family))

#Replace any family with order that end in _XXX with _XXXX, with variable1 = X
ALL_TAX.17<-ALL_TAX.16 %>% 
  mutate(Family = if_else(str_ends(Order, "_XXX"),mapply(paste0, Order, variable1), Family))

##replace families that are the same as order with order + _X
ALL_TAX.18<-ALL_TAX.17 %>%
  mutate(Family = ifelse(Family ==Order,mapply(paste0, Order, variable2), Family))

#Change all "_unclassified" to "_X" in Family
ALL_TAX.18$Family<-gsub("_unclassified", "_X", ALL_TAX.18$Family)

#Replace any order that are NA with class + _X
ALL_TAX.19<-ALL_TAX.18 %>% 
  mutate(Family = if_else(is.na(Family),mapply(paste0, Order, variable2),Family))

#Replace any Family with class that end in _uncultured bacterium with , with variable1 = X
ALL_TAX.19.1<-ALL_TAX.19 %>% 
  mutate(Family = if_else(str_detect(Order, 'uncultured'),mapply(paste0, Order, variable2), Family))

#Replace any amily with class that end in _uncultured bacterium with , with variable1 = X
ALL_TAX.19.2<-ALL_TAX.19.1%>% 
  mutate(Family = if_else(str_ends(Order, "_marine metagenome"),mapply(paste0, Order, variable2), Family))

ALL_TAX.19.3<-ALL_TAX.19.2 %>% 
  mutate(Family  = if_else(Family  == "Unknown_Family", mapply(paste0, Order, variable2), Family ))

ALL_TAX.19.4<-ALL_TAX.19.3 %>% 
  mutate(Family  = if_else(Family  == "Unknown Family", mapply(paste0, Order, variable2), Family ))

ALL_TAX.19.5<-ALL_TAX.19.4 %>% 
  mutate(Family  = if_else(Family  == "metagenome", mapply(paste0, Order, variable2), Family ))

ALL_TAX.19.6<-ALL_TAX.19.5%>% 
  mutate(Family = if_else(str_detect(Family, 'unidentified'),mapply(paste0, Order, variable2), Family))

ALL_TAX.19.7<-ALL_TAX.19.6 %>% 
  mutate(Family  = if_else(Family  == "hydrothermal vent metagenome", mapply(paste0, Order, sep = "_", Family), Family ))

ALL_TAX.19.8<-ALL_TAX.19.7 %>% 
  mutate(Family  = if_else(str_detect(Order, "Family"), mapply(paste0, Order, variable2), Family ))

ALL_TAX.19.9<-ALL_TAX.19.8 %>% 
  mutate(Family  = if_else(str_starts(Family, "Family"), mapply(paste0, Order, sep = "_", Family), Family ))

#Replace any Order that have uncultured with Class_uncultured
ALL_TAX.20<-ALL_TAX.19.9 %>% 
  mutate(Family  = if_else(str_detect(Family , "^uncultured"), mapply(paste0, Order, sep = "_", Family ), Family ))

ALL_TAX.20.1<-ALL_TAX.20 %>% 
  mutate(Family  = if_else(Family  == "marine metagenome", mapply(paste0, Order, sep = "_", Family ), Family ))

#Remove _ge from genera
ALL_TAX.21<-ALL_TAX.20.1 %>% 
  mutate_at("Genus", str_replace, "_ge", "")

#Replace any Genus with order that end in _X with _XX, with variable1 = X
ALL_TAX.22<-ALL_TAX.21 %>% 
  mutate(Genus = if_else(str_ends(Family, "_X"),mapply(paste0, Family, variable1), Genus))

#Replace any Genus with Family that end in _XX with _XXX, with variable1 = X
ALL_TAX.23<-ALL_TAX.22 %>% 
  mutate(Genus = if_else(str_ends(Family, "_XX"),mapply(paste0, Family, variable1), Genus))

#Replace any Genus with Family that end in _XXX with _XXXX, with variable1 = X
ALL_TAX.24<-ALL_TAX.23 %>% 
  mutate(Genus = if_else(str_ends(Family, "_XXX"),mapply(paste0, Family, variable1), Genus))

#Replace any Genus with Family that end in _XXXX with _XXXXX, with variable1 = X
ALL_TAX.25<-ALL_TAX.24 %>% 
  mutate(Genus = if_else(str_ends(Family, "_XXXX"),mapply(paste0, Family, variable1), Genus))

##replace genera that are the same as Family with Family + _X
ALL_TAX.26<-ALL_TAX.25 %>%
  mutate(Genus = ifelse(Genus ==Family,mapply(paste0, Family, variable2), Genus))

#Change all "_unclassified" to "_X" in Genus
ALL_TAX.26$Genus<-gsub("_unclassified", "_X", ALL_TAX.26$Genus)

#Replace any Family that are NA with class + _X
ALL_TAX.27<-ALL_TAX.26 %>% 
  mutate(Genus = if_else(is.na(Genus),mapply(paste0, Family, variable2),Genus))

#Replace any genus with class that end in _uncultured bacterium with , with variable1 = X
ALL_TAX.27.1<-ALL_TAX.27 %>% 
  mutate(Genus = if_else(str_detect(Family, 'uncultured'),mapply(paste0, Family, variable2), Genus))

#Replace any genus with class that end in _uncultured bacterium with , with variable1 = X
ALL_TAX.27.2<-ALL_TAX.27.1%>% 
  mutate(Genus = if_else(str_ends(Family, "_marine metagenome"),mapply(paste0, Family, variable2), Genus))

#Replace any genus that have uncultured with Class_uncultured
ALL_TAX.28<-ALL_TAX.27.2 %>% 
  mutate(Genus  = if_else(str_detect(Genus , "^uncultured"), mapply(paste0, Family, sep = "_", Genus ), Genus ))

ALL_TAX.28.1<-ALL_TAX.28 %>% 
  mutate(Genus  = if_else(Genus  == "marine metagenome", mapply(paste0, Family, sep = "_", Genus ), Genus ))

ALL_TAX.28.2<-ALL_TAX.28.1 %>% 
  mutate(Genus = if_else(str_detect(Family, "Family"), mapply(paste0, Family, variable2), Genus ))

ALL_TAX.28.3<-ALL_TAX.28.2 %>% 
  mutate(Genus  = if_else(str_starts(Genus, "Family"), mapply(paste0, Family, sep = "_", Genus), Genus ))

ALL_TAX.28.4<-ALL_TAX.28.3 %>% 
  mutate(Genus  = if_else(Genus  == "metagenome", mapply(paste0, Family, variable2), Genus ))

ALL_TAX.28.5<-ALL_TAX.28.4 %>% 
  mutate(Genus  = if_else(str_detect(Family, "hydrothermal vent metagenome"), mapply(paste0, Family, variable2), Genus ))

ALL_TAX.28.6<-ALL_TAX.28.5 %>% 
  mutate(Genus  = if_else(str_starts(Genus, "hydrothermal vent metagenome"), mapply(paste0, Family, sep = "_", Genus), Genus ))

ALL_TAX.28.7<-ALL_TAX.28.6 %>% 
  mutate(Genus  = if_else(Genus  == "Ambiguous_taxa", mapply(paste0, Family, sep = "_", Genus), Genus ))


#Remove _sp from species
ALL_TAX.29<-ALL_TAX.28.7 %>% 
  mutate_at("Species", str_replace, "_sp", "")

#Replace any Species with order that end in _X with _XX, with variable1 = X
ALL_TAX.30<-ALL_TAX.29 %>% 
  mutate(Species = if_else(str_ends(Genus, "_X"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XX with _XXX, with variable1 = X
ALL_TAX.31<-ALL_TAX.30 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XX"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XXX with _XXXX, with variable1 = X
ALL_TAX.32<-ALL_TAX.31 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XXX"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XXXX with _XXXXX, with variable1 = X
ALL_TAX.33<-ALL_TAX.32 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XXXX"),mapply(paste0, Genus, variable1), Species))

#Replace any Species with Genus that end in _XXXXX with _XXXXXX, with variable1 = X
ALL_TAX.34<-ALL_TAX.33 %>% 
  mutate(Species = if_else(str_ends(Genus, "_XXXXX"),mapply(paste0, Genus, variable1), Species))

##replace species that are the same as Genus with Genus + _X
ALL_TAX.35<-ALL_TAX.34 %>%
  mutate(Species = ifelse(Species ==Genus,mapply(paste0, Genus, variable2), Species))

#Change all "_unclassified" to "_X" in Species
ALL_TAX.35$Species<-gsub("_unclassified", "_X", ALL_TAX.35$Species)

#Replace any Species that are NA with class + _X
ALL_TAX.36<-ALL_TAX.35 %>% 
  mutate(Species = if_else(is.na(Species),mapply(paste0, Genus, variable2),Species))

#Replace any Species with class that end in _uncultured bacterium with , with variable1 = X
ALL_TAX.36.1<-ALL_TAX.36 %>% 
  mutate(Species = if_else(str_detect(Genus, 'uncultured'),mapply(paste0, Genus, variable2), Species))

#Replace any Species with class that end in _uncultured bacterium with , with variable1 = X
ALL_TAX.36.2<-ALL_TAX.36.1%>% 
  mutate(Species = if_else(str_ends(Genus, "_marine metagenome"),mapply(paste0, Genus, variable2), Species))

#Replace any Species that have uncultured with Class_uncultured
ALL_TAX.37<-ALL_TAX.36.2 %>% 
  mutate(Species  = if_else(str_detect(Species , "^uncultured"), mapply(paste0, Genus, sep = "_", Species ), Species ))

ALL_TAX.37.1<-ALL_TAX.37 %>% 
  mutate(Species  = if_else(str_ends(Species, 'metagenome'), mapply(paste0, Genus, sep = "_", Species ), Species ))

ALL_TAX.37.2<-ALL_TAX.37.1 %>% 
  mutate(Species = if_else(str_detect(Genus, "Family"), mapply(paste0, Genus, variable2), Species  ))

ALL_TAX.37.3<-ALL_TAX.37.2 %>% 
  mutate(Species  = if_else(str_starts(Species, "Family"), mapply(paste0, Genus, sep = "_", Species), Species ))

ALL_TAX.37.4<-ALL_TAX.37.3 %>% 
  mutate(Species = if_else(str_detect(Genus, "hydrothermal vent metagenome"), mapply(paste0, Genus, variable2), Species))

ALL_TAX.37.5<-ALL_TAX.37.4 %>% 
  mutate(Species  = if_else(str_starts(Species, "hydrothermal vent metagenome"), mapply(paste0, Genus, sep = "_", Species), Species ))

ALL_TAX.37.6<-ALL_TAX.37.5 %>% 
  mutate(Species  = if_else(str_detect(Genus, "Ambiguous_taxa"), mapply(paste0, Genus, variable2), Species ))

ALL_TAX.37.7<-ALL_TAX.37.6 %>% 
  mutate(Species  = if_else(Species == "Ambiguous_taxa", mapply(paste0, Genus, sep = "_", Species), Species ))

#Concatenate Genus+Species if species names starts lowercase except "uncultured"
ALL_TAX.38<-ALL_TAX.37.7 %>% 
  mutate(Species = if_else(str_detect(Species, '^[:lower:]'),mapply(paste0, Genus, sep = " ", Species),Species))


###replacing dupes
ALL_TAX.39<-ALL_TAX.38 %>% 
  mutate(Phylum = ifelse(Class== "Campylobacteria" & Phylum == "Epsilonbacteraeota","Campilobacterota",Phylum))

ALL_TAX.39.1<-ALL_TAX.39 %>% 
  mutate(Phylum = ifelse(Class== "Chlorobia" & Phylum == "Bacteroidetes","Chlorobi",Phylum))

ALL_TAX.39.2<-ALL_TAX.39.1 %>% 
  mutate(Phylum = ifelse(Class== "Ignavibacteria" & Phylum == "Ignavibacteriae","Bacteroidetes",Phylum))

ALL_TAX.39.3<-ALL_TAX.39.2 %>% 
  mutate(Phylum = ifelse(Class== "Rhodothermia" & Phylum == "Bacteroidetes","Rhodothermaeota",Phylum))

ALL_TAX.39.4<-ALL_TAX.39.3 %>% 
  mutate(Class = ifelse(Order== "Caldilineales" & Class  == "Anaerolineae","Caldilineae",Class ))

ALL_TAX.39.5<-ALL_TAX.39.4 %>% 
  mutate(Class = ifelse(Order== "Opitutales" & Class  == "Verrucomicrobiae","Opitutae",Class ))   

ALL_TAX.39.6<-ALL_TAX.39.5 %>% 
  mutate(Class = ifelse(Order== "Sphingobacteriales" & Class  == "Bacteroidia","Sphingobacteriia",Class ))   

ALL_TAX.39.7<-ALL_TAX.39.6 %>% 
  mutate(Class = ifelse(Order== "Chitinophagales" & Class  == "Bacteroidia","Chitinophagia",Class ))   

ALL_TAX.39.8<-ALL_TAX.39.7 %>% 
  mutate(Class = ifelse(Order== "Cytophagales" & Class  == "Bacteroidia","Cytophagia",Class ))

ALL_TAX.39.9<-ALL_TAX.39.8 %>% 
  mutate(Class = ifelse(Order== "Flavobacteriales" & Class  == "Flavobacteriia","Bacteroidia",Class ))

ALL_TAX.39.10<-ALL_TAX.39.9 %>% 
  mutate(Class = ifelse(Order== "Chlamydiales" & Class  == "Chlamydiia","Chlamydiae",Class ))

ALL_TAX.39.11<-ALL_TAX.39.10 %>% 
  mutate(Class = ifelse(Order== "Oligoflexales" & Class  == "Deltaproteobacteria","Oligoflexia",Class ))

ALL_TAX.39.12<-ALL_TAX.39.11 %>% 
  mutate(Class = ifelse(Order== "Balneolales" & Class  == "Rhodothermia","Balneolia",Class ))

ALL_TAX.39.13<-ALL_TAX.39.12 %>% 
  mutate(Class = ifelse(Order== "Blastocatellales" & Class  == "Blastocatellia","Blastocatellia_(Subgroup_4)",Class ))

ALL_TAX.39.14<-ALL_TAX.39.13 %>% 
  mutate(Class = ifelse(Class  == "Blastocatellia (Subgroup 4)","Blastocatellia_(Subgroup_4)",Class ))

ALL_TAX.39.15<-ALL_TAX.39.14 %>% 
  mutate(Order = ifelse(Family== "Campylobacteria_X" & Order  == "Blastocatellia","Blastocatellia_(Subgroup_4)",Order ))

ALL_TAX.39.16<-ALL_TAX.39.15 %>% 
  mutate(Phylum = ifelse(Class== "Balneolia" & Phylum == "Rhodothermaeota","Balneolaeota",Phylum))

ALL_TAX.39.17<-ALL_TAX.39.16 %>% 
  mutate(Order = ifelse(Family== "Nocardiaceae" & Order == "Corynebacteriales","Mycobacteriales",Order))

ALL_TAX.39.18<-ALL_TAX.39.17 %>% 
  mutate(Order = ifelse(Family== "Trueperaceae" & Order == "Deinococcales","Trueperales",Order))

ALL_TAX.39.19<-ALL_TAX.39.18 %>% 
  mutate(Order= ifelse(Family== "Acidaminococcaceae" & Order == "Selenomonadales","Acidaminococcales",Order))

ALL_TAX.39.20<-ALL_TAX.39.19 %>% 
  mutate(Order = ifelse(Family== "Ilumatobacteraceae" & Order == "Microtrichales","Acidimicrobiales",Order))

ALL_TAX.39.21<-ALL_TAX.39.20 %>% 
  mutate(Order = ifelse(Family== "Iamiaceae" & Order == "Microtrichales","Acidimicrobiales",Order))

ALL_TAX.39.22<-ALL_TAX.39.21 %>% 
  mutate(Order = ifelse(Family== "Bacteriovoracaceae" & Order == "Bdellovibrionales","Bacteriovoracales",Order))

ALL_TAX.39.22<-ALL_TAX.39.21 %>% 
  mutate(Order = ifelse(Family== "Lentimicrobiaceae" & Order == "Sphingobacteriales","Bacteroidales",Order))

ALL_TAX.39.23<-ALL_TAX.39.22 %>% 
  mutate(Order = ifelse(Family== "Leptospiraceae" & Order == "Spirochaetales","Leptospirales",Order))

ALL_TAX.39.24<-ALL_TAX.39.23 %>% 
  mutate(Order = ifelse(Family== "Rhodocyclaceae" & Order == "Betaproteobacteriales","Rhodocyclales",Order))

ALL_TAX.39.25<-ALL_TAX.39.24 %>% 
  mutate(Order = ifelse(Family== "Burkholderiaceae" & Order == "Betaproteobacteriales","Burkholderiales",Order))

ALL_TAX.39.26<-ALL_TAX.39.25 %>% 
  mutate(Order = ifelse(Family== "Methylophilaceae" & Order == "Betaproteobacteriales","Nitrosomonadales",Order))

ALL_TAX.39.27<-ALL_TAX.39.26 %>% 
  mutate(Order = ifelse(Family== "Gallionellaceae" & Order == "Betaproteobacteriales","Nitrosomonadales",Order))

ALL_TAX.39.28<-ALL_TAX.39.27 %>% 
  mutate(Order = ifelse(Family== "Nitrosomonadaceae" & Order == "Betaproteobacteriales","Nitrosomonadales",Order))

ALL_TAX.39.29<-ALL_TAX.39.28 %>% 
  mutate(Order = ifelse(Family== "Calditrichaceae" & Order == "Deferribacterales","Calditrichales",Order))

ALL_TAX.39.30<-ALL_TAX.39.29 %>% 
  mutate(Order = ifelse(Family== "Parvularculaceae" & Order == "Caulobacterales","Parvularculales",Order))

ALL_TAX.39.31<-ALL_TAX.39.30 %>% 
  mutate(Order = ifelse(Family== "Saprospiraceae" & Order == "Chitinophagales","Saprospirales",Order))

ALL_TAX.39.32<-ALL_TAX.39.31 %>% 
  mutate(Order = ifelse(Family== "Woeseiaceae" & Order == "Steroidobacterales","Chromatiales",Order))

ALL_TAX.39.33<-ALL_TAX.39.32 %>% 
  mutate(Order = ifelse(Family== "Thalassobaculaceae" & Order == "Thalassobaculales","Rhodospirillales",Order))

ALL_TAX.39.34<-ALL_TAX.39.33 %>% 
  mutate(Order = ifelse(Family== "Ectothiorhodospiraceae" & Order == "Ectothiorhodospirales","Chromatiales",Order))

ALL_TAX.39.35<-ALL_TAX.39.34 %>% 
  mutate(Order = ifelse(Family== "Enterobacteriaceae" & Order == "Enterobacteriales","Enterobacterales",Order))

ALL_TAX.39.36<-ALL_TAX.39.35 %>% 
  mutate(Order = ifelse(Family== "Francisellaceae" & Order == "Francisellales","Thiotrichales",Order))

ALL_TAX.39.37<-ALL_TAX.39.36 %>% 
  mutate(Order = ifelse(Family== "Corynebacteriaceae" & Order == "Corynebacteriales","Mycobacteriales",Order))

ALL_TAX.39.38<-ALL_TAX.39.37 %>% 
  mutate(Order = ifelse(Family== "Mycobacteriaceae" & Order == "Corynebacteriales","Mycobacteriales",Order))

ALL_TAX.39.39<-ALL_TAX.39.38 %>% 
  mutate(Order = ifelse(Family== "Gemmataceae" & Order == "Planctomycetales","Gemmatales",Order))

ALL_TAX.39.40<-ALL_TAX.39.39 %>% 
  mutate(Order = ifelse(Family== "Isosphaeraceae" & Order == "Planctomycetales","Isosphaerales",Order))

ALL_TAX.39.41<-ALL_TAX.39.40 %>% 
  mutate(Order = ifelse(Family== "Kiloniellaceae" & Order == "Rhodovibrionales","Kiloniellales",Order))

ALL_TAX.39.42<-ALL_TAX.39.41 %>% 
  mutate(Order = ifelse(Family== "Coxiellaceae" & Order == "Coxiellales","Legionellales",Order))

ALL_TAX.39.43<-ALL_TAX.39.42 %>% 
  mutate(Order = ifelse(Family== "Leptospiraceae" & Order == "Spirochaetales","Leptospirales",Order))

ALL_TAX.39.44<-ALL_TAX.39.43%>% 
  mutate(Order = ifelse(Family== "Steroidobacteraceae" & Order == "Steroidobacterales","Nevskiales",Order))

ALL_TAX.39.45<-ALL_TAX.39.44%>% 
  mutate(Order = ifelse(Family== "Veillonellaceae" & Order == "Selenomonadales","Veillonellales",Order))

ALL_TAX.39.46<-ALL_TAX.39.45%>% 
  mutate(Order = ifelse(Family== "Reyranellaceae" & Order == "Reyranellales","Rhodospirillales",Order))

ALL_TAX.39.47<-ALL_TAX.39.46%>% 
  mutate(Order = ifelse(Family== "Acetobacteraceae" & Order == "Acetobacterales","Rhodospirillales",Order))

ALL_TAX.39.48<-ALL_TAX.39.47%>% 
  mutate(Order = ifelse(Family== "Puniceicoccaceae" & Order == "Opitutales","Puniceicoccales",Order))

ALL_TAX.39.49<-ALL_TAX.39.48%>% 
  mutate(Order = ifelse(Family== "Acetobacteraceae" & Order == "Acetobacterales","Rhodospirillales",Order))

ALL_TAX.39.50<-ALL_TAX.39.49%>% 
  mutate(Order = ifelse(Family== "Geminicoccaceae" & Order == "Tistrellales","Rhodospirillales",Order))

ALL_TAX.39.51<-ALL_TAX.39.50%>% 
  mutate(Order = ifelse(Family== "Bacteriovoracaceae" & Order == "Bdellovibrionales","Bacteriovoracales",Order))

ALL_TAX.39.52<-ALL_TAX.39.51%>% 
  mutate(Order = ifelse(Family== "Parvibaculaceae" & Order == "Parvibaculales","Rhizobiales",Order))

ALL_TAX.39.53<-ALL_TAX.39.52%>% 
  mutate(Order = ifelse(Family== "Piscirickettsiaceae" & Order == "Piscirickettsiales","Thiotrichales",Order))

ALL_TAX.39.54<-ALL_TAX.39.53%>% 
  mutate(Family = ifelse(Genus== "Bryobacter" & Family == "Solibacteraceae (Subgroup 3)","Solibacteraceae_(Subgroup_3)",Family))

ALL_TAX.39.55<-ALL_TAX.39.54%>% 
  mutate(Family = ifelse(Genus== "Parabacteroides" & Family == "Porphyromonadaceae","Tannerellaceae",Family))

ALL_TAX.39.56<-ALL_TAX.39.55%>% 
  mutate(Family = ifelse(Genus== "Fulvivirga" & Family == "Cyclobacteriaceae","Fulvivirgaceae",Family))

ALL_TAX.39.57<-ALL_TAX.39.56%>% 
  mutate(Family = ifelse(Genus== "Marivirga" & Family == "Cyclobacteriaceae","Marivirgaceae",Family))

ALL_TAX.39.58<-ALL_TAX.39.57%>% 
  mutate(Family = ifelse(Genus== "Brumimicrobium" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

ALL_TAX.39.59<-ALL_TAX.39.58%>% 
  mutate(Family = ifelse(Genus== "Crocinitomix" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

ALL_TAX.39.60<-ALL_TAX.39.59%>% 
  mutate(Family = ifelse(Genus== "Leuconostoc" & Family == "Lactobacillaceae","Leuconostocaceae",Family))

ALL_TAX.39.61<-ALL_TAX.39.60%>% 
  mutate(Family = ifelse(Genus== "Mitsuokella" & Family == "Veillonellaceae","Selenomonadaceae",Family))

ALL_TAX.39.62<-ALL_TAX.39.61%>% 
  mutate(Family = ifelse(Genus== "Planctomicrobium" & Family == "Rubinisphaeraceae","Planctomycetaceae",Family))

ALL_TAX.39.63<-ALL_TAX.39.62%>% 
  mutate(Family = ifelse(Genus== "Methylobacterium" & Family == "Beijerinckiaceae","Methylobacteriaceae",Family))

ALL_TAX.39.64<-ALL_TAX.39.63%>% 
  mutate(Family = ifelse(Genus== "Nordella" & Family == "Rhizobiales Incertae Sedis","Rhizobiales_Incertae_Sedis",Family))

ALL_TAX.39.65<-ALL_TAX.39.64%>% 
  mutate(Family = ifelse(Genus== "Phreatobacter" & Family == "Rhizobiales_Incertae_Sedis","Phreatobacteraceae",Family))

ALL_TAX.39.66<-ALL_TAX.39.65%>% 
  mutate(Family = ifelse(Genus== "Phreatobacter" & Family == "Rhizobiales Incertae Sedis","Phreatobacteraceae",Family))

ALL_TAX.39.67<-ALL_TAX.39.66%>% 
  mutate(Family = ifelse(Genus== "Bradyrhizobium" & Family == "Xanthobacteraceae","Bradyrhizobiaceae",Family))

ALL_TAX.39.68<-ALL_TAX.39.67%>% 
  mutate(Family = ifelse(Genus== "Magnetospira" & Family == "Magnetospiraceae","Thalassospiraceae",Family))

ALL_TAX.39.69<-ALL_TAX.39.68%>% 
  mutate(Family = ifelse(Genus== "Marinobacter" & Family == "Marinobacteraceae","Alteromonadaceae",Family))

ALL_TAX.39.70<-ALL_TAX.39.69%>% 
  mutate(Family = ifelse(Genus== "Psychrobium" & Family == "Alteromonadales_incertae_sedis","Shewanellaceae",Family))

ALL_TAX.39.71<-ALL_TAX.39.70%>% 
  mutate(Family = ifelse(Genus== "Alcaligenes" & Family == "Burkholderiaceae","Alcaligenaceae",Family))

ALL_TAX.39.72<-ALL_TAX.39.71%>% 
  mutate(Family = ifelse(Genus== "Rhodoferax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.39.73<-ALL_TAX.39.72%>% 
  mutate(Family = ifelse(Genus== "Xylophilus" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.39.74<-ALL_TAX.39.73%>% 
  mutate(Family = ifelse(Genus== "Methyloversatilis" & Family == "Rhodocyclaceae","Sterolibacteriaceae",Family))

ALL_TAX.39.75<-ALL_TAX.39.74%>% 
  mutate(Family = ifelse(Genus== "Zoogloea" & Family == "Rhodocyclaceae","Zoogloeaceae",Family))

ALL_TAX.39.76<-ALL_TAX.39.75%>% 
  mutate(Family = ifelse(Genus== "Providencia" & Family == "Enterobacteriaceae","Morganellaceae",Family))

ALL_TAX.39.77<-ALL_TAX.39.76%>% 
  mutate(Family = ifelse(Genus== "Amphritea" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.39.78<-ALL_TAX.39.77%>% 
  mutate(Family = ifelse(Genus== "Profundimonas" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.39.79<-ALL_TAX.39.78%>% 
  mutate(Family = ifelse(Genus== "Oleibacter" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

ALL_TAX.39.80<-ALL_TAX.39.79%>% 
  mutate(Family = ifelse(Genus== "Nevskia" & Family == "Solimonadaceae","Nevskiaceae",Family))

ALL_TAX.39.81<-ALL_TAX.39.80%>% 
  mutate(Family = ifelse(Genus== "Luteolibacter" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

ALL_TAX.39.82<-ALL_TAX.39.81%>% 
  mutate(Family = ifelse(Genus== "Mitsuokella" & Family == "Veillonellaceae","Selenomonadaceae",Family))

ALL_TAX.39.83<-ALL_TAX.39.82%>% 
  mutate(Family = ifelse(Genus== "Selenomonas" & Family == "Veillonellaceae","Selenomonadaceae",Family))

ALL_TAX.39.84<-ALL_TAX.39.83%>% 
  mutate(Family = ifelse(Genus== "Acidovorax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.39.85<-ALL_TAX.39.84%>% 
  mutate(Family = ifelse(Genus== "Pseudacidovorax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.39.86<-ALL_TAX.39.85%>% 
  mutate(Family = ifelse(Genus== "Curvibacter" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.39.87<-ALL_TAX.39.86%>% 
  mutate(Family = ifelse(Genus== "Hydrogenophaga" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.39.88<-ALL_TAX.39.87%>% 
  mutate(Family = ifelse(Genus== "Thauera" & Family == "Rhodocyclaceae","Zoogloeaceae",Family))

ALL_TAX.39.89<-ALL_TAX.39.88%>% 
  mutate(Family = ifelse(Genus== "Neiella" & Family == "Alteromonadaceae","Alteromonadales_incertae_sedis",Family))

ALL_TAX.39.90<-ALL_TAX.39.89%>% 
  mutate(Family = ifelse(Genus== "Macellibacteroides" & Family == "Tannerellaceae","Porphyromonadaceae",Family))

ALL_TAX.39.91<-ALL_TAX.39.90%>% 
  mutate(Family = ifelse(Genus== "Salinirepens" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

ALL_TAX.39.92<-ALL_TAX.39.91%>% 
  mutate(Family = ifelse(Genus== "Dongia" & Family == "Dongiaceae","Rhodospirillaceae",Family))

ALL_TAX.39.93<-ALL_TAX.39.92%>% 
  mutate(Family = ifelse(Genus== "Novosphingobium" & Family == "Sphingomonadaceae","Erythrobacteraceae",Family))

ALL_TAX.39.94<-ALL_TAX.39.93%>% 
  mutate(Family = ifelse(Genus== "Sphaerochaeta" & Family == "Spirochaetaceae","Sphaerochaetaceae",Family))

ALL_TAX.39.95<-ALL_TAX.39.94%>% 
  mutate(Family = ifelse(Genus== "Thiomicrorhabdus" & Family == "Thiomicrospiraceae","Piscirickettsiaceae",Family))

ALL_TAX.39.96<-ALL_TAX.39.95%>% 
  mutate(Family = ifelse(Genus== "Methylophaga" & Family == "Methylophagaceae","Piscirickettsiaceae",Family))

ALL_TAX.39.97<-ALL_TAX.39.96%>% 
  mutate(Family = ifelse(Genus== "Pseudohongiella" & Family == "Pseudohongiellaceae","Oceanospirillales_incertae_sedis",Family))

ALL_TAX.39.98<-ALL_TAX.39.97%>% 
  mutate(Family = ifelse(Genus== "Reinekea" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

ALL_TAX.39.99<-ALL_TAX.39.98%>% 
  mutate(Family = ifelse(Genus== "Nitrincola" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.40<-ALL_TAX.39.99%>% 
  mutate(Family = ifelse(Genus== "Neptunomonas" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.40.1<-ALL_TAX.40%>% 
  mutate(Family = ifelse(Genus== "Neptuniibacter" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.40.2<-ALL_TAX.40.1%>% 
  mutate(Family = ifelse(Genus== "Marinobacterium" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.40.3<-ALL_TAX.40.2%>% 
  mutate(Family = ifelse(Genus== "Pantoea" & Family == "Enterobacteriaceae","Erwiniaceae",Family))

ALL_TAX.40.4<-ALL_TAX.40.3%>% 
  mutate(Family = ifelse(Genus== "Rheinheimera" & Family == "Alteromonadaceae","Chromatiaceae",Family))

ALL_TAX.40.4.1<-ALL_TAX.40.4%>% 
  mutate(Family = ifelse(Genus== "Desulfatiglans" & Family == "Desulfarculaceae","Desulfobacteraceae",Family))

ALL_TAX.40.5<-ALL_TAX.40.4.1%>% 
  mutate(Family = ifelse(Genus== "Ferribacterium" & Family == "Rhodocyclaceae","Azonexaceae",Family))

ALL_TAX.40.6<-ALL_TAX.40.5%>% 
  mutate(Family = ifelse(Genus== "Noviherbaspirillum" & Family == "Burkholderiaceae","Oxalobacteraceae",Family))

ALL_TAX.40.7<-ALL_TAX.40.6%>% 
  mutate(Family = ifelse(Genus== "Massilia" & Family == "Burkholderiaceae","Oxalobacteraceae",Family))

ALL_TAX.40.8<-ALL_TAX.40.7%>% 
  mutate(Family = ifelse(Genus== "Variovorax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.8.1<-ALL_TAX.40.8%>% 
  mutate(Family = ifelse(Genus== "Advenella" & Family == "Burkholderiaceae","Alcaligenaceae",Family))

ALL_TAX.40.9<-ALL_TAX.40.8.1%>% 
  mutate(Family = ifelse(Genus== "Altererythrobacter" & Family == "Sphingomonadaceae","Erythrobacteraceae",Family))

ALL_TAX.40.10<-ALL_TAX.40.9%>% 
  mutate(Family = ifelse(Genus== "Mesorhizobium" & Family == "Rhizobiaceae","Phyllobacteriaceae",Family))

ALL_TAX.40.11<-ALL_TAX.40.10%>% 
  mutate(Family = ifelse(Genus== "Ahrensia" & Family == "Rhizobiaceae","Ahrensiaceae",Family))

ALL_TAX.40.12<-ALL_TAX.40.11%>% 
  mutate(Family = ifelse(Genus== "Pseudahrensia" & Family == "Rhizobiaceae","Ahrensiaceae",Family))

ALL_TAX.40.13<-ALL_TAX.40.12%>% 
  mutate(Family = ifelse(Genus== "Rubinisphaera" & Family == "Rubinisphaeraceae","Planctomycetaceae",Family))

ALL_TAX.40.14<-ALL_TAX.40.13%>% 
  mutate(Family = ifelse(Genus== "Exiguobacterium" & Family == "Bacillales_Incertae Sedis XII","Bacillales_Incertae_Sedis_XII",Family))

ALL_TAX.40.15<-ALL_TAX.40.14%>% 
  mutate(Family = ifelse(Genus== "Oceanirhabdus" & Family == "Clostridiaceae 1","Clostridiaceae_1",Family))

ALL_TAX.40.16<-ALL_TAX.40.15%>% 
  mutate(Family = ifelse(Genus== "Bacillus" & Family == "Bacillaceae 1","Bacillaceae",Family))

ALL_TAX.40.17<-ALL_TAX.40.16%>% 
  mutate(Family = ifelse(Genus== "Bacillus" & Family == "Bacillaceae_1","Bacillaceae",Family))

ALL_TAX.40.18<-ALL_TAX.40.17%>% 
  mutate(Family = ifelse(Genus== "Roseibacillus" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

ALL_TAX.40.19<-ALL_TAX.40.18%>% 
  mutate(Family = ifelse(Genus== "Litoribacillus" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

ALL_TAX.40.20<-ALL_TAX.40.19%>% 
  mutate(Family = ifelse(Genus== "Wandonia" & Family == "Cryomorphaceae","Crocinitomicaceae",Family))

ALL_TAX.40.21<-ALL_TAX.40.20%>% 
  mutate(Family = ifelse(Genus== "Marinoscillum" & Family == "Cyclobacteriaceae","Reichenbachiellaceae",Family))

ALL_TAX.40.22<-ALL_TAX.40.21%>% 
  mutate(Family = ifelse(Genus== "Egicoccus" & Family == "Nitriliruptoraceae","Egicoccaceae",Family))

ALL_TAX.40.23<-ALL_TAX.40.22%>% 
  mutate(Family = ifelse(Genus== "Persicirhabdus" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

ALL_TAX.40.24<-ALL_TAX.40.23%>% 
  mutate(Family = ifelse(Genus== "Aestuariicella" & Family == "Cellvibrionaceae","Alteromonadales_incertae_sedis",Family))

ALL_TAX.40.25<-ALL_TAX.40.24%>% 
  mutate(Family = ifelse(Genus== "Afipia" & Family == "Xanthobacteraceae","Bradyrhizobiaceae",Family))

ALL_TAX.40.26<-ALL_TAX.40.25%>% 
  mutate(Family = ifelse(Genus== "Agarivorans" & Family == "Psychromonadaceae","Alteromonadaceae",Family))

ALL_TAX.40.27<-ALL_TAX.40.26%>% 
  mutate(Family = ifelse(Genus== "Alkanibacter" & Family == "Solimonadaceae","Nevskiaceae",Family))

ALL_TAX.40.28<-ALL_TAX.40.27%>% 
  mutate(Family = ifelse(Genus== "Anderseniella" & Family == "Rhizobiales Incertae Sedis","Parvibaculaceae",Family))

ALL_TAX.40.29<-ALL_TAX.40.28%>% 
  mutate(Family = ifelse(Genus== "Aquabacterium" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.30<-ALL_TAX.40.29%>% 
  mutate(Family = ifelse(Genus== "Aquincola" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.31<-ALL_TAX.40.30%>% 
  mutate(Family = ifelse(Genus== "Bythopirellula" & Family == "Pirellulaceae","Lacipirellulaceae",Family))

ALL_TAX.40.32<-ALL_TAX.40.31%>% 
  mutate(Family = ifelse(Genus== "Catellicoccus" & Family == "Carnobacteriaceae","Enterococcaceae",Family))

ALL_TAX.40.33<-ALL_TAX.40.32%>% 
  mutate(Family = ifelse(Genus== "Comamonas" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.34<-ALL_TAX.40.33%>% 
  mutate(Family = ifelse(Genus== "Corallomonas" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.40.35<-ALL_TAX.40.34%>% 
  mutate(Family = ifelse(Genus== "Dechloromonas" & Family == "Rhodocyclaceae","Azonexaceae",Family))

ALL_TAX.40.36<-ALL_TAX.40.35%>% 
  mutate(Family = ifelse(Genus== "Defluviicoccus" & Family == "Rhodopirillaceae","Geminicoccaceae",Family))

ALL_TAX.40.37<-ALL_TAX.40.36%>% 
  mutate(Family = ifelse(Genus== "Dyadobacter" & Family == "Cytophagaceae","Spirosomaceae",Family))

ALL_TAX.40.38<-ALL_TAX.40.37%>% 
  mutate(Family = ifelse(Genus== "Eilatimonas" & Family == "Kordiimonadales_Incertae_Sedis","Temperatibacteraceae",Family))

ALL_TAX.40.39<-ALL_TAX.40.38%>% 
  mutate(Family = ifelse(Genus== "Ekhidna" & Family == "Cyclobacteriaceae","Reichenbachiellaceae",Family))

ALL_TAX.40.40<-ALL_TAX.40.39%>% 
  mutate(Family = ifelse(Genus== "Erythrobacter" & Family == "Sphingomonadaceae","Erythrobacteraceae",Family))

ALL_TAX.40.41<-ALL_TAX.40.40%>% 
  mutate(Family = ifelse(Genus== "Fabibacter" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

ALL_TAX.40.42<-ALL_TAX.40.41%>% 
  mutate(Family = ifelse(Genus== "Ferrimonas" & Family == "Shewanellaceae","Ferrimonadaceae",Family))

ALL_TAX.40.43<-ALL_TAX.40.42%>% 
  mutate(Family = ifelse(Genus== "Fulvitalea" & Family == "Cyclobacteriaceae","Ferrimonadaceae",Family))

ALL_TAX.40.44<-ALL_TAX.40.43%>% 
  mutate(Family = ifelse(Genus== "Gimesia" & Family == "Gimesiaceae","Planctomycetaceae",Family))

ALL_TAX.40.45<-ALL_TAX.40.44%>% 
  mutate(Family = ifelse(Genus== "Granulosicoccus" & Family == "Thiohalorhabdaceae","Granulosicoccaceae",Family))

ALL_TAX.40.46<-ALL_TAX.40.45%>% 
  mutate(Family = ifelse(Genus== "Halobacteriovorax" & Family == "Bacteriovoracaceae","Halobacteriovoraceae",Family))

ALL_TAX.40.47<-ALL_TAX.40.46%>% 
  mutate(Family = ifelse(Genus== "Haloferula" & Family == "Rubritaleaceae","Verrucomicrobiaceae",Family))

ALL_TAX.40.48<-ALL_TAX.40.47%>% 
  mutate(Family = ifelse(Genus== "Laribacter" & Family == "Aquaspirillaceae","Chromobacteriaceae",Family))

ALL_TAX.40.49<-ALL_TAX.40.48%>% 
  mutate(Family = ifelse(Genus== "Lewinella" & Family == "Saprospiraceae","Lewinellaceae",Family))

ALL_TAX.40.50<-ALL_TAX.40.49%>% 
  mutate(Family = ifelse(Genus== "Limnohabitans" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.51<-ALL_TAX.40.50%>% 
  mutate(Family = ifelse(Genus== "Macromonas" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.52<-ALL_TAX.40.51%>% 
  mutate(Family = ifelse(Genus== "Malikia" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.53<-ALL_TAX.40.52%>% 
  mutate(Family = ifelse(Genus== "Marinomonas" & Family == "Marinomonadaceae","Oceanospirillaceae",Family))

ALL_TAX.40.54<-ALL_TAX.40.53%>% 
  mutate(Family = ifelse(Genus== "Methylibium" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.55<-ALL_TAX.40.54%>% 
  mutate(Family = ifelse(Genus== "Motiliproteus" & Family == "Nitrincolaceae","Oceanospirillaceae",Family))

ALL_TAX.40.56<-ALL_TAX.40.55%>% 
  mutate(Family = ifelse(Genus== "Nisaea" & Family == "Nisaeaceae","Thalassobaculaceae",Family))

ALL_TAX.40.57<-ALL_TAX.40.56%>% 
  mutate(Family = ifelse(Genus== "Oceaniserpentilla" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

ALL_TAX.40.58<-ALL_TAX.40.57%>% 
  mutate(Family = ifelse(Genus== "Oceanococcus" & Family == "Solimonadaceae","Ectothiorhodospiraceae",Family))

ALL_TAX.40.59<-ALL_TAX.40.58%>% 
  mutate(Family = ifelse(Genus== "Oceanospirillum" & Family == "Halomonadaceae","Oceanospirillaceae",Family))

ALL_TAX.40.60<-ALL_TAX.40.59%>% 
  mutate(Family = ifelse(Genus== "Odoribacter" & Family == "Marinifilaceae","Odoribacteraceae",Family))

ALL_TAX.40.61<-ALL_TAX.40.60%>% 
  mutate(Family = ifelse(Genus== "Oleispira" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

ALL_TAX.40.62<-ALL_TAX.40.61%>% 
  mutate(Family = ifelse(Genus== "Ottowia" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.63<-ALL_TAX.40.62%>% 
  mutate(Family = ifelse(Genus== "Pelagibius" & Family == "Kiloniellaceae","Rhodovibrionaceae",Family))

ALL_TAX.40.64<-ALL_TAX.40.63%>% 
  mutate(Family = ifelse(Genus== "Pelomonas" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.65<-ALL_TAX.40.64%>% 
  mutate(Family = ifelse(Genus== "Pelosinus" & Family == "Veillonellaceae","Sporomusaceae",Family))

ALL_TAX.40.66<-ALL_TAX.40.65%>% 
  mutate(Family = ifelse(Genus== "Persicobacter" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

ALL_TAX.40.67<-ALL_TAX.40.66%>% 
  mutate(Family = ifelse(Genus== "Portibacter" & Family == "Saprospiraceae","Lewinellaceae",Family))

ALL_TAX.40.68<-ALL_TAX.40.67%>% 
  mutate(Family = ifelse(Genus== "Pseudarcicella" & Family == "Cytophagaceae","Spirosomaceae",Family))

ALL_TAX.40.69<-ALL_TAX.40.68%>% 
  mutate(Family = ifelse(Genus== "Pseudobacteriovorax" & Family == "Oligoflexaceae","Pseudobacteriovoracaceae",Family))

ALL_TAX.40.70<-ALL_TAX.40.69%>% 
  mutate(Family = ifelse(Genus== "Reichenbachiella" & Family == "Cyclobacteriaceae","Reichenbachiellaceae",Family))

ALL_TAX.40.71<-ALL_TAX.40.70%>% 
  mutate(Family = ifelse(Genus== "Rubricoccus" & Family == "Rhodothermaceae","Rubricoccaceae",Family))

ALL_TAX.40.72<-ALL_TAX.40.71%>% 
  mutate(Family = ifelse(Genus== "Rubrivirga" & Family == "Rhodothermaceae","Rubricoccaceae",Family))

ALL_TAX.40.73<-ALL_TAX.40.72%>% 
  mutate(Family = ifelse(Genus== "Rubrivivax" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.74<-ALL_TAX.40.73%>% 
  mutate(Family = ifelse(Genus== "Simplicispira" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.75<-ALL_TAX.40.74%>% 
  mutate(Family = ifelse(Genus== "Spirosoma" & Family == "Cytophagaceae","Spirosomaceae",Family))

ALL_TAX.40.76<-ALL_TAX.40.75%>% 
  mutate(Family = ifelse(Genus== "Spongiimonas" & Family == "Flavobacteriaceae","Weeksellaceae",Family))

ALL_TAX.40.77<-ALL_TAX.40.76%>% 
  mutate(Family = ifelse(Genus== "Spongiispira" & Family == "Saccharospirillaceae","Oceanospirillaceae",Family))

ALL_TAX.40.78<-ALL_TAX.40.77%>% 
  mutate(Family = ifelse(Genus== "Tunicatimonas" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

ALL_TAX.40.79<-ALL_TAX.40.78%>% 
  mutate(Family = ifelse(Genus== "Tunicatimonas" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

ALL_TAX.40.80<-ALL_TAX.40.79%>% 
  mutate(Family = ifelse(Genus== "Ideonella" & Family == "Burkholderiaceae","Comamonadaceae",Family))

ALL_TAX.40.81<-ALL_TAX.40.80%>% 
  mutate(Class= ifelse(Order== "Bacteriovoracales" & Class == "Deltaproteobacteria","Oligoflexia",Class))

ALL_TAX.40.82<-ALL_TAX.40.81%>% 
  mutate(Class= ifelse(Order== "Bacteroidales" & Class == "Sphingobacteriia","Bacteroidia",Class))

ALL_TAX.40.83<-ALL_TAX.40.82%>% 
  mutate(Class= ifelse(Order== "Calditrichales" & Class == "Deferribacteres","Calditrichia",Class))

ALL_TAX.40.84<-ALL_TAX.40.83%>% 
  mutate(Class= ifelse(Order== "Leptospirales" & Class == "Leptospirae","Spirochaetia",Class))

ALL_TAX.40.85<-ALL_TAX.40.84%>% 
  mutate(Class= ifelse(Order== "Nitrosomonadales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

ALL_TAX.40.86<-ALL_TAX.40.85%>% 
  mutate(Class= ifelse(Order== "Rhodocyclales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

ALL_TAX.40.87<-ALL_TAX.40.86%>% 
  mutate(Class= ifelse(Order== "Saprospirales" & Class == "Chitinophagia","Saprospiria",Class))

ALL_TAX.40.88<-ALL_TAX.40.87%>% 
  mutate(Class= ifelse(Order== "Burkholderiales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

ALL_TAX.40.89<-ALL_TAX.40.88%>% 
  mutate(Order= ifelse(Family== "Alteromonadales_incertae_sedis" & Order == "Cellvibrionales","Alteromonadales",Order))

ALL_TAX.40.90<-ALL_TAX.40.89%>% 
  mutate(Order= ifelse(Family== "Chromatiaceae" & Order == "Alteromonadales","Chromatiales",Order))

ALL_TAX.40.91<-ALL_TAX.40.90%>% 
  mutate(Order= ifelse(Family== "Chromobacteriaceae" & Order == "Betaproteobacteriales","Neisseriales",Order))

ALL_TAX.40.92<-ALL_TAX.40.91%>% 
  mutate(Order= ifelse(Family== "Desulfobacteraceae" & Order == "Desulfarculales","Desulfobacterales",Order))

ALL_TAX.40.93<-ALL_TAX.40.92%>% 
  mutate(Order= ifelse(Family== "Ectothiorhodospiraceae" & Order == "Salinisphaerales","Chromatiales",Order))

ALL_TAX.40.94<-ALL_TAX.40.93%>% 
  mutate(Order= ifelse(Family== "Egicoccaceae" & Order == "Nitriliruptorales","Egicoccales",Order))

ALL_TAX.40.95<-ALL_TAX.40.94%>% 
  mutate(Order= ifelse(Family== "Ferrimonadaceae" & Order == "Cytophagales","Alteromonadales",Order))

ALL_TAX.40.96<-ALL_TAX.40.95%>% 
  mutate(Order= ifelse(Family== "Granulosicoccaceae" & Order == "Thiohalorhabdales","Chromatiales",Order))

ALL_TAX.40.97<-ALL_TAX.40.96%>% 
  mutate(Order= ifelse(Family== "Nevskiaceae" & Order == "Salinisphaerales","Nevskiales",Order))

ALL_TAX.40.98<-ALL_TAX.40.97%>% 
  mutate(Order= ifelse(Family== "Piscirickettsiaceae" & Order == "Nitrosococcales","Thiotrichales",Order))

ALL_TAX.40.99<-ALL_TAX.40.98%>% 
  mutate(Order= ifelse(Family== "Piscirickettsiaceae" & Order == "Thiomicrospirales","Thiotrichales",Order))

ALL_TAX.41<-ALL_TAX.40.99%>% 
  mutate(Order= ifelse(Family== "Rhodospirillaceae" & Order == "Dongiales","Rhodospirillales",Order))

ALL_TAX.41.1<-ALL_TAX.41%>% 
  mutate(Order= ifelse(Family== "Rhodovibrionaceae" & Order == "Kiloniellales","Rhodospirillales",Order))

ALL_TAX.41.2<-ALL_TAX.41.1%>% 
  mutate(Order= ifelse(Family== "Selenomonadaceae" & Order == "Veillonellales","Selenomonadales",Order))

ALL_TAX.41.3<-ALL_TAX.41.2%>% 
  mutate(Order= ifelse(Family== "Sporomusaceae" & Order == "Veillonellales","Selenomonadales",Order))

ALL_TAX.41.4<-ALL_TAX.41.3%>% 
  mutate(Order= ifelse(Family== "Sterolibacteriaceae" & Order == "Rhodocyclales","Nitrosomonadales",Order))

ALL_TAX.41.5<-ALL_TAX.41.4%>% 
  mutate(Order= ifelse(Family== "Thalassobaculaceae" & Order == "Thalassobaculales","Rhodospirillales",Order))

ALL_TAX.41.6<-ALL_TAX.41.5%>% 
  mutate(Family= ifelse(Family== "Rubrobacteriaceae" ,"Rubrobacteraceae",Family))

ALL_TAX.41.7<-ALL_TAX.41.6%>% 
  mutate(Genus= ifelse(Genus== "Prevotella 9" ,"Prevotella",Genus))

ALL_TAX.41.8<-ALL_TAX.41.7%>% 
  mutate(Family= ifelse(Genus== "Fulvitalea" & Family == "Ferrimonadaceae","Flammeovirgaceae",Family))

ALL_TAX.41.9<-ALL_TAX.41.8%>% 
  mutate(Order= ifelse(Family== "Flammeovirgaceae" & Order == "Alteromonadales","Cytophagales",Order))

ALL_TAX.41.10<-ALL_TAX.41.9%>% 
  mutate(Class= ifelse(Order== "Alteromonadales" & Class == "Cytophagia","Gammaproteobacteria",Class))

ALL_TAX.41.11<-ALL_TAX.41.10%>% 
  mutate(Class= ifelse(Order== "Neisseriales" & Class == "Gammaproteobacteria","Betaproteobacteria",Class))

ALL_TAX.41.12<-ALL_TAX.41.11%>% 
  mutate(Phylum= ifelse(Class== "Calditrichia" & Phylum == "Deferribacteres","Calditrichaeota",Phylum))

ALL_TAX.41.13<-ALL_TAX.41.12%>% 
  mutate(Phylum= ifelse(Class== "Oligoflexia" & Phylum == "Proteobacteria","Bdellovibrionota",Phylum))

ALL_TAX.41.14<-ALL_TAX.41.13%>% 
  mutate(Phylum= ifelse(Class== "Thermoleophilia" & Phylum == "Actinobacteriota","Actinobacteria",Phylum))

ALL_TAX.41.15<-ALL_TAX.41.14%>% 
  mutate(Phylum= ifelse(Class== "Chlamydiae" & Phylum == "Verrucomicrobia", "Chlamydiae",Phylum))

ALL_TAX.41.16<-ALL_TAX.41.15%>% 
  mutate(Class= ifelse(Class== "Chlamydiae" & Phylum == "Chlamydiae", "Chlamydiia",Class))

ALL_TAX.41.17<-ALL_TAX.41.16%>% 
  mutate(Phylum= ifelse(Class== "Acidimicrobiia" & Phylum == "Actinobacteriota", "Actinobacteria",Phylum))

ALL_TAX.41.18<-ALL_TAX.41.17%>% 
  mutate(Phylum= ifelse(Class== "Fibrobacteria" & Phylum == "Fibrobacterota", "Fibrobacteres",Phylum))

ALL_TAX.41.19<-ALL_TAX.41.18%>% 
  mutate(Phylum= ifelse(Class== "Ignavibacteria" & Phylum == "Bacteroidetes", "Ignavibacteriae",Phylum))

ALL_TAX.41.20<-ALL_TAX.41.19%>% 
  mutate(Phylum= ifelse(Class== "Chitinivibrionia" & Phylum == "Fibrobacterota", "Fibrobacteres",Phylum))

ALL_TAX.41.21<-ALL_TAX.41.20%>% 
  mutate(Phylum= ifelse(Class== "Kiritimatiellae" & Phylum == "Verrucomicrobia", "Kiritimatiellaeota",Phylum))

ALL_TAX.41.22<-ALL_TAX.41.21%>% 
  mutate(Phylum= ifelse(Class== "Omnitrophia" & Phylum == "Verrucomicrobia", "Omnitrophicaeota",Phylum))

ALL_TAX.41.23<-ALL_TAX.41.22%>% 
  mutate(Phylum= ifelse(Class== "Rubrobacteria" & Phylum == "Actinobacteriota", "Actinobacteria",Phylum))

ALL_TAX.41.24<-ALL_TAX.41.23%>% 
  mutate(Phylum= ifelse(Class== "Campylobacteria" & Phylum == "Epsilonbacteraeota", "Proteobacteria",Phylum))

ALL_TAX.41.25<-ALL_TAX.41.24%>% 
  mutate(Phylum= ifelse(Class== "Campylobacteria" & Phylum == "Campilobacterota", "Proteobacteria",Phylum))

ALL_TAX.41.26<-ALL_TAX.41.25%>% 
  mutate(Phylum= ifelse(Class== "Deinococci" & Phylum == "Deinococcota", "Deinococcus-Thermus",Phylum))

ALL_TAX.41.27<-ALL_TAX.41.26%>% 
  mutate(Phylum= ifelse(Class== "Lentisphaeria" & Phylum == "Verrucomicrobia", "Lentisphaerae",Phylum))

ALL_TAX.41.28<-ALL_TAX.41.27%>% 
  mutate(Phylum= ifelse(Class== "Lineage IIb" & Phylum == "Elusimicrobiota", "Elusimicrobia",Phylum))

ALL_TAX.41.29<-ALL_TAX.41.28%>% 
  mutate(Phylum= ifelse(Class== "Actinobacteria" & Phylum == "Actinobacteriota", "Actinobacteria",Phylum))

ALL_TAX.41.30<-ALL_TAX.41.29%>% 
  mutate(Phylum= ifelse(Class== "BD2-11 terrestrial group" & Phylum == "Gemmatimonadota", "Gemmatimonadetes",Phylum))

ALL_TAX.41.31<-ALL_TAX.41.30%>% 
  mutate(Phylum= ifelse(Class== "Campylobacteria" & Phylum == "Campylobacterota", "Proteobacteria",Phylum))

ALL_TAX.41.32<-ALL_TAX.41.31%>% 
  mutate(Class= ifelse(Order== "11-24" & Class == "Blastocatellia_(Subgroup_4)", "Blastocatellia",Class))

ALL_TAX.41.33<-ALL_TAX.41.32%>% 
  mutate(Class= ifelse(Order== "Acidimicrobiales" & Class == "Acidimicrobiia", "Actinobacteria",Class))

ALL_TAX.41.34<-ALL_TAX.41.33%>% 
  mutate(Class= ifelse(Order== "Bacteriovoracales" & Class == "Bdellovibrionia", "Oligoflexia",Class))

ALL_TAX.41.35<-ALL_TAX.41.34%>% 
  mutate(Class= ifelse(Order== "Bdellovibrionales" & Class == "Bdellovibrionia", "Deltaproteobacteria",Class))

ALL_TAX.41.36<-ALL_TAX.41.35%>% 
  mutate(Class= ifelse(Order== "Bradymonadales" & Class == "Desulfuromonadia", "Deltaproteobacteria",Class))

ALL_TAX.41.37<-ALL_TAX.41.36%>% 
  mutate(Class= ifelse(Order== "Caenarcaniphilales" & Class == "Vampirivibrionia", "Melainabacteria",Class))

ALL_TAX.41.38<-ALL_TAX.41.37%>% 
  mutate(Class= ifelse(Order== "Campylobacterales" & Class == "Campylobacteria", "Epsilonproteobacteria",Class))

ALL_TAX.41.39<-ALL_TAX.41.38%>% 
  mutate(Class= ifelse(Order== "Desulfobacterales" & Class == "Desulfobacteria", "Deltaproteobacteria",Class))

ALL_TAX.41.40<-ALL_TAX.41.39%>% 
  mutate(Class= ifelse(Order== "Desulfovibrionales" & Class == "Desulfovibrionia", "Deltaproteobacteria",Class))

ALL_TAX.41.41<-ALL_TAX.41.40%>% 
  mutate(Class= ifelse(Order== "Desulfuromonadales" & Class == "Desulfuromonadia", "Deltaproteobacteria",Class))

ALL_TAX.41.42<-ALL_TAX.41.41%>% 
  mutate(Class= ifelse(Order== "Gemmatales" & Class == "Planctomycetacia", "Planctomycetes",Class))

ALL_TAX.41.43<-ALL_TAX.41.42%>% 
  mutate(Class= ifelse(Order== "Leptolyngbyales" & Class == "Oxyphotobacteria", "Cyanobacteriia",Class))

ALL_TAX.41.44<-ALL_TAX.41.43%>% 
  mutate(Class= ifelse(Order== "Mycoplasmatales" & Class == "Bacilli", "Mollicutes",Class))

ALL_TAX.41.45<-ALL_TAX.41.44%>% 
  mutate(Class= ifelse(Order== "Myxococcales" & Class == "Myxococcia", "Deltaproteobacteria",Class))

ALL_TAX.41.46<-ALL_TAX.41.45%>% 
  mutate(Class= ifelse(Order== "Nitriliruptorales" & Class == "Actinobacteria", "Nitriliruptoria",Class))

ALL_TAX.41.47<-ALL_TAX.41.46%>% 
  mutate(Class= ifelse(Order== "Nitrospirales" & Class == "Nitrospiria", "Nitrospira",Class))

ALL_TAX.41.48<-ALL_TAX.41.47%>% 
  mutate(Class= ifelse(Order== "Oligosphaerales" & Class == "Lentisphaeria", "Oligosphaeria",Class))

ALL_TAX.41.49<-ALL_TAX.41.48%>% 
  mutate(Class= ifelse(Order== "P.palmC41" & Class == "Lentisphaeria", "Oligosphaeria",Class))

ALL_TAX.41.50<-ALL_TAX.41.49%>% 
  mutate(Class= ifelse(Order== "Phormidesmiales" & Class == "Oxyphotobacteria", "Cyanobacteriia",Class))

ALL_TAX.41.51<-ALL_TAX.41.50%>% 
  mutate(Class= ifelse(Order== "PB19" & Class == "Desulfuromonadia", "Deltaproteobacteria",Class))

ALL_TAX.41.52<-ALL_TAX.41.51%>% 
  mutate(Class= ifelse(Order== "Pirellulales" & Class == "Planctomycetacia", "Planctomycetia",Class))

ALL_TAX.41.53<-ALL_TAX.41.52%>% 
  mutate(Class= ifelse(Order== "Pirellulales" & Class == "Planctomycetes", "Planctomycetia",Class))

ALL_TAX.41.54<-ALL_TAX.41.53%>% 
  mutate(Class= ifelse(Order== "Planctomycetales" & Class == "Planctomycetacia", "Planctomycetia",Class))

ALL_TAX.41.55<-ALL_TAX.41.54%>% 
  mutate(Class= ifelse(Order== "Planctomycetales" & Class == "Planctomycetes", "Planctomycetia",Class))

ALL_TAX.41.56<-ALL_TAX.41.55%>% 
  mutate(Class= ifelse(Order== "Pseudanabaenales" & Class == "Oxyphotobacteria", "Cyanobacteriia",Class))

ALL_TAX.41.57<-ALL_TAX.41.56%>% 
  mutate(Class= ifelse(Order== "Rubrobacterales" & Class == "Rubrobacteria", "Actinobacteria",Class))

ALL_TAX.41.58<-ALL_TAX.41.57%>% 
  mutate(Class= ifelse(Order== "Solibacterales" & Class == "Acidobacteriia", "Acidobacteria",Class))

ALL_TAX.41.59<-ALL_TAX.41.58%>% 
  mutate(Class= ifelse(Order== "SS1-B-02-17" & Class == "Lentisphaeria", "Oligosphaeria",Class))

ALL_TAX.41.60<-ALL_TAX.41.59%>% 
  mutate(Class= ifelse(Order== "Synechococcales" & Class == "Oxyphotobacteria", "Cyanobacteriia",Class))

ALL_TAX.41.61<-ALL_TAX.41.60%>% 
  mutate(Class= ifelse(Order== "Thermosynechococcales" & Class == "Oxyphotobacteria", "Cyanobacteriia",Class))

ALL_TAX.41.62<-ALL_TAX.41.61%>% 
  mutate(Phylum= ifelse(Class== "Deltaproteobacteria" & Phylum == "Bdellovibrionota", "Proteobacteria",Phylum))

ALL_TAX.41.63<-ALL_TAX.41.62%>% 
  mutate(Phylum= ifelse(Class== "Deltaproteobacteria" & Phylum == "Desulfobacterota", "Proteobacteria",Phylum))

ALL_TAX.41.64<-ALL_TAX.41.63%>% 
  mutate(Phylum= ifelse(Class== "Deltaproteobacteria" & Phylum == "Myxococcota", "Proteobacteria",Phylum))

ALL_TAX.41.65<-ALL_TAX.41.64%>% 
  mutate(Phylum= ifelse(Class== "Mollicutes" & Phylum == "Firmicutes", "Tenericutes",Phylum))

ALL_TAX.41.66<-ALL_TAX.41.65%>% 
  mutate(Order= ifelse(Family== "Aeromonadaceae" & Order== "Enterobacterales", "Aeromonadales",Order))

ALL_TAX.41.67<-ALL_TAX.41.66%>% 
  mutate(Order= ifelse(Family== "Alteromonadaceae" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.41.68<-ALL_TAX.41.67%>% 
  mutate(Order= ifelse(Family== "Alteromonadaceae" & Order== "Pseudomonadales", "Alteromonadales",Order))

ALL_TAX.41.69<-ALL_TAX.41.68%>% 
  mutate(Order= ifelse(Family== "Alteromonadales_incertae_sedis" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.41.70<-ALL_TAX.41.69%>% 
  mutate(Order= ifelse(Family== "Azonexaceae" & Order== "Burkholderiales", "Rhodocyclales",Order))

ALL_TAX.41.71<-ALL_TAX.41.70%>% 
  mutate(Order= ifelse(Family== "BD2-7" & Order== "Pseudomonadales", "Cellvibrionales",Order))

ALL_TAX.41.72<-ALL_TAX.41.71%>% 
  mutate(Order= ifelse(Family== "Cellvibrionaceae" & Order== "Pseudomonadales", "Cellvibrionales",Order))

ALL_TAX.41.73<-ALL_TAX.41.72%>% 
  mutate(Order= ifelse(Family== "BIrii41" & Order== "Polyangiales", "Myxococcales",Order))

ALL_TAX.41.74<-ALL_TAX.41.73%>% 
  mutate(Order= ifelse(Family== "Chitinophagaceae" & Order== "Sphingobacteriales", "Chitinophagales",Order))

ALL_TAX.41.75<-ALL_TAX.41.74%>% 
  mutate(Order= ifelse(Family== "Colwelliaceae" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.41.76<-ALL_TAX.41.75%>% 
  mutate(Order= ifelse(Family== "Defluviitaleaceae" & Order== "Lachnospirales", "Clostridiales",Order))

ALL_TAX.41.77<-ALL_TAX.41.76%>% 
  mutate(Order= ifelse(Family== "Desulfobulbaceae" & Order== "Desulfobulbales", "Desulfobacterales",Order))

ALL_TAX.41.78<-ALL_TAX.41.77%>% 
  mutate(Order= ifelse(Family== "EC94" & Order== "Burkholderiales", "Betaproteobacteriales",Order))

ALL_TAX.41.79<-ALL_TAX.41.78%>% 
  mutate(Order= ifelse(Family== "Endozoicomonadaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.41.80<-ALL_TAX.41.79%>% 
  mutate(Order= ifelse(Family== "Ferrimonadaceae" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.41.81<-ALL_TAX.41.80%>% 
  mutate(Order= ifelse(Family== "Geobacteraceae" & Order== "Geobacterales", "Desulfuromonadales",Order))

ALL_TAX.41.82<-ALL_TAX.41.81%>% 
  mutate(Order= ifelse(Family== "Granulosicoccaceae" & Order== "Granulosicoccales", "Chromatiales",Order))

ALL_TAX.41.83<-ALL_TAX.41.82%>% 
  mutate(Order= ifelse(Family== "Gven-F17" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.41.84<-ALL_TAX.41.83%>% 
  mutate(Order= ifelse(Family== "Hahellaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.41.85<-ALL_TAX.41.84%>% 
  mutate(Order= ifelse(Family== "Haliangiaceae" & Order== "Haliangiales", "Myxococcales",Order))

ALL_TAX.41.86<-ALL_TAX.41.85 %>% 
  mutate(Order= ifelse(Family== "Halieaceae" & Order== "Pseudomonadales", "Cellvibrionales",Order))

ALL_TAX.41.87<-ALL_TAX.41.86 %>% 
  mutate(Order= ifelse(Family== "Halomonadaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.41.88<-ALL_TAX.41.87 %>% 
  mutate(Order= ifelse(Family== "Intrasporangiaceae" & Order== "Actinomycetales", "Micrococcales",Order))

ALL_TAX.41.89<-ALL_TAX.41.88 %>% 
  mutate(Order= ifelse(Family== "Kangiellaceae" & Order== "Enterobacterales", "Oceanospirillales",Order))

ALL_TAX.41.90<-ALL_TAX.41.89 %>% 
  mutate(Order= ifelse(Family== "Lachnospiraceae" & Order== "Lachnospirales", "Clostridiales",Order))

ALL_TAX.41.91<-ALL_TAX.41.90 %>% 
  mutate(Order= ifelse(Family== "Lewinellaceae" & Order== "Sphingobacteriales", "Saprospirales",Order))

ALL_TAX.41.92<-ALL_TAX.41.91 %>% 
  mutate(Order= ifelse(Family== "Litoricolaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.41.93<-ALL_TAX.41.92 %>% 
  mutate(Order= ifelse(Family== "Methylophilaceae" & Order== "Methylophilales", "Nitrosomonadales",Order))

ALL_TAX.41.94<-ALL_TAX.41.93 %>% 
  mutate(Order= ifelse(Family== "Microbacteriaceae" & Order== "Actinomycetales", "Micrococcales",Order))

ALL_TAX.41.95<-ALL_TAX.41.94 %>% 
  mutate(Order= ifelse(Family== "Micrococcaceae" & Order== "Actinomycetales", "Micrococcales",Order))

ALL_TAX.41.96<-ALL_TAX.41.95 %>% 
  mutate(Order= ifelse(Family== "Moritellaceae" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.41.97<-ALL_TAX.41.96 %>% 
  mutate(Order= ifelse(Family== "Mycobacteriaceae" & Order== "Actinomycetales", "Mycobacteriales",Order))

ALL_TAX.41.98<-ALL_TAX.41.97 %>% 
  mutate(Order= ifelse(Family== "Nannocystaceae" & Order== "Nannocystales", "Myxococcales",Order))

ALL_TAX.41.99<-ALL_TAX.41.98 %>% 
  mutate(Order= ifelse(Family== "Neisseriaceae" & Order== "Betaproteobacteriales", "Burkholderiales",Order))

ALL_TAX.42<-ALL_TAX.41.99 %>% 
  mutate(Order= ifelse(Family== "Nitrincolaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.42.1<-ALL_TAX.42 %>% 
  mutate(Order= ifelse(Family== "Nitrosomonadaceae" & Order== "Burkholderiales", "Nitrosomonadales",Order))

ALL_TAX.42.2<-ALL_TAX.42.1 %>% 
  mutate(Order= ifelse(Family== "Nocardiaceae" & Order== "Actinomycetales", "Mycobacteriales",Order))

ALL_TAX.42.3<-ALL_TAX.42.2 %>% 
  mutate(Order= ifelse(Family== "Oceanospirillaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.42.4<-ALL_TAX.42.3 %>% 
  mutate(Order= ifelse(Family== "Oceanospirillales_incertae_sedis" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.42.5<-ALL_TAX.42.4 %>% 
  mutate(Order= ifelse(Family== "Oleiphilaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.42.6<-ALL_TAX.42.5 %>% 
  mutate(Order= ifelse(Family== "Peptostreptococcaceae" & Order== "Peptostreptococcales-Tissierellales", "Clostridiales",Order))

ALL_TAX.42.7<-ALL_TAX.42.6 %>% 
  mutate(Order= ifelse(Family== "Porticoccaceae" & Order== "Pseudomonadales", "Cellvibrionales",Order))

ALL_TAX.42.8<-ALL_TAX.42.7 %>% 
  mutate(Order= ifelse(Family== "Pseudoalteromonadaceae" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.42.9<-ALL_TAX.42.8 %>% 
  mutate(Order= ifelse(Family== "Pseudobacteriovoracaceae" & Order== "Bdellovibrionales", "Oligoflexales",Order))

ALL_TAX.42.10<-ALL_TAX.42.9 %>% 
  mutate(Order= ifelse(Family== "Psychromonadaceae" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.42.11<-ALL_TAX.42.10 %>% 
  mutate(Order= ifelse(Family== "Rhodocyclaceae" & Order== "Burkholderiales", "Rhodocyclales",Order))

ALL_TAX.42.12<-ALL_TAX.42.11 %>% 
  mutate(Order= ifelse(Family== "Rhodothermaceae" & Order== "Sphingobacteriales", "Rhodothermales",Order))

ALL_TAX.42.13<-ALL_TAX.42.12 %>% 
  mutate(Order= ifelse(Family== "Rubricoccaceae" & Order== "Sphingobacteriales", "Rhodothermales",Order))

ALL_TAX.42.14<-ALL_TAX.42.13 %>% 
  mutate(Order= ifelse(Family== "Ruminococcaceae" & Order== "Clostridiales", "Oscillospirales",Order))

ALL_TAX.42.15<-ALL_TAX.42.14 %>% 
   mutate(Order= ifelse(Family== "Saccharospirillaceae" & Order== "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.42.16<-ALL_TAX.42.15 %>% 
  mutate(Order= ifelse(Family== "Salinisphaeraceae" & Order== "Salinisphaerales", "Nevskiales",Order))

ALL_TAX.42.17<-ALL_TAX.42.16 %>% 
  mutate(Order= ifelse(Family== "Sandaracinaceae" & Order== "Polyangiales", "Myxococcales",Order))

ALL_TAX.42.18<-ALL_TAX.42.17 %>% 
  mutate(Order= ifelse(Family== "Saprospiraceae" & Order== "Sphingobacteriales", "Saprospirales",Order))

ALL_TAX.42.19<-ALL_TAX.42.18 %>% 
  mutate(Order= ifelse(Family== "Shewanellaceae" & Order== "Enterobacterales", "Alteromonadales",Order))

ALL_TAX.42.20<-ALL_TAX.42.19 %>% 
  mutate(Order= ifelse(Family== "Spongiibacteraceae" & Order== "Pseudomonadales", "Cellvibrionales",Order))

ALL_TAX.42.21<-ALL_TAX.42.20 %>% 
  mutate(Order= ifelse(Family== "Staphylococcaceae" & Order== "Staphylococcales", "Bacillales",Order))

ALL_TAX.42.22<-ALL_TAX.42.21 %>% 
  mutate(Order= ifelse(Family== "Sterolibacteriaceae" & Order== "Burkholderiales", "Nitrosomonadales",Order))

ALL_TAX.42.23<-ALL_TAX.42.22 %>% 
  mutate(Order= ifelse(Family== "Succinivibrionaceae" & Order== "Enterobacterales", "Aeromonadales",Order))

ALL_TAX.42.24<-ALL_TAX.42.23 %>% 
  mutate(Order= ifelse(Family== "Thioalkalispiraceae" & Order== "Ectothiorhodospirales", "Chromatiales",Order))

ALL_TAX.42.25<-ALL_TAX.42.24 %>% 
  mutate(Order= ifelse(Family== "Thioglobaceae" & Order== "Pseudomonadales", "Thiomicrospirales",Order))

ALL_TAX.42.26<-ALL_TAX.42.25 %>% 
  mutate(Order= ifelse(Family== "Vibrionaceae" & Order== "Enterobacterales", "Vibrionales",Order))

ALL_TAX.42.27<-ALL_TAX.42.26 %>% 
  mutate(Order= ifelse(Family== "Xenococcaceae" & Order== "Nostocales", "Cyanobacteriales",Order))

ALL_TAX.42.28<-ALL_TAX.42.27 %>% 
  mutate(Order= ifelse(Family== "Zoogloeaceae" & Order== "Burkholderiales", "Rhodocyclales",Order))

ALL_TAX.42.29<-ALL_TAX.42.28 %>% 
  mutate(Order= ifelse(Order== "Betaproteobacteriales" , "Burkholderiales",Order))

ALL_TAX.42.30<-ALL_TAX.42.29%>% 
  mutate(Class= ifelse(Order== "Burkholderiales" & Class == "Betaproteobacteria", "Gammaproteobacteria",Class))

ALL_TAX.42.31<-ALL_TAX.42.30%>% 
  mutate(Class= ifelse(Order== "Chitinophagales" & Class == "Sphingobacteriia", "Chitinophagia",Class))

ALL_TAX.42.32<-ALL_TAX.42.31%>% 
  mutate(Class= ifelse(Order== "Cyanobacteriales" & Class == "Oxyphotobacteria", "Cyanobacteriia",Class))

ALL_TAX.42.33<-ALL_TAX.42.32%>% 
  mutate(Class= ifelse(Order== "Desulfobacterales" & Class == "Desulfobulbia", "Deltaproteobacteria",Class))

ALL_TAX.42.34<-ALL_TAX.42.33%>% 
  mutate(Class= ifelse(Order== "Desulfuromonadales" & Class == "Desulfuromonadia", "Deltaproteobacteria",Class))

ALL_TAX.42.35<-ALL_TAX.42.34%>% 
  mutate(Class= ifelse(Order== "Myxococcales" & Class == "Polyangia", "Deltaproteobacteria",Class))

ALL_TAX.42.36<-ALL_TAX.42.35%>% 
  mutate(Class= ifelse(Order== "Oligoflexales" & Class == "Deltaproteobacteria", "Oligoflexia",Class))

ALL_TAX.42.37<-ALL_TAX.42.36%>% 
  mutate(Class= ifelse(Order== "Rhodothermales" & Class == "Sphingobacteriia", "Rhodothermia",Class))

ALL_TAX.42.38<-ALL_TAX.42.37%>% 
  mutate(Class= ifelse(Order== "Saprospirales" & Class == "Sphingobacteriia", "Saprospiria",Class))

ALL_TAX.42.39<-ALL_TAX.42.38%>% 
  mutate(Phylum= ifelse(Class== "Oligoflexia" & Phylum == "Bdellovibrionota", "Proteobacteria",Phylum))

ALL_TAX.42.40<-ALL_TAX.42.39%>% 
  mutate(Phylum= ifelse(Class== "Rhodothermia" & Phylum == "Bacteroidetes", "Rhodothermaeota",Phylum))

ALL_TAX.42.41<-ALL_TAX.42.40%>% 
  mutate(Phylum= ifelse(Class== "Deltaproteobacteria" & Phylum == "Myxococcota", "Proteobacteria",Phylum))

ALL_TAX.42.42<-ALL_TAX.42.41%>% 
  mutate(Phylum= ifelse(Class== "Deltaproteobacteria" & Phylum == "Desulfobacterota", "Proteobacteria",Phylum))

ALL_TAX.42.43<-ALL_TAX.42.42%>% 
  mutate(Order= ifelse(Family== "Demequinaceae" & Order == "Actinomycetales", "Micrococcales",Order))

ALL_TAX.42.44<-ALL_TAX.42.43%>% 
  mutate(Order= ifelse(Family== "Methylophilaceae" & Order == "Burkholderiales", "Nitrosomonadales",Order))

ALL_TAX.42.45<-ALL_TAX.42.44%>% 
  mutate(Class= ifelse(Order== "Nitrosomonadales" & Class == "Gammaproteobacteria", "Betaproteobacteria",Class))

ALL_TAX.42.46<-ALL_TAX.42.45%>% 
  mutate(Class= ifelse(Order== "Nitrosomonadales" & Class == "Gammaproteobacteria", "Betaproteobacteria",Class))

ALL_TAX.42.47<-ALL_TAX.42.46%>% 
  mutate(Family = ifelse(Genus== "Adhaeribacter" & Family == "Cytophagaceae","Hymenobacteraceae",Family))

ALL_TAX.42.48<-ALL_TAX.42.47%>% 
  mutate(Family = ifelse(Genus== "Ahrensia" & Family == "Rhodobacteraceae","Ahrensiaceae",Family))

ALL_TAX.42.49<-ALL_TAX.42.48%>% 
  mutate(Family = ifelse(Genus== "Alcanivorax" & Family == "Alcanivoracaceae1","Alcanivoracaceae",Family))

ALL_TAX.42.50<-ALL_TAX.42.49%>% 
  mutate(Family = ifelse(Genus== "Anderseniella" & Family == "Rhodobiaceae","Parvibaculaceae",Family))

ALL_TAX.42.51<-ALL_TAX.42.50%>% 
  mutate(Family = ifelse(Genus== "Aquabacterium" & Family == "Burkholderiales_incertae_sedis","Comamonadaceae",Family))

ALL_TAX.42.52<-ALL_TAX.42.51%>% 
  mutate(Family = ifelse(Genus== "Aquicella" & Family == "Diplorickettsiaceae","Coxiellaceae",Family))

ALL_TAX.42.53<-ALL_TAX.42.52%>% 
  mutate(Family = ifelse(Genus== "Arcobacter" & Family == "Campylobacteraceae","Arcobacteraceae",Family))

ALL_TAX.42.54<-ALL_TAX.42.53%>% 
  mutate(Family = ifelse(Genus== "Aureibacter" & Family == "Cyclobacteriaceae","Flammeovirgaceae",Family))

ALL_TAX.42.55<-ALL_TAX.42.54%>% 
  mutate(Family = ifelse(Genus== "Balneola" & Family == "Chitinophagaceae","Balneolaceae",Family))

ALL_TAX.42.56<-ALL_TAX.42.55%>% 
  mutate(Family = ifelse(Genus== "Blastopirellula" & Family == "Planctomycetaceae","Pirellulaceae",Family))

ALL_TAX.42.57<-ALL_TAX.42.56%>% 
  mutate(Family = ifelse(Genus== "Bryobacter" & Family == "Solibacteraceae_(Subgroup_3)","Bryobacteraceae",Family))

ALL_TAX.42.58<-ALL_TAX.42.57%>% 
  mutate(Family = ifelse(Genus== "Chryseobacterium" & Family == "Flavobacteriaceae","Weeksellaceae",Family))

ALL_TAX.42.59<-ALL_TAX.42.58%>% 
  mutate(Family = ifelse(Genus== "Cloacibacterium" & Family == "Flavobacteriaceae","Weeksellaceae",Family))

ALL_TAX.42.60<-ALL_TAX.42.59%>% 
  mutate(Family = ifelse(Genus== "Dasania" & Family == "Pseudomonadales_incertae_sedis","Spongiibacteraceae",Family))

ALL_TAX.42.61<-ALL_TAX.42.60%>% 
  mutate(Family = ifelse(Genus== "Desulfatiglans" & Family == "Desulfatiglandaceae","Desulfobacteraceae",Family))

ALL_TAX.42.62<-ALL_TAX.42.61%>% 
  mutate(Family = ifelse(Genus== "Desulfocapsa" & Family == "Desulfocapsaceae","Desulfobulbaceae",Family))

ALL_TAX.42.63<-ALL_TAX.42.62%>% 
  mutate(Family = ifelse(Genus== "Desulfopila" & Family == "Desulfocapsaceae","Desulfobulbaceae",Family))

ALL_TAX.42.64<-ALL_TAX.42.63%>% 
  mutate(Family = ifelse(Genus== "Desulforhopalus" & Family == "Desulfocapsaceae","Desulfobulbaceae",Family))

ALL_TAX.42.65<-ALL_TAX.42.64%>% 
  mutate(Family = ifelse(Genus== "Desulfosarcina" & Family == "Desulfocapsaceae","Desulfobulbaceae",Family))

ALL_TAX.42.66<-ALL_TAX.42.65%>% 
  mutate(Family = ifelse(Genus== "Desulfotalea" & Family == "Desulfocapsaceae","Desulfobulbaceae",Family))

ALL_TAX.42.67<-ALL_TAX.42.66%>% 
  mutate(Family = ifelse(Genus== "Desulfuromusa" & Family == "Geopsychrobacteraceae","Desulfuromonadaceae",Family))

ALL_TAX.42.68<-ALL_TAX.42.67%>% 
  mutate(Family = ifelse(Genus== "Ekhidna" & Family == "Cytophagaceae","Reichenbachiellaceae",Family))

ALL_TAX.42.69<-ALL_TAX.42.68%>% 
  mutate(Family = ifelse(Genus== "Endozoicomonas" & Family == "Hahellaceae","Endozoicomonadaceae",Family))

ALL_TAX.42.70<-ALL_TAX.42.69%>% 
  mutate(Family = ifelse(Genus== "Fangia" & Family == "Francisellaceae","Thiotrichales_incertae_sedis",Family))

ALL_TAX.42.71<-ALL_TAX.42.70%>% 
  mutate(Family = ifelse(Genus== "Fluviicola" & Family == "Crocinitomicaceae","Cryomorphaceae",Family))

ALL_TAX.42.72<-ALL_TAX.42.71%>% 
  mutate(Family = ifelse(Genus== "Fulvivirga" & Family == "Flammeovirgaceae","Fulvivirgaceae",Family))

ALL_TAX.42.73<-ALL_TAX.42.72%>% 
  mutate(Family = ifelse(Genus== "Fusibacter" & Family == "Fusibacteraceae","Clostridiales_Incertae Sedis XII",Family))

ALL_TAX.42.74<-ALL_TAX.42.73%>% 
  mutate(Family = ifelse(Genus== "Gemmata" & Family == "Planctomycetaceae","Gemmataceae",Family))

ALL_TAX.42.75<-ALL_TAX.42.74%>% 
  mutate(Family = ifelse(Genus== "Haliea" & Family == "Alteromonadaceae","Halieaceae",Family))

ALL_TAX.42.76<-ALL_TAX.42.75%>% 
  mutate(Family = ifelse(Genus== "Hoeflea" & Family == "Phyllobacteriaceae","Rhizobiaceae",Family))

ALL_TAX.42.77<-ALL_TAX.42.76%>% 
  mutate(Family = ifelse(Genus== "Ideonella" & Family == "Burkholderiales_incertae_sedis","Comamonadaceae",Family))

ALL_TAX.42.78<-ALL_TAX.42.77%>% 
  mutate(Family = ifelse(Genus== "Ilumatobacter" & Family == "Acidimicrobiaceae","Ilumatobacteraceae",Family))

ALL_TAX.42.79<-ALL_TAX.42.78%>% 
  mutate(Family = ifelse(Genus== "Izimaplasma" & Family == "Izemoplasmataceae","Izimaplasmataceae",Family))

ALL_TAX.42.80<-ALL_TAX.42.79%>% 
  mutate(Family = ifelse(Genus== "Lentilitoribacter" & Family == "Phyllobacteriaceae","Rhizobiaceae",Family))

ALL_TAX.42.81<-ALL_TAX.42.80%>% 
  mutate(Family = ifelse(Genus== "Marinoscillum" & Family == "Flammeovirgaceae","Reichenbachiellaceae",Family))

ALL_TAX.42.82<-ALL_TAX.42.81%>% 
  mutate(Family = ifelse(Genus== "Maritalea" & Family == "Hyphomicrobiaceae","Devosiaceae",Family))

ALL_TAX.42.83<-ALL_TAX.42.82%>% 
  mutate(Family = ifelse(Genus== "MD3-55" & Family == "Fokiniaceae","Midichloriaceae",Family))

ALL_TAX.42.84<-ALL_TAX.42.83%>% 
  mutate(Family = ifelse(Genus== "Methylibium" & Family == "Comamonadaceae","Burkholderiales_incertae_sedis",Family))

ALL_TAX.42.85<-ALL_TAX.42.84%>% 
  mutate(Family = ifelse(Genus== "Nevskia" & Family == "Sinobacteraceae","Nevskiaceae",Family))

ALL_TAX.42.86<-ALL_TAX.42.85%>% 
  mutate(Family = ifelse(Genus== "Nisaea" & Family == "Rhodospirillaceae","Thalassobaculaceae",Family))

ALL_TAX.42.87<-ALL_TAX.42.86%>% 
  mutate(Family = ifelse(Genus== "Paludibacter" & Family == "Porphyromonadaceae","Paludibacteraceae",Family))

ALL_TAX.42.88<-ALL_TAX.42.87%>% 
  mutate(Family = ifelse(Genus== "Parvibaculum" & Family == "Parvibaculaceae","Rhodobiaceae",Family))

ALL_TAX.42.89<-ALL_TAX.42.88%>% 
  mutate(Family = ifelse(Genus== "Phreatobacter" & Family == "Rhizobiales_incertae_sedis","Phreatobacteraceae",Family))

ALL_TAX.42.90<-ALL_TAX.42.89%>% 
  mutate(Family = ifelse(Genus== "Phyllobacterium" & Family == "Rhizobiaceae","Phyllobacteriaceae",Family))

ALL_TAX.42.91<-ALL_TAX.42.90%>% 
  mutate(Family = ifelse(Genus== "Pseudahrensia" & Family == "Phyllobacteriaceae","Ahrensiaceae",Family))

ALL_TAX.42.92<-ALL_TAX.42.91%>% 
  mutate(Family = ifelse(Genus== "Pseudovibrio" & Family == "Rhodobacteraceae","Stappiaceae",Family))

ALL_TAX.42.93<-ALL_TAX.42.92%>% 
  mutate(Family = ifelse(Genus== "Reichenbachiella" & Family == "Flammeovirgaceae","Reichenbachiellaceae",Family))

ALL_TAX.42.94<-ALL_TAX.42.93%>% 
  mutate(Family = ifelse(Genus== "Rhizobacter" & Family == "Comamonadaceae","Pseudomonadaceae",Family))

ALL_TAX.42.95<-ALL_TAX.42.94%>% 
  mutate(Family = ifelse(Genus== "Rhodoligotrophos" & Family == "Rhodobiaceae","Parvibaculaceae",Family))

ALL_TAX.42.96<-ALL_TAX.42.95%>% 
  mutate(Family = ifelse(Genus== "Rhodopirellula" & Family == "Planctomycetaceae","Pirellulaceae",Family))

ALL_TAX.42.97<-ALL_TAX.42.96%>% 
  mutate(Family = ifelse(Genus== "RS62 marine group" & Family == "Comamonadaceae","Burkholderiaceae",Family))

ALL_TAX.42.98<-ALL_TAX.42.97%>% 
  mutate(Family = ifelse(Genus== "Saccharofermentans" & Family == "Hungateiclostridiaceae","Ruminococcaceae",Family))

ALL_TAX.42.99<-ALL_TAX.42.98%>% 
  mutate(Family = ifelse(Genus== "Saccharophagus" & Family == "Alteromonadaceae","Cellvibrionaceae",Family))

ALL_TAX.43<-ALL_TAX.42.99%>% 
  mutate(Family = ifelse(Genus== "Sulfuricurvum" & Family == "Sulfurimonadaceae","Helicobacteraceae",Family))

ALL_TAX.43.1<-ALL_TAX.43%>% 
  mutate(Family = ifelse(Genus== "Sulfuricurvum" & Family == "Thiovulaceae","Helicobacteraceae",Family))

ALL_TAX.43.2<-ALL_TAX.43.1%>% 
  mutate(Family = ifelse(Genus== "Sulfurimonas" & Family == "Thiovulaceae","Helicobacteraceae",Family))

ALL_TAX.43.3<-ALL_TAX.43.2%>% 
  mutate(Family = ifelse(Genus== "Sulfurimonas" & Family == "Sulfurimonadaceae","Helicobacteraceae",Family))

ALL_TAX.43.4<-ALL_TAX.43.3%>% 
  mutate(Family = ifelse(Genus== "Sulfurospirillum" & Family == "Sulfurospirillaceae","Campylobacteraceae",Family))

ALL_TAX.43.5<-ALL_TAX.43.4%>% 
  mutate(Family = ifelse(Genus== "Sulfurovum" & Family == "Sulfurovaceae","Helicobacteraceae",Family))

ALL_TAX.43.6<-ALL_TAX.43.5%>% 
  mutate(Family = ifelse(Genus== "Sva0081 sediment group" & Family == "Desulfosarcinaceae","Desulfobacteraceae",Family))

ALL_TAX.43.7<-ALL_TAX.43.6%>% 
  mutate(Family = ifelse(Genus== "Terasakiella" & Family == "Methylocystaceae","Terasakiellaceae",Family))

ALL_TAX.43.8<-ALL_TAX.43.7 %>% 
  mutate(Family = ifelse(Genus== "Teredinibacter" & Family == "Alteromonadales_incertae_sedis","Cellvibrionaceae",Family))

ALL_TAX.43.9<-ALL_TAX.43.8 %>% 
  mutate(Family = ifelse(Genus== "Thiovulum" & Family == "Thiovulaceae","Sulfurimonadaceae",Family))

ALL_TAX.43.10<-ALL_TAX.43.9 %>% 
  mutate(Family = ifelse(Genus== "Desulfosarcina" & Family == "Desulfosarcinaceae","Desulfobacteraceae",Family))

ALL_TAX.43.11<-ALL_TAX.43.10%>% 
  mutate(Order= ifelse(Family== "Ahrensiaceae" & Order == "Rhodobacterales", "Rhizobiales",Order))

ALL_TAX.43.12<-ALL_TAX.43.11%>% 
  mutate(Order= ifelse(Family== "Alcanivoracaceae" & Order == "Pseudomonadales", "Oceanospirillales",Order))

ALL_TAX.43.13<-ALL_TAX.43.12%>% 
  mutate(Order= ifelse(Family== "Balneolaceae" & Order == "Chitinophagales", "Balneolales",Order))

ALL_TAX.43.14<-ALL_TAX.43.13%>% 
  mutate(Order= ifelse(Family== "Bryobacteraceae" & Order == "Solibacterales", "Bryobacterales",Order))

ALL_TAX.43.15<-ALL_TAX.43.14%>% 
  mutate(Order= ifelse(Family== "Cellvibrionaceae" & Order == "Alteromonadales", "Cellvibrionales",Order))

ALL_TAX.43.16<-ALL_TAX.43.15%>% 
  mutate(Order= ifelse(Family== "Clostridiales_Incertae Sedis XII" & Order == "Peptostreptococcales-Tissierellales", "Clostridiales",Order))

ALL_TAX.43.17<-ALL_TAX.43.16%>% 
  mutate(Order= ifelse(Family== "Coxiellaceae" & Order == "Diplorickettsiales", "Legionellales",Order))

ALL_TAX.43.18<-ALL_TAX.43.17%>% 
  mutate(Order= ifelse(Family== "Desulfobacteraceae" & Order == "Desulfatiglandales", "Desulfobacterales",Order))

ALL_TAX.43.19<-ALL_TAX.43.18%>% 
  mutate(Order= ifelse(Family== "Desulfobulbaceae" & Order == "Desulfobulbales", "Desulfobacterales",Order))

ALL_TAX.43.20<-ALL_TAX.43.19%>% 
  mutate(Order= ifelse(Family== "Gemmataceae" & Order == "Planctomycetales", "Gemmatales",Order))

ALL_TAX.43.21<-ALL_TAX.43.20%>% 
  mutate(Order= ifelse(Family== "Halieaceae" & Order == "Alteromonadales", "Cellvibrionales",Order))

ALL_TAX.43.22<-ALL_TAX.43.21%>% 
  mutate(Order= ifelse(Family== "Izimaplasmataceae" & Order == "Izemoplasmatales", "Izimaplasmatales",Order))

ALL_TAX.43.23<-ALL_TAX.43.22%>% 
  mutate(Order= ifelse(Family== "Nevskiaceae" & Order == "Xanthomonadales", "Nevskiales",Order))

ALL_TAX.43.24<-ALL_TAX.43.23%>% 
  mutate(Order= ifelse(Family== "Pirellulaceae" & Order == "Planctomycetales", "Pirellulales",Order))

ALL_TAX.43.25<-ALL_TAX.43.24%>% 
  mutate(Order= ifelse(Family== "Pseudomonadaceae" & Order == "Burkholderiales", "Pseudomonadales",Order))

ALL_TAX.43.26<-ALL_TAX.43.25%>% 
  mutate(Order= ifelse(Family== "Spongiibacteraceae" & Order == "Pseudomonadales", "Cellvibrionales",Order))

ALL_TAX.43.27<-ALL_TAX.43.26%>% 
  mutate(Order= ifelse(Family== "Stappiaceae" & Order == "Rhodobacterales", "Rhizobiales",Order))

ALL_TAX.43.28<-ALL_TAX.43.27%>% 
  mutate(Order= ifelse(Family== "Terasakiellaceae" & Order == "Rhizobiales", "Rhodospirillales",Order))

ALL_TAX.43.29<-ALL_TAX.43.28%>% 
  mutate(Class= ifelse(Order== "Balneolales" & Class == "Chitinophagia", "Balneolia",Class))

ALL_TAX.43.30<-ALL_TAX.43.29%>% 
  mutate(Class= ifelse(Order== "Gemmatales" & Class == "Planctomycetes", "Planctomycetia",Class))

ALL_TAX.43.31<-ALL_TAX.43.30%>% 
  mutate(Class= ifelse(Order== "Izimaplasmatales" & Class == "Bacilli", "Mollicutes",Class))

ALL_TAX.43.32<-ALL_TAX.43.31%>% 
  mutate(Class= ifelse(Order== "Desulfobacterales" & Class == "Desulfobacteria", "Deltaproteobacteria",Class))

ALL_TAX.43.33<-ALL_TAX.43.32%>% 
  mutate(Class= ifelse(Order== "Desulfobacterales" & Class == "Desulfobulbia", "Deltaproteobacteria",Class))

ALL_TAX.43.34<-ALL_TAX.43.33%>% 
  mutate(Phylum= ifelse(Class== "Balneolia" & Phylum == "Bacteroidetes", "Balneolaeota",Phylum))

ALL_TAX.43.35<-ALL_TAX.43.34%>% 
  mutate(Phylum= ifelse(Class== "Planctomycetia" & Phylum == "Planctomycetes", "Planctomycetes",Phylum))

ALL_TAX.43.36<-ALL_TAX.43.35%>% 
  mutate(Phylum= ifelse(Class== "Mollicutes" & Phylum == "Firmicutes", "Tenericutes",Phylum))

ALL_TAX.43.37<-ALL_TAX.43.36%>% 
  mutate(Phylum= ifelse(Class== "Deltaproteobacteria" & Phylum == "Desulfobacterota", "Proteobacteria",Phylum))

ALL_TAX.43.38<-ALL_TAX.43.37%>% 
  mutate(Phylum= ifelse(Class== "Deltaproteobacteria" & Phylum == "Desulfobacterota", "Proteobacteria",Phylum))

ALL_TAX.43.39<-ALL_TAX.43.38%>% 
  mutate(Genus= ifelse(Species== "Cenchrus americanus" & Genus == "Ideonella", "Aquabacterium",Genus))

####Find duplicates####

tax<-ALL_TAX.43.39

#write.csv2(ALL_TAX,"Edited_ALL_TAX_table.csv" )

ALL_TAX.2<-as.matrix(tax)
otu<-read.csv2("ALL_OTU.csv", header=TRUE)
otu<-column_to_rownames(otu, "OTU")
mapfile = "Sample_data_all.csv"
map <- read.csv2(mapfile)
rownames(map) <- map$Sample_name
map<-map[,-1]              

ALL_TAX_final<-tax_table(ALL_TAX.2)
#write.csv2(ALL_TAX_final, "All_Final_TAX.csv")
otu_final<- otu_table(otu, taxa_are_rows = TRUE)
#write.csv2(otu_final, "All_Final_OTU.csv")

MAP<-sample_data(map)

data<-merge_phyloseq(ALL_TAX_final,otu_final,MAP)


data.1 <- subset_taxa(data, Kingdom == "Bacteria" & Class != "Chloroplast" & Order !="Chloroplast" &Family != "mitochondria" &  Phylum != "Cyanobacteria/Chloroplast")
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
dup_phyla

#Class duplicates
data_class_t<-as.data.frame(tax_table(data_class))
dup_class<- data_class_t %>% filter(duplicated(data_class_t$Class)| duplicated(Class, fromLast = TRUE))
dup_class

#Order duplicates
data_order_t<-as.data.frame(tax_table(data_order))
dup_order<- data_order_t %>% filter(duplicated(data_order_t$Order)| duplicated(Order, fromLast = TRUE))
dup_order

#Family duplicates
data_family_t<-as.data.frame(tax_table(data_family))
dup_Family<- data_family_t %>% filter(duplicated(data_family_t$Family)| duplicated(Family, fromLast = TRUE))
dup_Family

#Genus duplicates
data_genus_t<-as.data.frame(tax_table(data_genus))
dup_Genus<- data_genus_t %>% filter(duplicated(data_genus_t$Genus)| duplicated(Genus, fromLast = TRUE))
dup_Genus

#Species duplicates
data_species_t<-as.data.frame(tax_table(data_species))
dup_species<- data_species_t %>% filter(duplicated(data_species_t$Species)| duplicated(Species, fromLast = TRUE))
dup_species
