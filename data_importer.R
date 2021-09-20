data_importer <- function(datafile, taxfile) {
#read data
  dat <- read.table(datafile, header=TRUE)
 
  
#read tax file
  tax <- read.delim(taxfile, header=TRUE)
  tax <- subset(tax, select=-c(confidence))
  
#fix tax file so that the format fits qiime2 (remove D__, ...)  
  tax  <- separate(tax, taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";", fill="right")
  tax <- separate(tax, Kingdom, into = c(NA, "Kingdom"), sep="__", fill="right")
  tax <- separate(tax, Phylum, into = c(NA, "Phylum"), sep="__", fill="right")
  tax <- separate(tax, Class, into = c(NA, "Class"), sep="__", fill="right")
  tax <- separate(tax, Order, into = c(NA, "Order"), sep="__", fill="right")
  tax <- separate(tax, Family, into = c(NA, "Family"), sep="__", fill="right")
  tax <- separate(tax, Genus, into = c(NA, "Genus"), sep="__", fill="right")
  tax <- separate(tax, Species, into = c(NA, "Species"), sep="__", fill="right")

#concatenate all taxonomy into 1 coumn, separated by semi colons:
  tax$taxonomy<- paste(tax$Kingdom,tax$Phylum, tax$Class, tax$Order, tax$Family, tax$Genus, tax$Species, sep=";")
  
#remove Kingdom, ... columns:
  tax <- subset(tax, select=-c(Kingdom, Phylum, Class, Order, Family, Genus, Species))

#Merge tax and OTU table:
  data <- merge(dat, tax, by="OTUID" )


#Sum by tax name:
data <- ddply(data,"taxonomy",numcolwise(sum))

#repeat for all samples, then merge by taxonomy...sum by taxonomy again. split taxonomy names into columns again, 
#Ready to read into phyloseq.
}