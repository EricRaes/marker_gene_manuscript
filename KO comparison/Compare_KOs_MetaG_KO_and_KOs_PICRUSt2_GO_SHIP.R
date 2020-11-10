library(castor)
library(parallel)
library(stringi)
library(dplyr)
library(tidyr )
library(microbiome)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(broom)
library(ggpmisc)


############
#1# First select only the stations were there is both shotgun MGS and 16S rRNA gene sequence data
#2# Then find the overlap between the KOs from shotgun MGS and PICRUSt2 predictions (as done by Douglas et al. 2020; https://www.nature.com/articles/s41587-020-0548-6)
#3# Select your metabolic function of interest and all the corresponding KOs
#4# Make Spearman correlation plots between KOs from shotgun NGS and PICRUSt2 predictions for each station
############

#1# select only the samples where you have shotgun data and PICRUSt2 KO data
stations <- read.csv("./samples.csv",sep = ",", header=TRUE)
head(stations)

KO_PICRSUT2_All_samples <- read.csv("./KO_predicitions_PICRUSt2_97_percent_OTU_no_mitochondria_no_chloroplasts.csv",sep = ",", header=TRUE)
head(KO_PICRSUT2_All_samples)

KO_Table_long <- KO_PICRSUT2_All_samples %>% gather(ANDS, abundance, -c("KO"))
head(KO_Table_long)

Shotgun_stations_only<- merge(stations, KO_Table_long, by="ANDS",all.x = FALSE) 
head(Shotgun_stations_only)

Shotgun_stations_only_Matrix<-Shotgun_stations_only %>% spread(key = "ANDS", value = abundance)
head(Shotgun_stations_only_Matrix)

write.csv(Shotgun_stations_only_Matrix, "PICRUSt2_data_whith_Shotgun_stations_only.csv")


#2# find the overlap between the KOs from shotgun MGS and PICRUSt2 predictions;
#   As detailed by Douglas et al. 2020 - the script below are from https://github.com/gavinmdouglas/picrust2_manuscript
#   https://www.nature.com/articles/s41587-020-0548-6

# overlap KOs between MEtaG and PICRUSt2
possible_picrust2_kos <- read.table("./PICRUSt2_ko.txt", header=F, stringsAsFactors = FALSE)$V1
possible_mgs_kos <- read.table("./humann2_ko.txt", header=F, stringsAsFactors = FALSE)$V1

# Identify subset of KOs that could have been output by all approaches.
overlapping_possible_kos <- possible_picrust2_kos[which(possible_picrust2_kos %in% possible_mgs_kos)]
overlapping_possible_kos <-data.frame(overlapping_possible_kos)
colnames(overlapping_possible_kos)[1] <- "KO"
write.csv(overlapping_possible_kos, "overlapping_possible_kos.csv")

# read in MetaG KOs + PICRUSt2 KOs + Overlap file
picrust2_ko_nsti2 <- read.csv("PICRUSt2_data_whith_Shotgun_stations_only.csv", header=T)
mgs_ko <- read.csv("KO_profiled_from_shotgun_Metagenome_data.csv", header=T)
overlap <- read.csv( "overlapping_possible_kos.csv")

# Merge MetaG KOs with overlap KO file; then also merge PICRUSt2 KO prediction file
# MAKE SURE you open the files and check headings e.g., KO and/or function and sample id column names are all the same.
data_merge_overlap_MetaG <- merge(mgs_ko,overlap , by = "KO", all.y=TRUE)
data_merge_overlap_MetaG_Picrust2 <- merge(picrust2_ko_nsti2, data_merge_overlap_MetaG, by = "KO", all.y=TRUE)
d <- as.data.frame(data_merge_overlap_MetaG_Picrust2)
d[is.na(d)] <- 0 #convert NAs to 0
write.csv(d, "data_merge_overlap_MetaG_Picrust2.csv", row.names=FALSE)


#3# now do Spearman correlation for all teh KOs at each station where we have shotgun MGS and 16S rRNA gene data
## The below code will subset one of the 11 shotgun samples and then subet the KO and pathway table to a specific function of interest
## Apologies in advance for the redundancy in the code. I am aware that teh code can be opimised with a loop

overlap_MetaG_Picr_KO <- read.csv( "data_merge_overlap_MetaG_Picrust2.csv",header=T )
overlap_MetaG_Picr_KO_sample <- select(overlap_MetaG_Picr_KO,"KO", "G34369.x","G34369.y")

#Sample_ids 
## G34369	#G34393	#G34402	#G34415	#G34421	#G34430	#G34450	#G34458	#G34473
## G34481	#G34499
## just copy paste the sample ids in the 'select' line above but leave kepe the .x and .y
## .x means KOs predicted by PICRUSt2
## .y means KOs profiled from shotgun metagenomes 

## now read in the KOs profiled from shotgun metagenomes data as a phyloseq object
## Three files
##  1) KO abundance file
##  2) KO pathways, e.g., similar to a taxonomy file
##  3) metadata file
otu_mat <- "./Shotgun_Stations_GO_SHIP_KOs.csv"
tax_mat <- "./KO_numbers_with_pathways_names.csv"
samples_df <- "./META_GO_SHIP.csv"
file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

KO_MGS_GO_SHIP<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

# Now subset to the specific function of interest and all the KOs that cover this function 
KO_MGS_GO_SHIP_Function = subset_taxa(KO_MGS_GO_SHIP, TAX2==c(" Metabolism of cofactors and vitamins"))

# the other functions presented in the manuscript are shown below and can be copy pasted in the subset_taxa function above.
## TAX2==c(" Lipid metabolism"))
## TAX2==c(" Amino acid metabolism"))
## TAX2==c(" Carbohydrate metabolism"))
## TAX2==c(" Biosynthesis of other secondary metabolites"))
## TAX2==c(" Metabolism of cofactors and vitamins"))
## TAX2==c(" Nucleotide metabolism"))
## TAX2==c(" Cell growth and death"))
## TAX2==c(" Energy metabolism"))
# Other functions
##carbom_Functions = subset_taxa(carbom, TAX3a==c("Fermentation"))
##carbom_Functions = subset_taxa(carbom, TAX3a==c("Methane metabolism"))
##carbom_Functions = subset_taxa(carbom, TAX3a==c("Sulfur metabolism"))
##carbom_Functions = subset_taxa(carbom, TAX3a==c("Carbon fixation pathways in prokaryotes"))
##carbom_Functions = subset_taxa(carbom, TAX3a=="Carbon fixation in photosynthetic organisms"|
##                                       TAX3a=="Photosynthesis - antenna proteins"|
##                                       TAX3a=="Photosynthesis"|
##                                       TAX3a=="Photosynthesis proteins"|
##                                       TAX3a=="Carbon fixation pathways in prokaryotes")


# now extract the subsetted (to your function of interest)'KO table' using 'tax_table'
GO_SHIP_Function_selected<-tax_table(KO_MGS_GO_SHIP_Function)
write.csv(GO_SHIP_Function_selected,"GO_SHIP_Function_selected.csv")

# load the selected function file of interest and merge with the overlap file created in #2#
subset_function_get_KOs<-read.csv("GO_SHIP_Function_selected.csv")
merge_KOs_Picurst2_Shotgun_subset<-merge(subset_function_get_KOs, overlap_MetaG_Picr_KO_sample, by="KO", all.x = TRUE )
d <- as.data.frame(merge_KOs_Picurst2_Shotgun_subset)
d[is.na(d)] <- 0 #convert NAs to 0
d <- d [-c(2:229)] #get rid of extra TAX columns

write.csv(d, "d.csv")
d1<- read.csv("d.csv",header=T)
d1<-d1[,-1]

d_long <- d1 %>% gather(Station, abundance, -c("KO"))
head(d_long)

d_long_split<-d_long %>% separate(Station, c("Station", "Method"))
head(d_long_split)
dlong_PICrust<-subset(d_long_split, Method=="x") # add method x inferred KOs from PCIRUSt2
dlong_Shotgun<-subset(d_long_split, Method=="y") # add method y profiled KOs from shotgun MGS data

KO_abundance_PI_SHOT<-merge(dlong_PICrust,dlong_Shotgun, by=c("KO","Station")) # merge

#4# make Spearman Correlation plots for each function and for each sample
## repeat the script above and change the sample names and subset to each specific function
## Not the fastest way but it does the job for 11 samples

p<-ggscatter(KO_abundance_PI_SHOT, x = "abundance.x", y = "abundance.y", 
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regression line
          add.params = list(color = "red", fill = "lightgray"), # Customize reg. line
          xlab = "KOs_PICRUSt2", ylab = "KOs_MetaG")+
          ggtitle("sample 34369 \nFunctional Pathway of interest")+
          stat_cor(method = "spearman", label.sep = sprintf(", n = %s, ",
          sum(complete.cases(KO_abundance_PI_SHOT[c("abundance.x","abundance.y")]))))
p          

## plot the Spearman correlation coeficients

Corr<-read.csv("ALL_SPEARMAN_CORR_samples.csv", header=TRUE)
colnames(Corr)

dat.m <- reshape2::melt(Corr,id.vars='Sample_id', measure.vars=c("Amino.acid.metabolism","Sulfur.Metabolism",
"Carbohydrate.metabolism",
"Carbon.fixation.pathways.in.prokaryotes", "Carbon.fixation.pathways","Cell.growth.and.death",
"Energy.metabolism", "Lipid.metabolism", "metabolism",
"Methane.metabolism", "Nucleotide.metabolism"))

dat.m$variable = factor(dat.m$variable, levels=c("Carbon.fixation.pathways.in.prokaryotes", 
                                                 "Carbon.fixation.pathways",
                                                 "Lipid.metabolism",
                                                 "Cell.growth.and.death",
                                                 "Methane.metabolism",
                                                 "Energy.metabolism", 
                                                 "Carbohydrate.metabolism",
                                                 "Sulfur.Metabolism",
                                                 "Amino.acid.metabolism", 
                                                 "Nucleotide.metabolism",
                                                 "Metabolism.of.cofactors.and.vitamins",
                                                 "metabolism"),
                                       labels=c( "Carbon fixation pathways \n(prokaryotes)", 
                                                 "Carbon fixation pathways",
                                                 "Lipid metabolism",
                                                 "Cell growth and death",
                                                "Methane metabolism",
                                                 "Energy metabolism", 
                                                "Carbohydrate metabolism",
                                                "Sulfur Metabolism",
                                                 "Amino acid metabolism", 
                                                "Nucleotide metabolism",
                                                 "Cofactor and vitamin metabolism",
                                                "All Metabolism Pathways"))

# Plot the data
p <- ggplot(dat.m, aes(x=value, y=variable))
p+geom_boxplot(fill = "grey80")+
  geom_point(color="black", size=0.4, alpha=0.9) +
  scale_x_continuous( limits=c(0, 1), breaks=seq(0,1,0.2)) +
    theme_bw()+ ylab("")+xlab("Spearman coefficient ")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=16, vjust = 0.3),
        axis.title.x = element_text(size=15, vjust = 0.9),
              axis.title.y = element_text(size=14, vjust = 0.9),
              axis.text.y = element_text(size=14, vjust = 0.9),
              legend.position = "none",
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Metabolic Functions")


#####