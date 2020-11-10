
#load libraries
library("phyloseq"); packageVersion("phyloseq")
library(tidyr )
library(ggpubr)
library(vegan)
library(scales)
library(grid)
library(reshape2)
library(ggplot2)
library(vegan)
library(data.table)
library(plyr)
library(dplyr)
library(microbiome)


#####load PICRUSt2 output into a Phyloseq object
#1. Pathway abundance file
#2. MetaCyc pathway file;
#To link the metaCyc pathway identifiers to their respective functions at the secondary superclass level
#using MetaCyc pathway hierarchy system go to https://metacyc.org/ and create a SmartTable (see SmartTables menu)
#by importing those pathways from a file. Then click the menu "Add Transform Column" and choose the transform "Ontology -- Parents of Class".
#Then click on the header of that new column, and again select that same transform to create another column that goes up one level in the ontology.
#To further narrow the Parents of Class pathways down, click on the header of that new column, and again select that same transform to create another column that goes up one level in the ontology."
#3. Metadata file

otu_mat <- "./Metacyc_Pathways_Abundances_PICRUSt2.csv"
tax_mat <- "./Metacyc_Pathways_PICRUSt2.csv"
samples_df <- "./Metadata_GO_SHIP_P15S_PICRUSt2.csv"

file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

GO_SHIP<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

########## filter pathways with less than 10 occurences and rarify to min(sample_sums)
GO_SHIP_10 = filter_taxa(GO_SHIP, function(x) sum(x) > 10, TRUE)
GO_SHIP_rar <-rarefy_even_depth(GO_SHIP_10, sample.size = min(sample_sums(GO_SHIP_10)),
                              rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)
sample_sums(GO_SHIP_rar)
taxa_names(GO_SHIP_rar)

########## Plot beta-diveristy for the functional community

my.cols <- c("blue", "orange", "darkgreen", "red")
sample_data(GO_SHIP_rar)$Water.mass = factor(sample_data(GO_SHIP_rar)$Water.mass, levels = c("Southern Ocean","STF","SPSG","PED"))

logt=transform_sample_counts(GO_SHIP_rar, function(x) sqrt(x))
out.pcoa.logt <- ordinate(logt, method = "CAP", distance = "bray", formula=~Water.mass)
evals <- out.pcoa.logt$values$Eigenvalues
p3<-plot_ordination(logt, out.pcoa.logt, type = "samples", 
                    color = "Water.mass",shape = "Depth")
p3 + scale_colour_manual(values=my.cols) + theme_bw()+ geom_point(size=3)+
  ggtitle("Beta diversity - 3 depths within in the Mixed Layer Depth") +
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  

######### Plot alpha-diveristy for the functional community

p<-plot_richness(GO_SHIP_rar,x="Latitude..decimal.degrees.", color="Water.mass",
                 measures=c("Chao1"))
p + theme_bw()+geom_point(shape = 1,size = 1.7,colour = "black")+
  scale_color_manual(values=c("red", "blue", "darkgreen", "orange"))+
  ggtitle("Alpha diversity - 3 depths within in the Mixed Layer Depth") +
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name="Latitude", breaks=seq(-70, 5, 10))+
  scale_y_continuous(name="Alpha Diversity Measure", limits=c(280, 380),breaks=seq(280,380,40))+
  theme(axis.text=element_text(size=12),
  axis.title=element_text(size=13,face="bold"))+
  theme(panel.grid.major = element_blank())+stat_smooth(method = "lm", formula = y ~ x, se = FALSE)


######### Plot latitudinal changes

#convert to relative abudances
GO_SHIP_rel = transform_sample_counts(GO_SHIP_rar, function(x) x / sum(x) )

#now subset to a specific secondary superclass e.g., Amino-Acid-Biosynthesis
GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Amino-Acid-Biosynthesis")) #1

#then go through the list below
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Nucleotide-Biosynthesis")) #2
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Energy-Metabolism"))       #3
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Lipid-Biosynthesis"))      #4
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Carbohydrates-Biosynthesis")) #5
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Cell-Structure/ Cell-Wall-Biosynthesis")) #6
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Cofactor-Biosynthesis"))   #7
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Secondary-Metabolite-Biosynthesis")) #8
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Vitamin-Biosynthesis"))     #9
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("Fermentation"))       #10
#GO_SHIP_rel_MetaCyc = subset_taxa(GO_SHIP_rel, Pathway==c("CO2-Fixation")) #11

#now glomerate 
glom <- tax_glom(GO_SHIP_rel_MetaCyc, taxrank = 'Pathway')
glom # should list MetaCyc pathway identifiers as secondary superclass levels 
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Pathway <- as.character(data_glom$Pathway) #convert to character
data_glom$Water.mass <- factor(data_glom$Water.mass, levels = c("Southern Ocean","STF","SPSG","PED"))
abundance_percent<-data_glom$Abundance*100 #convert to percentages
#write.csv(data_glom, "xxx_pathway_MetaCyc.csv") #export glomerated pathway file if you want

# Plot data as points
p<-ggplot(data_glom, aes(x=Latitude..decimal.degrees. , y=abundance_percent, fill=Water.mass)) 
p <- p +geom_point(shape = 21,size = 2,colour = "black")+
  scale_fill_manual(values=c("blue", "orange", "green", "red")) +
  ggtitle("Latitudinal change of a functional pathway of interest") +
  scale_x_continuous(name="Latitude", limits=c(-70, 10), breaks=seq(-70,5,10)) +
  scale_y_continuous(name="Relative Abundance (%)", labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'))+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

########### Wilcoxon Statistics and averages 
compare_means(Abundance ~ Water.mass,  data = data_glom, p.adjust.method = "bonferroni")
ddply(data_glom, ~Water.mass, plyr:::summarise, mean = mean(Abundance), sd = sd(Abundance),
      max = max(Abundance), min = min(Abundance))

# Plot data as boxplots with jitter
p<-ggplot(data_glom, aes(x=Latitude..decimal.degrees. , y=abundance_percent, fill=Water.mass))
p + geom_boxplot(notch=FALSE, outlier.shape=NA)+
  geom_jitter(shape = 1,colour="black", size = 1, width = 0.05, height=0.05)+
  scale_fill_manual(values=c("blue", "orange", "green", "red")) +
  ggtitle("Latitudinal change of a functional pathway of interest") +
  scale_x_continuous(name="Latitude", limits=c(-70, 10), breaks=seq(-70,5,10)) +
  scale_y_continuous(name="Relative Abundance (%)",labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'))+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Plot Cell-Structure/ Cell-Wall-Biosynthesis vs temperature 
p<-ggplot(data_glom, aes(x=T , y=abundance_percent, fill=Water.mass)) #Latitude..decimal.degrees.
p <- p +geom_point(shape = 21,size = 2,colour = "black")+#geom_box(fill=Water.mass)
  scale_fill_manual(values=c("blue", "orange", "green", "red")) +
  ggtitle("Cell-Structure/ Cell-Wall-Biosynthesis") +
    scale_x_continuous(limits=c(-2, 32), breaks=seq(-2,32,5))+
  scale_y_continuous(name="Relative Abundance (%)", labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'))+
    theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p + labs(x = expression(paste("Temperature ", (C^o))))


################ INDISPECIES
 
glom <- tax_glom(GO_SHIP_rar, taxrank = 'Pathway') # Glomerate to 'Pathway'
glom # should list MetaCyc pathway identifiers as secondary superclass levels 

Y <- veganotu(glom)
write.csv(Y, "Indispecies_Pathways_glom.csv")

dat <- data.frame(sample_data(glom))
write.csv(dat, "Indispecies_Pathways_glom_Metadata.csv")

#then open file "Indispecies_Pathways_glom.csv" and add factor 'water.mass' from the metadata file
#this will file will now be the input file for the indicspecies analysis.

install.packages("indicspecies")
library(indicspecies)

Indic_Spec = read.csv("Indispecies_Pathways_glom.csv", header= TRUE)

abund = Indic_Spec[,3:ncol(Indic_Spec)]
#abund = pc_all[,3:ncol(pc_all)]
Water.mass = Indic_Spec$Water.mass
#Water.mass = pc_all$Water.mass

Indic_Spec_run = multipatt(abund, Water.mass, func = "r.g", control = how(nperm=9999))
summary(Indic_Spec_run)


############### Boosted Regression Tree Modelling (BRT)
#install.packages("dismo")
library(dismo)
#install.packages("gbm")
library(gbm)

#make phyloseq object from pathways; rarefy; get relative abundances; then glom each pathway and make a new file for that
# The relative abudance for each pathwway would be saved as a .csv file and can then be used to
# to predict/estimate which parameters contribute to a certain pathway
#get n= as for this data set we will only have prob ~n=75 from the sfc were we have PON and POC data

GO_SHIP_rar_rel = transform_sample_counts(GO_SHIP_rar, function(x) x / sum(x) )
GO_SHIP_rar_rel_depth_1 = subset_samples(GO_SHIP_rar_rel, Depth_BIN==1)

BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("CO2-Fixation"))
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Amino-Acid-Biosynthesis")) #1
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Nucleotide-Biosynthesis")) #2
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Energy-Metabolism"))       #3
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Lipid-Biosynthesis"))      #4
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Carbohydrates-Biosynthesis")) #5
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Cell-Structure/ Cell-Wall-Biosynthesis")) #6
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Cofactor-Biosynthesis"))   #7
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Secondary-Metabolite-Biosynthesis")) #8
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Vitamin-Biosynthesis"))     #9
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Sulfur-Metabolism")) #10
#BRT_GO_SHIP = subset_taxa(GO_SHIP_rar_rel_depth_1, Pathway==c("Fermentation"))  #11

glom <- tax_glom(BRT_GO_SHIP, taxrank = 'Pathway')
glom 
data_glom<- psmelt(glom) # create dataframe from phyloseq object
write.csv(data_glom, "PATHWAY_for_BRT.csv")

Pathways<-read.csv(file="PATHWAY_for_BRT.csv",header = TRUE)

## for Pathways data remove NAs
Pathways_cols_1_59 <- subset(Pathways, select = c(1:21,27,29,38,41:42,45:59))
Pathways_no_NAs <- na.omit(Pathways_cols_1_59)

# Run BRT for Pathways
Pathways_no_NAs_BRT <- gbm.step(data=Pathways_no_NAs, gbm.x = 14:35, gbm.y = 4, family = "gaussian",
                                tree.complexity = 10, learning.rate = 0.001, bag.fraction = 0.5)

output_Pathways_1<-summary(Pathways_no_NAs_BRT)
summary(Pathways_no_NAs_BRT)
names(Pathways_no_NAs_BRT)
write.csv(Pathways_no_NAs_BRT, file ="output_Pathways_Picrust2.3.0b_specific_pathway.csv")

## for primary productivity (P.P.) data remove NAs
carbom_PP= meta(carbom_3)
Pathways_cols_1_59 <- subset(carbom_PP, select = c(1:21,27,29,38,41:42,45:59))
Pathways_no_NAs_PP <- na.omit(Pathways_cols_1_59)

# Run BRT for primary productivity (P.P.)
Pathways_no_NAs_BRT_PP <- gbm.step(data=Pathways_no_NAs_PP, gbm.x = 10:31, gbm.y = 34, family = "gaussian",
                                tree.complexity = 10, learning.rate = 0.005, bag.fraction = 0.5,
                                tolerance.method = "fixed", tolerance = 0.001)

output_Pathways_1<-summary(Pathways_no_NAs_BRT_PP)
summary(Pathways_no_NAs_BRT_PP)
names(Pathways_no_NAs_BRT_PP)

########
