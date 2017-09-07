# Visual genetic characterization of NAM population for manuscript
# February 20, 2017
# Alex Ollhoff

library(ggplot2)
library(dplyr)
library(reshape2)

# Pairwise genetic distance between Ras and donor parents
#Distance <- read.csv("~/Documents/PhD/NAM/Population_description/Genetic_distance/NAM_parents_with_genetic_distance.csv", header = T)
Distance <- read.csv("~/Desktop/data_for_thesis/NAM_parents_with_genetic_distance.csv", header = T)
head(Distance)

Distance$Pop_location_sorted = factor(Distance$Pop_location, levels = c('Central European', 'Coastal Mediterranean', 'East African', 'Asian', 'Admixed', 'Un'))

# Sort by proportion disimilarity
Distance_by_diff <- Distance[order(Distance[,11]),]

# Add new sort column
Distance_by_diff$Sort_by_diff <- c(1:88)

# Plot
ggplot() +
  labs(x = "Donor Parents", y = "Proportion of Variants\nDiffering from Rasmusson") +
  theme(strip.background = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(), legend.position = "none", axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
  geom_bar(data = Distance_by_diff, aes(x = as.numeric(Sort_by_diff), y = as.numeric(Proportion_disimilarity_toRas), fill = factor(Pop_location)), na.rm = T, stat = "identity") +
  scale_fill_manual(name = "NAM Parent\nSubpopulations", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red", "grey55")) +
  #facet_grid( Pop_location_sorted ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_continuous(breaks = c(0, 0.010, 0.020, 0.030, 0.040), labels = c("0", "0.01", "0.02", "0.03", "0.04"), limits = c(0, 0.042), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), labels = Distance_by_diff$NAM_names) +
  coord_flip() 
# Save as portrait 3x7.5 device size in Population_description as "Pairwise_distance_by either genetic/pheno/overalldist"

# Quantify trend
Distance_few <- select(Distance, c(line_name, Proportion_disimilarity_toRas, Pop_location))

gendist <- lm(Proportion_disimilarity_toRas ~ factor(Pop_location), data = Distance_few)
anova(gendist)
summary(gendist)
####################################################################################################################################
####################################################################################################################################
###>>>*** Segregation distortion ***>>>###
rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape2)
Seg_dist <- read.table("/Users/agonzale/Documents/SmithLab/NAM/Analysis/SNPdensity_and_segDist/88Fam_genmap_segDist_binned.txt", header = T)

# Separate chromsome number from bin label
Chrom_num <- gsub("_[0-9]*","", row.names(Seg_dist))
Seg_dist$chrom <- Chrom_num

# Add bin name as a column 
Seg_dist$bin <- row.names(Seg_dist)

Melted_dist <- melt(Seg_dist, id.vars = c("bin", "chrom"))

# Replace names of chromsomes
Melted_chrom <- Melted_dist %>%
  mutate(chrom = replace(chrom, chrom=="CHR1", "1H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR2", "2H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR3", "3H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR4", "4H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR5", "5H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR6", "6H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR7", "7H")) %>%
  as.data.frame()

# Add ID column
Melted_chrom$ID <- 1:nrow(Melted_chrom)

# Sort and then group families by subpop
# Import family information, with RILs and parents to know which color to paint them
family_info <- read.csv("/Users/agonzale/Dropbox/GITHUB/BRIDG6/Datasets/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA")

# Separate family in populations
Admix<-subset(family_info, family_info$Pop_location == "Admixed")
coastal_m<-subset(family_info, family_info$Pop_location == "Coastal Mediterranean")
Asian<-subset(family_info, family_info$Pop_location == "Asian")
central_e<-subset(family_info, family_info$Pop_location == "Central European")
east_af<-subset(family_info, family_info$Pop_location == "East African")
unassigned<-subset(family_info, family_info$Pop_location == "Un")

#get identifier for individuals in a family
List_indiv_admix<-c(as.character(Admix$Family), as.character(Admix$NAM_name))
List_indiv_coastal_m <-c(as.character(coastal_m$Family), as.character(coastal_m$NAM_name))
List_indiv_Asian <-c(as.character(Asian$Family), as.character(Asian$NAM_name))
List_indiv_central_e <-c(as.character(central_e$Family), as.character(central_e$NAM_name))
List_indiv_east_af <-c(as.character(east_af$Family), as.character(east_af$NAM_name))
List_indiv_unknown <-c(as.character(unassigned$Family), as.character(unassigned$NAM_name))

# Get all families in each subpopulation
admix_indiv<-NULL
for (p in 1:length(List_indiv_admix)){
  admix_indiv <-c(admix_indiv, (Melted_chrom$ID)[grep(List_indiv_admix[p], (Melted_chrom$variable))])
}

coastalm_indiv<-NULL
for (p in 1:length(List_indiv_coastal_m)){
  coastalm_indiv <-c(coastalm_indiv ,(Melted_chrom$ID)[grep(List_indiv_coastal_m[p], (Melted_chrom$variable))])
}

asian_indiv<-NULL
for (p in 1:length(List_indiv_Asian)){
  asian_indiv <-c(asian_indiv ,(Melted_chrom$ID)[grep(List_indiv_Asian[p], (Melted_chrom$variable))])
}

centrale_indiv<-NULL
for (p in 1:length(List_indiv_central_e)){
  centrale_indiv <-c(centrale_indiv ,(Melted_chrom$ID)[grep(List_indiv_central_e[p], (Melted_chrom$variable))])
}

eastaf_indiv<-NULL
for (p in 1:length(List_indiv_east_af)){
  eastaf_indiv <-c(eastaf_indiv ,(Melted_chrom$ID)[grep(List_indiv_east_af[p], (Melted_chrom$variable))])
}
unknown_indiv<-NULL
for (p in 1:length(List_indiv_unknown)){
  unknown_indiv <-c(unknown_indiv ,(Melted_chrom$ID)[grep(List_indiv_unknown[p], (Melted_chrom$variable))])
}

# Get binned allele proportions for families in each subpopulation
admixed_fams <- Melted_chrom[admix_indiv,c(1:5)]
asian_fams <- Melted_chrom[asian_indiv,c(1:5)]
coastalmed_fams <- Melted_chrom[coastalm_indiv,c(1:5)]
euro_fams <- Melted_chrom[centrale_indiv,c(1:5)]
eastaf_fams <- Melted_chrom[eastaf_indiv,c(1:5)]
un_fams <- Melted_chrom[unknown_indiv,c(1:5)]

# Add Subpop column to each
admixed_fams$Pop_location <- "Admixed"
asian_fams$Pop_location <- "Asian"
coastalmed_fams$Pop_location <- "Coastal Med."
euro_fams$Pop_location <- "C. European"
eastaf_fams$Pop_location <- "E. African"
un_fams$Pop_location <- "Un"

# Put them all together
Fam_seg_dist <- rbind(admixed_fams, asian_fams, coastalmed_fams, euro_fams, eastaf_fams, un_fams)

# Reorder factors in new dataframe
#Fam_seg_dist$Pop_location_sort = factor(Fam_seg_dist$Pop_location, levels = c('Central European', 'Coastal Mediterranean', 'East African', 'Asian', 'Admixed', 'Un'))
Fam_seg_dist$Pop_location_sort = factor(Fam_seg_dist$Pop_location, levels = c('C. European', 'Coastal Med.', 'E. African', 'Asian', 'Admixed', 'Un'))

# Reorder parent data
#family_info_sorted <- arrange(family_info, man_sort)

# Reorder families within factor
# Family_info should be sorted in the order you want to families in this case by man_sort
Arranged_seg_dist <- Fam_seg_dist %>%
  mutate(variable = factor(variable, levels = family_info$NAM_name)) %>%
  arrange(variable)

# Plot heat map of segregation distortion
ggplot(data = Arranged_seg_dist, aes(x = bin, y=factor(variable), fill=as.numeric(value))) +
  labs(x = "Chromosome", y = "BRIDG Family") +
  facet_grid( Pop_location_sort ~ chrom, scales = "free", space = "free", switch = "both") +
  #scale_fill_gradient2(high = "purple", low = "green", mid = "white", na.value = NA, midpoint = 0.50, guide = "colorbar", "Freuency of\nRasmusson Allele", limits = c(0.00, 1.00)) +
  scale_fill_distiller(palette = "PRGn", limits = c(0.00, 1.00), "Proportion\nRasmusson\nAllele", guide = "colorbar", labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), na.value = "white") +
  theme(axis.text.y = element_blank(), strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_rect(color = "grey"), 
        panel.grid = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "right") +
  geom_tile() 
# Export as device size landscape 7x12 Manuscripts/Figures/Seg_dist_bypopthenpheno

###########
#visualize one chromosome for one population at the time
#CHR5_OUTLIERS<-which(Arranged_seg_dist$chrom == "5H" & Arranged_seg_dist$Pop_location == "Coastal Med." & Arranged_seg_dist$value >0.75 )
CHR5_OUTLIERS<-which(Arranged_seg_dist$chrom == "3H" & Arranged_seg_dist$Pop_location == "E. African" )
SEG_SEGMENTS<-Arranged_seg_dist[CHR5_OUTLIERS,]
table(SEG_SEGMENTS[,1])

# get the unmelted data and visualize it in excell
Seg_dist_chr5<-Seg_dist[grep("CHR3",row.names(Seg_dist)),]
write.table(Seg_dist_chr5,"~/Desktop/Seg_dist_chr3_eafric.xls",quote=F,row.names=T,col.names=T,sep="\t")
# get only segments 
ggplot(data = SEG_SEGMENTS, aes(x = bin, y=factor(variable), fill=as.numeric(value))) +
  labs(x = "Chromosome", y = "BRIDG Family") +
  facet_grid( Pop_location_sort ~ chrom, scales = "free", space = "free", switch = "both") +
  #scale_fill_gradient2(high = "purple", low = "green", mid = "white", na.value = NA, midpoint = 0.50, guide = "colorbar", "Freuency of\nRasmusson Allele", limits = c(0.00, 1.00)) +
  scale_fill_distiller(palette = "PRGn", limits = c(0.00, 1.00), "Proportion\nRasmusson\nAllele", guide = "colorbar", labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), na.value = "white") +
  theme(axis.text.y = element_blank(), strip.background = element_blank(), axis.text.x = element_text(face="bold", color="#993333", size=6, angle=90), axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_rect(color = "grey"), 
        panel.grid = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "right") +
  geom_tile() 

#############
############
#PRGn is good
#RdGy
#RdYlBu shows that many sites are in equilibrium

# Estimate segregation distortion using R qtl
# This file contains 
genos_forSD <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/genos_qtl_SD.csv", header = T)

library(qtl)

BRIDG <- read.cross(format="csv", dir="~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/", file = "genos_qtl_SD.csv", na.strings="NA",
        genotypes=c("0","1","2"), alleles=c("0", "2"), estimate.map=F, error.prob=0.0001)

SD <- geno.table(BRIDG, chr = 1)
SD[SD$P.value < 0.01,]
SD_out <- geno.table(BRIDG, scanone.output=TRUE)
plot(SD_out)
plot(SD_out, lod=2)

data(listeria)
geno.table(listeria)
geno.table(listeria, chr=13)
gt <- geno.table(listeria)
gt[gt$P.value < 0.01,]
out <- geno.table(listeria, scanone.output=TRUE)
plot(out)
plot(out, lod=2)



# NSGC and NAM donor PCA
NSGC_6 <- read.csv("~/Documents/PhD/NAM/Parents/NSGC6_NAM_subpop_filled_nodups.csv", header=T, sep = ",", na.strings="NA")
# PArents with PC1 and PC2, subpopulation, and haplotype
NSGC_haplo <- read.csv("~/Documents/PhD/NAM/Parents/NSGC_parents_coords.csv", header=T, sep = ",", na.strings="NA")
NSGC_haplo <- as.character(NSGC_haplo$HaplotypeClasses)
# Merge so we have all NSGC, too
NSGC_parents <- merge(NSGC_haplo, NSGC_6, by.y = "Synonyms", by.x = "synonym", all.x = T, all.y = T)

# Reorder subpopulation factor for legend
NSGC_6$Subpop_sorted = factor(NSGC_6$Pop_location, levels = c("Rasmusson", "Central European", "Coastal Mediterranean", "East African", "Asian", "Admixed", "Unassigned", "NSGC Core"))

# PC1 and PC2 by region
ggplot(NSGC_6, aes(x = -(NSGC_6$V1), y = (NSGC_6$V2))) + 
  labs(x = "PC1 (15%)", y = "PC2 (8%)") +
  scale_color_manual(name = "BRIDG Parent\nSubpopulations", values = c("black", "olivedrab3", "turquoise1", "red", "violet", "goldenrod1", "gray35", "gray75")) +
  scale_size_manual(name = "BRIDG Parent\nSubpopulations", values = c(3, 1, 1, 1, 1, 1, 1, 1, .5)) +
  scale_shape_manual(name = "BRIDG Parent\nSubpopulations", values = (c(1, 1, 1, 1, 1, 1, 1, 1))) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = "grey"), axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10), legend.key = element_blank(), legend.background = element_rect(color = "grey")) +
  geom_point(aes(color = factor(Subpop_sorted), size = factor(Subpop_sorted), shape = factor(Subpop_sorted)), alpha = 0.75, stroke = 1.5, position = position_jitter(width = 0.008, height = 0.008), na.rm = T) +
  scale_x_continuous(limits = c(-0.052, 0.066), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.073, 0.066), expand = c(0,0)) 

# no color
ggplot(NSGC_6, aes(x = -(NSGC_6$V1), y = (NSGC_6$V2))) + 
  labs(x = "PC1 (15%)", y = "PC2 (8%)") +
  scale_color_manual(name = "BRIDG Parent\nSubpopulations", values = "gray75") +
  scale_size_manual(name = "BRIDG Parent\nSubpopulations", values = 1) +
  scale_shape_manual(name = "BRIDG Parent\nSubpopulations", values = 1) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = "grey"), axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10), legend.key = element_blank(), legend.background = element_rect(color = "grey")) +
  geom_point(na.rm = T, color = "gray75", size = 0.5, shape = 1, stroke = 1.5, alpha = 0.75) +
  scale_x_continuous(limits = c(-0.052, 0.066), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.073, 0.066), expand = c(0,0)) 
# Export in deivce size portrait 5x3 Population_description/PCA/NSGC_donor_PCA
# Principle component analysis of the NAM founders and 1000 6-row National Small Grains Collection core accessions based on 1027 BOPA SNP markers.
# Eigenvalue 1 explains approximately 15% of total variation and eigenvalue 2 explains approximately 8% of total variation. 

# Colored by Ppd-H1 haplotype
NSGC_parents$Subpop_sorted = factor(NSGC_parents$Pop_location, levels = c("Rasmusson", "Central European", "Coastal Mediterranean", "East African", "Asian", "Admixed", "Unassigned", "NSGC Core", "NA"))

# Get only columns of necessary data
NSGC_6_few <- select(NSGC_6, c(Synonyms, V1, V2, Pop_location))
NSGC_haplo_few <- select(NSGC_haplo, c(synonym, HaplotypeClasses))

NSGC_parents <- merge(NSGC_6_few, NSGC_haplo_few, by.x = "Synonyms", by.y = "synonym", all.x = T)
write.csv(NSGC_parents, "~/Documents/PhD/NAM/Parents/NSGC_parents_coords.csv") 
NSGC_parents <- read.csv("~/Documents/PhD/NAM/Parents/NSGC_parents_coords_filledhaps.csv", header= T)

ggplot(NSGC_parents, aes(x = -(NSGC_parents$V1), y = (NSGC_parents$V2))) + 
  labs(x = "PC1 (15%)", y = "PC2 (8%)") +
  scale_color_manual(name = "BRIDG Parent\nSubpopulations", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red",  "gray85", "black", "grey75")) +
  #scale_size_manual(name = "BRIDG Parent\nSubpopulations", values = c(3, 1, 1, 1, 1, 1, 1, 1, .5)) +
  scale_shape_manual(name = "Ppd-H1\nHaplotype",  na.value = "o", values = (c("1", "2", "3", "4", "5", "6", "7", "8", "9", "NA"))) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), 
        panel.grid = element_blank(), panel.background = element_rect(color = "grey"), 
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.key = element_blank(), legend.background = element_rect(color = "grey")) +
  geom_point(aes(color = factor(Pop_location), shape = factor(HaplotypeClasses)), 
        size = 4, alpha = 0.9, stroke = 1.5, position = position_jitter(width = 0.01, height = 0.01), na.rm = F) 
  scale_x_continuous(limits = c(-0.052, 0.066), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.073, 0.066), expand = c(0,0)) 
# NAM donor and progeny PCA

rm(list=ls())

# Import family information, with RILs and parents to know which color to paint them
family_info<-read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks.csv", header=T)
head(family_info)

# Separate family in populations
Admix<-subset(family_info, family_info$Pop_location == "Admixed")
coastal_m<-subset(family_info, family_info$Pop_location == "Coastal Mediterranean")
Asian<-subset(family_info, family_info$Pop_location == "Asian")
central_e<-subset(family_info, family_info$Pop_location == "Central European")
east_af<-subset(family_info, family_info$Pop_location == "East African")
unassigned<-subset(family_info, family_info$Pop_location == "Un")

#get identifier for individuals in a family
List_indiv_admix<-c(as.character(Admix$Family), as.character(Admix$NAM_name))
List_indiv_coastal_m <-c(as.character(coastal_m$Family), as.character(coastal_m$NAM_name))
List_indiv_Asian <-c(as.character(Asian$Family), as.character(Asian$NAM_name))
List_indiv_central_e <-c(as.character(central_e$Family), as.character(central_e$NAM_name))
List_indiv_east_af <-c(as.character(east_af$Family), as.character(east_af$NAM_name))
List_indiv_unknown <-c(as.character(unassigned$Family), as.character(unassigned$NAM_name))

#Calculate percentage of variance explained
ANALYSIS<-c("NAM_RILs")
MAIN<-c("NAM RIL and Parents")
DIR<-"20missing"

  EVE<-read.table(paste("~/Documents/PhD/NAM/Population_description/PCA/fromNAM80_afterQC_RIL_parents/20missing/output/NAM_RILs.eval",sep=""))
  pc1<-round((EVE[1,]/sum(EVE)*100),2) #for PC1
  pc2<-round((EVE[2,]/sum(EVE)*100),2)#for PC2
 
  
  TOTAL<-sum(EVE)
  EXPLAIN<-(EVE[,1]/TOTAL)*100
  plot(EXPLAIN)
  DATA<-read.table(paste("~/Documents/PhD/NAM/Population_description/PCA/fromNAM80_afterQC_RIL_parents/20missing/output/NAM_RILs.pca.evec",sep=""),header=F,row.names=1)
  
  dim(DATA)
  admix_indiv<-NULL
  for (p in 1:length(List_indiv_admix)){
    admix_indiv <-c(admix_indiv ,row.names(DATA)[grep(List_indiv_admix[p], row.names(DATA))])
  }
  
  coastalm_indiv<-NULL
  for (p in 1:length(List_indiv_coastal_m)){
    coastalm_indiv <-c(coastalm_indiv ,row.names(DATA)[grep(List_indiv_coastal_m[p], row.names(DATA))])
  }
  
  asian_indiv<-NULL
  for (p in 1:length(List_indiv_Asian)){
    asian_indiv <-c(asian_indiv ,row.names(DATA)[grep(List_indiv_Asian[p], row.names(DATA))])
  }
  
  centrale_indiv<-NULL
  for (p in 1:length(List_indiv_central_e)){
    centrale_indiv <-c(centrale_indiv ,row.names(DATA)[grep(List_indiv_central_e[p], row.names(DATA))])
  }
  
  eastaf_indiv<-NULL
  for (p in 1:length(List_indiv_east_af)){
    eastaf_indiv <-c(eastaf_indiv ,row.names(DATA)[grep(List_indiv_east_af[p], row.names(DATA))])
  }
  unknown_indiv<-NULL
  for (p in 1:length(List_indiv_unknown)){
    unknown_indiv <-c(unknown_indiv ,row.names(DATA)[grep(List_indiv_unknown[p], row.names(DATA))])
  }
  
# Get NAM parents in subpopulations
admixed_parents <- as.data.frame(cbind(DATA[admix_indiv[grep("PI|CIho", admix_indiv)],1],DATA[admix_indiv[grep("PI|CIho", admix_indiv)],2]))
asian_parents <- as.data.frame(cbind(DATA[asian_indiv[grep("PI|CIho", asian_indiv)],1],DATA[asian_indiv[grep("PI|CIho", asian_indiv)],2]))
coastalmed_parents <- as.data.frame(cbind(DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],1],DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],2]))
eastaf_parents <- as.data.frame(cbind(DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],1],DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],2]))
un_parents <- as.data.frame(cbind(DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],1],DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],2]))
euro_parents <- as.data.frame(cbind(DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],1],DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],2]))
Ras <- as.data.frame(cbind(DATA[grep("Ras",row.names(DATA)),1],DATA[grep("Ras",row.names(DATA)),2]))

# Get PC1 and PC2 for RILs in each subpopulation
admixed_RILs <- DATA[admix_indiv,c(1,2)]
asian_RILs <- DATA[asian_indiv,c(1,2)]
coastalmed_RILs <- DATA[coastalm_indiv,c(1,2)]
euro_RILs <- DATA[centrale_indiv,c(1,2)]
eastaf_RILs <- DATA[eastaf_indiv,c(1,2)]
un_RILs <- DATA[unknown_indiv,c(1,2)]

# Plot PC1 and PC2 for NAM parents and RILs colored by subpopulation
ggplot() +
  labs(x = "PC1 (4.16%)", y = "PC2 (2.76%)") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = "grey"), axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
  geom_point(data = un_RILs, aes(x = as.numeric(V2), y = as.numeric(-V3)), color = "gray35", size = .5, shape = 1, alpha = 0.5, position = position_jitter(width = 0.005, height = 0.005)) +
  geom_point(data = admixed_RILs, aes(x = as.numeric(V2), y = as.numeric(-V3)), color = "goldenrod1", size = .5, shape = 1, alpha = 0.5, position = position_jitter(width = 0.005, height = 0.005)) +
  geom_point(data = asian_RILs, aes(x = as.numeric(V2), y = as.numeric(-V3)), color = "violet", size = .5, shape = 1, alpha = 0.5, position = position_jitter(width = 0.005, height = 0.005)) +
  geom_point(data = coastalmed_RILs, aes(x = as.numeric(V2), y = as.numeric(-V3)), color = "turquoise1", size = .5, alpha = 0.5, shape = 1, position = position_jitter(width = 0.005, height = 0.005)) +
  geom_point(data = euro_RILs, aes(x = as.numeric(V2), y = as.numeric(-V3)), color = "olivedrab3", size = .5, shape = 1, alpha = 0.5, position = position_jitter(width = 0.005, height = 0.005)) +
  geom_point(data = eastaf_RILs, aes(x = as.numeric(V2), y = as.numeric(-V3)), color = "red", size = .5, shape = 1, alpha = 0.5, position = position_jitter(width = 0.005, height = 0.005)) +
  geom_point(data = admixed_parents, aes(x = as.numeric(V1), y = as.numeric(-V2)), color = "goldenrod1", size = 1, shape = 1, stroke = 1.5, position = position_jitter(width = 0.004, height = 0.004)) +
  geom_point(data = asian_parents, aes(x = as.numeric(V1), y = as.numeric(-V2)), color = "violet", size = 1, shape = 1, stroke = 1.5, position = position_jitter(width = 0.004, height = 0.004)) +
  geom_point(data = coastalmed_parents, aes(x = as.numeric(V1), y = as.numeric(-V2)), color = "turquoise1", size = 1, stroke = 1.5, shape = 1, position = position_jitter(width = 0.0004, height = 0.004)) +
  geom_point(data = euro_parents, aes(x = as.numeric(V1), y = as.numeric(-V2)), color = "olivedrab3", size = 1, stroke = 1.5, shape = 1, position = position_jitter(width = 0.004, height = 0.004)) +
  geom_point(data = eastaf_parents, aes(x = as.numeric(V1), y = as.numeric(-V2)), color = "red", size = 1, stroke = 1.5, shape = 1, position = position_jitter(width = 0.004, height = 0.004)) +
  geom_point(data = un_parents, aes(x = as.numeric(V1), y = as.numeric(-V2)), color = "gray35", size = 1, stroke = 1.5, shape = 1, position = position_jitter(width = 0.004, height = 0.004)) +
  geom_point(data = Ras, aes(x = as.numeric(V1), y = as.numeric(-V2)), color = "black", size = 3, shape = 1, stroke = 1.5) +
  scale_x_continuous(limits = c(-0.021, 0.085), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.085, 0.02), expand = c(0,0)) 
# save as device size portrait 3x3 in Population_description/PCA/Donor_progeny_PCA
dev.off()

# Calculate statistics on segregation distortion
install.packages("qtl")
library(qtl)
library(NAM)

# Get genotypes
genos <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/genos_80miss_byfam.csv", header = T)
phenos <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/phenos_80miss_byfam.csv", header = T)

# Get the number of SNPs per chromosome
chr = data.frame(table(gsub("._.+$", "",colnames(genos))))[,2]
fam = phenos$family

my_FST = Fst ( gen = genos, fam = fam )
plot(my_FST)
my_FST_data <- as.data.frame.array(my_FST)


Melted_chrom <- Melted_dist %>%
  mutate(chrom = replace(chrom, chrom=="CHR1", "1H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR2", "2H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR3", "3H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR4", "4H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR5", "5H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR6", "6H")) %>%
  mutate(chrom = replace(chrom, chrom=="CHR7", "7H")) %>%
  as.data.frame()
