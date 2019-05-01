#Step3_PlotSegDist.R

# Author: Alex Ollhoff
# Description: Plot Rasmusson contribution to homozygous genotypes
################################################################################################################
rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape2)
Seg_dist <- read.table("~/Documents/Manuscript/88Fam_genmap_segDist_binned_AA.txt", header = T)

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
family_info <- read.csv("~/Documents/Manuscript/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA")

# Separate family in populations
Admix<-subset(family_info, family_info$Pop_location == "Admixed")
coastal_m<-subset(family_info, family_info$Pop_location == "Coastal Mediterranean")
Asian<-subset(family_info, family_info$Pop_location == "Asian")
central_e<-subset(family_info, family_info$Pop_location == "Central European")
east_af<-subset(family_info, family_info$Pop_location == "East African")
unassigned<-subset(family_info, family_info$Pop_location == "Un")

# Get identifier for individuals in a family
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

# Reorder families within factor
# Family_info should be sorted in the order you want to families, in this case by man_sort
Arranged_seg_dist <- Fam_seg_dist %>%
  mutate(variable = factor(variable, levels = family_info$NAM_name)) %>%
  arrange(variable)

# Plot heat map of segregation distortion
pdf("~/Documents/Manuscript/Figure2_ParentalContributions.pdf",width=10,height=7)
ggplot(data = Arranged_seg_dist, aes(x = bin, y=factor(variable), fill=as.numeric(value))) +
  labs(x = "Chromosome", y = "BRIDG6 Family") +
  facet_grid( Pop_location_sort ~ chrom, scales = "free", space = "free", switch = "both") +
  #scale_fill_gradient2(high = "purple", low = "green", mid = "white", na.value = NA, midpoint = 0.50, guide = "colorbar", "Freuency of\nRasmusson Allele", limits = c(0.00, 1.00)) +
  scale_fill_distiller(palette = "PRGn", limits = c(0.00, 1.00), "Proportion\nRasmusson\nAllele", guide = "colorbar", labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), na.value = "white") +
  theme(axis.text.y = element_blank(), strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_rect(color = "grey"), 
        panel.grid = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "right") +
  geom_tile() 
dev.off()

###########
# Visualize one chromosome for one population at a time
CHR5_OUTLIERS<-which(Arranged_seg_dist$chrom == "5H" & Arranged_seg_dist$Pop_location == "Coastal Med."  )
#CHR5_OUTLIERS<-which(Arranged_seg_dist$chrom == "3H" & Arranged_seg_dist$Pop_location == "E. African" )
#CHR5_OUTLIERS<-which(Arranged_seg_dist$chrom == "1H" )
SEG_SEGMENTS<-Arranged_seg_dist[CHR5_OUTLIERS,]


# Get the unmelted data and visualize it in excel
Seg_dist_chr5<-Seg_dist[grep("CHR5",row.names(Seg_dist)),]
write.table(Seg_dist_chr5,"~/Desktop/Seg_dist_chr3_eafric.xls",quote=F,row.names=T,col.names=T,sep="\t")

HR645S<-SEG_SEGMENTS[which(SEG_SEGMENTS$variable == "HR645S"),]
# Get only segments 
ggplot(data = SEG_SEGMENTS, aes(x = bin, y=factor(variable), fill=as.numeric(value))) +
  labs(x = "Chromosome", y = "BRIDG Family") +
  facet_grid( Pop_location_sort ~ chrom, scales = "free", space = "free", switch = "both") +
  #scale_fill_gradient2(high = "purple", low = "green", mid = "white", na.value = NA, midpoint = 0.50, guide = "colorbar", "Freuency of\nRasmusson Allele", limits = c(0.00, 1.00)) +
  scale_fill_distiller(palette = "PRGn", limits = c(0.00, 1.00), "Proportion\nRasmusson\nAllele", guide = "colorbar", labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), na.value = "white") +
  theme(axis.text.y = element_text(face="bold", color="#993333", size=6), strip.background = element_blank(), axis.text.x = element_text(face="bold", color="#993333", size=6, angle=90), axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_rect(color = "grey"), 
        panel.grid = element_blank(), text = element_text(size = 16), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "right") +
  geom_tile() 
