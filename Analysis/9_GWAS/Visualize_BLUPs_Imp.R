# Visualize mapping results for manuscript
# 80miss_byfam 
# March 4, 2017
# Alex Ollhoff
#######################################################################################################################################
rm(list=ls())
library(ggplot2)
library(dplyr)

# RIL phenotypic data
BRIDG_heading <- read.csv("~/Documents/PhD/NAM/BLUPs/fam_fixed/fam_BLUEs_nob_no19_famcorrect_forboxplot.csv", header=T, sep = ",", na.strings="NA")
# Parent phenotypes
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/BRIDG_parents_computing.csv", header = T, na.strings = "NA")
# Significnace values
SNP_info_complete <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/80miss_byfam_SNP_info.csv", header = T)
# Get GWAS results 
Allele_effects_binned <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/results_for_doesnaughtcompute.csv", header = T, check.names = F)


#################*** Format results for visualization ***##################

# Load marker data 
RESULTS <- read.csv("~/Dropbox/GITHUB/BRIDG6/Datasets/GWAS_Results/my_gwas_PolyTest_out_80miss_BLUE.csv", header = T,row.names=1)

# Load marker names
NAMES <- read.csv("~/Dropbox/GITHUB/BRIDG6/Datasets/GWAS_Results/my_gwas_SNPs_out_80miss_BLUE.csv", header = T)

# Calculate the 5% false discovery rate threshold
THR = -log10(0.05 / ( nrow(RESULTS) * (1-0.05)))

# See how many SNPs are above the significance threshold
Results_sig <- filter(RESULTS, lod > THR) # 100 significant markers
w <- length(Results_sig)

# Add SNP names to mapping results
row.names(RESULTS)<-as.character(NAMES[,2])

# Get SNP names, separate the position and chromosome number
SNP_original_names <- as.data.frame(gsub("X", "", row.names(RESULTS)))
Physical_positions<- apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][2])
Physical_positions<-as.data.frame(Physical_positions)
Physical_positions<-cbind(SNP_original_names[,1],Physical_positions)
ADD_CHR<-c(312837513,393532674,394310633,355061206,380865482,294822070,325797516)

for (i in 1:7){
  CHR1<-Physical_positions [grep(paste(i,"H1", sep=""),Physical_positions[,1]),]
  CHR2<-Physical_positions [grep(paste(i,"H2", sep=""),Physical_positions[,1]),]
  
  CHRboth<-rbind(CHR1,CHR2)
  
  CHR2_corrected<-as.numeric(as.character(CHR2[,2])) + ADD_CHR[i]
  CHRall_pos<-c(as.numeric(as.character(CHR1[,2])), CHR2_corrected)
  
  Positions_corrected<-cbind(CHRboth, CHRall_pos)
  assign(paste("CHR_cor_",i,sep=""), Positions_corrected)
  
}

NEW_positions<-rbind(CHR_cor_1, CHR_cor_2, CHR_cor_3, CHR_cor_4, CHR_cor_5, CHR_cor_6, CHR_cor_7)


# Add column with chromosome names
CHROM<-apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][1])
Chrom_chr<-sub("H1|H2","H",CHROM)

# Add new column to data frame 
Results_complete<- cbind(NEW_positions[,3],Chrom_chr,RESULTS)
names(Results_complete)[1]<-"Cumulative_bp"


# UPDATE to corresponding mapping directory 
write.table(Results_complete, "~/Desktop/80miss_byfam_SNP_info.csv",quote=F,row.names=T,col.names=T,sep="\t")
SNP_info_complete <- read.table("~/Desktop/80miss_byfam_SNP_info.csv", header = T,row.names=1)

###<<<*** Manhattan plot ***>>>###
THR =  -log10(0.05 / ( nrow(RESULTS) * (1-0.05))) 

#THR =  -log10(0.05 / ( ncol(SNP_info_complete) * (1-0.05)))
w = which((SNP_info_complete$pval) > (-log10(THR)) )

# Simple manhattan, manuscript version
ggplot(SNP_info_complete, aes(Cumulative_bp, pval)) +
  labs(x = "Chromosome (base pairs)", y = "-log(p-value)") +
  theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), 
        panel.background = element_rect(color = "grey"), axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
  geom_point(size = 1, alpha = 1/2) +
  geom_hline(yintercept = THR) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 50), position = "right") +
  scale_x_continuous(expand = c(0,0)) +
  facet_grid( ~ Chrom_chr, scales = "free_x", switch = "x") 
# Export as 3x12 landscape deivce size BLUPs_GxE/Simple_man_hd
# This warning message indicates that there are points above the limit on the y axis. Change using scale_y_continuous limits = c(0, my_new_limit)
    #Warning message:
    #Removed 15 rows containing missing values (geom_point). 


###*** Allele effects heat map ***###
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)

# Heatmap of allele effects 
Fam_QTL_vis <- read.table("~/Desktop/80miss_byfam_SNP_info.csv", header = T,row.names=1)
parent_heading <- read.csv("~/Dropbox/GITHUB/BRIDG6/Datasets/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA")
Fam_QTL_vis[1:10,1:10]

# add column with SNP names
SNP<-sub("X","",row.names(Fam_QTL_vis))

Fam_QTL_vis<-cbind(SNP,Fam_QTL_vis)
# Get useful columns
Fam_QTL_vis_few <- select(Fam_QTL_vis, c(SNP, Chrom_chr, Cumulative_bp, lod, (dim(Fam_QTL_vis)[2]-87) : dim(Fam_QTL_vis)[2]))

CHROM_LENGTH<-c(558535432, 768075024, 699711114, 647060158, 289164678, 583380513, 657224000)

#work one chromsoome at a time
for (n in 1:7){
  print(n)
  # Create an empty vector where empty bins will go
  BIN_empy<-NULL
  CHR<-Fam_QTL_vis_few[grep(paste(n,"H",sep=""), Fam_QTL_vis_few$Chrom_chr),]
  BIN<-as.data.frame(rep(NA,dim(CHR)[1]))
  colnames(BIN)<-"BIN"
  CHR<-cbind(CHR,BIN)
  
  #add a row of bin number 
  BIN_START<- 1
  Bp_start<-1
  
  # while the next bin is less than 1M to the end of the chromosome
  WindowSize<-20000000
  while(Bp_start < CHROM_LENGTH[n] ){
    START<-Bp_start
    END<-Bp_start + (WindowSize -1)
    
    SNPsInBin<-which(CHR$Cumulative_bp >= START & CHR$Cumulative_bp <= END)
    # if there are not SNPs in a bin then ad the bin number to the emtpy vector, otherwise add it to
    # the new column in the DATA matrix
    if (length(SNPsInBin) >0){
      CHR[SNPsInBin,dim(CHR)[2]]<-paste("",BIN_START,sep="")
    }else{BIN_empy<-c(BIN_empy,paste("",BIN_START,sep=""))}
    
    Bp_start<-END+1
    BIN_START<-BIN_START + 1
  }
  
  if (length(BIN_empy)>0){
    BIN_extra<-as.data.frame(BIN_empy)
    #matrix to add Bins without information
    BIM_matrix<-matrix(0,ncol=(dim(CHR)[2] -1), nrow=length(BIN_empy))
    BIM_full_matrix<-cbind(BIM_matrix,BIN_extra)
    
    colnames(BIM_full_matrix)<-colnames(CHR)
    
    # Combine bins in SNPs with those that didn't have any SNP
    CHR_complete<-rbind(CHR,BIM_full_matrix)
    
  }else{CHR_complete<-CHR}
  
  assign(paste("CHR_",n,"",sep=""),CHR_complete)
  
}

DATA_w_bins<-rbind(CHR_1,CHR_2,CHR_3,CHR_4,CHR_5,CHR_6,CHR_7)

# Bin marker effects and get the minimum, maximum, and absolute max for each bin
Binned_means_max <- DATA_w_bins %>% 
  group_by(Chrom_chr, BIN) %>%
  summarise_each(funs(max), allele.eff.founder.1:allele.eff.founder.88)

Binned_means_min <- DATA_w_bins %>% 
  group_by(Chrom_chr, BIN) %>%
  summarise_each(funs(min), allele.eff.founder.1:allele.eff.founder.88)

Binned_means_abs <- DATA_w_bins %>% 
  group_by(Chrom_chr, BIN) %>%
  summarise_each(funs(.[which.max(abs(.))]), allele.eff.founder.1:allele.eff.founder.88)

# Make these values numeric so they plot nice
Binned_means_abs$BIN <- as.numeric(Binned_means_abs$BIN)

# Sort data by chromsome and then bin
#Binned_means_sorted <- Binned_means_abs %>%
  #arrange(Chrom_chr, BIN)

# Melt into long format for heatmap
melted_QTL_try <- melt(Binned_means_abs, id.vars = c("Chrom_chr", "BIN"))

# Order by family then chrom then bins on each chromsome
melted_QTL_try_ordered <- melted_QTL_try %>%
  arrange(variable, Chrom_chr, BIN)

nbins <- nrow(melted_QTL_try_ordered)/88

# Get family number
Fam_num <- apply(melted_QTL_try_ordered, 1, function(x) strsplit(as.character(x), ".")[[1]][1])
Fam_num<-sub("allele.eff.founder.","",  melted_QTL_try_ordered$variable)
head(Fam_num)
melted_QTL_withfam <- cbind(melted_QTL_try_ordered, Fam_num)

# Reorder families by subpopulation and then from earliest parent hd to latest))
parent_heading_sorted <- arrange(parent_heading, man_sort_reverse)

#melted_QTL_withfam$fam_man <- factor(melted_QTL_withfam$Fam_num, levels = c(parent_heading$man_sort_reverse))
arranged_effects <- melted_QTL_withfam %>%
  mutate(Fam_num = factor(Fam_num, levels = parent_heading_sorted$new_fam)) %>%
  arrange(Fam_num)

# Add subpopulation
# The number of rows in the rep = (the number of families in the subpop)*(the number of bins)
subpopulation_levels <- c(rep("Central European", (17*nbins)), rep("Coastal Mediterranean", (19*nbins)), rep("East African", (8*nbins)), rep("Asian", (19*nbins)), rep("Admixed", (23*nbins)), rep("Unassigned", (2*nbins)))
arranged_effects$subpop <- subpopulation_levels
head(arranged_effects)

# Reorder subpopulations
arranged_effects$subpop_sort = factor(arranged_effects$subpop, levels = c('Central European', 'Coastal Mediterranean', 'East African', 'Asian', 'Admixed', 'Un'))

# Remove chromsomes NA
arranged_effects_noNA <- filter(arranged_effects, Chrom_chr != "NA")

# plot
ggplot(data = arranged_effects_noNA, aes(x = BIN, y=factor(Fam_num), fill=value)) +
  labs(x = "Chromosome", y = "Family and\nSubpopulation") +
  facet_grid(subpop_sort ~ Chrom_chr, scales = "free", space = "free_y", switch = "both") +
  scale_fill_gradient2(high = "red3", low = "navy", mid = "white", na.value = "grey", midpoint = 0, guide = "colorbar", "Allele\nEffect\n(Days)", limits = c(-4, 3)) +
  geom_tile() +
  theme(axis.title.y = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(), 
        axis.text.x = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), text = element_text(size = 14),
        panel.background = element_rect(color = "white"), panel.grid = element_blank(), legend.text = element_text(size = 14)) +
  scale_x_continuous(expand = c(0,0)) 
# Export as device size landscape 7x12 "Heat_hd" in 80miss_byfam
##################################################################################################################################
############################### Revised upto this point ============#################################################################
##################################################################################################################################

# Add cM positions for markers that are in the genetic map Ana Poets made
SNP_info_complete_QTL <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/80miss_byfam_SNP_info_QTL_byhand_nogenes.csv", header = T)
cM_positions <- as.data.frame(read.table("~/Documents/PhD/NAM/NAM_mapping/Genotypes/Composite_genMap_complete.txt", header = TRUE))
colnames(cM_positions)[1] <- "SNP"
markers_cM_QTL <- merge(SNP_info_complete_QTL, cM_positions, by.x = "SNP", by.y = "SNP", all = T)
head(markers_cM_QTL)

write.csv(markers_cM_QTL, "~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/SNP_info_cM.csv")
QTL_cM <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/SNP_info_cM.csv", header = T)
QTL_cM <- as.data.frame(QTL_cM)

# Look at cM region covered by each QTL
QTL_cM <- QTL_cM %>%
  group_by(QTL) %>%
  summarise(Total_SNPs = n(), 
            Minimum_bp = min(as.numeric(as.character(Cumulative_bp))/1000000, na.rm = T),
            Maximum_bp = max(as.numeric(as.character(Cumulative_bp))/1000000, na.rm = T),
            QTL_mbp = (Maximum_bp - Minimum_bp),
            Minimum_cM = min(Cumulative_cM, na.rm = T), 
            Maximum_cM = max(Cumulative_cM, na.rm = T), 
            QTL_cM = Maximum_cM - Minimum_cM, na.rm = T)

# Import genotypes 
genos_for_vis <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/genos_GxE.csv", header = T)
genos_for_vis[1:10, 1:10]

# Pull out family from each individual line name
Family_num <- as.data.frame(gsub("HR", "", genos_for_vis$X))
Family_num_only <- as.data.frame(gsub("S[0-9][0-9][0-9]", "", Family_num[,1]))

# Add family to original genotype frame
genos_with_fam <- cbind((genos_for_vis), Family_num_only)

# Rename column
colnames(genos_with_fam)[7775] <- "family"

# Remove first column
genos_with_fam_nonames <- genos_with_fam[,-1]
genos_with_fam_nonames[1:10, 7770:7774]

# take only a few columns to test with
genos_sub <- select(genos_with_fam_nonames, X7H2_330413962, family)

# 

# Count the number of families segregating for each marker
#is.polymorphic <- function(x) length(unique(x)) == 0
#fam_per_marker <- genos_sub %>% group_by(family) %>% summarize_at(vars(starts_with("X")), is.polymorphic) %>% summarize_at(vars(starts_with("X")), sum)
#write.csv(fam_per_marker, "~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/Seg_fams_per_marker.csv")


# Count the number of markers segregating in each family
is.polymorphic <- function(x) length(unique(x)) == 0
seg_per_fam <- genos_sub %>% group_by(family) %>% summarize_at(vars(starts_with("X")), is.polymorphic) 
seg_per_fam_t <- t(seg_per_fam)
seg_per_fam_sum <- summarize(seg_per_fam_t)
write.csv(seg_per_fam, "~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/Seg_markers_per_fam.csv")


###*** Phenotypic distribution ***###

library(dplyr)
library(ggplot2)

# Box plot of family heading date distributions and parent mean
NAM_heading <- read.csv("~/Documents/PhD/NAM/BLUPs/fam_fixed/fam_BLUEs_nob_no19_famcorrect_forboxplot.csv", header=T, sep = ",", na.strings="NA", stringsAsFactors = F)
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA", stringsAsFactors = F)
head(parent_heading)

summary(NAM_heading)

# Add population mean to BLUEs
NAM_heading$hd <- (NAM_heading$BLUE + 56.44459)

# Look at distribution of BLUEs
ggplot() +
  labs(x = "Days to Heading", y = "Frequency") +
  theme(strip.background = element_blank(), axis.ticks = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), panel.background = element_rect(color = "grey"), legend.position = "bottom",
        axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  #facet_grid( subpop_sort ~ ., scales = "free_y", space = "free_y", switch = "y") +
  #geom_violin(data = NAM_heading_sorted, aes(x = factor(man_sort_reverse), y = hd), na.rm = T, draw_quantiles = 0.5) +
  #geom_point(data = parent_heading_sorted, aes(x = factor(man_sort_reverse), y =  Mean, color = parent_heading_sorted$Mean), shape = 21, fill = "black", stroke = 1, size = 1, na.rm = T) +
  #scale_color_distiller(palette = "RdBu", na.value = NA, guide = FALSE) +
  #scale_color_gradient2(high = "red3", low = "navy", mid = "white", na.value = NA, midpoint = 52.5, breaks = c(40.8, 50, 60, 71), limits = c(40, 72), "Days to\nHeading") +
  scale_x_continuous(limits = c(40, 73), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_color_gradient2(high = "red3", low = "navy", mid = "white", na.value = NA, midpoint = 52.5, breaks = c(40.8, 50, 60, 71), limits = c(40, 72), "Days to\nHeading") +
  geom_histogram(data = NAM_heading, aes(BLUE), bins =100) +
  geom_vline(xintercept = 52.5) 


# Separate family in populations
Admix<-subset(parent_heading, parent_heading$Pop_location == "Admixed")
coastal_m<-subset(parent_heading, parent_heading$Pop_location == "Coastal Mediterranean")
Asian<-subset(parent_heading, parent_heading$Pop_location == "Asian")
central_e<-subset(parent_heading, parent_heading$Pop_location == "Central European")
east_af<-subset(parent_heading, parent_heading$Pop_location == "East African")
unassigned<-subset(parent_heading, parent_heading$Pop_location == "Un")

#get identifier for individuals in a family
List_indiv_admix<-c(as.character(Admix$new_fam))
List_indiv_coastal_m <-c(as.character(coastal_m$new_fam))
List_indiv_Asian <-c(as.character(Asian$new_fam))
List_indiv_central_e <-c(as.character(central_e$new_fam))
List_indiv_east_af <-c(as.character(east_af$new_fam))
List_indiv_unknown <-c(as.character(unassigned$new_fam))

# Add ID column to BLUPs
NAM_heading$ID <- c(1:nrow(NAM_heading))

# Get all families in each subpopulation
admix_indiv<-NULL
for (p in 1:length(List_indiv_admix)){
  admix_indiv <-c(admix_indiv, (NAM_heading$ID)[grep(List_indiv_admix[p], (NAM_heading$family))])
}

coastalm_indiv<-NULL
for (p in 1:length(List_indiv_coastal_m)){
  coastalm_indiv <-c(coastalm_indiv ,(NAM_heading$ID)[grep(List_indiv_coastal_m[p], (NAM_heading$family))])
}

asian_indiv<-NULL
for (p in 1:length(List_indiv_Asian)){
  asian_indiv <-c(asian_indiv ,(NAM_heading$ID)[grep(List_indiv_Asian[p], (NAM_heading$family))])
}

centrale_indiv<-NULL
for (p in 1:length(List_indiv_central_e)){
  centrale_indiv <-c(centrale_indiv ,(NAM_heading$ID)[grep(List_indiv_central_e[p], (NAM_heading$family))])
}

eastaf_indiv<-NULL
for (p in 1:length(List_indiv_east_af)){
  eastaf_indiv <-c(eastaf_indiv ,(NAM_heading$ID)[grep(List_indiv_east_af[p], (NAM_heading$family))])
}

unknown_indiv<-NULL
for (p in 1:length(List_indiv_unknown)){
  unknown_indiv <-c(unknown_indiv ,(NAM_heading$ID)[grep(List_indiv_unknown[p], (NAM_heading$family))])
}

# Get families in each subpopulation
admixed_fams <- NAM_heading[admix_indiv,c(1:9)]
asian_fams <- NAM_heading[asian_indiv,c(1:9)]
coastalmed_fams <- NAM_heading[coastalm_indiv,c(1:9)]
euro_fams <- NAM_heading[centrale_indiv,c(1:9)]
eastaf_fams <- NAM_heading[eastaf_indiv,c(1:9)]
un_fams <- NAM_heading[unknown_indiv,c(1:9)]

# Add Subpop column to each
admixed_fams$Pop_location <- "Admixed"
asian_fams$Pop_location <- "Asian"
coastalmed_fams$Pop_location <- "Coastal Med."
euro_fams$Pop_location <- "C. European"
eastaf_fams$Pop_location <- "E. African"
un_fams$Pop_location <- "Un"

# Put them all together
NAM_heading_pops <- rbind(admixed_fams, asian_fams, coastalmed_fams, euro_fams, eastaf_fams, un_fams)

# Rename subpopulations that have too long of labels
NAM_heading_renamed <- NAM_heading_pops %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Central European", "C. European")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Coastal Mediterranean", "Coastal Med.")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="East African", "E. African")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Asian", "Asian")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Admixed", "Admixed")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Un", "Un")) %>%
  as.data.frame()

# Rename subpopulations that have too long of labels
parent_heading_renamed <- parent_heading %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Central European", "C. European")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Coastal Mediterranean", "Coastal Med.")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="East African", "E. African")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Asian", "Asian")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Admixed", "Admixed")) %>%
  mutate(Pop_location = replace(Pop_location, Pop_location=="Un", "Un")) %>%
  as.data.frame()

# Assign new order to subpopulation factor
parent_heading_renamed$subpop_sort = factor(parent_heading_renamed$Pop_location, levels = c('C. European', 'Coastal Med.', 'E. African', 'Asian', 'Admixed', 'Un'))
NAM_heading_renamed$subpop_sort = factor(NAM_heading_renamed$Pop_location, levels = c('C. European', 'Coastal Med.', 'E. African', 'Asian', 'Admixed', 'Un'))

# Reorder families by subpopulation and then from earliest parent hd to latest))
parent_heading_sorted <- arrange(parent_heading_renamed, man_sort_reverse)

# Add a new sort column to each family so they can be sorted by subpopulation and then family mean
NAM_heading_renamed$man_sort <- factor(NAM_heading_renamed$family, levels = c(parent_heading$man_sort_reverse))
NAM_heading_sorted <- NAM_heading_renamed %>%
  mutate(man_sort_reverse = factor(family, levels = parent_heading_sorted$new_fam)) %>%
  arrange(man_sort_reverse)

# Plot
ggplot() +
  labs(x = "Exotic Parent Subpopulations", y = "Days to Heading") +
  theme(strip.background = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), panel.background = element_rect(color = "grey"), legend.position = "bottom",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  facet_grid( subpop_sort ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_violin(data = NAM_heading_sorted, aes(x = factor(man_sort_reverse), y = hd), na.rm = T, draw_quantiles = 0.5) +
  geom_point(data = parent_heading_sorted, aes(x = factor(man_sort_reverse), y =  Mean, color = parent_heading_sorted$Mean), shape = 21, fill = "black", stroke = 1, size = 1, na.rm = T) +
  #scale_color_distiller(palette = "RdBu", na.value = NA, guide = FALSE) +
  scale_color_gradient2(high = "red3", low = "navy", mid = "white", na.value = NA, midpoint = 52.5, breaks = c(40.8, 50, 60, 71), limits = c(40, 72), "Days to\nHeading") +
  ylim(40, 72) +
  geom_hline(yintercept = 52.5) +
  coord_flip() 
# Save in device portrait 3x6 "NAM_pheno_dist" in Population_description

ggplot() +
  labs(x = "Exotic Parent Subpopulations", y = "Days to Heading") +
  theme(strip.background = element_blank(), axis.ticks = element_blank(), panel.grid.minor = element_blank(), 
       panel.background = element_rect(color = "grey"), legend.position = "bottom",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  #facet_grid( subpop_sort ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_boxplot(data = NAM_heading_sorted, aes(x = factor(family), y = hd), na.rm = T, draw_quantiles = 0.5) +
  #geom_point(data = parent_heading_sorted, aes(x = factor(man_sort_reverse), y =  Mean, color = parent_heading_sorted$Mean), shape = 21, fill = "black", stroke = 1, size = 1, na.rm = T) +
  #scale_color_distiller(palette = "RdBu", na.value = NA, guide = FALSE) +
  scale_color_gradient2(high = "red3", low = "navy", mid = "white", na.value = NA, midpoint = 52.5, breaks = c(40.8, 50, 60, 71), limits = c(40, 72), "Days to\nHeading") +
  ylim(40, 72) +
  geom_hline(yintercept = 52.5) +
  coord_flip() 

# Get only Parent mean HD
parent_hd <- select(parent_heading, c(man_sort, subpop_sort, Mean, fam_mean))
head(parent_hd)
parent_melted <- melt(parent_hd, id.vars = c("man_sort", "subpop_sort"))
head(parent_melted)
# Make heat map of parent mean and family mean
ggplot(data = parent_melted, aes(x = variable, y = -man_sort, fill = value)) +
  labs(x = "SNP Marker and\nCandidate Gene", y = "Family and\nSubpopulation") +
  facet_grid( subpop_sort ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_distiller(palette = "RdBu", "Allele\nEffect\n(Days)", na.value = NA, limits = c(40, 70)) +
  #scale_color_gradient2(high = "red3", low = "navy", mid = "white", na.value = "white", limits = c(40, 78)) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +
  geom_tile() 


# Get Parent data
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA", stringsAsFactors = F)
head(parent_heading)

# Get all effects for major QTL



# Get NSGC data with lat and long
NSGC <- read.csv("~/Documents/PhD/NAM/Parents/NSGC_info_from_maria.csv", header = T, na.strings = "NA", stringsAsFactors = F)
head(NSGC)

# Merge
NSGC_parents <- merge(parent_heading, NSGC, by.x = "synonym", by.y = "Accession")
head(NSGC_parents)

################ Map with colored points #####################
GPS <- parent_phenos <- read.csv("~/Documents/PhD/NAM/Parents/GPS_coords_parents_only.csv", header=T, sep = ",", na.strings="NA") 

mp <- NULL
mapWorld <- borders("world", colour="gray78", fill="gray50") # create a layer of borders
mp <- ggplot() +   mapWorld
mp <- mp+ geom_point(aes(x=NSGC_parents$Longitude, y=NSGC_parents$Latitude) ,color=factor(NSGC_parents$color), size=5) +
  labs(title = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), title=element_text(size = 30)) 

mp <- ggplot() +   mapWorld
mp <- mp+ geom_point(aes(x=NSGC_parents$LONGITUDE, y=NSGC_parents$LATITUDE) ,color=factor(NSGC_parents$color), size=2) +
  labs(title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), title=element_text(size = 30)) 

mp

mp <- ggplot() +   mapWorld
mp <- mp + geom_point(data = NSGC_parents, aes(x=LONGITUDE, y=LATITUDE, color=X2H1_27204805), size=2) +
  scale_colour_distiller(palette = "RdBu", "Allele\nEffect\n(Days)", na.value = NA) +
  labs(x = "Latitude", y = "Longitude") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), title=element_text(size = 10)) 

mp

mp <- ggplot() +   mapWorld
mp <- mp + 
  geom_point(data = NSGC_parents, aes(x=LONGITUDE, y=LATITUDE, color=X2H1_27204805), size=2) +
  #geom_point(data = NSGC_parents, aes(x=LONGITUDE+4, y=LATITUDE, color=X7H2_320469033), size=2) +
  #geom_point(data = NSGC_parents, aes(x=LONGITUDE+2, y=LATITUDE, color=X7H1_39192808), size=2) +
  scale_colour_distiller(palette = "RdBu", "Allele\nEffect\n(Days)", na.value = NA) +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), title=element_text(size = 10)) 

mp

# Correlation between family means and sum of allele effects for significant markers
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA", stringsAsFactors = F)

# Get family and heading date
NAM_heading <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/phenos_MAF05.csv", header = T)

# Get mean heading date for each family
fam_heading_means <- NAM_heading %>%
  group_by(family) %>%
  summarise(fam_mean = mean(BLUE))

# Get allele effects
SNP_info_complete_QTL <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/SNP_info_QTL.csv", header = T)

# Melt 
SNP_info_select <- select(SNP_info_complete_QTL, c(SNP, eff.1:eff.88))
colnames(SNP_info_select) <- c("SNP", 1:88) # This will give columns the same # as their family
SNP_info_melted <- melt(SNP_info_select, id.vars = "SNP")
colnames(SNP_info_melted) <- c("SNP", "family", "effects")
SNP_info_melted$effects <- as.numeric(SNP_info_melted$effects)

# Sum allele effects for each family
fam_sum_effects <- SNP_info_melted %>%
  group_by(family) %>%
  summarise(effect_sum = sum(effects, na.rm = T))

fam_sum_effects$family <- as.numeric(fam_sum_effects$family)

# Join Family means and allele effect sums
Mean_sum_all <- full_join(fam_heading_means, fam_sum_effects, by = "family")

# Plot relationship between family mean and sum of effects for ALL markers
summary(Mean_sum_all)
lm_all <- lm(fam_mean ~ effect_sum, data = Mean_sum_all)
summary(lm_all)
anova(lm_all)

ggplot(Mean_sum_all, aes(effect_sum, fam_mean)) +
  labs(x = "Sum of Allele Effects", y = "Family Mean Days to Flowering") +
  theme(axis.ticks = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = "grey"), axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 10)) +
  #scale_color_manual(name = "NAM Parent\nSubpopulations", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red", "gray75")) +
  geom_abline(aes(intercept = 53.56, slope = 0.036)) +
  geom_text(x = 50, y = 50, label = "y = 0.04x + 53.56") +
  geom_text(x = 50, y = 49.3, label = "R^2 == 0.50", parse = TRUE) +
  scale_y_continuous(limits = c(49, 59)) +
  scale_x_continuous(limits = c(-82, 100)) +
  geom_point(size = 1)
# Export in device landscape 3x3 80miss_MAF05/Mean_sum_all

# Get only significant markers
SNP_info_complete_QTL_sig <- filter(SNP_info_complete_QTL, lod > 2.4)

# Melt 
SNP_info_select_sig <- select(SNP_info_complete_QTL_sig, c(SNP, eff.1:eff.88))
colnames(SNP_info_select_sig) <- c("SNP", 1:88) # This will give columns the same # as their family
SNP_info_melted_sig <- melt(SNP_info_select_sig, id.vars = "SNP")
colnames(SNP_info_melted_sig) <- c("SNP", "family", "effects")
SNP_info_melted_sig$effects <- as.numeric(SNP_info_melted_sig$effects)

# Sum allele effects for each family
fam_sum_effects_sig <- SNP_info_melted_sig %>%
  group_by(family) %>%
  summarise(effect_sum = sum(effects, na.rm = T))

fam_sum_effects_sig$family <- as.numeric(fam_sum_effects_sig$family)

# Join Family means and allele effect sums
Mean_sum_sig <- full_join(fam_heading_means, fam_sum_effects_sig, by = "family")

# Plot relationship between family mean and sum of effects for ALL markers
summary(Mean_sum_sig)
lm_sig <- lm(fam_mean ~ effect_sum, data = Mean_sum_sig)
summary(lm_sig)
anova(lm_sig)

ggplot(Mean_sum_sig, aes(effect_sum, fam_mean)) +
  labs(x = "Sum of Significant Allele Effects", y = "Family Mean Days to Flowering") +
  theme(axis.ticks = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = "grey"), axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
  #scale_color_manual(name = "NAM Parent\nSubpopulations", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red", "gray75")) +
  geom_abline(aes(intercept = 53.58, slope = 0.036)) +
  geom_text(x = 50, y = 50, label = "y = 0.04x 53.58") +
  geom_text(x = 50, y = 49.3, label = "R^2 == 0.32", parse = TRUE) +
  scale_y_continuous(limits = c(49, 59)) +
  scale_x_continuous(limits = c(-82, 100)) +
  geom_point(size = 1)
# Export in device landscape 3x3 80miss_MAF05/Mean_sum_sig

#Get allele effects with fam subpop
Fam_QTL_vis <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/80miss_byfam_SNP_info.csv", header = T, check.names = F)
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA")
head(Fam_QTL_vis)
# Get useful columns
Fam_QTL_vis_few <- select(Fam_QTL_vis, c(SNP, Chrom_chr, Cumulative_bp, pval, 19:106))

Fam_QTL_melted <- melt(Fam_QTL_vis_few, id.vars = c("SNP", "Chrom_chr", "Cumulative_bp", "pval"))

# Get family number
Fam_num <- apply(Fam_QTL_melted, 1, function(x) strsplit(as.character(x), ".")[[1]][1])
Fam_num<-sub("eff.","",  Fam_QTL_melted$variable)
head(Fam_num)
melted_QTL_withfam <- cbind(Fam_QTL_melted, Fam_num)

# Reorder families by subpopulation and then from earliest parent hd to latest))
parent_heading_sorted <- arrange(parent_heading, man_sort_reverse)

#melted_QTL_withfam$fam_man <- factor(melted_QTL_withfam$Fam_num, levels = c(parent_heading$man_sort_reverse))
arranged_effects <- melted_QTL_withfam %>%
  mutate(Fam_num = factor(Fam_num, levels = parent_heading_sorted$new_fam)) %>%
  arrange(Fam_num)

nbins <- nrow(melted_QTL_withfam)/88

# Add subpopulation
# The number of rows in the rep is the number of families in the subpop*the number of bins
subpopulation_levels <- c(rep("Central European", (17*nbins)), rep("Coastal Mediterranean", (19*nbins)), rep("East African", (8*nbins)), rep("Asian", (19*nbins)), rep("Admixed", (23*nbins)), rep("Unassigned", (2*nbins)))
#below is wrong to test
#subpopulation_levels <- c(rep("Central European", (6*nbins)), rep("Coastal Mediterranean", (20*nbins)), rep("East African", (15*nbins)), rep("Asian", (22*nbins)), rep("Admixed", (15*nbins)), rep("Unassigned", (10*nbins)))

# Add subpopulation data
arranged_effects$subpop <- subpopulation_levels
head(arranged_effects)

# Reorder subpopulations
arranged_effects$subpop_sort = factor(arranged_effects$subpop, levels = c('Central European', 'Coastal Mediterranean', 'East African', 'Asian', 'Admixed', 'Un'))

# Do an ANOVA to see if allele effects for significant SNP vary by subpopulation
subpop_lm <- lm(value ~ factor(subpop), data = arranged_effects)
anova(subpop_lm)
summary(subpop_lm)

fam_lm <- lm(value ~ factor(Fam_num), data = arranged_effects)
summary(fam_lm)

# Remove allele effects for non significant markers
arranged_eff_sig <- filter(arranged_effects, pval > 3.308)

subpop_sig_lm <- lm(value ~ factor(subpop), data = arranged_eff_sig)
anova(subpop_sig_lm)
summary(subpop_sig_lm)

fam_sig_lm <- lm(value ~ factor(Fam_num), data = arranged_eff_sig)
summary(fam_sig_lm)

# Summarize familes BLUEs
library(dplyr)
library(ggplot2)
library(reshape2)

y <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/phenos_GxE.csv", header = T)

# Add a column with heading date in days 
y$DAP <- y$BLUE + 56.44459

BLUE_summary <- y %>% 
  group_by(family) %>%
  summarise(Minimum = min(DAP), 
            Maximum = max(DAP), 
            Mean = mean(DAP), 
            Range = (max(DAP) - min(DAP)), 
            Variance = var(DAP))

# Get the minimum and maximum mean in any family
all_min_Mean <- min(BLUE_summary$Mean) #48.94306
all_max_Mean <- max(BLUE_summary$Mean) #58.37186
  
# Get the minimum and maximum range in any family
all_min_range <- min(BLUE_summary$Range) #6.277296
all_max_range <- max(BLUE_summary$Range) #24.09969

# Get the minimum and maximum variance in any family
all_min_Variance <- min(BLUE_summary$Variance) #2.157817
all_max_Variance <- max(BLUE_summary$Variance) #29.32443

# Add subpopulations
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA", stringsAsFactors = F)

# Get family number and subpop
Parent_subpop <- select(parent_heading, c(new_fam, Pop_location))

# Rename the columns
colnames(Parent_subpop)[1] <- "family"

# Add subpops to summary
Summary_subpop <- full_join(Parent_subpop, BLUE_summary, by = "family")

# See if subpops have significant differences in heading date
heading_lm <- lm(Variance ~ Pop_location, data = Summary_subpop)
anova(heading_lm)
summary(heading_lm)
mean(Summary_subpop$Variance)

Sum_subpop <- Summary_subpop %>%
  group_by(Pop_location) %>%
    summarise(variance_sub = mean(Variance), 
              range_sub = mean(Range), 
              mean_sub = mean(Mean))
mean(Summary_subpop$Range)


