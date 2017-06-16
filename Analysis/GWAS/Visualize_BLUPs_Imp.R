# Visualize mapping results for manuscript
# 80miss_byfam 
# March 4, 2017
# Alex Ollhoff

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
RESULTS <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Imp_GxE/my_gwas_PolyTest_out_80miss_BLUE.csv", header = T)

# Load marker names
NAMES <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Imp_GxE/my_gwas_SNPs_out_80miss_BLUE.csv", header = T)

# Calculate the 5% false discovery rate threshold
THR = -log10(0.05 / ( ncol(RESULTS) * (1-0.05)))

# See how many SNPs are above the significance threshold
Results_sig <- filter(RESULTS, lod > THR) # 100 significant markers
w <- length(Results_sig)

# Add SNP names to mapping results
row.names(RESULTS)<-as.character(NAMES[,2])

# Get SNP names, separate the position and chromosome number
SNP_original_names <- as.data.frame(gsub("X", "", row.names(RESULTS)))
Physical_positions<- apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][2])
Chrom <- apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][1])
Chrom_chr<-gsub("H1|H2","H",Chrom)

SNPs_info <- cbind((SNP_original_names), Physical_positions, Chrom_chr)
SNP_info_complete <- cbind(SNPs_info, RESULTS)
row.names(SNP_info_complete) <- as.character(SNP_original_names[,1])

# Store the length of each chromsome part from ChromPartLength.txt
Chr1_L <- 312837513
Chr2_L <- 393532674
Chr3_L <- 394310633
Chr4_L <- 355061206
Chr5_L <- 380865482
Chr6_L <- 294822070
Chr7_L <- 325797516

# Get the first part of each chromosome
Chr1_1 <- grep("1H1", SNP_info_complete[,1])
Chr2_1 <- grep("2H1", SNP_info_complete[,1])
Chr3_1 <- grep("3H1", SNP_info_complete[,1])
Chr4_1 <- grep("4H1", SNP_info_complete[,1])
Chr5_1 <- grep("5H1", SNP_info_complete[,1])
Chr6_1 <- grep("6H1", SNP_info_complete[,1])
Chr7_1 <- grep("7H1", SNP_info_complete[,1])

# Get each chromosome part 2  
Chr1_2 <- grep("1H2", SNP_info_complete[,1])
Chr2_2 <- grep("2H2", SNP_info_complete[,1])
Chr3_2 <- grep("3H2", SNP_info_complete[,1])
Chr4_2 <- grep("4H2", SNP_info_complete[,1])
Chr5_2 <- grep("5H2", SNP_info_complete[,1])
Chr6_2 <- grep("6H2", SNP_info_complete[,1])
Chr7_2 <- grep("7H2", SNP_info_complete[,1])

# Add the end of each chromsome to the beginning
NewPositionsCHR1_2<-as.numeric(as.character(SNP_info_complete[Chr1_2,2])) + Chr1_L
NewPositionsCHR2_2<-as.numeric(as.character(SNP_info_complete[Chr2_2,2])) + Chr2_L
NewPositionsCHR3_2<-as.numeric(as.character(SNP_info_complete[Chr3_2,2])) + Chr3_L
NewPositionsCHR4_2<-as.numeric(as.character(SNP_info_complete[Chr4_2,2])) + Chr4_L
NewPositionsCHR5_2<-as.numeric(as.character(SNP_info_complete[Chr5_2,2])) + Chr5_L
NewPositionsCHR6_2<-as.numeric(as.character(SNP_info_complete[Chr6_2,2])) + Chr6_L
NewPositionsCHR7_2<-as.numeric(as.character(SNP_info_complete[Chr7_2,2])) + Chr7_L

# Add all chromsome parts together 
Cumulative_SNP_info_chr1 <- c(as.character(SNP_info_complete[Chr1_1,2]), NewPositionsCHR1_2) 
Cumulative_SNP_info_chr2 <- c(as.character(SNP_info_complete[Chr2_1,2]), NewPositionsCHR2_2) 
Cumulative_SNP_info_chr3 <- c(as.character(SNP_info_complete[Chr3_1,2]), NewPositionsCHR3_2) 
Cumulative_SNP_info_chr4 <- c(as.character(SNP_info_complete[Chr4_1,2]), NewPositionsCHR4_2) 
Cumulative_SNP_info_chr5 <- c(as.character(SNP_info_complete[Chr5_1,2]), NewPositionsCHR5_2) 
Cumulative_SNP_info_chr6 <- c(as.character(SNP_info_complete[Chr6_1,2]), NewPositionsCHR6_2) 
Cumulative_SNP_info_chr7 <- c(as.character(SNP_info_complete[Chr7_1,2]), NewPositionsCHR7_2) 

# get "positions" of UN markers
UNnames<-SNP_original_names[grep("UN", SNP_original_names[,1]),]
UNpositions<-gsub("UN_","", UNnames)

# Join all new cumulative bp positions
Cumulative_allCHR<-c(Cumulative_SNP_info_chr1, Cumulative_SNP_info_chr2, Cumulative_SNP_info_chr3, Cumulative_SNP_info_chr4, Cumulative_SNP_info_chr5, Cumulative_SNP_info_chr6, Cumulative_SNP_info_chr7, UNpositions)

# Add new column to data frame 
SNP_info_complete$Cumulative_bp <- Cumulative_allCHR

# Rename column
colnames(SNP_info_complete)[1] <- "SNP"

# UPDATE to corresponding mapping directory 
write.csv(SNP_info_complete, "~/Documents/PhD/NAM/NAM_mapping/Imp_GxE/80miss_byfam_SNP_info.csv")
SNP_info_complete <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Imp_GxE/80miss_byfam_SNP_info.csv", header = T)

###<<<*** Manhattan plot ***>>>###
THR =  -log10(0.05 / ( ncol(RESULTS) * (1-0.05))) #*** Ana, these THR are different but I can't figure out why. They should be identical'

THR =  -log10(0.05 / ( ncol(SNP_info_complete) * (1-0.05)))
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
library(reshape2)
library(ggplot2)
library(dplyr)

# Heatmap of allele effects 
Fam_QTL_vis <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Imp_GxE/80miss_byfam_SNP_info.csv", header = T, check.names = F)
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA")
head(Fam_QTL_vis)
# Get useful columns
Fam_QTL_vis_few <- select(Fam_QTL_vis, c(SNP, Chrom_chr, Cumulative_bp, lod, 19:107))

CHROM_LENGTH<-c(558535432, 768075024, 699711114, 647060158, 289164678, 583380513, 657224000)

#work one chromsoome at a time
for (n in 1:7){
  # Create an empty vector where empty bins will go
  BIN_empy<-NULL
  CHR<-Fam_QTL_vis_few[grep(paste(n,"H",sep=""), Fam_QTL_vis_few[,1]),]
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

# Look at allele effects
Fam_QTL_vis <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/80miss_byfam_SNP_info_QTL_byhand_nogenes.csv", header = T)

# Filter for SNPs and effects 
Few <- select(Fam_QTL_vis, c(QTL, pval, eff.1:eff.88))
#Few <- select(Fam_QTL_vis, c(SNP, eff.1:eff.88))

# Get allele effects for QTL2.1 
Effects_1.1 <- filter(Few, QTL == "QTL1.1" & pval >= 3.308)
Effects_2.1 <- filter(Few, QTL == "QTL2.1" & pval >= 3.308)
Effects_2.2 <- filter(Few, QTL == "QTL2.2" & pval >= 3.308)
Effects_2.3 <- filter(Few, QTL == "QTL2.3" & pval >= 3.308)
Effects_2.4 <- filter(Few, QTL == "QTL2.4" & pval >= 3.308)
Effects_3.1 <- filter(Few, QTL == "QTL3.1" & pval >= 3.308)
Effects_4.1 <- filter(Few, QTL == "QTL4.1" & pval >= 3.308)
Effects_4.2 <- filter(Few, QTL == "QTL4.2" & pval >= 3.308)
Effects_5.1 <- filter(Few, QTL == "QTL5.1" & pval >= 3.308)
Effects_5.2 <- filter(Few, QTL == "QTL5.2" & pval >= 3.308)
Effects_6.1 <- filter(Few, QTL == "QTL6.1" & pval >= 3.308)
Effects_7.1 <- filter(Few, QTL == "QTL7.1" & pval >= 3.308)
Effects_7.2 <- filter(Few, QTL == "QTL7.2" & pval >= 3.308)
Effects_7.3 <- filter(Few, QTL == "QTL7.3" & pval >= 3.308)

#Effects <- filter(Few, SNP == "2H1_28501792")

# Get largest allele effect in each family 
QTL1.1 = Effects_1.1 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL2.1 = Effects_2.1 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL2.2 = Effects_2.2 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL2.3 = Effects_2.3 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL2.4 = Effects_2.4 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL3.1 = Effects_3.1 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL4.1 = Effects_4.1 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL4.2 = Effects_4.1 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL5.1 = Effects_5.1 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL5.2 = Effects_5.2 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL7.1 = Effects_7.1 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL7.2 = Effects_7.2 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)
QTL7.3 = Effects_7.3 %>%
  summarise_each(funs(.[which.max(abs(.))]), eff.1:eff.88)

# Melt into long format
#melted_1.1 <- melt(Effects_max, id.vars = c("QTL", "pval"))
melted_1.1 <- melt(QTL1.1)
melted_2.1 <- melt(QTL2.1)
melted_2.2 <- melt(QTL2.2)
melted_2.3 <- melt(QTL2.3)
melted_2.4 <- melt(QTL2.4)
melted_3.1 <- melt(QTL3.1)
melted_4.1 <- melt(QTL4.1)
melted_4.2 <- melt(QTL4.2)
melted_5.1 <- melt(QTL5.1)
melted_7.1 <- melt(QTL7.1)
melted_7.2 <- melt(QTL7.2)
melted_7.3 <- melt(QTL7.3)

# Get allele effects for all QTL
QTL_effects <- cbind( melted_1.1, melted_2.1$value, melted_2.2$value, melted_2.3$value, melted_2.4$value, melted_3.1$value, 
                      melted_4.1$value, melted_5.1$value, melted_4.2$value, melted_7.1$value, melted_7.2$value, melted_7.3$value)
# melted_1.1, melted_2.3$value, melted_2.4$value, melted_5.1$value,melted_7.2$value,  melted_4.1$value,  

# rename columns
colnames(QTL_effects) <- c("family", "QTL1.1", "QTL2.1", "QTL2.2", "QTL2.3", "QTL2.4", "QTL3.1",
                           "QTL4.1", "QTL5.1", "QTL4.2", "QTL7.1",  "QTL7.2", "QTL7.3")
#"QTL1.1","QTL2.3", "QTL2.4", "QTL5.1",  "QTL7.2", "QTL4.1",   

# Rename families
parent_heading <- read.csv("~/Documents/PhD/NAM/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA", stringsAsFactors = F)

# Get family number and subpop
Parent_subpop <- select(parent_heading, c(new_fam, Pop_location))

# Rename the columns
colnames(Parent_subpop)[1] <- "family"

QTL_effects$family <- c(1:88)

# Get subpopulations
Effect_subpop <- full_join(Parent_subpop, QTL_effects, by = "family")

# melt
melted_eff <- melt(Effect_subpop, id.vars = c("family", "Pop_location"))

# Set values of zero to NA
melted_eff[melted_eff == 0] <- NA

ggplot(melted_eff, aes(value)) +
  labs(x = "Allele Effect", y = "Frequency") +
  facet_grid(variable ~ ., switch = "y", scales = "free") +
  scale_fill_manual(name = "BRIDG6 Diverse\nParent Subpopulation", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red", "grey55")) +
  theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks = element_blank(), panel.border = element_blank(), text = element_text(size = 14),
        legend.text = element_text(size = 14), strip.placement = "outside") +
  scale_x_continuous(expand = c(0,0)) +
  geom_histogram(binwidth = .2, aes(fill = factor(Pop_location)), na.rm = T)

lm_1.1 <- lm(QTL1.1 ~ Pop_location, data = Effect_subpop)
anova(lm_1.1)
summary(lm_1.1)

lm_2.1 <- lm(QTL2.1 ~ Pop_location, data = Effect_subpop)
anova(lm_2.1)
summary(lm_2.1)

Estimate Std. Error t value Pr(>|t|)   
(Intercept)                        -0.7151     0.3565  -2.006   0.0482 * 
  Pop_locationAsian                   0.4843     0.5301   0.914   0.3636   
Pop_locationCentral European        1.3617     0.5469   2.490   0.0148 * 
  Pop_locationCoastal Mediterranean  -0.6734     0.5301  -1.270   0.2075   
Pop_locationEast African            2.2412     0.7018   3.193   0.0020 **
  Pop_locationUn                      0.3317     1.2606   0.263   0.7931  

lm_2.2 <- lm(QTL2.2 ~ Pop_location, data = Effect_subpop)
anova(lm_2.2)
summary(lm_2.2)

lm_2.3 <- lm(QTL2.3 ~ Pop_location, data = Effect_subpop)
anova(lm_2.3)
summary(lm_2.3)

lm_2.4 <- lm(QTL2.4 ~ Pop_location, data = Effect_subpop)
anova(lm_2.4)
summary(lm_2.4)

lm_3.1 <- lm(QTL3.1 ~ Pop_location, data = Effect_subpop)
anova(lm_3.1)
summary(lm_3.1)

lm_4.1 <- lm(QTL4.1  ~ Pop_location, data = Effect_subpop)
anova(lm_4.1 )
summary(lm_4.1 )

lm_5.1 <- lm(QTL5.1 ~ Pop_location, data = Effect_subpop)
anova(lm_5.1)
summary(lm_5.1)

lm_7.1 <- lm(QTL7.1  ~ Pop_location, data = Effect_subpop)
anova(lm_7.1 )
summary(lm_7.1 )

(Intercept)                        0.23331    0.19218   1.214 0.228211    
Pop_locationAsian                 -1.13652    0.28573  -3.978 0.000149 ***
  Pop_locationCentral European       0.19230    0.29479   0.652 0.516007    
Pop_locationCoastal Mediterranean -0.32226    0.28573  -1.128 0.262673    
Pop_locationEast African          -0.07339    0.37830  -0.194 0.846665    
Pop_locationUn                    -1.22736    0.67945  -1.806 0.074524 .  


lm_7.2 <- lm(QTL7.2 ~ Pop_location, data = Effect_subpop)
anova(lm_7.2)
summary(lm_7.2)

lm_7.3 <- lm(QTL7.3  ~ Pop_location, data = Effect_subpop)
anova(lm_7.3)
summary(lm_7.3)

Estimate Std. Error t value Pr(>|t|)    
(Intercept)                        -0.5136     0.1568  -3.276  0.00154 ** 
  Pop_locationAsian                   0.6773     0.2331   2.906  0.00471 ** 
  Pop_locationCentral European        0.6992     0.2405   2.907  0.00469 ** 
  Pop_locationCoastal Mediterranean   0.2172     0.2331   0.932  0.35424    
Pop_locationEast African            1.5228     0.3086   4.934 4.15e-06 ***
  Pop_locationUn                      0.3465     0.5543   0.625  0.53364   

# Look at map of all NSGC colored by my colors

################ Map of NSGC #####################
GPS <- read.csv("~/Documents/PhD/NAM/Parents/NSGC_info_from_Maria.csv", header=T, sep = ",", na.strings="NA") 

mp <- NULL
mapWorld <- borders("world", colour="gray78", fill="gray50") # create a layer of borders
mp <- ggplot() +   mapWorld
mp <- mp + 
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.title = element_text(size = 12), axis.ticks = element_blank()) +
  geom_jitter(data = GPS, aes(x=LONGITUDE, y=LATITUDE) ,color=factor(GPS$color), size=1.5, width = 0.1, show.legend = T) 
  
mp


mp <- ggplot() +   mapWorld
mp <- mp+ geom_point(aes(x=GPS$LONGITUDE, y=GPS$LATITUDE) ,color=factor(GPS$color), size=2) +
  labs(x = "Latitude", y = "Longitude") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), title=element_text(size = 30)) 

mp

mp <- ggplot() +   mapWorld
mp <- mp + geom_point(data = GPS, aes(x=LONGITUDE, y=LATITUDE, color=color), size=2) +
  scale_colour_distiller(palette = "RdBu", "Allele\nEffect\n(Days)", na.value = NA) +
  labs(x = "Latitude", y = "Longitude") +
  scale_color_discrete(name = "BRIDG6 Diverse\nParent Subpopulation", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), title=element_text(size = 10)) 



###*** Annotations ***### This section no longer works
# Get mapping results
SNP_info_complete <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/80miss_byfam_SNP_info.csv", header = T, stringsAsFactors = F)
head(SNP_info_complete_few)
# Filter for significant markers
#SNP_info_complete_sig <- filter(SNP_info_complete, lod > 3.275)

# Filter for marker info
SNP_info_complete_few <- select(SNP_info_complete, c(SNP, Chrom_chr, Cumulative_bp))

SNP_info_complete_few$Platform <- "QTL"

# Remove last row
SNP_info_complete_few <- SNP_info_complete_few[-36932,]

# Get marker annoataions
annotated_markers <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Flowering_description/FT_annotations_select.csv", header = T, stringsAsFactors = F)

Annotated_few <- select(annotated_markers, c(Abrevation, Gene_Ch, Gene_mid))
head(Annotated_few)

# Filter annotations for genes on chromsomes, not Unk
Annotated_few <- filter(Annotated_few, Gene_Ch != "UN")

# Rename SNP positions 
colnames(Annotated_few)[1] <- "SNP"
colnames(Annotated_few)[2] <- "Chrom_chr"
colnames(Annotated_few)[3] <- "Cumulative_bp"
Annotated_few$Platform <- "Annotations"

#Annotated_few$Cumulative_bp <- as.character(Annotated_few$Cumulative_bp)

# Combine SNPs
SNP_annotation <- rbind.data.frame(SNP_info_complete_few, Annotated_few)
head(Annotated_few)

#SNP_annotation$Cumulative_bp <- as.numeric(SNP_annotation$Cumulative_bp)

# Reorder facet levels
#SNP_annotation$platform_sorted = factor(SNP_annotation$Platform, levels = c('Annotations', 'QTL'))
library(ggplot2)
# install.packages("devtools")
# devtools::install_github("slowkow/ggrepel")
library(ggrepel)
SNP_annotation <- as.data.frame(SNP_annotation)

#write.csv(SNP_annotation, "~/Documents/PhD/NAM/NAM_mapping/Flowering_description/Annotations_GBS_combined.csv")

# Plot annotations, GBS markers
ggplot(data = SNP_annotation, aes(as.numeric(Cumulative_bp), stat = "identity") +
         labs(x = "Chromosome (base pair)") +
         theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(color = "grey")) +
         scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
         #scale_x_continuous(expand = c(0,0)) +
         geom_histogram(na.rm = TRUE, show.legend = F, fill = "white") +
         geom_text(data = Annotated_few, na.rm = TRUE, aes(x = Cumulative_bp, y = 1, label = SNP), angle = 270, size = 3, hjust = 0) +
         facet_grid( ~ Chrom_chr, scales = "free_x", space = "free_x", switch = "x") 
       #geom_text_repel(data = annotated_markers, na.rm = TRUE, aes(x = Cumulative_bp, y = 1, label = SNP), angle = 270, size = 4) 
       
# export in 2x10 landscape device size as "Annotations" in 80miss_byfam
       
       
       ###*** Marker density across the genome ***### This section is no longer relevant
       # Color genes in annotation that are within 5cM of a QTL
       # Import all genotypes
       SNP_info_complete <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/80miss_byfam_SNP_info.csv", header = TRUE, stringsAsFactors = F)
       as.data.frame(SNP_info_complete)
       
       # Add platform
       SNP_info_complete$Platform <- "GBS"
       
       # Import 9K
       old_SNPs <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Genotypes/BOPA_9K_phys_2_7_17.csv", header = TRUE, stringsAsFactors = FALSE)
       head(old_SNPs)
       
       # Rename SNP positions 
       colnames(old_SNPs)[3] <- "Cumulative_bp"
       
       # Get only 9K SNPs
       NineK <- filter(old_SNPs, Platform == "9K")
       
       # Combine SNPs
       GBS_9K_BOPA <- intersect(colnames(SNP_info_complete), colnames(NineK))
       GBS_9K_BOPA_SNPs <- rbind(subset(SNP_info_complete, select = GBS_9K_BOPA), subset(NineK, select = GBS_9K_BOPA))
       
       head(GBS_9K_BOPA_SNPs)
       
       # Filter for only SNPs on chromsomses
       chroms <- c("1H", "2H", "3H", "4H","5H", "6H","7H")
       
       SNPs_chrom <- filter(GBS_9K_BOPA_SNPs, Chrom_chr %in% chroms)
       
       # Save
       write.csv(SNPs_chrom, "~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/GBS_9K.csv")
       
       # Reorder platforms 
       SNPs_chrom$Platform_sorted = factor(SNPs_chrom$Platform, levels = c('GBS', '9K'))
       head(SNPs_chrom)
       # Plot SNPs as a percent of total markers on each chromsome per SNP platform
       ggplot(SNPs_chrom, aes(as.numeric(Cumulative_bp))) +
         scale_fill_manual(values = c("black", "white"), "Genotyping\nPlatform") +
         labs(x = "Chromosome (base pairs)", y = "Marker Density", title = element_blank()) +    
         theme(legend.position = c(0.95,0.8), strip.background = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         facet_grid( ~ Chrom_chr, scales = "free", space = "free", switch = "x") +
         scale_y_continuous(expand = c(0,0)) +
         scale_x_continuous(expand = c(0,0)) +
         geom_density(aes(fill = factor(Platform)), alpha = 2/3, na.rm = T, adjust = 1/8) # Export in device size 3x8 "SNP_density"
       dev.off()
       
       
