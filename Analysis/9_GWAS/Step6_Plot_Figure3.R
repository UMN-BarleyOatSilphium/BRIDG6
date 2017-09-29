# Author: Alex Ollhoff & Ana Poets
# Description: Visualize mapping resutls- Figure 3.
#######################################################################################################################################
rm(list=ls())
library(ggplot2)
library(dplyr)

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

