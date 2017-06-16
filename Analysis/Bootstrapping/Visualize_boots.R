# Visualize results of GWAS in bootstrapped samples of the NAM population
# Alex Ollhoff
# February 17, 2017

library(reshape2)
library(ggplot2)
library(dplyr)

# Import the number of markers found significant in each sample size

Sample_250 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_250/Number_sig_SNPs.txt", header = F)
Sample_500 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_500/Number_sig_SNPs.txt", header = T)
Sample_1000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_1000/Number_sig_SNPs.txt", header = T)
Sample_2000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_2000/Number_sig_SNPs.txt", header = T)
Sample_3000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_3000/Number_sig_SNPs.txt", header = F)
Sample_4000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_4000/Number_sig_SNPs.txt", header = F)
Sample_5000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping//100_repsSample_5000/Number_sig_SNPs.txt", header = F)

# Join them together into a table
Boot_all <- cbind(Sample_250, Sample_500, Sample_1000, Sample_2000)

# Give each column a header
names(Boot_all) <- c("Sample_250", "Sample_500", "Sample_1000", "Sample_2000")
, "Sample_3000", "Sample_4000", "Sample_5000"
# Melt data into long format
melted_Boot <- melt(Boot_all)
head(melted_Boot)

# get correlation between sample size and # significant markers 
summary(Boot_all)
size_QTL <- lm(value ~ variable, data = melted_Boot)
summary(size_QTL)
anova(size_QTL)

# Plot the number of significant markers in samples of each size with confidence interval 
ggplot(melted_Boot, aes(variable, value)) +
  labs(x = "Sample Size", y = "Number of Significant Markers") +
  theme(legend.background = element_rect(color = "white"), axis.ticks = element_blank(), panel.background = element_rect(color = "grey"), panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c(250, 500, 1000, 2000, 3000, 4000, 5000)) +
  ylim(0, 300) +
  geom_boxplot() +
  geom_point(size = 1) +
  geom_hline(yintercept = 293) 
# Export in device size 3x2 landscape to /PhD/Manuscripts/Figures/Sample_box

# Import the number of markers found significant in each sample size

Sample_5 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/5_fams/Number_sig_SNPs.txt", header = F)
Sample_10 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/10_fams/Number_sig_SNPs.txt", header = F)
Sample_25 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/25_fams/Number_sig_SNPs.txt", header = F)
Sample_50 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/50_fams/Number_sig_SNPs.txt", header = F)
Sample_75 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/75_fams/Number_sig_SNPs.txt", header = F)

# Join them together into a table
Boot_all <- cbind(Sample_5, Sample_10, Sample_25, Sample_50, Sample_75)

# Give each column a header
names(Boot_all) <- c("Sample_5", "Sample_10", "Sample_25", "Sample_50", "Sample_75")

# Melt data into long format
melted_Boot <- melt(Boot_all)
head(melted_Boot)

# get correlation between sample size and # significant markers 
summary(Boot_all)
size_QTL <- lm(value ~ variable, data = melted_Boot)
summary(size_QTL)
anova(size_QTL)

# Plot the number of significant markers in samples of each size with confidence interval 
ggplot(melted_Boot, aes(variable, value)) +
  labs(x = "Sample Size", y = "Number of Significant Markers") +
  theme(legend.background = element_rect(color = "white"), axis.ticks = element_blank(), panel.background = element_rect(color = "grey"), panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c(5, 10, 25, 50, 75)) +
  ylim(0, 300) +
  geom_boxplot() +
  geom_point(size = 1) +
  geom_hline(yintercept = 293) 
# Export in device size 3x2 landscape to /PhD/Manuscripts/Figures/Sample_box

# Get all mapping results 
all_results <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/80miss_byfam_SNP_info.csv", header = T)

# Define QTL region with list of all markers, check significant marker list against that
Hv_ELF3_SNPs <- all_results$SNP = c(1H2_243985899:2H1_8905766)
Hv_ELF3_SNPs <- grep(1H2_243985899:2H1_8905766, all_results$SNP)

Hv_ELF3 <- filter(all_results, SNP == c(1H2_243985899:2H1_8905766))

# Make file of markers that are in QTL with each QTL named based on whether there is a nearby gene
# Filter bootstrap results for markers that are in QTL 

# Import SNP significance counts
count_SNP <- read.csv 
# For each sample size count number of QTL that are represented
# Plot the number of QTL detected by at least 1 of the same markers in each bootstrap sample


# Import which markers were significant and the number of times they were found to be
Sig_250 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_250/Count_sig_per_SNP.txt", header = T)
Sig_500 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_500/Count_sig_per_SNP.txt", header = T)
Sig_1000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_1000/Count_sig_per_SNP.txt", header = T)
Sig_2000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_2000/Count_sig_per_SNP.txt", header = T)
Sig_3000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_3000/Count_sig_per_SNP.txt", header = T)
Sig_4000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_4000/Count_sig_per_SNP.txt", header = T)
Sig_5000 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/Sample_5000/Count_sig_per_SNP.txt", header = T)

All_SNPs <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/my_gwas_PolyTest_out_80miss_BLUE.csv", header = T)
SNP_list <- select(All_SNPs, c(X))
SNP_list$count_all <- 1

# Rename columns
colnames(SNP_list) <- c("SNP", "count_all")
colnames(Sig_250) <- c("SNP", "count_250")
colnames(Sig_500) <- c("SNP", "count_500")
colnames(Sig_1000) <- c("SNP", "count_1000")
colnames(Sig_2000) <- c("SNP", "count_2000")
colnames(Sig_3000) <- c("SNP", "count_3000")
colnames(Sig_4000) <- c("SNP", "count_4000")
colnames(Sig_5000) <- c("SNP", "count_5000")

# Join results in a table, set uncommon mSNp significance counts to NA
Sig_sum_all <- full_join(SNP_list, Sig_250, by = "SNP")
Sig_sum_1 <- full_join(Sig_sum_all, Sig_500, by = "SNP")
Sig_sum_2 <- full_join(Sig_sum_1, Sig_1000, by = "SNP")
Sig_sum_3 <- full_join(Sig_sum_2, Sig_2000, by = "SNP")
Sig_sum_4 <- full_join(Sig_sum_3, Sig_3000, by = "SNP")
Sig_sum_5 <- full_join(Sig_sum_4, Sig_4000, by = "SNP")
Sig_sum_all <- full_join(Sig_sum_5, Sig_5000, by = "SNP")

# Get SNP names, cumulative bp position, and chromosome
SNP_original_names <- as.data.frame(gsub("X", "", Sig_sum_all$SNP))
Physical_positions<- apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][2])
Chrom <- apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][1])
Chrom_chr<-gsub("H1|H2","H",Chrom)

SNPs_info <- cbind((SNP_original_names), Physical_positions, Chrom_chr)

head(SNPs_info)

SNP_info_complete <- cbind(SNPs_info, Sig_sum_all)
SNP_info_complete[1:10,1:6]

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

# get positions of UN markers
UNnames<-SNP_original_names[grep("UN", SNP_original_names[,1]),]
UNpositions<-gsub("UN_","", UNnames)

# Join all new cumulative bp positions
Cumulative_allCHR<-c(Cumulative_SNP_info_chr1, Cumulative_SNP_info_chr2, Cumulative_SNP_info_chr3, Cumulative_SNP_info_chr4, Cumulative_SNP_info_chr5, Cumulative_SNP_info_chr6, Cumulative_SNP_info_chr7, UNpositions)

# Add new column to data frame 
SNP_info_complete$Cumulative_bp <- Cumulative_allCHR

# Rename columns
colnames(SNP_info_complete)[1] <- "SNP"
head(SNP_info_complete)

# UPDATE to corresponding mapping directory 
write.csv(SNP_info_complete, "~/Documents/PhD/NAM/NAM_mapping/bootstrapping/Boot_SNP_info.csv")
Boot_info_complete <- read.csv("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/Boot_SNP_info.csv", header = T)

# Melt data into long format
Boots_melted <- melt(Boot_info_complete, id.vars = c("X", "SNP", "Physical_positions", "Chrom_chr", "SNP.1", "Cumulative_bp"))

# Plot each sample size with frequency of detection on the y axis
ggplot(Boots_melted, aes(x = as.numeric(Cumulative_bp), y = value)) +
  facet_grid( variable ~ Chrom_chr, scales = "free", space = "free", switch = "both" ) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(color = "grey")) +
  scale_x_continuous(expand = c(0,0)) +
  geom_point(size = 0.3, alpha = 0.5)

### Numbers of families
#Import which markers were significant and the number of times they were found to be
Sig_5 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/5_fams/Count_sig_per_SNP.txt", header = T)
Sig_10 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/10_fams/Count_sig_per_SNP.txt", header = T)
Sig_25 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/25_fams/Count_sig_per_SNP.txt", header = T)
Sig_50 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/50_fams/Count_sig_per_SNP.txt", header = T)
Sig_75 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/75_fams/Count_sig_per_SNP.txt", header = T)

All_SNPs <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/my_gwas_PolyTest_out_80miss_BLUE.csv", header = T)
SNP_list <- select(All_SNPs, c(X))
SNP_list$count_all <- 1

# Rename columns
colnames(SNP_list) <- c("SNP", "count_all")
colnames(Sig_5) <- c("SNP", "count_5")
colnames(Sig_10) <- c("SNP", "count_10")
colnames(Sig_25) <- c("SNP", "count_25")
colnames(Sig_50) <- c("SNP", "count_50")
colnames(Sig_75) <- c("SNP", "count_75")

# Join results in a table, set uncommon mSNp significance counts to NA
Sig_sum_all <- full_join(SNP_list, Sig_5, by = "SNP")
Sig_sum_1 <- full_join(Sig_sum_all, Sig_10, by = "SNP")
Sig_sum_2 <- full_join(Sig_sum_1, Sig_25, by = "SNP")
Sig_sum_3 <- full_join(Sig_sum_2, Sig_50, by = "SNP")
Sig_sum_4 <- full_join(Sig_sum_3, Sig_75, by = "SNP")

# Get SNP names, cumulative bp position, and chromosome
SNP_original_names <- as.data.frame(gsub("X", "", Sig_sum_4$SNP))
Physical_positions<- apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][2])
Chrom <- apply(SNP_original_names,1, function(x) strsplit(as.character(x), "_")[[1]][1])
Chrom_chr<-gsub("H1|H2","H",Chrom)

SNPs_info <- cbind((SNP_original_names), Physical_positions, Chrom_chr)

head(SNPs_info)

SNP_info_complete <- cbind(SNPs_info, Sig_sum_4)
SNP_info_complete[1:10,1:10]

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

# get positions of UN markers
UNnames<-SNP_original_names[grep("UN", SNP_original_names[,1]),]
UNpositions<-gsub("UN_","", UNnames)

# Join all new cumulative bp positions
Cumulative_allCHR<-c(Cumulative_SNP_info_chr1, Cumulative_SNP_info_chr2, Cumulative_SNP_info_chr3, Cumulative_SNP_info_chr4, Cumulative_SNP_info_chr5, Cumulative_SNP_info_chr6, Cumulative_SNP_info_chr7, UNpositions)

# Add new column to data frame 
SNP_info_complete$Cumulative_bp <- Cumulative_allCHR

# Rename columns
colnames(SNP_info_complete)[1] <- "SNP"
head(SNP_info_complete)

# UPDATE to corresponding mapping directory 
write.csv(SNP_info_complete, "~/Documents/PhD/NAM/NAM_mapping/bootstrapping/Boot_SNP_info_fams.csv")
Boot_info_complete <- read.csv("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/Boot_SNP_info_fams.csv", header = T)

# Melt data into long format
Boots_melted <- melt(Boot_info_complete, id.vars = c("X", "SNP", "Physical_positions", "Chrom_chr", "SNP.1", "Cumulative_bp"))

# Plot each sample size with frequency of detection on the y axis
ggplot(Boots_melted, aes(x = as.numeric(Cumulative_bp), y = value)) +
  facet_grid( variable ~ Chrom_chr, scales = "free", space = "free", switch = "both" ) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(color = "grey")) +
  scale_x_continuous(expand = c(0,0)) +
  geom_point(size = 0.3, alpha = 0.5)

### Plot detection of QTL ###

# Get mapping results
# This file was made manually by grouping markers within 5mbp of any significant marker
SNP_info_complete_QTL <- read.csv("~/Documents/PhD/NAM/NAM_mapping/BLUPs_GxE/80miss_byfam_SNP_info_QTL_byhand_nogenes.csv", header = T)
head(SNP_info_complete_QTL)

# Merge QTL and bootstraps
### Numbers of families
#Import which markers were significant and the number of times they were found to be
Sig_5 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/5_fams/Count_sig_per_SNP.txt", header = T)
Sig_10 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/10_fams/Count_sig_per_SNP.txt", header = T)
Sig_25 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/25_fams/Count_sig_per_SNP.txt", header = T)
Sig_50 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/50_fams/Count_sig_per_SNP.txt", header = T)
Sig_75 <- read.table("~/Documents/PhD/NAM/NAM_mapping/bootstrapping/100_reps/75_fams/Count_sig_per_SNP.txt", header = T)

SNP_list <- select(SNP_info_complete_QTL, c(QTL, X))
SNP_list$count_all <- 1

# Rename columns
colnames(SNP_list) <- c("QTL", "SNP", "count_all")
colnames(Sig_5) <- c("SNP", "count_5")
colnames(Sig_10) <- c("SNP", "count_10")
colnames(Sig_25) <- c("SNP", "count_25")
colnames(Sig_50) <- c("SNP", "count_50")
colnames(Sig_75) <- c("SNP", "count_75")

# Join results in a table, set uncommon mSNp significance counts to NA
Sig_sum_all <- full_join(SNP_list, Sig_5, by = "SNP")
Sig_sum_1 <- full_join(Sig_sum_all, Sig_10, by = "SNP")
Sig_sum_2 <- full_join(Sig_sum_1, Sig_25, by = "SNP")
Sig_sum_3 <- full_join(Sig_sum_2, Sig_50, by = "SNP")
Sig_sum_4 <- full_join(Sig_sum_3, Sig_75, by = "SNP")

# Replace all NA with O
Sig_for_counts <- as.data.frame(Sig_sum_4)
Sig_for_counts[is.na(Sig_for_counts)] <- 0


# Get markers in each QTL
QTL_count <- Sig_for_counts %>%
  group_by(QTL) %>%
  summarise(QTL_n = n(), 
            QTL_5 = max(count_5), 
            QTL_10 = max(count_10), 
            QTL_25 = max(count_25), 
            QTL_50 = max(count_50), 
            QTL_75 = max(count_75))
            
# Look at random samples
# Rename columns
colnames(SNP_list) <- c("QTL", "SNP", "count_all")
colnames(Sig_250) <- c("SNP", "count_250")
colnames(Sig_500) <- c("SNP", "count_500")
colnames(Sig_1000) <- c("SNP", "count_1000")
colnames(Sig_2000) <- c("SNP", "count_2000")
colnames(Sig_3000) <- c("SNP", "count_3000")
colnames(Sig_4000) <- c("SNP", "count_4000")
colnames(Sig_5000) <- c("SNP", "count_5000")

# Join results in a table, set uncommon mSNp significance counts to NA
Ran <- full_join(SNP_list, Sig_250, by = "SNP")
Sig_sum_1 <- full_join(Ran, Sig_500, by = "SNP")
Sig_sum_2 <- full_join(Sig_sum_1, Sig_1000, by = "SNP")
Sig_sum_3 <- full_join(Sig_sum_2, Sig_2000, by = "SNP")
Sig_sum_4 <- full_join(Sig_sum_3, Sig_3000, by = "SNP")
Sig_sum_5 <- full_join(Sig_sum_4, Sig_4000, by = "SNP")
Sig_sum_ran <- full_join(Sig_sum_5, Sig_5000, by = "SNP")

# Replace all NA with O
Sig_ran_counts <- as.data.frame(Sig_sum_ran)
Sig_ran_counts[is.na(Sig_ran_counts)] <- 0

# Get markers in each QTL
QTL_count_ran <- Sig_ran_counts %>%
  group_by(QTL) %>%
  summarise(QTL_n = n(), 
            QTL_250 = max(count_250), 
            QTL_500 = max(count_500), 
            QTL_1000 = max(count_1000), 
            QTL_2000 = max(count_2000), 
            QTL_3000 = max(count_3000), 
            QTL_4000 = max(count_4000), 
            QTL_5000 = max(count_5000))

