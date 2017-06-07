# Author: Alex Ollhoff
# Description: Prep genos and phenos for mapping
# genos:
  # It has been filtered using the following criteria:
  # 1. Remove SNPs that were NA or Hete in Rasmusson
  # 2. Call genotypes respect to Rasmusson (0,2,1 and NA)
  # 3. Remove SNPs missing more than 80% of data
  # 4. Remove SNPs in the Unknown chromosome
  # 5. Remove SNPs in >=80% LD 
  # 6. Remove SNPs with MAF <0.05
  # The filtered imputed data for FILLIN has 6846 SNPs for 6033 samples.
# phenos: BLUPs calculated across 3 environments with fam as fixed, line and env as random, GxE
#######################################################################################################################################
rm(list=ls())
# load packages
library(NAM)
library(ggplot2)
library(dplyr)

# Format raw data files for processing
# Import genotypes. SNPs are in rows, and samples are in columns. Genotypes have been imputed. Genotypes coded based on Rasmusson = 2, diverse parent = 0, Heterozygous=1, missing=NA
gen_raw <- read.table("~/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/Final_filt_Fillin_noMis_noNAHHras_rasBased_MAF5_noLD.txt", header = T)
# Turn table to have SNPs in columns and samples in rows
gen_raw_t <- t(gen_raw)
# Save genotypes for later use
write.csv(gen_raw_t, "~/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/Genotypes/Final_filt_Fillin_noMis_noNAHHras_rasBased_MAF5_noLD_t.csv")

# import files for processing
# I manually added the population mean to the adjusted mean for each RIL (DAP: 56.44459, height: 70.11302757)
# Import phenotypes
y <- read.csv("~/Dropbox/GITHUB/BRIDG6/Datasets/Phenotypes/DAP_BLUPs_nob_no2_no19.csv", header=T)

# Import genotypes. SNPs are in columns and samples in rows.
gen <- read.csv("~/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/Genotypes/Final_filt_Fillin_noMis_noNAHHras_rasBased_MAF5_noLD_t.csv", header = T)

# See what I'm working with
class(y)
class(gen)
dim(y)
dim(gen)

# Remove Ras from the first row
gen_noRas <- gen[-1,]

# Filter for rows present in both genotypes and phenotypes
gen_1<-gen_noRas[(gen_noRas[,1] %in% y[,1] ),]

# Remove rows where all genos are NA in phenos
gen_2 <- gen_1[rowSums(is.na(gen_1)) != ncol(gen_1),]

print (paste("There were", dim(gen_1)[1] - dim(gen_2)[1], "samples with all SNPs missing", sep=" "))

# filter for rows present in both in phenos
y_1 <- y[match(gen_2$X, y$line_name, nomatch = NA, incomparables = F),]
y_2 <- y_1[rowSums(is.na(y_1)) != ncol(y_1),]

# Sort genotypes to match samples in phenotypes
gen_2_or<-gen_2[match(y_2[,1] , gen_2[,1]),]

# Confirm that samples are in the same order in the genotypes and phenotypes sets
if (identical(as.character(gen_2_or[,1]), as.character(y_2[,1])) == FALSE) stop("Phenotypes and Genotypes samples do not match")
# View(gen_2)
# View(y_2)

# Remove line names
gen_naked <- select(gen_2_or, -(X))

# Remove markers with MAF below 0.05 after intersect
adjusted_genotypes = snpQC( gen=gen_naked, MAF=0.05, impute=FALSE)
rownames(adjusted_genotypes) = y_2$line_name
dim(adjusted_genotypes)

# 5141 individuals with 8101 markers

# Look at the data a little before saving it for mapping

# Count the number of individuals with each genotype at each marker
na_count <-sapply(adjusted_genotypes, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count$name<-rownames(na_count)

# Count hets
het_count <-sapply(adjusted_genotypes, function(y) sum(length(which(y == 1))))
het_count <- data.frame(het_count)
het_count$name<-rownames(het_count)

# Count alt alleles
alt_count <-sapply(adjusted_genotypes, function(y) sum(length(which(y == 0))))
alt_count <- data.frame(alt_count)
alt_count$name<-rownames(alt_count)

ggplot(alt_count, aes(alt_count)) +
  geom_histogram(bins = 50)

# Count Ras alleles
Ras_count <-sapply(adjusted_genotypes, function(y) sum(length(which(y == 2))))
Ras_count <- data.frame(Ras_count)
Ras_count$name<-rownames(Ras_count)

ggplot(Ras_count, aes(Ras_count)) +
  geom_histogram(bins = 50)

# Combine columns
#geno_freq <- cbind(alt_count, Ras_count, het_count, na_count)
geno_freq <- cbind(alt_count, het_count,Ras_count )

# remove marker names
geno_frequency <- geno_freq[,-c(2,4,6)]
geno_frequency$name<-rownames(geno_freq)

# Proportions by individuals

library(reshape2)
melted_freq <- melt(geno_frequency, id.vars = "name")

# Plot the frequency of each genotype class at each marker in the BRIDG6
ggplot(melted_freq, aes(value, fill = factor(variable))) +
  labs(x = "Number of Observations in the BRIDG6", y = "Number of Markers") +
  scale_fill_manual(name="Genotype",values= c("Blue", "yellow", "red"), labels = c( "Donor Parent", "Heterozygous","Rasmusson")) +
  geom_histogram(bins = 50, alpha = 0.6)



# Get the number of markers on each chromosome
chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
class(chr) # should be integer

# Make file for MSI
write.csv(adjusted_genotypes, "~/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/genos_fillin.csv")
write.csv(y_2, "~/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/phenos_fillin.csv")
