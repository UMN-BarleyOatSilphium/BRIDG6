# Author: Alex Ollhoff
# Description: Prep genos and phenos for mapping
# genos: 
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
# Import genotypes. SNPs are in rows, and samples are in collumns. Genotypes have been imputed.
gen_raw <- read.table("~/Documents/PhD/NAM/NAM_mapping/Genotypes/Final_filt_Fillin_noMis_noNAHHras_rasBased_MAF5_noLD.txt", header = T)
# Turn table to have SNPs in columns and samples in rows
gen_raw_t <- t(gen_raw)
# Save genotypes for later use
write.csv(gen_raw_t, "~/Documents/PhD/NAM/NAM_mapping/Genotypes/Final_filt_Fillin_noMis_noNAHHras_rasBased_MAF5_noLD_t.csv")

# import files for processing
# I manually added the population mean to the adjusted mean for each RIL (DAP: 56.44459, height: 70.11302757)
# Import phenotypes
y <- read.csv("~/Desktop/BRIDG6/Phenotypes/GxE/DAP_BLUPs_nob_no2_no19.csv", header=T)

# Import genotypes. SNPs are in columns and samples in rows.
gen <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Genotypes/Final_filt_Fillin_noMis_noNAHHras_rasBased_MAF5_noLD_t.csv", header = T)

# Recode genotypes so Ras = 2, diverse parent = 0
#Ras_based_recode <- function(dat){
#  dat[which(dat == "0")]<-"3"
#  dat[which(dat == "2")]<-"4"
#  dat[which(dat == "1")]<-"1"
#  dat[which(dat == "NA")]<-"NA"
#  return(dat)
#}

#gen_recode <- as.data.frame(t(apply(gen, 1, Ras_based_recode)))

gen_recode<- gen

# See what I'm working with
class(y)
class(gen)
dim(y)
dim(gen)

# Remove Ras from the first row
gen_noRas <- gen[-1,]

# Filter for rows present in both genotypes and phenotypes
gen_1 <- gen_noRas[match(y$line_name, gen_noRas$X, nomatch = NA, incomparables = F),]

# Remove rows where all genos are NA in phenos
gen_2 <- gen_1[rowSums(is.na(gen_1)) != ncol(gen_1),]

# filter for rows present in both in phenos
y_1 <- y[match(gen_2$X, y$line_name, nomatch = NA, incomparables = F),]
y_2 <- y_1[rowSums(is.na(y_1)) != ncol(y_1),]

# View(gen_2)
# View(y_2)

# Remove line names
gen_naked <- select(gen_2, -(X))

# Remove markers with MAF below 0.05 after intersect
adjusted_genotypes = snpQC( gen=gen_naked, MAF=0.05, impute=FALSE)
rownames(adjusted_genotypes) = y_2$line_name
dim(adjusted_genotypes)

# 5141 individuals with 6758 markers

# Look at the data a little before saving it for mapping

# Count the number of individuals with each genotype at each marker
na_count <-sapply(adjusted_genotypes, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count$name<-rownames(na_count)

# Count NAs for each individual
genos_for_viz_t <- t(adjusted_genotypes)
na_count_indiv <-sapply(genos_for_viz_t, function(y) sum(length(which(is.na(y)))))
na_count_indiv <- data.frame(na_count_indiv)
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
geno_freq <- cbind(alt_count, Ras_count, het_count, na_count)

# remove marker names
geno_freq_1 <- geno_freq[,-2]
geno_freq_1 <- geno_freq_1[,-3]
geno_freq_1 <- geno_freq_1[,-4]
geno_freq_1 <- geno_freq_1[,-5]
geno_freq_1$name<-rownames(geno_freq)

geno_frequency <- geno_freq_1[,-6]
library(reshape2)
melted_freq <- melt(geno_frequency, id.vars = "name")

# Plot the number of markers at any given frequency in the BRIDG6
ggplot(melted_freq, aes(value, fill = factor(variable))) +
  labs(x = "Number of Observations in the BRIDG6", y = "Number of Markers") +
  scale_fill_manual(name="Genotype",values= c("grey25", "yellow", "orange", "grey75"), labels = c("Rasmusson", "Donor Parent", "Heterozygous", "Missing")) +
  geom_histogram(bins = 200, alpha = 0.6)

# Get the number of markers on each chromosome
chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
class(chr) # should be integer

# Make file for MSI
write.csv(adjusted_genotypes, "~/Documents/PhD/NAM/NAM_mapping/Imp_GxE/genos_fillin.csv")
write.csv(y_2, "~/Documents/PhD/NAM/NAM_mapping/Imp_GxE/phenos_fillin.csv")
