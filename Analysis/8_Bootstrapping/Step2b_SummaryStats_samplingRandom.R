# Author: Ana Poets
# Description: Calculate the families and number of individuals per family 
# represented every iteration, during the random individual selection strategy.
# Sampling random number across the BRIDG6 (264,528,968,2012)
# Information used to add family size information to table S4.
##########################################################################################################################
rm(list=ls())

# load package
require(NAM)
#library(dplyr)

# import files for processing. Phenotype and Genotypes have been filtered after imputation of genotypes.
# only samples with phenotypic and genotypic data are included. SNPs in LD >0.8 have been removed, and genotypes are 
# coded based on Rasmusson allele.
# Data set was generated using Step1_Prep_Imp_GxE.R

# import phenotypes
y <- read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Input/phenos_LDKNNI.csv", header=T)

# import genotypes
gen_noRas <- read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Input/genos_LDKNNI.csv", header = T)

#LD <- read.table("~/Documents/PhD/NAM/NAM_mapping/Genotypes/List_SNPtoKeep_noLD80.txt")

# See what I'm working with
class(y)
class(gen_noRas)
dim(y)
dim(gen_noRas)


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
if (identical(as.character(gen_2[,1]), as.character(y_2[,2])) == FALSE)stop("Error! sample order in phen and gen is not the same")
gen_naked <- gen_2[,-c(1)]

# Marker QC + imputation
adjusted_genotypes = snpQC( gen=gen_naked, MAF=0.05, impute=FALSE)
rownames(adjusted_genotypes) = y_2$line_name
dim(adjusted_genotypes)

# 8,021 markers

# getting chr
chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
class(chr)

SAMPLE_NAMES<-y_2[,2]


##### for bootstrapping even number of individuals per family
TotalNumber<-c(264,528,968,2024)

#Summary of number of individuals per family captured randomly
SUMMARY_COUNT_IndFam<-matrix(NA,ncol=5,nrow=4)
colnames(SUMMARY_COUNT_IndFam)<-c("min","median","mean","max","Fam>25")
row.names(SUMMARY_COUNT_IndFam)<-TotalNumber

for (t in 1:length(TotalNumber)){
  TOTALrandomSamples<- TotalNumber[t]
  
  #Variables to store output
  Family_counts<-NULL
  # individuals per family summary
  SUMM_IND_counts<-matrix(NA,ncol=5,nrow=100)
  colnames(SUMM_IND_counts)<-c("min","median","mean","max","Fam>25")
  pdf(paste("~/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forBoothstrapping/GWAS_bootstrap/Plots/hist_",TOTALrandomSamples,".pdf",sep=""),width=20,height=20)
  par(mfrow=c(10,10))
  # for i number of runs
  for (i in 1:100){
    a<-i
    
    #select subsamples regardless of family
    SamplesTotal <- SAMPLE_NAMES[ sample(1:length(SAMPLE_NAMES), TOTALrandomSamples, replace = F)]
    # how many families are there?
    SamplesTotal <-as.data.frame(SamplesTotal)
    Sam_fam <-  apply(SamplesTotal,1, function(x) strsplit(as.character(x), "S")[[1]][1])
    Family_counts<-c(Family_counts,length(names(table(Sam_fam))))
    
    # Summary statistics
    SUMM_IND_counts[i,1]<- min(as.data.frame(table(Sam_fam))[,2])
    SUMM_IND_counts[i,2]<-median(as.data.frame(table(Sam_fam))[,2])
    SUMM_IND_counts[i,3]<-mean(as.data.frame(table(Sam_fam))[,2])
    SUMM_IND_counts[i,4]<-max(as.data.frame(table(Sam_fam))[,2])
    # how many families had more than 25 individuals
    SUMM_IND_counts[i,5]<-length(which(as.data.frame(table(Sam_fam))[,2] >= 25))
    
    # PLOT
    hist(table(Sam_fam), main=as.character(TOTALrandomSamples))
  }
  dev.off()
  assign(paste("Family_counts_",TOTALrandomSamples,sep=""),Family_counts)
  
  # Summary across 100 runs
  SUMMARY_total<-apply(SUMM_IND_counts,2,mean)
  SUMMARY_COUNT_IndFam[t,]<-SUMMARY_total
}

# Number of families represented each ran
summary(Family_counts_264)
summary(Family_counts_528)
summary(Family_counts_968)
summary(Family_counts_2024)
