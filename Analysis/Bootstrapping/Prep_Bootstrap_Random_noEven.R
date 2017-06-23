# Author: Alex Ollhoff & Ana Poets
# Description: Bootstrap sampling mapping panels for the BRIDG6, Sampling random number of individuals from each family
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
#y <- read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Input/phenos_LDKNNI.csv", header=T)
y <- read.csv("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/Data/phenos_LDKNNI.csv", header=T)

# import genotypes
#gen_noRas <- read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Input/genos_LDKNNI.csv", header = T)
gen_noRas <- read.csv("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/Data/genos_LDKNNI.csv", header = T)

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

# Make file for MSI
#write.csv(adjusted_genotypes, "~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forBoothstrapping/geno_pheno/genos_MAF05.csv")
#write.csv(y_2, "~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forBoothstrapping/geno_pheno/phenos_MAF05.csv")

#write.csv(adjusted_genotypes, "/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/geno_pheno/genos_MAF05.csv")
#write.csv(y_2, "/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/geno_pheno/phenos_MAF05.csv")

### Prep loops for mapping below, then copy and edit in map_bootstrap script ###
# Store a list of samples
SAMPLE_NAMES<-y_2[,2]

############
# Import list of families
#Families<-read.table("~/Dropbox/SmithLab/NAM/Data/Alex/Pop_structure_propDiss_JULY2016_parentLessMissingData_grep.txt",header=T)
Families<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/Data/Pop_structure_propDiss_JULY2016_parentLessMissingData_grep.txt",header=T)

# Empty matrix to store results
RESULTS_SIGNIF<-NULL
SNPs_significant<-NULL
##### for bootstrapping even number of individuals per family
#TotalNumber<-c(264,528,968,2012)
TOTALrandomSamples<- 264

# Create a directory to put output files
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam/Random_",TOTALrandomSamples,sep=""))
# Create a directory for plots
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam/Random_",TOTALrandomSamples,"/Plot", sep=""))

# for i number of runs
for (i in 1:100){
  a<-i
 
  #select subsamples regardless of family
  SamplesTotal <- SAMPLE_NAMES[ sample(1:length(SAMPLE_NAMES), TOTALrandomSamples, replace = F)]

  # Filter for phenotypes present in the sample
  Sampled_phen<-y_2[(y_2[,2] %in% SamplesTotal),]
  
  # Filter for genotypes present in the sample
  Sampled_gen<-gen_2[(gen_2[,1] %in% Sampled_phen$line_name),]
  
  # Strip off first column that contains the sample names
  gen_naked <- Sampled_gen[,-1]
  
  # QC for MAF of 1/250 (0.00047)
  adjusted_genotypes <- snpQC(gen = gen_naked, MAF = 0.003, impute = FALSE)
  rownames(adjusted_genotypes) <- Sampled_phen$line_name
  head(Sampled_phen)
  dim(adjusted_genotypes)
  
  # Get chrom 
  chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
  class(chr) # should be integer
  
  # Need to write to 10 different csv... 
  my_gwas = gwas2(Sampled_phen$DAP_BLUPs,adjusted_genotypes,Sampled_phen$family,chr)
  number_of_markers = ncol(adjusted_genotypes)
  
  #= Saving plots ====
  #significant threshold
  #THR = -log10(0.05 / ( ncol(RESULTS) * (1-0.05)))
  # (log(number_of_markers*0.05))
  
  THR<- -log10(0.05 / ( dim(my_gwas$PolyTest)[1] * (1-0.05)))
  pdf(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam/Random_",TOTALrandomSamples,"/Plot/GWAS_",i,".pdf",sep=""),width=7,height=5)
  
  plot( my_gwas, FDR = 0.05)
  abline(h=THR,col="red")
  dev.off()
  
  #= Assign output to a variable and safe also to a file
  RESULTS<-my_gwas$PolyTest
  row.names(RESULTS)<-sub("X","",as.character(my_gwas$SNPs))
  if (length(which(RESULTS$pval > THR)) >0){
    Significant_Results<-RESULTS[which(RESULTS$pval >THR),]
    SNPs_significant<-c(SNPs_significant,row.names(Significant_Results))
    RESULTS_SIGNIF<-c(RESULTS_SIGNIF,length(which(RESULTS$pval > THR))) # change significance threshold
  }
  
  row.names(my_gwas$PolyTest)<-sub("X","",my_gwas$SNPs)
  assign(paste("NAM_",i,sep=""),my_gwas$PolyTest)
  write.table(my_gwas$PolyTest, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam/Random_",TOTALrandomSamples,"/PolyTest_RUN_",i,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
  write.table(my_gwas$SNPs, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam/Random_",TOTALrandomSamples,"/SNPs_RUN_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  
}

# Write table of the number of significant SNPs 
write.table(as.data.frame(RESULTS_SIGNIF), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam/Random_",TOTALrandomSamples,"/Number_sig_SNPs.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
# The number of times that each SNP was found significant 
write.table(as.data.frame(table(SNPs_significant)), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam/Random_",TOTALrandomSamples,"/Count_sig_per_SNP.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
