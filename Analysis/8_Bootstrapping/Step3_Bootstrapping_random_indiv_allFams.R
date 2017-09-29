# Author: Alex Ollhoff & Ana Poets
# Description: Bootstrap sampling mapping panels for the BRIDG6
#     Sampling even number of individuals from each family (264,528,968,2012)
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
#gen_noLD <- read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Input/genos_LDKNNI.csv", header = T)
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

############
# Import list of families
Families<-read.table("~/GITHUB/BRIDG6/Datasets/Intermidiate_dataset/Pop_structure_propDiss_JULY2016_parentLessMissingData_grep.txt",header=T)

#Count how many individuals per family we have
Count_ind_fam<-matrix(NA,ncol=2,nrow=dim(Families)[1])
for(x in 1:dim(Families)[1]){
  Fam<-grep(Families[x,"Family_grep"],SAMPLE_NAMES)
  Count_ind_fam[x,1]<-as.character(Families[x,"Family_grep"])
  Count_ind_fam[x,2]<-length(Fam)
}
Count_ind_fam<-as.data.frame(Count_ind_fam)
# separate the families by those that have more than 25 individuals and those with less than 25 individuals
# if the family has <25 individuals, all the individuals in the family will be used every time

FamLess25<-Count_ind_fam[which(as.numeric(as.character(Count_ind_fam[,2])) <= 25),]
FamMore25<-Count_ind_fam[which(as.numeric(as.character(Count_ind_fam[,2])) > 25),]

##############
# Empty matrix to store results
RESULTS_SIGNIF<-NULL
SNPs_significant<-NULL
##### for bootstrapping even number of individuals per family
TotalNumber<-c(264,528,968,2024)
for ( t in 1:length(TotalNumber)){
TOTALrandomSamples<- TotalNumber[t]

# Create a directory to put output files
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT/Random_",TOTALrandomSamples,sep=""))
# Create a directory for plots
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT/Random_",TOTALrandomSamples,"/Plot", sep=""))
# matrix of how many samples to sample per family in each escenario
SUB_Sampling<-data.frame(c(264,528,968,2024),c(3,6,11,23))
colnames(SUB_Sampling)<-c("TotalSize","IndPerFam")
# for i number of runs
for (i in 1:100){
  a<-i
  SamplesTotal<-NULL
  
  # if total random samples <= 968
  if (TOTALrandomSamples <=968){
    for (h in 1:dim(Families)[1]){
      Fam_all<-SAMPLE_NAMES[grep(Families[h,"Family_grep"],SAMPLE_NAMES)]
      #select subsamples per family
      Subsample <- Fam_all[ sample(1:length(Fam_all), SUB_Sampling[which(SUB_Sampling[,1] ==TOTALrandomSamples ),2] , replace = F)]
      SamplesTotal<-c(SamplesTotal,as.character(Subsample))
    }
  }
  
  # if total random samples > 968 then for family HR639S use all 11 samples , while for other families get 23 samples
  if (TOTALrandomSamples >968){
    
    for (m in 1:dim(FamLess25)[1]){
      SamplesLess25<-SAMPLE_NAMES[grep(FamLess25[m,1], SAMPLE_NAMES)]
      SamplesTotal<-c(SamplesTotal,as.character(SamplesLess25))
    }
    
    for (n in 1:dim(FamMore25)[1]){
      SamplesMore25<-SAMPLE_NAMES[grep(FamMore25[n,1], SAMPLE_NAMES)]
      # Select X individuals randomly from a family 
      Subsample20 <- SamplesMore25[ sample(1:length(SamplesMore25), SUB_Sampling[which(SUB_Sampling[,1] ==TOTALrandomSamples ),2], replace = F)]
      SamplesTotal<-c(SamplesTotal,as.character(Subsample20))
    }
  }  
  
  
  
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
  pdf(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT/Random_",TOTALrandomSamples,"/Plot/GWAS_",i,".pdf",sep=""),width=7,height=5)
  
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
  write.table(my_gwas$PolyTest, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT/Random_",TOTALrandomSamples,"/PolyTest_RUN_",i,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
  write.table(my_gwas$SNPs, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT/Random_",TOTALrandomSamples,"/SNPs_RUN_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  
}

# Write table of the number of significant SNPs 
write.table(as.data.frame(RESULTS_SIGNIF), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT/Random_",TOTALrandomSamples,"/Number_sig_SNPs.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
# The number of times that each SNP was found significant 
write.table(as.data.frame(table(SNPs_significant)), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT/Random_",TOTALrandomSamples,"/Count_sig_per_SNP.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
}