# Author: Ana Poets
# Description: Bootstrap sampling mapping panels for the BRIDG6, sampling X number of families. 
# 1) Select even number of individuals per family.
# 2) Select all the individuals for the same families
###################################################################################################################################################
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


##############

##### for bootstrapping even number of individuals per family, where family number varies.
#TotalNumber<-c(5,10,25,50,75)
for (t in 1:length(TotalNumber)){
TOTALrandomFamilies<- TotalNumber[t]

# Create a directory to put output files from even number of individuals per family
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd/Random_",TOTALrandomFamilies,sep=""))
# Create a directory for plots
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd/Random_",TOTALrandomFamilies,"/Plot", sep=""))

# Create a directory to put files from all individuals in a subset of families
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd/Random_",TOTALrandomFamilies,sep=""))
# Create a directory for plots
dir.create(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd/Random_",TOTALrandomFamilies,"/Plot", sep=""))


# matrix of how many samples to sample per family in each escenario
# add a column to indicate the number of individuals we want to subsample from each families, if the family
# has less than 25 individuals then take all the individuals, otherwise take 25 individuals.
Count_ind_fam_Sampling<-cbind(Count_ind_fam,rep(NA,dim(Count_ind_fam)[1]))
Count_ind_fam_Sampling[which(as.numeric(as.character(Count_ind_fam_Sampling[,2])) <=20),3]<-as.numeric(as.character(Count_ind_fam_Sampling[which(as.numeric(as.character(Count_ind_fam_Sampling[,2])) <=20),2]))
Count_ind_fam_Sampling[which(as.numeric(as.character(Count_ind_fam_Sampling[,2])) >20),3]<-25
colnames(Count_ind_fam_Sampling)<-c("Family_grep","NumberInd","SampleSize")
Count_ind_fam_Sampling<-as.data.frame(Count_ind_fam_Sampling)


# Empty matrix to store results for even number of individuals per family
RESULTS_SIGNIF<-NULL
SNPs_significant<-NULL

# Empty matrix to store results for ALL number of individuals per family
RESULTS_SIGNIF_all<-NULL
SNPs_significant_all<-NULL

# for i number of runs
for (i in 1:100){
  a<-i
  SamplesTotal<-NULL
  SamplesTotal_all<-NULL
    #Select families to work with
    FamiliesSelected <- Count_ind_fam_Sampling[ sample(1:dim(Count_ind_fam_Sampling)[1],TOTALrandomFamilies , replace = F),1]
    
    ##### 1) select an even set of individuals from each family
      for (x in 1:length(FamiliesSelected)){
        #get X individuals per families selected
        Fam_all<-SAMPLE_NAMES[grep(FamiliesSelected[x],SAMPLE_NAMES)]
        # select even number of individuals per family
        IndFam<-Fam_all[sample(c(1:length(Fam_all)), as.numeric(as.character(Count_ind_fam_Sampling[grep(FamiliesSelected[x],Count_ind_fam_Sampling[,1]),3])), replace=F )]
        SamplesTotal<-c(SamplesTotal,as.character(IndFam))
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
    #head(Sampled_phen)
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
    pdf(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd/Random_",TOTALrandomFamilies,"/Plot/GWAS_",i,".pdf",sep=""),width=7,height=5)
      plot( my_gwas)
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
    write.table(my_gwas$PolyTest, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd/Random_",TOTALrandomFamilies,"/PolyTest_RUN_",i,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
    write.table(my_gwas$SNPs, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd/Random_",TOTALrandomFamilies,"/SNPs_RUN_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
 
    #############################################################################
    ##### 2) select an even set of individuals from each family
    for (y in 1:length(FamiliesSelected)){
      #get all individuals per families selected
      Fam_all_ind<-SAMPLE_NAMES[grep(FamiliesSelected[y],SAMPLE_NAMES)]
      # select all the individuals per family
      SamplesTotal_all<-c(SamplesTotal_all,as.character(Fam_all_ind))
    }
    
    
    # Filter for phenotypes present in the sample
    Sampled_phen_all<-y_2[(y_2[,2] %in% SamplesTotal_all),]
    
    # Filter for genotypes present in the sample
    Sampled_gen_all<-gen_2[(gen_2[,1] %in% Sampled_phen_all$line_name),]
    
    # Strip off first column that contains the sample names
    gen_naked_all <- Sampled_gen_all[,-1]
    
    # QC for MAF of 1/250 (0.00047)
    adjusted_genotypes_all <- snpQC(gen = gen_naked_all, MAF = 0.003, impute = FALSE)
    rownames(adjusted_genotypes_all) <- Sampled_phen_all$line_name
    #head(Sampled_phen_all)
    dim(adjusted_genotypes_all)
    
    # Get chrom 
    chr_all = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes_all))))[,2]
    class(chr_all) # should be integer
    
    # Need to write to 10 different csv... 
    my_gwas_all = gwas2(Sampled_phen_all$DAP_BLUPs,adjusted_genotypes_all,Sampled_phen_all$family,chr_all)
    number_of_markers_all = ncol(adjusted_genotypes_all)
    
    #= Saving plots ====
    #significant threshold
    #THR = -log10(0.05 / ( ncol(RESULTS) * (1-0.05)))
    # (log(number_of_markers*0.05))
    
    THR_all<- -log10(0.05 / ( dim(my_gwas_all$PolyTest)[1] * (1-0.05)))
    pdf(paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd/Random_",TOTALrandomFamilies,"/Plot/GWAS_",i,".pdf",sep=""),width=7,height=5)
    plot( my_gwas_all)
    abline(h=THR_all,col="red")
    dev.off()
    
    #= Assign output to a variable and safe also to a file
    RESULTS_all<-my_gwas_all$PolyTest
    row.names(RESULTS_all)<-sub("X","",as.character(my_gwas_all$SNPs))
    if (length(which(RESULTS_all$pval > THR_all)) >0){
      Significant_Results_all<-RESULTS_all[which(RESULTS_all$pval >THR_all),]
      SNPs_significant_all<-c(SNPs_significant_all,row.names(Significant_Results_all))
      RESULTS_SIGNIF_all<-c(RESULTS_SIGNIF_all,length(which(RESULTS_all$pval > THR_all))) # change significance threshold
    }
    
    row.names(my_gwas_all$PolyTest)<-sub("X","",my_gwas_all$SNPs)
    assign(paste("NAM_all_",i,sep=""),my_gwas_all$PolyTest)
    write.table(my_gwas_all$PolyTest, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd/Random_",TOTALrandomFamilies,"/PolyTest_RUN_",i,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
    write.table(my_gwas_all$SNPs, paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd/Random_",TOTALrandomFamilies,"/SNPs_RUN_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  }

# Write table of the number of significant SNPs 
write.table(as.data.frame(RESULTS_SIGNIF), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd/Random_",TOTALrandomFamilies,"/Number_sig_SNPs.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
# The number of times that each SNP was found significant 
write.table(as.data.frame(table(SNPs_significant)), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd/Random_",TOTALrandomFamilies,"/Count_sig_per_SNP.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

# Write table of the number of significant SNPs 
write.table(as.data.frame(RESULTS_SIGNIF_all), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd/Random_",TOTALrandomFamilies,"/Number_sig_SNPs.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
# The number of times that each SNP was found significant 
write.table(as.data.frame(table(SNPs_significant_all)), paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd/Random_",TOTALrandomFamilies,"/Count_sig_per_SNP.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
}