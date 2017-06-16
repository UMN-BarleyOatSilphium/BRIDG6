# Plot PCA of NAM donors and NSGC 
# Plot PCA of NAM donors and RILs
# Calcuate Fst between NSGC and NAM donors

# install.packages("adegenet")
# install.packages("hierfstat")

library(dplyr)
library(adegenet)
library(ggplot2)
library(hierfstat)

###### calculate Fst between NAM donor parents and 6 row NSGC accessions

# Get NSGC genos and NAM donor parent genos
NSGC6_markers <- read.table("~/Documents/PhD/NAM/Parents/NSGC_6_AA_BB_genos/NSGC_6_genotype.hmp.txt", header=T, na.strings = "NA", row.names = 1)
NSGC6_markers_clean <- NSGC6_markers[,-c(1:3)]

NAM_markers <- read.table("~/Documents/PhD/NAM/Parents/NAM_donor_AA_BB_genos/NAM_donor_genotype.hmp.txt", header=T, na.strings="NA", row.names = 1)
NAM_markers_clean <- NAM_markers[,-c(1:3)]

#NSGC6_NAM_markers <- merge(NSGC6_markers_clean, NAM_markers_clean, by.x = "rs", by.y = "rs", nomatch = NA, all.x = F)
#as.data.frame(NSGC6_NAM_markers)

# filter NSGC markers for markers present in the NAM donors
NSGC6_matched <- NSGC6_markers_clean[match(NAM_markers_clean$rs, NSGC6_markers_clean$rs, nomatch = NA, incomparables = F),]  
NAM_matched <- NAM_markers_clean[match(NSGC6_matched$rs, NAM_markers_clean$rs, nomatch = NA, incomparables = F),]  


# Recode genotypes for use in Hierfstat
# AA/1 = 1
# BB/-1 = 2
# AB/0 = NA
# missing/NA = NA
# Set Rassmusson allele as "2",  Donor Parent as "0" and heterozygotes = "NA".
Fst_based<-function(dat){
  dat<-(dat)
  #AA
  dat[which(dat == "1")] <-"1"
  #BB
  dat[which(dat == "-1")] <-"2"
  #Hetes
  dat[which(dat == "0")] <- "NA"
  return(dat)
}

NSGC6_markers_new<-NSGC6_matched[,-1]
Fst_based_NSGC6<-(apply(NSGC6_markers_new, 1, Fst_based))
Fst_based_NSGC6 <- as.data.frame(Fst_based_NSGC6)

NAM_markers_new<-NAM_matched[,-1]
Fst_based_NAM<-(apply(NAM_markers_new, 1, Fst_based))
Fst_based_NAM <- as.data.frame(Fst_based_NAM)

write.csv(Fst_based_NAM, "~/Documents/PhD/NAM/Population_description/Fst_autocorrelation/NAM_markers.csv")
Fst_based_NAM <- read.csv("~/Documents/PhD/NAM/Population_description/Fst_autocorrelation/NAM_markers.csv", header = T)


## ==== ================Applying QC filterings =======================================================

#Thresholds
#Hetes=0.20
Missingness=0.97 # No partition can be missing any marker
#Missing_samples=0.5
Maf=0.001 # at least 1 individual have to have the minor allele

#==== Functions to remove Missingness and LOW  MAF SNPs =======================

MISSING<-function(dat){
  missing<-length(dat[is.na(dat)])/(length(dat))
  return(missing)
}

MAF<-function(dat){
  missing<-length(dat[is.na(dat)])
  
  #Count the alleles
  AllelesPresent<-table(dat)
  if (length(AllelesPresent) == 1){maf<-0}else{
    hete<-which(names(AllelesPresent) == "NA")
    
    Homo<-which(names(AllelesPresent) != "NA")
    
    #look for alleles different from Hete
    if (length(hete) >0){
      hete_count<-AllelesPresent[hete]			
    }else{hete_count<-0}
    
    if (length(Homo) ==1){AlleleA<-AllelesPresent[Homo]; AlleleB<-0}
    if(length(Homo) == 2){AlleleA<-AllelesPresent[Homo[1]]; AlleleB<-AllelesPresent[Homo[2]]}
    MAJOR<-((max(AlleleA,AlleleB))*2) + hete_count[1]
    MINOR<-((min(AlleleA,AlleleB))*2) + hete_count[1]
    
    maf<-(MINOR)/(2*(length(dat)-missing))
  }
  return(maf)
}


# Remove SNPs with more than 10% missingness in NSGC
miss_count<-apply(Fst_based_NSGC6,1, MISSING)

RILs_ras_missFiltered<-Fst_based_NSGC6[-c(which(miss_count > Missingness )),] # 24768SNP  6161
RILS_filtered_missing <- as.data.frame(RILs_ras_missFiltered)

# Make sure SNPs are still segreagting in the NSGC
maf_count<-apply(RILS_filtered_missing, 1, MAF)
#RILs_ras_miss_mafFiltered<-RILS_filtered_missing[-c(which(maf_count < Maf )),] # 6825SNPs 6161samples
summary(maf_count)

#write.csv(RILS_filtered_missing, "~/Documents/PhD/NAM/Population_description/Fst_autocorrelation/NSGC6_markers.csv")
NSGC_markers <- read.csv("~/Documents/PhD/NAM/Population_description/Fst_autocorrelation/NSGC6_markers.csv", header = T)
colnames(NSGC_markers) <- c("X", 1:6648)

# NAM donors all have more stringent filtering for both individuals and markers, no need to do step above

# Randomly sample 90 individuals with replacement 1000 times from the NSGC genos 
rownames(Fst_based_NAM) <- c()
colnames(Fst_based_NAM) <- c("X", 1:6648)
NAM_for_Fst <- Fst_based_NAM[,-1]

NAM_for_Fst_noNA <- NAM_for_Fst[rowSums(is.na(NAM_for_Fst)) != ncol(NAM_for_Fst),]

# Calculate Fst for each sample with the NAM donors
# Make a file to store Sample statistics
SAMPLE_NAMES<-as.character(NSGC_markers$X)
########################################################

# Function to count percentage of missing values
MISSING<-function(dat){
  missing<-length(which(is.na(dat)))/length(dat)
  return(missing)
}

# Get NSGC genos and NAM donor parent genos
NSGC6_markers <- read.table("~/Documents/PhD/NAM/Parents/NSGC_6_AA_BB_genos/NSGC_6_genotype.hmp.txt", header=T, na.strings = "NA")
NSGC6_markers_clean <- NSGC6_markers[,-c(2:4)]

NAM_markers <- read.table("~/Documents/PhD/NAM/Parents/NAM_donor_AA_BB_genos/NAM_donor_genotype.hmp.txt", header=T, na.strings="NA")
NAM_markers_clean <- NAM_markers[,-c(2:4)]
# Make a file to store Sample statistics
SAMPLE_NAMES<-as.character(colnames(NSGC6_markers_clean)[-1])
SNP_names<-NSGC6_markers_clean[,1]
# Empty matrix to store results
Fst_values<-NULL
Average_Fst<-NULL
# for b bootstrappings do:

# matrix to retain all Fst values for the X reps
REP<-2
#Fst_inReps<-matrix(NA, nrow=(dim(NSGC6_markers_clean)[2] -1), ncol=REP)
for (i in 1:REP){
  x<-i
  # Select 20 individuals randomly from the population
  Sample_NSGC <- SAMPLE_NAMES[sample(1:length(SAMPLE_NAMES), 90, replace = F)]
  
  # Filter for genotypes present in the sample
  Sample_NSGC_genos<-NSGC6_markers_clean[,(colnames(NSGC6_markers_clean) %in% Sample_NSGC)]
  
  # Name columns by number
  #colnames(Sample_NSGC_genos) <- c("X", 1:6648)
  
  # Calculate percentage of missing values per SNP
  Sample_NSGC_genos_missing<-apply(Sample_NSGC_genos,1,MISSING)
  # Filter for markers with all NA
  if (length(which(Sample_NSGC_genos_missing == 1)) >0){
    Sample_NSGC_genos_noNA<-Sample_NSGC_genos[-c(which(Sample_NSGC_genos_missing == 1)),]
    SNPnames_noNA<-SNP_names[-which(Sample_NSGC_genos_missing == 1)]
    Sample_NSGC_genos_noNA_names<-cbind(SNP_names,Sample_NSGC_genos_noNA)
  }else{
    Sample_NSGC_genos_noNA<-Sample_NSGC_genos
    SNPnames_noNA<-SNP_names
    Sample_NSGC_genos_noNA_names<-cbind(SNP_names,Sample_NSGC_genos_noNA)
  }
  
  
  # Add NAM donor genos and NSGC sample genos
  common_row <- intersect(NAM_markers_clean[,1], Sample_NSGC_genos_noNA_names[,1])
  
  # get only common markers from the NSGC
  NAM_genotypes<-NAM_markers_clean[(NAM_markers_clean[,1] %in% as.character(common_row)),]
  NSGC_genotypes<-Sample_NSGC_genos_noNA_names[(Sample_NSGC_genos_noNA_names[,1] %in% as.character(common_row)),]
  
  #sort NSGC genotypes according to NAM parents SNP's names
  NSGC_genotypes_or<-NSGC_genotypes[match(NAM_genotypes[,1],NSGC_genotypes[,1]),]
  
  if (identical(as.character(NSGC_genotypes_or[,1]),as.character(NAM_genotypes[,1])) == FALSE)stop("Error! SNP names in NAM parents and NSGC are not in the same order")
  Sample_combined <- cbind(NAM_genotypes,NSGC_genotypes_or[,-1])
  
  gen<-t(Sample_combined[,-1])
  
 
  # Create levels for NAM parents and NSGC sample
  levels = c(rep("1", 90), rep("2", 90))
  
  # Calculate Fst between NAM parents and each NSGC sample
  Fst<-NULL
  
  for (j in 1:dim(gen)[2]) {
    
    Fst[j]<-varcomp(data.frame(levels,as.numeric(gen[,j])),diploid=FALSE)$F[1,1]   
    
  }

  # Store Fst value for each SNP in a dataframe
  Fst_results<-as.data.frame(Fst)
  
  # Combine dataframes based on the SNPs present in each sample
  rownames(Fst_results)<-Sample_combined[,1]
  colnames(Fst_results)<-paste("Fst_",i,sep="")
  assign(paste("Fst_run_",i,sep=""),Fst_results)

}

Fst_results_all <- cbind("Fst_run_",i,"")


write.csv(Fst_results_all, "~/Documents/PhD/NAM/Population_description/Fst_autocorrelation/Fst_NAM_NSGC_",i,".csv", sep=""), quote=F,row.names=F,col.names=T, sep = "\t")

  # Take mean of each comparison
  
write.table(as.data.frame(RESULTS_SIGNIF), "~/NAM/mapping/10_fams/OUTPUT/Number_sig_SNPs.txt",quote=F,row.names=F,col.names=F,sep="\t")

  # Assign output to a variable and then save to a file
  RESULTS<-my_gwas$PolyTest
  if (length(which(RESULTS$pval >(log(number_of_markers*0.05)))) >0){
    Significant_Results<-RESULTS[which(RESULTS$pval >(log(number_of_markers*0.05))),]
    SNPs_significant<-c(SNPs_significant,row.names(Significant_Results))
    RESULTS_SIGNIF<-c(RESULTS_SIGNIF,length(which(RESULTS$pval > (log(number_of_markers*0.05))))) # significance threshold is function of the number of markers and FDR
    
  }
  assign(paste("~/NAM/mapping/bootstrap250_80miss/OUTPUT/NAM_",i,sep=""),my_gwas$PolyTest)
  write.table(my_gwas$PolyTest, paste("~/NAM/mapping/bootstrap250_80miss/OUTPUT/NAM_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  write.table(my_gwas$SNPs, paste("~/NAM/mapping/bootstrap250_80miss/OUTPUT/NAM_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  



