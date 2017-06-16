# NAM of bootstrap samples of different sizes in MSI
# Alex Ollhoff
# February 1, 2017

# Map across all families using environment BLUPs and genotypes filtered by family with up to 80% missingness
# Alex Ollhoff
# February 1, 2017

# pbs file 
# vim NAM_mapping_bootstrap250_80missR.pbs

#!/bin/bash
#PBS    -l      walltime=80:00:00,nodes=2:ppn=2,mem=990gb       
#PBS    -o      /home/smithkp/ollho007/NAM/mapping/bootstrap250_80miss/mapping_output
#PBS    -e      /home/smithkp/ollho007/NAM/mapping/bootstrap250_80miss/mapping_error_output
#PBS    -N      Bootstrap250_80miss
#PBS    -M      ollho007@umn.edu
#PBS    -m      abe 
#PBS    -r      n
#PBS    -q      ram1t
cd /home/smithkp/ollho007/NAM/mapping/bootstrap250_80miss
module load R/3.3.1
R CMD BATCH --slave /home/smithkp/ollho007/NAM/mapping/bootstrap250_80miss/NAM_mapping_bootstrap250_80miss.R

##########################################################################
# vim NAM_mapping_bootstrap250_80miss.R

# Map across all families using environment BLUPs and genotypes filtered by family with up to 80% missingness
# Alex Ollhoff
# February 1, 2017

library(NAM)

# Load genotypic data
# ref allele = 2 = major allele
# 0 = minor allele, this allows minor alleles to have different effects if stratification is provided
gen = read.csv("~/NAM/mapping/bootstrap250_80miss/genos_80miss_forbootstrap.csv", header = T)
gen = data.frame(gen)

# Load phenotypic and family data
# phe1 = pheno in col 3, phe2 = fam in col 2
y = read.csv("~/NAM/mapping/bootstrap250_80miss/phenos_80miss_forbootstrap.csv", header = TRUE)[,-1]
y = data.frame(y)

# Make a file to store Sample statistics
SAMPLE_NAMES<-c(1:88)

# Empty matrix to store results
RESULTS_SIGNIF<-NULL
SNPs_significant<-NULL
# for b bootstrappings do:
for (i in 1:10){
  # Select 20 individuals randomly from the population
  Sample200 <- SAMPLE_NAMES[ sample(1:length(SAMPLE_NAMES), 5, replace = F)]
  
  # Filter for phenotypes present in the sample
  Sampled_phen<-y_2[(y_2[,3] %in% Sample200),]
  
  # Filter for genotypes present in the sample
  Sampled_gen<-gen_2[(gen_2[,1] %in% Sampled_phen$line_name),]
  
  # Strip off first column 
  gen_naked <- Sampled_gen[,-1]
  
  # QC for MAF of 1/250 (0.00047)
  adjusted_genotypes <- snpQC(gen = gen_naked, MAF = 0.00047, impute = FALSE)
  rownames(adjusted_genotypes) <- Sampled_phen$line_name
  head(Sampled_phen)
  dim(adjusted_genotypes)
  
  # Get chrom 
  chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
  class(chr) # should be integer
  
  # Need to write to 10 different csv... 
  my_gwas = gwas2(Sampled_phen$BLUE,adjusted_genotypes,Sampled_phen$family,chr)
  number_of_markers = nrow(adjusted_genotypes)
  
  #= Saving plots ====
  pdf(paste("~/NAM/mapping/bootstrap250_80miss/OUTPUT/GWAS_",i,".pdf",sep=""),width=7,height=5)
  plot( my_gwas, FDR = 0.05)
  dev.off()

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
  
}

# Write table of the number of significant SNPs 
write.table(as.data.frame(RESULTS_SIGNIF), "~/NAM/mapping/bootstrap250_80miss/OUTPUT/Number_sig_SNPs.txt",quote=F,row.names=F,col.names=F,sep="\t")

# The number of times that each SNP was found significant 
write.table(as.data.frame(table(SNPs_significant)), "~/NAM/mapping/bootstrap250_80miss/OUTPUT/Count_sig_per_SNP.txt",quote=F,row.names=F,col.names=F,sep="\t")
