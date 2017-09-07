# Author: Ana Poets
# Description: Get a summary of significant QTL in each Boothstrapping
# 1) by random families (individuals in families)
# 2) by a subset of individuals (individuals per family)
###################################################################################
rm(list=ls())

# 1. Boothstrapping by subsampling families
DIR<-c("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesAllInd","/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_familiesEvenInd")
SUBDIR<-c("Random_10" , "Random_25" , "Random_5" , "Random_50" , "Random_75")

for (D in 1:2){
  for (S in 1:length(SUBDIR)){
    print(paste("D=",D, " S=", S))
    
# Empty matrix to store results
RESULTS_SIGNIF<-NULL
SNPs_significant<-NULL
for (i in 1:100){
RESULTS<-read.table(paste(DIR[D],"/",SUBDIR[S],"/PolyTest_RUN_",i,".txt",sep=""),header=T)

THR<- -log10(0.05 / ( dim(RESULTS)[1] * (1-0.05)))

#= Assign output to a variable and safe also to a file
if (length(which(RESULTS$pval > THR)) >0){
  Significant_Results<-RESULTS[which(RESULTS$pval >THR),]
  SNPs_significant<-c(SNPs_significant,row.names(Significant_Results))
  RESULTS_SIGNIF<-c(RESULTS_SIGNIF,length(which(RESULTS$pval > THR))) # change significance threshold
}
}
# Write table of the number of significant SNPs 
write.table(as.data.frame(RESULTS_SIGNIF), paste(DIR[D],"/",SUBDIR[S],"/Number_sig_SNPs.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
# The number of times that each SNP was found significant 
write.table(as.data.frame(table(SNPs_significant)), paste(DIR[D],"/",SUBDIR[S],"/Count_sig_per_SNP.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
  
  rm(RESULTS_SIGNIF,SNPs_significant)
  }
}

rm(DIR,SUBDIR)

## 2. RANDOM INDIVIDUALS
DIR<-c("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT_randomNoFam","/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/OUTPUT")
SUBDIR<-c("Random_2012" , "Random_264",  "Random_528" , "Random_968")

for (D in 1:2){
  for (S in 1:length(SUBDIR)){
    print(paste("D=",D, " S=", S))
    # Empty matrix to store results
    RESULTS_SIGNIF<-NULL
    SNPs_significant<-NULL
    for (i in 1:100){
      RESULTS<-read.table(paste(DIR[D],"/",SUBDIR[S],"/PolyTest_RUN_",i,".txt",sep=""),header=T)
      
      THR<- -log10(0.05 / ( dim(RESULTS)[1] * (1-0.05)))
      
      #= Assign output to a variable and safe also to a file
      if (length(which(RESULTS$pval > THR)) >0){
        Significant_Results<-RESULTS[which(RESULTS$pval >THR),]
        SNPs_significant<-c(SNPs_significant,row.names(Significant_Results))
        RESULTS_SIGNIF<-c(RESULTS_SIGNIF,length(which(RESULTS$pval > THR))) # change significance threshold
      }
    }
    # Write table of the number of significant SNPs 
    write.table(as.data.frame(RESULTS_SIGNIF), paste(DIR[D],"/",SUBDIR[S],"/Number_sig_SNPs.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
    # The number of times that each SNP was found significant 
    write.table(as.data.frame(table(SNPs_significant)), paste(DIR[D],"/",SUBDIR[S],"/Count_sig_per_SNP.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
    rm(RESULTS_SIGNIF,SNPs_significant)
    
     }
}