# Author: Ana Poets
# Description: Identify how many of the QTL detected with the whone NAM set, are detected
# during subsampling
##########################################################################################

rm(list=ls())

# Import QTL assigned from GWAS using the whole NAM set

#NAM<-read.csv("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/Alex_QTLassignation_3mbp/80miss_byfam_SNP_info_withQTL.csv",header=T,sep=",")
NAM<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/GenomeWide_QTLassigned_5e+06.xls",header=T)

# List of the number of QTL assigned
QTL_all<-names(table(NAM[,2]))

# Directory with boothstrapping summaries
#DIR_SUMMARY<-c("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forBoothstrapping/GWAS_bootstrap")
DIR_SUMMARY<-c("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap")

# 1. Boothstrapping by subsampling families
Subsam1<-c("OUTPUT_familiesAllInd",	"OUTPUT_familiesEvenInd")
  SAMPLING1<-c("Random_5",	"Random_10",	"Random_25",	"Random_50",	"Random_75")

## 2. Boothstrapping by random individuals
Subsam2<-c("OUTPUT_randomNoFam",	"OUTPUT_random_wFam")
  SAMPLING2<-c("Random_264",	"Random_528",	"Random_968",	"Random_2024")

## For each subsampling (100 runs) make a table of the number of times a QTL was identified

# For SUBSAMPLING
for (S in 1:2) {
  
  # For All or subset 
  for (x in 1:2){
    
    # Create a table to store the results of how many QTL were identified in each Subsampling strategy
    TABLE<-matrix(0,nrow=length(get(paste("SAMPLING",S, sep=""))),ncol=length(QTL_all))
    row.names(TABLE)<-get(paste("SAMPLING",S, sep=""))
    colnames(TABLE)<-QTL_all
  # for SAMPLING by fam or by individuals
    for (y in 1:length(get(paste("SAMPLING",S,sep="")))){
      for (r in 1:100){
      print(paste(DIR_SUMMARY,"/",get(paste("Subsam",S,sep=""))[x],"/",get(paste("SAMPLING",S,sep=""))[y],"/","PolyTest_RUN_",r,".txt",sep=""))
      
      RESULTS<-read.table(paste(DIR_SUMMARY,"/",get(paste("Subsam",S,sep=""))[x],"/",get(paste("SAMPLING",S,sep=""))[y],"/","PolyTest_RUN_",r,".txt",sep=""))
      # set threshold of significanse
      THR<- -log10(0.05 / ( dim(RESULTS)[1] * (1-0.05)))
      if (length(which(RESULTS$pval > THR)) >0){
        Significant_Results<-row.names(RESULTS[which(RESULTS$pval >THR),])
        
      # Check if SNPs in QTL for the whole NAM where captured. If at least one SNP in a QTL is significant in the boothstrapping
      # count the QTL as present.
       
      for (Q in 1:length(QTL_all)){
        # get SNP names
        SNPinQTL_NAM<-row.names(NAM[grep(QTL_all[Q],NAM[,2]),])
        #if there is at least one SNP present in QTL count the QTL.
        SNPsigSubset<-intersect(Significant_Results,SNPinQTL_NAM)
        if(length(SNPsigSubset)>0){
          TABLE[y,Q]<-TABLE[y,Q]+1
        }else{TABLE[y,Q]<-TABLE[y,Q]}
      }
      }
      }
      write.table(TABLE,paste("/home/smithkp/agonzale/Projects/NAM/Analysis/GWAS_bootstrap/TABLES/",paste(get(paste("Subsam",S,sep=""))[x],"_Table",sep=""),".xls",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
    assign(paste(get(paste("Subsam",S,sep=""))[x],"_Table",sep=""),TABLE)
    }
  }
}
