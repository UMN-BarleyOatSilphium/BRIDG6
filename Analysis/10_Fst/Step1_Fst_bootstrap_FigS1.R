# Author: Ana Poets
# Description: Estimate Fst between NSGC and Donor Parents using 100 bootstraps of 88 samples in each partition 
#############################################################################################################


rm(list=ls())

library(hierfstat)

### Functions
# Change genotypes to 1 allele A,2 allele B, NA for hete 
FstBased<-function(dat){
  dat<-(dat)
  #AA
  dat[which(dat == "1")] <-"1"
  #BB
  dat[which(dat == "-1")] <-"2"
  #Hetes
  dat[which(dat == "0")] <- "NA"
  return(dat)
}

# identify amount of missing values
MISSING<-function(dat){
  missing<-length(dat[is.na(dat)])/(length(dat))
  return(missing)
}

# calculate MAF
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
### Analysis
# matrix ,rows are SNPs, columns are samples
Donor<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/Data_SFS/NAM_donor_genotype.hmp.txt",header=T,row.names=1)
nsgc<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/Data_SFS/NSGC_6_genotype.hmp.txt",header=T,row.names=1)

# get markers shared
SNPshared<-intersect(row.names(Donor),row.names(nsgc))

donor_sh<-Donor[(row.names(Donor) %in% SNPshared), -c(1:3)]
nsgc_sh<-nsgc[(row.names(nsgc) %in% SNPshared), -c(1:3)]

for (i in 1:100){
  print(i)
  # Select 90 individuals from each partition
  Donor_sampleNames<-colnames(donor_sh)
  nsgc_sh_sampleNames<-sample(colnames(nsgc_sh), 90, replace=F)
  
  #get genotypes for only those samples
  donor_sh_forFST<-donor_sh[,(colnames(donor_sh) %in% Donor_sampleNames)]
  nsgc_sh_forFST<-nsgc_sh[,(colnames(nsgc_sh) %in% nsgc_sh_sampleNames)]
  
  # Change genotypes 
  donor_sh_forFST_gen<-apply(donor_sh_forFST,2,FstBased)
  nsgc_sh_forFST_gen<-apply(nsgc_sh_forFST,2,FstBased)
  
  # find markers missing 100% in either partition
  donor_miss<-apply(donor_sh_forFST_gen,1,MISSING)
  nsgc_miss<-apply(nsgc_sh_forFST_gen,1,MISSING)
  
  # make a list of all markers that are missing in either population and remove them from both populations
  SNPsRemove<-NULL
  if(length(which(donor_miss == 1)) >0 | length(which(nsgc_miss == 1))>0){
    if(length(which(donor_miss == 1)) > 0){
      SNPsRemove<-c(SNPsRemove,row.names(donor_sh_forFST_gen)[which(donor_miss == 1)])
    }
    if(length(which(donor_miss == 1)) > 0){
      SNPsRemove<-c(SNPsRemove,row.names(nsgc_sh_forFST_gen)[which(nsgc_miss == 1)])
    }
  }
  if(length(SNPsRemove)>0){
    Donor_complete<-donor_sh_forFST_gen[!(row.names(donor_sh_forFST_gen) %in%SNPsRemove),]
    nsgc_complete<-nsgc_sh_forFST_gen[!(row.names(nsgc_sh_forFST_gen) %in%SNPsRemove),]
  }else{Donor_complete<-donor_sh_forFST_gen;nsgc_complete<-nsgc_sh_forFST_gen}
  
  # turn tables to have SNPs in columns and samples in rows
  t_donor<-t(Donor_complete)
  t_nsgc<-t(nsgc_complete)
  
  # combine tables and remove SNPs that are monomorphic
  Data_all<-rbind(t_donor,t_nsgc)
  
  # find and remove monomorphic SNPs
  
  Mono_snp<-apply(Data_all,2,MAF)
  
  if (length(which(Mono_snp >0.99)) >0){
    Data_all_ready<-Data_all[,-c(which(Mono_snp >0.99))]
  }else{Data_all_ready<-Data_all}

  ## Calculate Fst
  gen<-Data_all_ready
  
  
  # Create levels for NAM parents (donor on Top) and NSGC sample
  levels = c(rep("1", 90), rep("2", 90))
  
  # Calculate Fst between NAM parents and each NSGC sample
  Fst<-NULL
  
  for (j in 1:dim(gen)[2]) {
    
    Fst[j]<-varcomp(data.frame(levels,as.numeric(gen[,j])),diploid=FALSE)$F[1,1]   
    
  }
  #add SNP name
  FstResults<-as.data.frame(Fst)
  row.names(FstResults)<-colnames(gen)
  names(Fst_results)[1]<-paste("Fst_",i,sep="")
  assign(paste("Fst_",i,sep=""),FstResults)
  write.table(FstResults,paste("~/Desktop/Fst_temp/Fst_",i,".txt",sep=""),quote=F,row.names=T,col.names=T,sep="\t")
}

#Combine the results from the 100 funs
RESULTS<-read.table(paste("~/Desktop/Fst_temp/Fst_",i,".txt",sep=""),header=T,row.names=1)
for (i in 2:100){
  MY_TABLE<-read.table(paste("~/Desktop/Fst_temp/Fst_",i,".txt",sep=""),header=T,row.names=1)
  RESULTS<-cbind(RESULTS,MY_TABLE)
}

#Get mean Fst
MEAN_Fst<-as.data.frame(apply(RESULTS,1,mean))
MEAN_Fst<-cbind(row.names(RESULTS),MEAN_Fst)

# Sort SNPs by physical position
SNPpositions<-read.csv("~/Dropbox/GITHUB/BRIDG6/Datasets/NSGC/iSelect/iSelect_physicalPosition.txt",skip=7,sep="\t")
head(SNPpositions)

# Find positions for SNPs in Fst comparison
SNPpos_sh<-intersect(SNPpositions$ID, MEAN_Fst[,1]) # 328 SNPs don't have known physical position

SNPpositions_sh<-SNPpositions[(SNPpositions$ID %in% SNPpos_sh),]
MEAN_Fst_sh<-MEAN_Fst[((MEAN_Fst[,1]) %in% SNPpos_sh),]

SNPpos_sh_or<-SNPpositions_sh[match(MEAN_Fst_sh[,1],SNPpositions_sh$ID),]

if (identical(as.character(SNPpos_sh_or$ID), as.character(MEAN_Fst_sh[,1])) == FALSE)stop("SNPs are in different order")

FST_positions<-cbind(SNPpos_sh_or[,1:3],MEAN_Fst_sh)

# Sort SNPs by chromosome and by position
FST_positions_or<-FST_positions[order(FST_positions$X.CHROM , FST_positions$POS),]

# Get cumulative positions genome-wide
Physical_positions<-FST_positions_or[,1:2]

#total bp prior to each chromosome
ADD_CHR<-c(0,558535432,1326610456,2026321570,2673381728,3343411888,3926792401,4584016401)

for (i in 1:8){
  if(i <8){
  CHR<-Physical_positions [grep(paste("chr",i, sep=""),Physical_positions[,1]),]
  CHR_corrected<-as.numeric(as.character(CHR[,2])) + ADD_CHR[i]
  }
  if(i==8){
    CHR<-Physical_positions [grep("chrUn",Physical_positions[,1]),]
    CHR_corrected<-as.numeric(as.character(CHR[,2])) + ADD_CHR[i]
  }
  assign(paste("CHR_cor_",i,sep=""), CHR_corrected)
}

NEW_positions<-c(CHR_cor_1, CHR_cor_2, CHR_cor_3, CHR_cor_4, CHR_cor_5, CHR_cor_6, CHR_cor_7,CHR_cor_8)

Fst_ordered<-cbind(as.character(FST_positions_or$X.CHROM),as.data.frame(NEW_positions), as.data.frame(FST_positions_or[,5]))
Fst_ordered<-as.data.frame(Fst_ordered)
colnames(Fst_ordered)<-c("CHR","Cumul_pos","Mean_Fst")

### ========== PLOT mean FST accross 100 runs ===============================================

pdf("~/Documents/SmithLab/NAM/Analysis/Fst/Plot/Fst_mean_100_NSGC_DP.pdf",width=8,height=6)
plot(Fst_ordered$Cumul_pos/1000000, Fst_ordered$Mean_Fst, ylab="Fst", xlab="Physical Position (Mbp)",cex=0.6)
# color each chromosome
COLOR<-c("red","blue","green","magenta","orange","cyan","light green")
for(x in 1:7){
  CHR<-Fst_ordered[grep(paste("chr",x,"H",sep=""), Fst_ordered[,1]),]
  points(CHR$Cumul_pos/1000000, CHR$Mean_Fst, cex=0.6, col=COLOR[x])
}

legend("topright",pch=c(1,1,1,1,1,1,1,1),col=c("red","blue","green","magenta","orange","cyan","light green","black"),c("1H","2H","3H","4H","5H","6H","7H","Unknown"),cex=0.6)
dev.off()
