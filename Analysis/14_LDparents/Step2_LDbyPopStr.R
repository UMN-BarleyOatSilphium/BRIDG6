# Author: Ana Poets
# Description: LD estimation in BRIDG6 parents, using GBS data prior interpolation of markers
##################################################################################################
rm(list=ls())
library(gtools)

# Load genotypic data (same used for GWAS)
# ref allele = 2 = major allele
# 0 = minor allele, this allows minor alleles to have different effects if stratification is provided
# SNPs are in rows. Samples are in rows(5,141)
gen = read.table("~/Dropbox/SmithLab/NAM/Analysis/LD_BRIDG6parents/data/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_noSegDist_rmIndmisHete_monoSNP_wPARENTS.txt",header=T)
parentsALL<-gen[,1:92]

popStr<-read.csv("~/Dropbox/SmithLab/NAM/Data/Alex/Pop_structure92.csv", sep=",")

admixed<-popStr[which(popStr[,2] == "Admixed"),1]
asian<-popStr[which(popStr[,2] == "Asian"),1]
centralEu<-popStr[which(popStr[,2] == "Central European"),1]
coastalM<-popStr[which(popStr[,2] == "Coastal Mediterranean"),1]
african<-popStr[which(popStr[,2] == "East African"),1]


# Convert To Rasmusson = 2 , donor=0, hete = HH=1
ConvertRass<-function(dat){
  Ras<-dat[1]
  dat<-dat[-1]
  dat[which(dat == "HH")]<-1
  dat[which(dat == as.character(Ras))]<-2
  dat[which(dat != 2 & dat != 1)]<-0
  return(dat)
}

parentsALL_Ras<-t(as.data.frame(apply(parentsALL,1,ConvertRass)))
## Remove SNPs missing more than 80% 

MISSING<-function(dat){
  missing<-length(which(is.na(dat)))/length(dat)
  return(missing)
}

## Remove SNPs with MAF <10%
MAF<-function(dat){
      allele1<-length(which(dat == 0))
      allele2<-length(which(dat == 2))
      freqAl<-(min(allele1,allele2))/(allele1+allele2)
      return(freqAl)
}

# Function to convert genotypes to PLINK format
# Convert genotypes to A/C: 2 =AA; 0=CC;1=AC, NA=0
# Convert genotypes to A/C: 2 =AA; 0=CC;1=AC, NA=0
CONVERTgenoPlink<-function(dat){
  dat[which(dat == "2")]<- paste("1","1",sep="\t")
  dat[which(dat == "0")]<-paste("2","2",sep="\t")
  dat[which(dat == "1")]<-paste("1","2",sep="\t")
  dat[which(is.na(dat))]<-paste("0","0",sep="\t")
  return(dat)
}

pops<-c("admixed","asian","centralEu","coastalM","african")


for (p in 1:length(pops)){
  POP<-get(pops[p])
  parents<-parentsALL_Ras[,(colnames(parentsALL_Ras) %in% POP)]
  
  missing_snp<-apply(parents,1,MISSING) #14478 SNPs removed due to 10%missing

  if (length(which(missing_snp >0.1)) >0){
    DATA_ms<-parents[-(which(missing_snp >0.1)),]
  }else{DATA_ms<-parents}

  print (paste("SNP removed for missing >10% = ", length(which(missing_snp >0.1)) ,sep=" ")) #0

  # Remove SNP with MAF <10%
  mafparents<-apply(DATA_ms,1,MAF)

  if (length(which(mafparents <0.1)) >0){
    DATA_ms<-DATA_ms[-(which(mafparents <0.1)),]
  }else{DATA_ms<-DATA_ms}  
  
  print (paste("SNP removed for MAF <10% = ", length(which(mafparents <0.1)) ,sep=" ")) #0
  
########## Remove SNPs in UNknown chromosome ########################
  DATA_NA<-DATA_ms[grep('UN',row.names(DATA_ms)),]
  
  #write.table(row.names(DATA_NA),"/Users/Mia/Dropbox/SmithLab/NAM/Analysis/LD_BRIDG6parents/data/List_UNknown_snps.txt",quote=F,row.names=F,col.names=F,sep="\t")
  DATA_noUN_remain<-DATA_ms[-grep('UN',row.names(DATA_ms)),] #160376 SNP remain
  dir.create(paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/LD_BRIDG6parents/r2_plink/populations/",pops[p],sep=""))
  dir.create(paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/LD_BRIDG6parents/r2_plink/populations/",pops[p],"/input",sep=""))
  dir.create(paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/LD_BRIDG6parents/r2_plink/populations/",pops[p],"/output",sep=""))
    # Create new map and ped files to run r2 in PLINK
    #Length of part1 of each chromosome
    Length_p1<-c(312837513,393532674,394310633,355061206,380865482,294822070,325797516)
    # Separate the data by chromosomes, then create the PED and MAP files for each chromosome
    for (n in 1:7){
      DATA_chr<-DATA_noUN_remain[grep(paste(n,"H",sep=""),row.names(DATA_noUN_remain)),]
      
      #### Create MAP file
      ####a file with marker information, including chromosome number, the SNP name, the position in morgans on the chromosome, and the position in base pairs on the chromosome , with the extension .map. Again, this file has NO HEADER so looks like this:
      
      SNPs<-as.data.frame(row.names(DATA_chr))
      SNP_positions <-  as.data.frame(apply(SNPs,1, function(x) strsplit(as.character(x), "_")[[1]][2]))
      Positions<-cbind(SNPs[,1], SNP_positions)
      #get cumulative bp count within the chromosome
      PART1<-Positions[grep("H1",(Positions[,1])),2]
      PART2<-Positions[grep("H2",(Positions[,1])),2]
      
      #correct positions
      PART2_cum<-as.numeric(as.character(PART2)) + Length_p1[n]
      Total_positions<-c(as.numeric(as.character(PART1)), as.numeric(as.character(PART2_cum)))
      print(Total_positions[length(Total_positions)])
      #we don't have morgan positions, so use the physical position instead.
      #MAPfile<-cbind(rep(n,dim(DATA_chr)[1]), row.names(DATA_chr), as.numeric(as.character(Total_positions)), as.numeric(as.character(Total_positions)))
      MAPfile<-cbind(rep(n,dim(DATA_chr)[1]), row.names(DATA_chr), NA,as.numeric(as.character(Total_positions)))
      write.table(MAPfile,paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/LD_BRIDG6parents/r2_plink/populations/",pops[p],"/input/CHR_",n,".map",sep=""), quote=F,row.names=F,col.names=F,sep="\t")
      #Convert genotypes to A/C. and Turns table to have samples in rows
      t_DATAchr_noUN_genotypes<-as.data.frame(apply(DATA_chr,1, CONVERTgenoPlink))
      
      # Get Family names
      Sample_ID<-row.names(t_DATAchr_noUN_genotypes)
      
      #remove Rasmusson
      #Progeny_ID<-as.data.frame(Sample_ID[-1])
      #Family_name <-  as.data.frame(apply(Progeny_ID,1, function(x) strsplit(as.character(x), "S")[[1]][1]))
      #Family_ID<-c("Ras_consensus", as.character(Family_name[,1]))
      
      # Ped format
      PEDfile<-cbind(as.data.frame(Sample_ID), as.data.frame(Sample_ID), rep(0, dim(t_DATAchr_noUN_genotypes)[1]), rep(0, dim(t_DATAchr_noUN_genotypes)[1]) , rep("other", dim(t_DATAchr_noUN_genotypes)[1]), rep(0, dim(t_DATAchr_noUN_genotypes)[1]), t_DATAchr_noUN_genotypes)
      
      write.table(PEDfile,paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/LD_BRIDG6parents/r2_plink/populations/",pops[p],"/input/CHR_",n,".ped",sep=""), quote=F,row.names=F,col.names=F,sep="\t")
      
    }
}
    # Set --ld-window-kb to the length of the longest chromosome in kb. I am not sure what this is --ld-window, but I am told to set it to some stupidly large number.
    #plink --bfile "snp-thin" --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out "snp-thin"
    #largest chromosome 767855.1 kb


    #plink --file CHR_7 --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 767855.1 --out "snp-thin_chr7"
#for i in 1 2 3 4 5 6 7
#do
#plink --file CHR_${i} --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 767855.1 --out "snp-thin_chr${i}"
#done

###############################################################################################################
###### Plot R2 vs Distance
###########
pops<-c("admixed","asian","centralEu","coastalM","african")

for (p in 1:length(pops)){
    LDresults<-data.frame()
    for (chr in 1:7){
      output<-read.table(paste("/home/iaa/agonzale/LD/",pops[p],"/input/output/","snp-thin_chr",chr,".ld",sep=""),header=T)
      
      
      distance<-abs(output[,2] - output[,5])
      
      info<-cbind(distance, output[,7])
      LDresults<-rbind(LDresults,info)
    }
    
    
    # get average every 0.5 Mbp
    
    SortedResults<-LDresults[order(LDresults[,1]),]
    SortedResults_Mbp<-cbind(((SortedResults[,1])/1000000),SortedResults[,2])
    
    lengthwindows=0.5
    start= 0
    end = 0.5
    AvePoints<-round(max(SortedResults_Mbp[,1])/end)
    Average<-matrix(NA,ncol=2,nrow=AvePoints)
    
    for (i in 1:AvePoints){
      aveWindows<-mean(SortedResults_Mbp[which(as.numeric(SortedResults_Mbp[,1]) >= start & as.numeric(SortedResults_Mbp[,1]) <= end),2])
      Average[i,1]<-start
      Average[i,2]<-aveWindows
      start<-end+1
      end<-start + lengthwindows
    }
    
    #pdf("~/Desktop/LD_bridgParents.pdf",width=7,height=5)
    pdf(paste("/home/iaa/agonzale/LD/Plot_pops/LD_",pops[p],"_noThin.pdf",sep=""),width=7,height=5)
    
    plot(LDresults[,1]/1000000,LDresults[,2],cex=0.2,col="light gray",xlab="Distance (Mbp)",ylab="r^2")
    
    points(Average[,1],Average[,2],cex=0.2,col="red")
    lines(Average[,1],Average[,2],col="red")
    dev.off()
}