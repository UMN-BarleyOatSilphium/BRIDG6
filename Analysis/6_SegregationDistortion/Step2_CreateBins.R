# Author: Ana M Poets
# Description: Create bins to plot the Excess contribution of donor parent at each family using X2. Expectation 1:1 ratio Donor:Common parent genotype frequency at F5. Df = 2 -1 = 1 (two classes allowed AA and BB)
################################################################################################################################################

rm(list=ls())
library(stringr)
#===Bin chromosomes and estimate mean allele frequency of donor parent ===================================================================
# =========== Calculate the proportion of Rasmusson allele (1- donor allele freq ) in windows of ... SNPs =====================================
# Import table of frequency of donor allele
DATA<-read.table("~/Documents/SmithLab/NAM/Analysis/SNPdensity_and_segDist/SegDist_091817/Output/88Fam_genmap_SNPdonorgenotypeFreqAA.txt", header=T, row.names=1 )

# Separate by chromosomes

CHR1<-DATA[grep("1H", row.names(DATA)),]
CHR2<-DATA[grep("2H", row.names(DATA)),]
CHR3<-DATA[grep("3H", row.names(DATA)),]
CHR4<-DATA[grep("4H", row.names(DATA)),]
CHR5<-DATA[grep("5H", row.names(DATA)),]
CHR6<-DATA[grep("6H", row.names(DATA)),]
CHR7<-DATA[grep("7H", row.names(DATA)),]
CHR_UN<-DATA[grep("UN", row.names(DATA)),]


# Find the positions of SNPs at the begining of 100 snp segments
RESULTS_mean<-matrix(NA, ncol=dim(DATA)[2])
colnames(RESULTS_mean)<-colnames(DATA)

CHRnames<-c("CHR1", "CHR2", "CHR3", "CHR4", "CHR5", "CHR6", "CHR7") # ignore the UNknown chrosmosome since we cannot group it by physical proximity


# Function to calculate the mean, setting it to NA if there are not values in that segment
GetMean<-function(dat){
  if (length(which(is.na(dat))) == length(dat)){mean_value<-NA}
  mean_value<-mean(dat[which(!is.na(dat))])
  return(mean_value)
}

for (n in 1:length(CHRnames)){
  print (n)
  # Devide each chromosome in 200 parts and get the average frequency of donor allele for that segment
  CHR_X<-get(CHRnames[n])
  # number of SNPs per segment to get 200 segments
  SNPsPerSegment<-round(dim(CHR_X)[1]/200, digits=0)
  
  # Total posible segments
  dim(CHR_X)[1]/SNPsPerSegment
  TotalSegments<-round(dim(CHR_X)[1]/SNPsPerSegment)
  
  EndPosition<-seq(TotalSegments , dim(CHR_X)[1], by= SNPsPerSegment)
  
  # if the last position is less than the SNPsperSegment from the end of the chromosome, change the position to be the end of the chromosome
  if ((dim(CHR_X)[1]  - EndPosition[length(EndPosition)]) <SNPsPerSegment) {EndPosition[length(EndPosition)]<-(dim(CHR_X)[1])}
  
  
  start<-1
  for (s in 1: length(EndPosition)){
    #print (s)
    start<-start
    end<-EndPosition[s]
    
    SEGMENT<-CHR_X[start:end,]
    #MEAN<-as.data.frame(t(apply(SEGMENT, 2, summary)[4,]))
    MEAN<-as.data.frame(t(apply(SEGMENT , 2, GetMean)))
    
    #row.names(MEAN)<-paste(CHRnames[n], "_", s, sep="")
    row.names(MEAN)<-paste(CHRnames[n], "_", str_pad(s, 4, pad = "0"), sep="")
    # add values to our Results table
    #if (s==1){RESULTS_mean[s,]<-MEAN}else{RESULTS_mean<-rbind(RESULTS_mean, (MEAN))}
    RESULTS_mean<-rbind(RESULTS_mean, (MEAN))
    start<-(end+1)	
  }
  
}

#
RESULTS_mean<-RESULTS_mean[-1,]

write.table(RESULTS_mean, "~/Documents/SmithLab/NAM/Analysis/SNPdensity_and_segDist/SegDist_091817/Output/88Fam_genmap_segDist_binned_AA.txt", quote=F,row.names=T,col.names=T,sep="\t")
