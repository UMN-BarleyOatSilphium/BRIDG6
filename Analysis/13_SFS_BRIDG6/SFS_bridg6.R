#Author: Ana Poets
# Description: Estimate SFS in the RILs BRIDG6 and within each family
#################################################################################################
rm(list=ls())

# Load genotypic data (same used for GWAS)
# ref allele = 2 = major allele
# 0 = minor allele, this allows minor alleles to have different effects if stratification is provided
# SNPs are in columns 8,101. Samples are in rows(5,141)
gen = read.csv("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/LDKNNI/Input/genos_LDKNNI.csv")
row.names(gen)<-gen[,1]
gen<-gen[,-1]

# Assume samples are haploid (low heterozygosity)
MinorAlleleCount<-function(dat){
  #alleleHete <-length(which(dat == 1))
  #ras allele
  allele2<-(length(which(dat == '2')))
  #donor allele
  allele0<-length(which(dat == '0'))
  minAlle_count <-min(allele2,allele0)
  return(minAlle_count)
}

MinorCount <-apply(gen,2,MinorAlleleCount)

tableMinCount<-as.data.frame(table(MinorCount))

write.table(tableMinCount,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/SFS/SFS_count_bridg6_hap.txt",quote=F,row.names=F,col.names=T,sep="\t")

# Get the SFS for each population in the BRIDG6 RILs
PopulationsInfo<-read.csv("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/LDKNNI/Input/Family_summary.csv",header=T,sep=",")
# Plot SFS across all RILs

admixed<-PopulationsInfo[which(PopulationsInfo$Pop_assignment == "Admixed"),10]
asian<-PopulationsInfo[which(PopulationsInfo$Pop_assignment == "Asian"),10]
mediterranean<-PopulationsInfo[which(PopulationsInfo$Pop_assignment == "Coastal Mediterranean"),10]
eastAfrican<-PopulationsInfo[which(PopulationsInfo$Pop_assignment == "East African"),10]
centralEuropean<-PopulationsInfo[which(PopulationsInfo$Pop_assignment == "Central European"),10]

pops<-c("admixed","asian","mediterranean","eastAfrican","centralEuropean")

SummaryTable<-matrix(NA,ncol=2,nrow=5)

for (p in 1:length(pops)){
  subpop<-get(pops[p])
  pattern<-paste(subpop,collapse="|")
  temp<-gen[grep(pattern,row.names(gen)),]
  SummaryTable[p,1]<-pops[p]
  SummaryTable[p,2]<-dim(temp)[1]
  MinorCount <-apply(temp,2,MinorAlleleCount)
  
  tableMinCount<-as.data.frame(table(MinorCount))
  
  write.table(tableMinCount,paste("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/SFS/SFS_count_",pops[p],"_hap.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  
}

colnames(SummaryTable)<-c("Populations","NumberIndiv")
write.table(SummaryTable,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/SFS/SummaryTable.txt",quote=F,row.names=F,col.names=T,sep="\t")


#######################################################
#### Plots
########################################################
tableMinCount<-read.table(paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/SFS_Bridg6/output/SFS_count_","bridg6","_hap.txt",sep=""),header=T)
# Remove the 0 category
if (which(tableMinCount[,1] == 0)){
  tableMinCount<-tableMinCount[-which(tableMinCount[,1] == 0),]
}
# The whole population has 5,141 samples, since we threated the data as haploid there are classes=5,141
# for windows of 100 count (classes)
i=1 
n = 50

Results<-data.frame("Start"=NA,"End"= NA,"freq"=NA)
#5141/2 = max MAF
while (i < round(5141/2)){
x<-tableMinCount[which(tableMinCount[,1] >=i & tableMinCount[,1] < (i+n)),]
if (dim(x)[1] == 0){
  freqX = 0
}else{
freqX<-(sum(x[,2])/sum(tableMinCount[,2]))*100} # percentage of SNPs that are present at the freq of each class

temp<-data.frame("Start"=i,"End"= (i+n),"freq"=freqX)

Results<-rbind(Results,temp)
i= (i+n)

}

Results<-Results[-1,]



# Replace freq=0 with 
pdf("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/SFS_Bridg6/Plots/BRIDGrils_SFS.pdf",width=7, height=5)
plt<-barplot(Results[,3],Results[,2],xaxt="n",space=1, ylim=c(0,15), width=1, ylab="Percentage SNPs", xlab = "Minor Allele Count", main ="SFS BRIDG6 RILs")
text(plt, par("usr")[3], labels = paste(seq(1,round(5200/2),by=50),seq(50,round(5200/2),by=50),sep="-"), srt = 90, adj = c(1.1,0.5), xpd = TRUE, cex=0.6) 
dev.off()

## Plot the SFS in each of the 5 populations

#par(mfrow=c(3,2))
PopName = c("Admixed","Asian","Mediterranean","East African","Central European")
popsAll<-c("admixed","asian","mediterranean","eastAfrican","centralEuropean")
SummaryTable<-read.table("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/SFS_Bridg6/output/SummaryTable.txt",header=T)

for (p in 1:length(popsAll)) {
tableMinCount<-read.table(paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/SFS_Bridg6/output/SFS_count_",popsAll[p],"_hap.txt",sep=""),header=T)
TotalRILs<-SummaryTable[which(SummaryTable[,1] == popsAll[p]),2]

# Remove the 0 category
if (which(tableMinCount[,1] == 0)){
  tableMinCount<-tableMinCount[-which(tableMinCount[,1] == 0),]
}
  temp<-cbind(tableMinCount,tableMinCount[,1]/1103)
# The whole population has 5,141 samples, since we threated the data as haploid there are classes=5,141
# for windows of 100 count (classes)
i=1 
n = 50

Results<-data.frame("Start"=NA,"End"= NA,"freq"=NA)
#5141/2 = max MAF
while (i < round(TotalRILs/2)){
  x<-tableMinCount[which(tableMinCount[,1] >=i & tableMinCount[,1] < (i+n)),]
  if (dim(x)[1] == 0){
    freqX = 0
  }else{
    freqX<-(sum(x[,2])/sum(tableMinCount[,2]))*100} # percentage of SNPs that are present at the freq of each class
  
  temp<-data.frame("Start"=i,"End"= (i+n),"freq"=freqX)
  
  Results<-rbind(Results,temp)
  i= (i+n)
  
}

Results<-Results[-1,]


# Replace freq=0 with 
pdf(paste("/Users/Mia/Dropbox/SmithLab/NAM/Analysis/SFS_Bridg6/Plots/SFS_",popsAll[p],".pdf",sep=""),width=7,height=5)
plt<-barplot(Results[,3],Results[,2],xaxt="n",space=1, ylim=c(0,40), width=1, ylab="Percentage SNPs", xlab = "Minor Allele Count", main =paste("SFS ",PopName[p]," (",TotalRILs," RILs)",sep=""))
#text(plt, par("usr")[3], labels = seq(50,mround(TotalRILs,50),by=50), srt = 90, adj = c(1.1,0.5), xpd = TRUE, cex=0.8) 
text(plt, par("usr")[3], labels =paste(Results[,1],(Results[,2]-1),sep="-"), srt = 90, adj = c(1.1,0.5), xpd = TRUE, cex=0.8) 
dev.off()

}

