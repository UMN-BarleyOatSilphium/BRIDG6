# Author: Ana Poets
# Description: All the SNPs in the GBS data for the parents, when at least one parent has the SNP. 
# Calculate pairwise percentage genetic distance between the parents, using the function dist.gene from the ape package
#=========================================================================================================================
rm(list =ls ())

library(ape)

# Import Parents genotypes 92 parents
PARENTS <-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/MNS_July2016_Parents_DP5_GQ30.recode_transformedHETE_toNA_HH_hmp.txt", header=T,row.names=1)

data<-PARENTS[,-c(1:3)]

### Find for each pair of DonorParent and other Donor Parent which SNPs are present and homozygous in at least one parent.

DifferencesBtwParents<-matrix(NA, nrow=(dim(data)[2]), ncol=(dim(data)[2]))
row.names(DifferencesBtwParents)<-names(data)
colnames(DifferencesBtwParents)<-names(data)

# matrix to save percent differences using heterozygous
DifferencesBtwParents_hete<-matrix(NA, nrow=(dim(data)[2]), ncol=(dim(data)[2]))
row.names(DifferencesBtwParents_hete)<-names(data)
colnames(DifferencesBtwParents_hete)<-names(data)

# matrix to track the number of SNPs in each comparison
SNPinvolved<-matrix(NA, nrow=(dim(data)[2]), ncol=(dim(data)[2]))
row.names(SNPinvolved)<-names(data)
colnames(SNPinvolved)<-names(data)

for ( i in 1:(dim(data)[2])){
	
	#compare each parent to all other parents
	for (j in 1:(dim(data)[2])){
		Parents_pair<-cbind(as.character(data[,i]), as.character(data[,j]))
		row.names(Parents_pair)<-row.names(data)
		Parents_pair<-as.data.frame(Parents_pair)
			
		#Remove SNPs that are missing in both or either parents
		BOTH_NA<-which(is.na(Parents_pair[,1]) | is.na(Parents_pair[,2]))
		if (length(BOTH_NA) >0){
		Parents_present<-Parents_pair[-c(BOTH_NA),]}else{Parents_present<-Parents_pair}
		
		# Find which SNPs are HETE in either parent or both parents
		BOTH_HH<-which(Parents_present[,1] == 'HH' | Parents_present[,2] == 'HH')
		if(length(BOTH_HH) >0){
		Parents_noBothHete<-Parents_present [-c(BOTH_HH),]}else{Parents_noBothHete<-Parents_present }
		
		
		GoodSNPs<-row.names(Parents_noBothHete)

		SNPinvolved[i,j]<-length(GoodSNPs)
		###=========convert to numeric to use correlation function ============
		DIFFERENCES<-length(which(Parents_noBothHete[,1] != Parents_noBothHete[,2]))
		DIFFERENCES<-dist.gene(t(Parents_noBothHete), method="percentage", pairwise.deletion = T, variance = FALSE)
		# amount differences based on sites compared
		DifferencesBtwParents[i,j]<- round(DIFFERENCES, digits=3)
		
		# Make comparisons using all the Hete classes too
		DIFFERENCES_hete<-dist.gene(t(Parents_present), method="percentage", pairwise.deletion = T, variance = FALSE)
		# amount differences based on sites compared
		DifferencesBtwParents_hete[i,j]<- round(DIFFERENCES_hete, digits=3)

		}
}

write.table(DifferencesBtwParents,"/home/smithkp/agonzale/Projects/NAM/Analysis/GeneticDistanceParents/Pairwise_genDistance_92parents_gbs_onlyHomo.txt",quote=F,row.names=T,col.names=T,sep="\t")

write.table(DifferencesBtwParents_hete,"/home/smithkp/agonzale/Projects/NAM/Analysis/GeneticDistanceParents/Pairwise_genDistance_92parents_gbs_includingHete.txt",quote=F,row.names=T,col.names=T,sep="\t")

write.table(SNPinvolved,"/home/smithkp/agonzale/Projects/NAM/Analysis/GeneticDistanceParents/Pairwise_genDistance_92parents_gbs_onlyHomo_SNPcount.txt",quote=F,row.names=T,col.names=T,sep="\t")