#Author: Ana M. Poets
# Description: Get consensus Rasmusson call from the 4 Ras_ copies of GBS data
###################################################################################################

### === Function to identify a consensus calls for Rasmusson ===========================
CONSENSUS<-function(dat){
		IS_NA<-length(dat[is.na(dat)])
		TOTAL<-length(dat)
		if (IS_NA == TOTAL) {CALL <-NA}
			if ((IS_NA) == 0) {NO_NA<-dat}else {NO_NA<- dat[-which(is.na(dat))] }
		
			ALLELES<-table(t(NO_NA))
			if	(length(ALLELES) == 1){CALL<-names(ALLELES)}
			if (length(ALLELES) > 1) {CALL<-NA}			
		return(CALL)
}



## =========== INPUT FILTER RAS.VCF in NUCLEOTIDE HAPMAP===============================

DATA<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Ras_DP5_GQ30.recode_transformedHETE_toNA_nuclHH_hmp.txt",header= T, row.names= 1)

# Remove SNPs information, Rasmusson, and Ras_1

DATA_ready<-DATA[,-c(1:4,6)]

Ras_gbs<-(as.data.frame(apply(DATA_ready,1, CONSENSUS)))
names(Ras_gbs)<-"Ras_consensus"

write.table(Ras_gbs,"/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/Ras_consensus_gbs_gbsnames.txt", quote=F, row.names=T, col.names=T, sep="\t")
# Convert GBS data SNPs calls to long names "chr1H_part1" as in the Exome capture datasets.

GBS_1<-gsub("1H1","chr1H_part1", row.names(Ras_gbs))
GBS_1<-gsub("2H1","chr2H_part1", GBS_1)
GBS_1<-gsub("3H1","chr3H_part1", GBS_1)
GBS_1<-gsub("4H1","chr4H_part1", GBS_1)
GBS_1<-gsub("5H1","chr5H_part1", GBS_1)
GBS_1<-gsub("6H1","chr6H_part1", GBS_1)
GBS_1<-gsub("7H1","chr7H_part1", GBS_1)

GBS_1<-gsub("1H2","chr1H_part2", GBS_1)
GBS_1<-gsub("2H2","chr2H_part2", GBS_1)
GBS_1<-gsub("3H2","chr3H_part2", GBS_1)
GBS_1<-gsub("4H2","chr4H_part2", GBS_1)
GBS_1<-gsub("5H2","chr5H_part2", GBS_1)
GBS_1<-gsub("6H2","chr6H_part2", GBS_1)
GBS_1<-gsub("7H2","chr7H_part2", GBS_1)

GBS_1<-gsub("UN","chrUn", GBS_1)

#Replace new SNP names in GBS dataset
row.names(Ras_gbs)<- GBS_1

write.table(Ras_gbs,"/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/Ras_consensus_gbs_exomenames.txt", quote=F, row.names=T, col.names=T, sep="\t")