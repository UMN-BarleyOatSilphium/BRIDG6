#Author: Ana Poets
# Description: Using the Consensus (gbs and exome) rasmusson calls convert the SNP names to the format used by priyanka
#=========================================================================================
rm(list=ls())

#import conesensus calls
Rasmusson<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome.txt",header=F,row.names=1)


	GBS_1<-gsub("chr1H_part1", "1H1", row.names(Rasmusson))
	GBS_1<-gsub("chr2H_part1","2H1", GBS_1)
	GBS_1<-gsub("chr3H_part1", "3H1",GBS_1)
	GBS_1<-gsub("chr4H_part1", "4H1",GBS_1)
	GBS_1<-gsub("chr5H_part1","5H1", GBS_1)
	GBS_1<-gsub("chr6H_part1", "6H1",GBS_1)
	GBS_1<-gsub("chr7H_part1", "7H1",GBS_1)
	
	GBS_1<-gsub("chr1H_part2", "1H2",GBS_1)
	GBS_1<-gsub("chr2H_part2", "2H2",GBS_1)
	GBS_1<-gsub("chr3H_part2", "3H2",GBS_1)
	GBS_1<-gsub("chr4H_part2", "4H2",GBS_1)
	GBS_1<-gsub("chr5H_part2", "5H2",GBS_1)
	GBS_1<-gsub("chr6H_part2", "6H2",GBS_1)
	GBS_1<-gsub("chr7H_part2", "7H2",GBS_1)
	
	SNPnames_asGBS<-gsub("chrUn", "UN",GBS_1)
	row.names(Rasmusson)<-SNPnames_asGBS
	
	write.table(Rasmusson,"/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt",quote=F,row.names=T,col.names=F,sep="\t")