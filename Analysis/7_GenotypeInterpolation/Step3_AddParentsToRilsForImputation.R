# Author: Ana M Poets
# Combine RIL with Parents for interpolation
########################################################################################################################

rm(list=ls())
Family<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/List_Fam_afterQC.txt",header=F)

RILs_ras_complete_RasBased <-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_rasBased_rm100SNPclose_naExcessHH_naMAF_NAonlyHeteHomo.txt",header=T,row.names=1) # 168672 SNPs   6060samples including ras

RASMUSSON<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt",header=F)

PARENTS<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/MNS_July2016_Parents_DP5_GQ30.recode_transformedHETE_toNA_HH_hmp.txt",header=T)


# Get only markers that are in the RILS
RAS_ril<-RASMUSSON[(RASMUSSON[,1] %in% row.names(RILs_ras_complete_RasBased) ),] #168672 snps
PARENTS_ril<-PARENTS[(PARENTS[,1] %in% row.names(RILs_ras_complete_RasBased) ), -c(2:4)] # 168672 snps

# Combine the parents and transform calls to Rasmusson

# order SNPs according to parents
RAS_ril_or_par<-RAS_ril[match(PARENTS_ril[,1], RAS_ril[,1]),]

identical(as.character(RAS_ril_or_par[,1]), as.character(PARENTS_ril[,1]))

ALL_parents<-cbind(RAS_ril_or_par, PARENTS_ril)
row.names(ALL_parents)<-ALL_parents[,1]

ALL_parents<-ALL_parents[,-c(1,3)]

names(ALL_parents)[1]<-"Ras_consensus"


# Rass = 2, donor parent = 0, hete =1
ConvertToRas<-function(dat){
	dat[which(dat == dat[1])]<-"2"
	dat[which(dat != dat[1] & dat != 'HH' & !is.na(dat))]<-"0"
	dat[which(dat == 'HH')]<-"1"
	return(dat)
}

ALL_parents_RasBased<-apply(ALL_parents, 1, ConvertToRas)
ALL_parents_RasBased<-as.data.frame(ALL_parents_RasBased)
t_ALL_parents_RasBased<-t(ALL_parents_RasBased)
t_ALL_parents_RasBased[1:16,1:6]


# Combine parents to RILs
t_ALL_parents_RasBased_or<-t_ALL_parents_RasBased[match(row.names(RILs_ras_complete_RasBased),row.names(t_ALL_parents_RasBased)),]

identical(as.character(row.names(t_ALL_parents_RasBased_or)), as.character(row.names(RILs_ras_complete_RasBased)))

Whole_NAM<-cbind(t_ALL_parents_RasBased, RILs_ras_complete_RasBased[,-c(1)])

write.table(Whole_NAM,"/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/Imputation/ForPri_imp/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_rasBased_rm100SNPclose_naExcessHH_naMAF_NAonlyHeteHomo_PARENTS.txt",quote=F,row.names=T,col.names=T,sep="\t")


# Use this file as input for TASSEL LD-kNNI