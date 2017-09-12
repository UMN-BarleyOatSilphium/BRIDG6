# Author: Ana Poets
# Description: Combine genotypes from RILs and parents after removing SNPs closer than 100 bp
###########################################################################################################################
rm(list=ls())

# Import RILs and SNPs that passed QC:6059 Samples and  173261 SNPs 
QCrils<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_rasBased_naExcessHH_naMAF_NAonlyHeteHomo.txt", header=T,row.names=1)
# Import NAM 80% missing values and 6336 individuals
RILs_all<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp.txt", header = T ,row.names=1)

# Import Consensus Rasmusson
RAS<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt", header=F,row.names=1)

# Import Parents genotypes
Parents<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/MNS_July2016_Parents_DP5_GQ30.recode_transformedHETE_toNA_HH_hmp.txt", header=T,row.names=1)

### === Get Nucleotides for genotypes that passed QC in RILs
RILs<-RILs_all[(row.names(RILs_all) %in% row.names(QCrils)), (colnames(RILs_all) %in% colnames(QCrils))]
### ======== Add parental genotypes =========================================================================================
# identify the SNPs present in RILs and Parents
SHARED_ras_rils <-(intersect(row.names(RAS), row.names(RILs))) 

SHARED_all<-intersect(SHARED_ras_rils, row.names(Parents)) #173261 shared SNPs
# Get the SNPs in the RILs, Parenst and Rasmusson that are shared
RILsh<-subset(RILs, (row.names(RILs) %in% SHARED_all))

RAS_sh<-subset(RAS, (row.names(RAS) %in% SHARED_all))

# order all files according to the RILs

RAS_sh<-as.data.frame(cbind(row.names(RAS_sh), as.data.frame(RAS_sh[,1])))
RAS_sh_or<-RAS_sh[match(row.names(RILsh),(RAS_sh[,1])),]
colnames(RAS_sh_or)<-c("SNP","Ras_consensus")
if (identical(as.character(RAS_sh_or[,1]), as.character(row.names(RILsh))) == FALSE)stop("Rasmusson is in different order than SNPs in RILs")

# donor parents
Parents_sh<-subset(Parents, (row.names(Parents) %in% SHARED_all))

# order all files according to the RILs
Parents_sh_or<-Parents_sh[match(row.names(RILsh),row.names(Parents_sh)),]
if (identical(as.character(RAS_sh_or[,1]), as.character(row.names(Parents_sh_or))) == FALSE)stop("Rasmusson is in different order than SNPs in donor parents")

## Combine parents with RILs
NAM_complete<-cbind(as.data.frame(RAS_sh_or[,2]),  as.data.frame(Parents_sh_or [,-c(1:3)]), as.data.frame(RILsh))
names(NAM_complete)[1]<-"Ras_consensus"

# bring parents to the first columns
NAM_complete_or<-cbind(as.data.frame(NAM_complete [,grep("Ras_|PI|CIho",colnames(NAM_complete))]), as.data.frame(NAM_complete [,-(grep("Ras_|PI|CIho",colnames(NAM_complete)))]))
write.table(NAM_complete_or,"/home/smithkp/agonzale/Projects/NAM/DATA/HAPMAP_RIL_parents_Ras_80mis/HAPMAP_RIL_parents_Ras_80mis_20het_hmp.txt",quote=F,row.names=T,col.names=T,sep="\t")







