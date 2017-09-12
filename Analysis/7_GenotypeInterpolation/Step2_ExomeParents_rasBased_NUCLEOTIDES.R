# Author: Ana Poets
# Description:
# 1. Get exome capture data for only parents. 
# 2. Get GBS data in nucleotides for all parents
# 3. Add PI328632 exome, to the GBS data ready for imputation
################################################################################################################################
rm(list=ls())

# The parents are divided in two files, so import both files
EXOME1<-read.table("/home/smithkp/agonzale/Projects/NAM/DATA/ExomeCapture/Barley_NAM_Regionlast26_parents.txt",header=T)
EXOME2<-read.table("/home/smithkp/agonzale/Projects/NAM/DATA/ExomeCapture/Barley_GATK_Calling_parents.txt",header=T)

# Check for duplicate SNPs, these come from Lin's way of handleling the large file
which(table(EXOME1[,1]) >1) # EXOME1 doesn't have duplicated SNPs
which(table(EXOME2[,1]) >1) # has two duplicated SNPs 1H2_229894927,3H1_37439055; Remove SNPs in positions 508451 and 2375353

EXOME2_uniq<-EXOME2[-c(508451,2375353),]

# Set SNPs as row names
row.names(EXOME1)<-EXOME1[,1]
row.names(EXOME2_uniq)<-EXOME2_uniq[,1]

# Import Rasmusson consensus 
Rasmusson<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt",header=F)


# Import the table with the list of parents in the 88 families
Fam_table<-read.table("/home/smithkp/agonzale/Projects/NAM/DATA/Pop_structure_propDiss_JULY2016_parentLessMissingData.txt",header=T)

# Get parents from exom capture data (there are 74 in this data set, we don't take the duplicate parents)
Exome1_parents<-EXOME1[,(colnames(EXOME1) %in% Fam_table$Parent_lessMissingData)]
Exome2_parents<-EXOME2_uniq[,(colnames(EXOME2_uniq) %in% Fam_table$Parent_lessMissingData)]

write.table(Exome1_parents,"/home/smithkp/agonzale/TEMP/Exome1_parents_NUC.txt",quote=F,row.names=T,col.names=T,sep="\t")
write.table(Exome2_parents,"/home/smithkp/agonzale/TEMP/Exome2_parents_NUC.txt",quote=F,row.names=T,col.names=T,sep="\t")

# Sort Exome by Rasmusson SNP order
Exome1_parents_orRas<-Exome1_parents[match(Rasmusson[,1],row.names(Exome1_parents)),]
Exome2_parents_orRas<-Exome2_parents[match(Rasmusson[,1],row.names(Exome2_parents)),]

# add Rasmusson SNP names
Exome1_par_orRas_nameRas<-cbind(Rasmusson[,1],Exome1_parents_orRas)
Exome2_par_orRas_nameRas<-cbind(Rasmusson[,1],Exome2_parents_orRas)

# Make sure that SNP names match among the Exome and Ras files and cbind them.
identical(as.character(Exome1_par_orRas_nameRas[,1]), as.character(Exome2_par_orRas_nameRas[,1]))
identical(as.character(Exome1_par_orRas_nameRas[,1]), as.character(Rasmusson[,1])) 

Exome_all_parents<-cbind(Rasmusson[,2],Exome1_par_orRas_nameRas[,-c(1)], Exome2_par_orRas_nameRas[,-c(1)])
row.names(Exome_all_parents)<-Rasmusson[,1]
names(Exome_all_parents)[1]<-"Ras_consensus"

# Remove rows where all samples are NA
Exome_all_parents_noNAall<-Exome_all_parents[rowSums(is.na(Exome_all_parents))!= dim(Exome_all_parents)[2], ] # 155910 SNPs removed. There are 4601059 SNPs and 75 Samples

write.table(Exome_all_parents_noNAall,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Exome_PARENTS_Nucleotides_NOrasBased.txt",quote=F,row.names=T,col.names=T,sep="\t")

##################################################################################################################################
##########=== Get Exome capture data for PI328632 based on Rasmusson ========================
Exome_all_parents_noNAall<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Exome_PARENTS_Nucleotides_NOrasBased.txt",header=T,row.names=1)

Ras_PI328632_exome<-Exome_all_parents_noNAall[,grep("Ras|PI328632", colnames(Exome_all_parents_noNAall))]
# Remove markers that have HH in Ras

Ras_PI328632_exome_noHH_noNA<-Ras_PI328632_exome[-c(which(Ras_PI328632_exome[,1] == "HH")),]

# Get the markers used in the file ready for imputation
# import data for imputation 
PARENTS_RILSforimputation<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_rmIndmisHete_monoSNP_wPARENTS.txt",header=T,row.names=1)

Ras_PI328632_exome_noHH_inImpu<-Ras_PI328632_exome_noHH_noNA[(row.names(Ras_PI328632_exome_noHH_noNA) %in% row.names(PARENTS_RILSforimputation)),]

# order SNPs in this parent as the imputation data, and merge data sets
Ras_PI328632_exome_noHH_inImpu_rasBased_or<-Ras_PI328632_exome_noHH_inImpu[match(row.names(PARENTS_RILSforimputation),row.names(Ras_PI328632_exome_noHH_inImpu)),]
if(identical(as.character(row.names(Ras_PI328632_exome_noHH_inImpu_rasBased_or)), as.character(row.names(PARENTS_RILSforimputation))) == FALSE)stop("Error! Imputation-ready data and PPI328632 have diff SNP order")

Imputation_RIL_Parents_ready<-cbind(as.data.frame(PARENTS_RILSforimputation[,1]),as.data.frame(Ras_PI328632_exome_noHH_inImpu_rasBased_or[,2]), as.data.frame(PARENTS_RILSforimputation[,-c(1)]))
names(Imputation_RIL_Parents_ready)[2]<-"PI328632"
names(Imputation_RIL_Parents_ready)[1]<-"Ras_consensus"

write.table(Imputation_RIL_Parents_ready,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_noRasHH_NOrasBased_Polymorphic_rm100SNPclose_naExcessHH_rmIndmisHete_monoSNP_wPARENTS_PI328632_NUC.txt",quote=F,row.names=T,col.names=T,sep="\t")
##################################################################################################################################
###### === Get GBS data for all parents, remove parent for HR655 PI328632 ===============================================

# Import Rasmusson consensus 
Rasmusson<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt",header=F)

# Import GBS parents
GBS_parents<-read.table("~/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/MNS_July2016_Parents_DP5_GQ30.recode_transformedHETE_toNA_HH_hmp.txt",header=T,row.names=1)

# Import the table with the list of parents in the 88 families
Fam_table<-read.table("/home/smithkp/agonzale/Projects/NAM/DATA/Pop_structure_propDiss_JULY2016_parentLessMissingData.txt",header=T)

GBS_parents_noHR655<-GBS_parents[,-c(which(colnames(GBS_parents) == Fam_table[which(Fam_table$Family == "HR655S"),6]))]

# Get unique parents for the 87 families remaining families
GBS_parents_noHR655_87fam<-GBS_parents_noHR655[,(colnames(GBS_parents_noHR655) %in%Fam_table$Parent_lessMissingData )]

# Add Rasmusson
# Get SNPs in GBS only

Ras_inGBS<-Rasmusson[(Rasmusson[,1] %in% row.names(GBS_parents_noHR655_87fam) ),]

Ras_inGBS_or<-Ras_inGBS[match(row.names(GBS_parents_noHR655_87fam),Ras_inGBS[,1]),]

if (identical(as.character(Ras_inGBS_or[,1]), as.character(row.names(GBS_parents_noHR655_87fam))) == FALSE)stop("Error! Cannot merge GBS parents to Ras. SNPs dont match")

GBS_DP_Ras<-cbind(Ras_inGBS_or[,2],GBS_parents_noHR655_87fam)

names(GBS_DP_Ras)[1]<-"Ras_consensus"

write.table(GBS_DP_Ras,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/GBS_PARENTS_Nucleotides_NOrasBased.txt",quote=F,row.names=T,col.names=T,sep="\t")



