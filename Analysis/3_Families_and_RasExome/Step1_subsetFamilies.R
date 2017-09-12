# Author: Ana Poets
# Description: This code takes the HAPMAP file after basic vcf filtering and devides it in families. Adding a consensus copy of Rasmusson. It selects only those SNPs that 
# are present in all the parents pairs.
#########################################################################################################################################

rm(list=ls())

# List of all PI/CIho parents and their families. This list has the correct assignations for off-type parents.
Families<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Parent_fam_list_fixedOffType.txt")


# Filtered hapmap DP5, GQ30 and transformedHomoHete. Coded by nucleotide and HH is used for hetes.
#HAPMAP<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/MNS_July2016_NAM6ROW_DP5_GQ30_mis50.recode_transformedHETE_toNA_nuclHH_hmp.txt",header=T, row.names=1)

#HAPMAP without converting homo to hete missing
HAPMAP <-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/MNS_July2016_NAM6ROW_DP5_GQ30_mis50.recode_hmp.txt",header=T,row.names=1)

## Import the exome capture data just to get the parent missing and Rasmusson
#Exome1<-read.table("~/Dropbox/SmithLab/NAM/Analysis/Exome_parents/HAPMAP_nucleotides/Barley_GATK_Calling_NAM_genomewide_Parents_DP5_miss1_hmp.txt",header=T)
# Select HR655 parent PI328632 that in this file is 1F8
#Exome1_parent<-Exome1[,c(1,which(names(Exome1) == "X1F8"))]
#colnames(Exome1_parent)<-c("rs","PI328632")
#write.table(Exome1_parent,"~/Documents/smithlab/NAM/Analysis/NAM_parents_Pri/Exome_PI328632.txt",quote=F,row.names=F,col.names=T,sep="\t")

HR655_pi<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/Exome_PI328632.txt",header=F,row.names=1)

# Import rasmusson consensus GBS exome capture datasets 
Rasmusson<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt",header=F,row.names=1)

#Use only the calls among Ras_2 - Ras_5

#Rasmusson<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/Ras_consensus_gbs_gbsnames.txt",header=T)

#Converrt HR655_pi SNP names to GBS names
DataConvert<-c("HR655_pi")
for (j in 1:length(DataConvert)){
	DATA<-as.data.frame(get(DataConvert[j]))
	GBS_1<-gsub("chr1H_part1", "1H1", row.names(DATA))
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
	row.names(DATA)<-SNPnames_asGBS
	assign(DataConvert[j], DATA)
}

# Find SNPs shared between Rasmusson and HAPMAP
SharedSNPs<-intersect(row.names(Rasmusson) ,row.names(HAPMAP)) 

# Select the markers that are shared from HR655_pi and Rasmusson, so they are ready to be merged with the families divisions
HR655_pi_shared<-subset(HR655_pi,(row.names(HR655_pi) %in% SharedSNPs ))
Rasmusson_shared<-subset(Rasmusson,(row.names(Rasmusson) %in% SharedSNPs))
###
# Filtered hapmap DP5, GQ30 and transformedHomoHete. Coded by nucleotide and HH is used for hetes. Select SNPs that are in Rasmusson_shared

HAPMAP_shared<-subset(HAPMAP,(row.names(HAPMAP) %in% SharedSNPs))


#Get the SNPs from Rasmusson that match the markers in the HAPMAP
#Rasmusson_sh_hpm<-subset(Rasmusson_shared,(row.names(Rasmusson_shared) %in% row.names(HAPMAP_shared)))

#order Rasmusson  in the same SNP order as the HAPMAP
Ras_temp<-as.data.frame(cbind(row.names(Rasmusson_shared), as.character(Rasmusson_shared[,1])))
Rasmusson_sh_hpm_or<-Ras_temp[match(row.names(HAPMAP_shared),(Ras_temp[,1])),]

if (identical(as.character((Rasmusson_sh_hpm_or[,1])), as.character(row.names(HAPMAP_shared))) == FALSE)stop("Error: Rasmusson and HAPMAP SNPs are in different order")


#===============  Create HAPMAP files per each family =========================================================

#If the family is HR655 get the PI parent HR655_pi_shared, not from the family files.

# Get a vector with all the column names in the HAPMAP
NAMnames<-(names(HAPMAP_shared))
for (i in 1:(dim(Families)[1])){
	if (Families[i,2] != "HR655S"){
		#Select the columns for all samples including parents for each family
		FamIndiv<-c(grep(Families[i,1], NAMnames),grep(Families[i,2], NAMnames))
		
		#Form HAPMAP per family with all the SNPs in Rasmusson and the family. Some SNPs will be all NA for the family except for at least one of the parents
		HAPMAP_family<-cbind(as.data.frame(HAPMAP_shared[,1:3]), as.character(Rasmusson_sh_hpm_or[,2]), HAPMAP_shared[,c(FamIndiv)])
		names(HAPMAP_family)[4]<-"Ras_consensus"
		write.table(HAPMAP_family, paste("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Analysis/Hapmap_mis50/",Families[i,2],"_shrSNP_hmp.txt",sep =""), quote=F, row.names=T, col.names=T, sep="\t")
	}	
	if (Families[i,2] == "HR655S"){
		#Select the columns for all samples Excluding the parent for this family
		FamIndiv<-c(grep(Families[i,2], NAMnames))
		 # Get the HAPMAP only for the RILs with SNP information
		HAPMAP_family<-HAPMAP_shared[,c(1:3,FamIndiv)]
		# Select the SNPs in the HAPMAP present in the parent HR655_pi
		Shared_hmp_pi<-intersect(row.names(HAPMAP_family),row.names(HR655_pi_shared))
		# Select the SNPs that are present also in Rasmusson
		Shared_hmp_pi_ras<-intersect(Shared_hmp_pi, Rasmusson_sh_hpm_or[,1])
		
		#Get HAPMAP, Ras and PI data for those markers
		HAPMAP_family_sh<-HAPMAP_family[(row.names(HAPMAP_family) %in% Shared_hmp_pi_ras),]
		Ras_sh<-Rasmusson_sh_hpm_or[(Rasmusson_sh_hpm_or[,1] %in% Shared_hmp_pi_ras ),]
		
		#make sure SNPs are in the same order
		if (identical(row.names(HAPMAP_family_sh), as.character(Ras_sh[,1])) == FALSE)stop("Error! in family HR655 rasmusson and hampap SNPs are in different order")
		
		PI_sh<-subset(HR655_pi_shared, (row.names(HR655_pi_shared) %in% Shared_hmp_pi_ras ))
		# order PI SNPs according to HAPMAP
		PI_temp<-as.data.frame(cbind(row.names(PI_sh), as.character(PI_sh[,1])))
		PI_sh_hpm_or<-PI_temp[match(row.names(HAPMAP_family_sh),(PI_temp[,1])),]
	
		if (identical(row.names(HAPMAP_family_sh), as.character(PI_sh_hpm_or[,1])) == FALSE)stop("Error! in family HR655 PI and hampap SNPs are in different order")

		HAPMAP_family_HR655<-cbind(as.data.frame(HAPMAP_family_sh[,1:3]), as.character(Ras_sh[,2]), as.character(PI_sh_hpm_or[,2]) , as.data.frame(HAPMAP_family_sh[,-c(1:3)]))
		names(HAPMAP_family_HR655)[4]<-"Ras_consensus"
		names(HAPMAP_family_HR655)[5]<-"PI328632"
		write.table(HAPMAP_family_HR655, paste("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Analysis/Hapmap_mis50/",Families[i,2],"_shrSNP_hmp.txt",sep =""), quote=F, row.names=T, col.names=T, sep="\t")
	}
	
}

