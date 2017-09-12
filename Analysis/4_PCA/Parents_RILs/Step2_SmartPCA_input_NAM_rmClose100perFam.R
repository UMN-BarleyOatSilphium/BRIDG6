# Author: Ana Poets
# Description: Estimate population structure among NAM in the file used originaly with 10% missing that ploted three families outliers
###########################################################################################################


#Create input files for smpartPCA. Using the ancestrymap format:
# Use minor allele as the reference allele
# Set gender to "U" Unknown 
# Set population group as 1 for all, to have an independent result from Structure. We don't have case/controls
# The genotype file uses "9" for missing data and "0" if reference allele (minor allele) is not present


rm(list=ls())

# matrix ,rows are SNPs columns are samples
# input file comes from Step1_CombineGentoyes_NAM_RIL_DonorP_Ras.R
# This files has 80% missing in the RILs
genotype_ALL <-read.table("/home/smithkp/agonzale/Projects/NAM/DATA/HAPMAP_RIL_parents_Ras_80mis/HAPMAP_RIL_parents_Ras_80mis_20het_hmp.txt",header=T,row.names=1)

# Remove SNPs missing in more than 50% of the data set
MISSING<-function(dat){
	miss<-length(which(is.na(dat)))/length(dat)
	return(miss)
}
missing_snp<-apply(genotype_ALL, 1, MISSING)

if (length(which(missing_snp >0.5))>0){
	genotype<-genotype_ALL[-c(which(missing_snp >0.5)),]
}else{genotype<-genotype_ALL}

# Remove all SNPs for which Rasmusson is NA

if (length(which(is.na(genotype[,1]))) >0){
genotype<-genotype[-c(which(is.na(genotype[,1]))),]
} else{genotype<-genotype}

if (length(which(genotype[,1] == "HH"))>0){
genotype<-genotype[-c(which(genotype[,1] == "HH")),]
}else{genotype<-genotype}


# Convert genotypes respect to Rasmusson (common parent). Rasmusson is in the first column
	# Set Rassmusson allele as "AA" and Donor Parent as "BB", heterozygotes = "AB"	
	Ras_based<-function(dat){
		dat<-(dat)
		dat[which(dat ==(dat[1]))]<-"XX"
		dat[which(dat !=(dat[1]) & dat !=("HH"))]<-"YY"
		dat[which(dat == "HH")]<-"ZZ"
		return(dat)
		}

	CHANGE<-function(dat){
		dat<-(dat)
		dat[which(dat == "XX")]<-"AA"
		dat[which(dat == "YY")]<-"BB"
		dat[which(dat == "ZZ")]<-"AB"
		return(dat)
		}

## Divide the genotype file into 7 smaller files to process the next two steps
genotype_1<-genotype[1:21657,]
genotype_2<-genotype[21658:43314,]
genotype_3<-genotype[43315:64971,]
genotype_4<-genotype[64972:86628,]
genotype_5<-genotype[86629:108285,]
genotype_6<-genotype[108286:129942,]
genotype_7<-genotype[129943:(dim(genotype)[1]),]

for (n in 1:7){
DATA<-get(paste("genotype_",n,sep=""))
DATA_ras<-as.data.frame(t(apply(DATA,1, Ras_based)))
All_NAM_rasbased<-as.data.frame(t(apply(DATA_ras,1,CHANGE)))
assign(paste("All_NAM_rasbased_",n,sep=""), All_NAM_rasbased)

if (n==1){
write.table(All_NAM_rasbased,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/TEMP/All_NAM_rasbased_1.txt",quote=F,row.names=T,col.names=T,sep="\t")
}

if(n!=1){
	write.table(All_NAM_rasbased,paste("/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/TEMP/All_NAM_rasbased_",n,".txt",sep=""),quote=F,row.names=T,col.names=F, sep="\t")
}
}

All_NAM_rasbased_ALL<-rbind(All_NAM_rasbased_1,All_NAM_rasbased_2,All_NAM_rasbased_3,All_NAM_rasbased_4,All_NAM_rasbased_5,All_NAM_rasbased_6,All_NAM_rasbased_7)

write.table(All_NAM_rasbased_ALL,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/ALL_NAMrasbased.txt",quote=F,row.names=T,col.names=T,sep="\t")

## Import the All_NAM_rasbased_ALL

All_NAM_rasbased_ALL<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/All_NAM_rasbased_ALL.txt", header=T, row.names=1)

##################################
#Allow 20% missingness in SNPs 
MISSING_ALLOWED<-0.1
All_NAM_missing_snp<-apply(All_NAM_rasbased_ALL, 1, MISSING)

if (length(which(All_NAM_missing_snp > MISSING_ALLOWED))>0){
	All_NAM_rasbased_ALL_1 <-All_NAM_rasbased_ALL[-c(which(All_NAM_missing_snp > MISSING_ALLOWED)),]
}else{All_NAM_rasbased_ALL_1<-All_NAM_rasbased_ALL}

All_NAM_rasbased_ALL <-All_NAM_rasbased_ALL_1 #26445 SNPs for 6152 individuals
#################################

#Remove monomorphic markers
mono<-function(dat){
	AlleleA<-length(which(dat == 'AA'))
	AlleleB<-length(which(dat == 'BB'))
	if(AlleleA == 0 | AlleleB == 0){monomorphic<-'yes'}else{monomorphic<-'no'}
	return(monomorphic)
}

MONO<-apply(All_NAM_rasbased_ALL,1, mono)

monoSites<-which(MONO == "yes")

if(length(monoSites) >0){
	All_NAM_rasbased<-All_NAM_rasbased_ALL[-c(monoSites),]
}else{All_NAM_rasbased<-All_NAM_rasbased_ALL}

dim(All_NAM_rasbased) #5828 SNP

# Remove parent PI328632 that we know is the wrong sample. For analysis we used the genotypes from the exome capture data
All_NAM_rasbased <-All_NAM_rasbased [,-grep("PI328632", colnames(All_NAM_rasbased))]


write.table(All_NAM_rasbased,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/10missing/data/All_NAM_rasbased_genotypes.txt",quote=F,row.names=T,col.names=T,sep="\t")

########################################################################################################################
# ============= Identify which SNPs on the per family basis , are closer than 100 bp. Make a list of them and remove them from the big file =========================
# Import families info
Fam_info<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/10missing/data/Pop_structure_propDiss_JULY2016.txt",header=T)

#function to obtain only polymorphic markers
#Remove monomorphic markers
mono<-function(dat){
	AlleleA<-length(which(dat == 'AA'))
	AlleleB<-length(which(dat == 'BB'))
	if(AlleleA == 0 | AlleleB == 0){monomorphic<-'yes'}else{monomorphic<-'no'}
	return(monomorphic)
}

# Create a variable to retain the SNPs closer than 100bp in each family
SNPclose_allFam<-NULL

for (f in 1:dim(Fam_info)[1]){
	if (Fam_info[f,3] != "HR648_HR649_HR650" & Fam_info[f,3] !=  "HR632_HR651"){
		#Select all the individuals in a family
		Fam_genotypes<-All_NAM_rasbased [,grep(Fam_info[f,3], colnames(All_NAM_rasbased))]
		}
	
	if (Fam_info[f,3] != "HR648_HR649_HR650"){
		#Select all the individuals in a family
		Fam_genotypes<-All_NAM_rasbased [,grep("HR648|HR649|HR650", colnames(All_NAM_rasbased))]
	}

	if (Fam_info[f,3] != "HR632_HR651"){
		#Select all the individuals in a family
		Fam_genotypes<-All_NAM_rasbased [,grep("HR632|HR651", colnames(All_NAM_rasbased))]
	}

		#Remove non-segregating markers
		mono_obs<-apply(Fam_genotypes, 1, mono)
		if(length(which(mono_obs == "yes")) >0){
			Fam_genotypes_poly<-Fam_genotypes[-c(which(mono_obs == "yes")),]
		}else{Fam_genotypes_poly<-Fam_genotypes}
		
		# Using only the polymorphic markers identify those closer than 100bp
		# =========== Remove SNPs that are separated by 100 bp or less =======================================================
				  Chromosome_parts<-c("1H1" , "1H2" ,"2H1" ,"2H2" , "3H1", "3H2" , "4H1" , "4H2" , "5H1" , "5H2" , "6H1" , "6H2" ,"7H1" ,"7H2")
				  for (q in 1:length(Chromosome_parts)){
				    chr_part<-Chromosome_parts[q]
				    # Separate each CHR part
				    chrom<-Fam_genotypes_poly[grep(chr_part, row.names(Fam_genotypes_poly)),]
				    
				    POSITIONS<-sub(paste(chr_part,"_",sep=""),"",row.names(chrom))
				    
				    #starting from the second position, get the difference between x2-x1
				    
				    BasePairs<-NA
				    for (p in 2:length(POSITIONS)){
				      BasePairs<-c( BasePairs,(as.numeric(POSITIONS[p])-as.numeric(POSITIONS[p-1])))
				    }
				    
				    if (length(which(BasePairs <=100)) >0){
				    # find which SNPs are less than 100 bp appart. 
				    FirstItem<-which(BasePairs <=100)
				    #Since the comparision was done x2-x1, find the position of Xi-1, this is the pair that is 100 bp from the FirstItem 
				    SecondItem<-(FirstItem -1)
				    
				    CLOSE_SNPs<-sort(c(FirstItem, SecondItem))
				    
				    # Get the name of SNPs that have to be removed
				    SNP_name_remove<-row.names(chrom)[CLOSE_SNPs]
				    
				    }else{SNP_name_remove<-NA}
				    assign(paste("SNPfor_", q, sep=""),SNP_name_remove )
				    
				  }
				  
				  # Get all the SNPs that have to be remove # removes 27,674 SNPs
				  SNP_removeALL<-c(SNPfor_1, SNPfor_2, SNPfor_3, SNPfor_4 ,SNPfor_5, SNPfor_6, SNPfor_7, SNPfor_8, SNPfor_9 ,SNPfor_10, SNPfor_11, SNPfor_12, SNPfor_13, SNPfor_14)
				
				List_uniq_SNPremove<-unique(SNP_removeALL)
				SNPclose_allFam<-c(SNPclose_allFam, List_uniq_SNPremove)
	
}

# Remove SNPs that are 100 bp close in any given family #497 SNPs total
All_NAM_rasbased_rm<-All_NAM_rasbased[!(row.names(All_NAM_rasbased) %in% SNPclose_allFam),]
###################################################################################################################################################
#Find Minor allele at each SNP
minorAllele<-function(dat){
	AlleleA<-length(which(dat == 'AA'))
	AlleleB<-length(which(dat == 'BB'))
	minor_allele<-if(AlleleA <=AlleleB) {'A'} else {'B'}
	return(minor_allele)
}

MinorAllele_ref<-apply(All_NAM_rasbased_rm,1, minorAllele)
MinorAllele_ref<-as.data.frame(MinorAllele_ref)

#Make a list of Major allele
MajorAllele<-function(dat){
	major_allele <-if (dat == 'A') {'B'} else {'A'}
	return (major_allele)
	}
Major_Allele_ref<-apply(MinorAllele_ref,1, MajorAllele)


# Set SNP information, then MAF, then genotypes
GENOTYPE_READY<-cbind(as.data.frame(MinorAllele_ref),as.data.frame(All_NAM_rasbased_rm))


##1. Create SNP.snp file

#Since we don't know the centimorgan position of the SNP set them to a vector 1:length(SNP)
#MORGANS<-c(0: ((dim(GENOTYPE_READY)[1]) -1))
MORGANS<-rep(NA, (dim(GENOTYPE_READY)[1]))
# Get SNP information Chromosome and position
CHROMSOME_INFO<-c(rep("1", length(grep("1H", row.names(GENOTYPE_READY)))),  rep("2", length(grep("2H", row.names(GENOTYPE_READY)))) , rep("3", length(grep("3H", row.names(GENOTYPE_READY)))) , rep("4", length(grep("4H", row.names(GENOTYPE_READY)))), rep("5", length(grep("5H", row.names(GENOTYPE_READY)))) , rep("6", length(grep("6H", row.names(GENOTYPE_READY)))), rep("7", length(grep("7H", row.names(GENOTYPE_READY)))), rep("8", length(grep("UN", row.names(GENOTYPE_READY)))) )

# == Calculate cumulative physical position ====
  #get physical positions
  Positions_phy_mbp <-as.data.frame(as.character(row.names(GENOTYPE_READY)))
  Physical_positions <-  as.data.frame(apply(Positions_phy_mbp,1, function(x) strsplit(as.character(x), "_")[[1]][2]))
  
  Positions_phy_mbp<-as.data.frame(Physical_positions)
  row.names(Positions_phy_mbp)<-row.names(GENOTYPE_READY)
  
  head(Positions_phy_mbp)
  
  CHR1_1<-grep("1H1",row.names(Positions_phy_mbp))
  CHR1_2<-grep("1H2",row.names(Positions_phy_mbp))
  
  CHR2_1<-grep("2H1",row.names(Positions_phy_mbp))
  CHR2_2<-grep("2H2",row.names(Positions_phy_mbp))
  
  CHR3_1<-grep("3H1",row.names(Positions_phy_mbp))
  CHR3_2<-grep("3H2",row.names(Positions_phy_mbp))
  
  CHR4_1<-grep("4H1",row.names(Positions_phy_mbp))
  CHR4_2<-grep("4H2",row.names(Positions_phy_mbp))
  
  CHR5_1<-grep("5H1",row.names(Positions_phy_mbp))
  CHR5_2<-grep("5H2",row.names(Positions_phy_mbp))
  
  CHR6_1<-grep("6H1",row.names(Positions_phy_mbp))
  CHR6_2<-grep("6H2",row.names(Positions_phy_mbp))
  
  CHR7_1<-grep("7H1",row.names(Positions_phy_mbp))
  CHR7_2<-grep("7H2",row.names(Positions_phy_mbp))
  
  CHR_UN<-grep("UN",row.names(Positions_phy_mbp))
  
  List_chr_parts<-c("CHR1_1","CHR1_2","CHR2_1", "CHR2_2", "CHR3_1", "CHR3_2", "CHR4_1", "CHR4_2", "CHR5_1", "CHR5_2", "CHR6_1", "CHR6_2", "CHR7_1", "CHR7_2","CHR_UN")
   
  #Calculate joint physical positions for both parts of a chromosome, then calculate Cumulative positons genome wide
  Positions_original<-Positions_phy_mbp[,1]

  # Number of bp to add to each part to have a cumulative physical position, including the UNKnown CHR
  ADD_CHR<-c(0,312837513,558535432,952068106,1326610456,1720921089,2026321570,2381382776,2673381728,3054247210,3343411888,3638233958,3926792401,4252589917, 4584016401)
  
  #Calculate joint physical positions for both parts of a chromosome, then calculate Cumulative positons genome wide
  Positions_original<-Physical_positions[,1]
  
  New_order_phy<-NULL
  for (C in 1:15){
	CHR_pos<-as.numeric(as.character(Physical_positions[get(List_chr_parts[C]),1]))
	NEW_chr_pos<-CHR_pos + ADD_CHR[C]
	New_order_phy <-c(New_order_phy ,NEW_chr_pos)
    }  
	Positions_phy_mbp_cumul<-cbind(Positions_phy_mbp, New_order_phy)


## === put SNP information data together ===============================================================================================
SNP_file<-cbind(row.names(GENOTYPE_READY),as.data.frame(CHROMSOME_INFO), as.data.frame(MORGANS),as.data.frame(Positions_phy_mbp_cumul[,2]), as.data.frame(GENOTYPE_READY[,1]), as.data.frame(Major_Allele_ref))
colnames(SNP_file)<-c("SNP_name","Chromosome","cM","Position","Reference_mino","Reference_major")

#write.table(SNP_file,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/input/NAM_RILs.snp",quote=F,row.names=F,col.names=F,sep="\t")
#when only 20 % missinges is allowed
write.table(SNP_file,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/10missing/input_rm100/NAM_RILs.snp",quote=F,row.names=F,col.names=F,sep="\t")


## 2. Create Sample.ind
#Remove parents
#Parents_position<-grep("Ras_|PI|CIho", colnames(GENOTYPE_READY))
#GENOTYPE_READY <-GENOTYPE_READY[,-c(Parents_position)]

# Get only genotypes names without SNP information
Genotypes_names<-colnames(GENOTYPE_READY)[-1]
Samples_ind<-cbind(as.data.frame(Genotypes_names), as.data.frame(rep('U',length(Genotypes_names))), as.data.frame(rep(1,length(Genotypes_names))))
colnames(Samples_ind)<-c("SampleID","Gender","pop_group")

#write.table(Samples_ind,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/input/NAM_RILs.ind",quote=F,row.names=F,col.names=F,sep="\t")
#when only 20 % missinges is allowed
write.table(Samples_ind,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/10missing/input_rm100/NAM_RILs.ind",quote=F,row.names=F,col.names=F,sep="\t")

## 3. Create Genotype.eigenstratgeno

COUNT_ALLELE<-function(dat){
	#remove maf
	GENOTYPES<-dat[-c(1)]
	MAF<-dat[1]
	
	AA<-which(GENOTYPES == 'AA')
	BB<-which(GENOTYPES == 'BB')
	AB<-which(GENOTYPES == 'AB')
	Missing<- which(is.na(GENOTYPES) == 'TRUE')
	
	if(MAF == 'A'){
		GENOTYPES[AA]<-'2'
		GENOTYPES[BB]<-'0'
	}
	if(MAF == 'B'){
		GENOTYPES[AA]<-'0'
		GENOTYPES[BB]<-'2'
	}
	
	#All AB are 1
	GENOTYPES[AB]<-1
	
	#All missing values are = 9
	GENOTYPES[Missing]<-'9'
	return(GENOTYPES)
}

RESULT<-as.data.frame(apply(GENOTYPE_READY,1, COUNT_ALLELE))

#Turn table to have SNPs in rows and samples in columns

Genotype_counts<-as.data.frame(t(RESULT))

write.table(Genotype_counts,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/10missing/input_rm100/NAM_RILs.eigenstratgeno",quote=F,row.names=F,col.names=F,sep="")
#when only 20 % missinges is allowed
#write.table(Genotype_counts,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/20missing/input/NAM_RILs.eigenstratgeno",quote=F,row.names=F,col.names=F,sep="")

print(paste("SNPs total:",dim(Genotype_counts)[1]))

#Now use the files as input files for smartPCA in the command line in MSI. package part of Eigensatat.Using a par file containing the parameters used. All the files should be in the "bin/" folder 
# of EIGENSATATS.
# [agonzale@login03] (~/Software/EIG-6.1.4/bin) $ ./smartpca -p NAM_all_smartPCA.par

#Use Plot_smartPCA_output.R to process output