# Author: Ana Poets
# Description: Estimate population structure between the NAM parents and the NSGC accessions
###########################################################################################################


#Create input files for smpartPCA. Using the ancestrymap format:
# Use minor allele as the reference allele
# Set gender to "U" Unknown 
# Set population group as 1 for all, to have an independent result from Structure. We don't have case/controls
# The genotype file uses "9" for missing data and "0" if reference allele (minor allele) is not present


rm(list=ls())

# Genotypes downlowaded from T3 in A/B calls
Genotypes_ALL<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/download_BAQY/snpfile.txt",header=T,row.names=1)

# List of donor parents to be found among the genotypes. Rasmusson es M109
DPlist<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/download_BAQY/DPnames.txt")

if(FALSE){# matrix ,rows are SNPs, columns are samples
Donor<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/Data_SFS/NAM_donor_genotype.hmp.txt",header=T,row.names=1)
nsgc<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/Data_SFS/NSGC_6_genotype.hmp.txt",header=T,row.names=1)

# get markers shared
SNPshared<-intersect(row.names(Donor),row.names(nsgc))

donor_sh<-Donor[(row.names(Donor) %in% SNPshared), -c(1:3)]
nsgc_sh<-nsgc[(row.names(nsgc) %in% SNPshared), -c(1:3)]

# Transform numeric genotypes to A/B
genoNotation<-function(dat){
  dat[which(dat == "-1")]<-"BB"
  dat[which(dat == "1")]<-"AA"
  dat[which(dat == "0")]<-"AB"
  return(dat)
}

donor_sh_nuc<-apply(donor_sh,2,genoNotation)
nsgc_sh_nuc<-apply(nsgc_sh,2,genoNotation)

#combine parents and nsgc genotypes
genotype_ALL<-cbind(donor_sh_nuc,nsgc_sh_nuc)
}
#turn table to have SNPs in rows and samples in columns
genotype_ALL<-as.data.frame(t(Genotypes_ALL))
# Remove SNPs missing in more than 50% of the data set
MISSING<-function(dat){
	miss<-length(which(is.na(dat)))/length(dat)
	return(miss)
}
missing_snp<-apply(genotype_ALL, 1, MISSING)

if (length(which(missing_snp >0.5))>0){
	genotype<-genotype_ALL[-c(which(missing_snp >0.5)),]
}else{genotype<-genotype_ALL}



write.table(genotype,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/NSCG_Diverse/NSGC_Diverse_genotypesBOPA.txt",quote=F,row.names=T,col.names=T,sep="\t")
#write.table(genotype,"~/Desktop/forPCA/NSGC_Diverse_genotypesBOPA.txt",quote=F,row.names=T,col.names=T,sep="\t")

## Import the genotypes

DATA<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/NSCG_Diverse/NSGC_Diverse_genotypesBOPA.txt", header=T, row.names=1)
#DATA<-read.table("~/Desktop/forPCA/NSGC_Diverse_genotypesBOPA.txt", header=T, row.names=1)

##################################
#Allow 20% missingness in SNPs 
MISSING_ALLOWED<-0.1
DATA_missing_snp<-apply(DATA, 1, MISSING)

if (length(which(DATA_missing_snp > MISSING_ALLOWED))>0){
	All_parentsNSGC <-DATA[-c(which(DATA_missing_snp > MISSING_ALLOWED)),]
}else{All_parentsNSGC<-DATA}

#################################

#Remove monomorphic markers
mono<-function(dat){
	AlleleA<-length(which(dat == 'AA'))
	AlleleB<-length(which(dat == 'BB'))
	if(AlleleA == 0 | AlleleB == 0){monomorphic<-'yes'}else{monomorphic<-'no'}
	return(monomorphic)
}

MONO<-apply(All_parentsNSGC,1, mono)

monoSites<-which(MONO == "yes")

if(length(monoSites) >0){
  All_parentsNSGC_poly<-All_parentsNSGC[-c(monoSites),]
}else{All_parentsNSGC_poly<-All_parentsNSGC}

dim(All_parentsNSGC_poly)


###################################################################################################################################################
#Find Minor allele at each SNP
minorAllele<-function(dat){
	AlleleA<-length(which(dat == 'AA'))
	AlleleB<-length(which(dat == 'BB'))
	minor_allele<-if(AlleleA <=AlleleB) {'A'} else {'B'}
	return(minor_allele)
}

MinorAllele_ref<-apply(All_parentsNSGC_poly,1, minorAllele)
MinorAllele_ref<-as.data.frame(MinorAllele_ref)

#Make a list of Major allele
MajorAllele<-function(dat){
	major_allele <-if (dat == 'A') {'B'} else {'A'}
	return (major_allele)
	}
Major_Allele_ref<-apply(MinorAllele_ref,1, MajorAllele)


# Set SNP information, then MAF, then genotypes
GENOTYPE_READY<-cbind(as.data.frame(MinorAllele_ref),as.data.frame(All_parentsNSGC_poly))



##1. Create SNP.snp file

#Since we don't know the centimorgan position of the SNP set them to a vector 1:length(SNP)
#MORGANS<-c(0: ((dim(GENOTYPE_READY)[1]) -1))
MORGANS<-rep(NA, (dim(GENOTYPE_READY)[1]))
# Get SNP information Chromosome and position, since we don't know the true position or all the SNPs, set them to 1
CHROMSOME_INFO<-rep("1", dim(GENOTYPE_READY)[1])

## === put SNP information data together ===============================================================================================
SNP_file<-cbind(row.names(GENOTYPE_READY),as.data.frame(CHROMSOME_INFO), as.data.frame(MORGANS),rep(NA,dim(GENOTYPE_READY)[1]), as.data.frame(GENOTYPE_READY[,1]), as.data.frame(Major_Allele_ref))
colnames(SNP_file)<-c("SNP_name","Chromosome","cM","Position","Reference_mino","Reference_major")

#write.table(SNP_file,"/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/input/NAM_RILs.snp",quote=F,row.names=F,col.names=F,sep="\t")
#when only 20 % missinges is allowed

#write.table(SNP_file,"~/Documents/SmithLab/NAM/Analysis/smartPCA/NSGC_diverseParents/NSGC_DP.snp",quote=F,row.names=F,col.names=F,sep="\t")
write.table(SNP_file,"~/Desktop/forPCA/NSGC_DP.snp",quote=F,row.names=F,col.names=F,sep="\t")


## 2. Create Sample.ind
# Get only genotypes names without SNP information
Genotypes_names<-colnames(GENOTYPE_READY)[-1]
Samples_ind<-cbind(as.data.frame(Genotypes_names), as.data.frame(rep('U',length(Genotypes_names))), as.data.frame(rep(1,length(Genotypes_names))))
colnames(Samples_ind)<-c("SampleID","Gender","pop_group")

#when only 20 % missinges is allowed
#write.table(Samples_ind,"~/Documents/SmithLab/NAM/Analysis/smartPCA/NSGC_diverseParents/NSGC_DP.ind",quote=F,row.names=F,col.names=F,sep="\t")
write.table(Samples_ind,"~/Desktop/forPCA/NSGC_DP.ind",quote=F,row.names=F,col.names=F,sep="\t")
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

#write.table(Genotype_counts,"~/Documents/SmithLab/NAM/Analysis/smartPCA/NSGC_diverseParents/NSGC_DP.eigenstratgeno",quote=F,row.names=F,col.names=F,sep="")
write.table(Genotype_counts,"/Users/agonzale/Desktop/forPCA/NSGC_DP.eigenstratgeno",quote=F,row.names=F,col.names=F,sep="")

#when only 20 % missinges is allowed

print(paste("SNPs total:",dim(Genotype_counts)[1]))

#Now use the files as input files for smartPCA in the command line in MSI. package part of Eigensatat.Using a par file containing the parameters used. All the files should be in the "bin/" folder 
# of EIGENSATATS.
# [agonzale@login03] (~/Software/EIG-6.1.4/bin) $ ./smartpca -p NAM_all_smartPCA.par

#Use Plot_smartPCA_output.R to process output