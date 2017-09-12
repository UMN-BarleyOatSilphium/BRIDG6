# Author: Ana Poets
# Description: QC per family using the 80% missing data file. Using Nucleotide calls.
# 0. Remove SNPs that have HH in Ras
# 1. Remove monomorphic SNPs
# 2. Remove SNPs closer than 100 bp
# 3. Convert to NA sites cousing bias in segregation distortion
# 4. Convert to NA sites that lead to >20% hete per SNP
# 5. Remove samples missing >90% SNPs
# 6. Remove Samples >20% heterozygosity
# 7. Remove monomorphic SNPs again
# 8. Remove samples missing >90% again
#######################################################################################################################################
rm(list=ls())

HAPMAP<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp.txt",header=T,row.names=1)

Rasmusson<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt",header=F)


# Add Rasmusson to the file
library(gtools)

MISSING<-function(dat){
  missing<-length(dat[is.na(dat)])/(length(dat))
  return(missing)
}

# Subset Rasmusson SNPs that are in the NAM
Rasmusson_sh<-Rasmusson[(Rasmusson[,1] %in% row.names(HAPMAP) ),]

# order Rasmusson SNPs as in HAPMAP
Rasmusson_sh_or<-Rasmusson_sh[match(row.names(HAPMAP), Rasmusson_sh[,1]),]

if ( identical(as.character(Rasmusson_sh_or[,1]), row.names(HAPMAP)) == FALSE) stop ("Rasmusson and RILs are not in the same order")

# Combine RILs and Rasmusson genotypes
RILs_ras<-cbind(Rasmusson_sh_or ,HAPMAP)

RILs_ras<-RILs_ras[,-c(3:5)]

names(RILs_ras)[2]<-"Ras_consensus"

#write.table(RILs_ras,"/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/NAM_MNS_July2016_DP5_GQ30_mis80_5775ind.recode_transformedHETE_toNA_HH_hmp_withRAS.txt", quote=F,row.names=F,col.names=T,sep="\t")

### Remove SNPs that have HH in Rasmusson

RILs_ras_complete<-RILs_ras [-c(which(RILs_ras$Ras_consensus == "HH" )),]

write.table(RILs_ras_complete,"/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/KeepNucleotides/NAM_MNS_July2016_DP5_GQ30_mis80_5775ind.recode_transformedHETE_toNA_HH_hmp_RasHHremoved.txt", quote=F,row.names=F,col.names=T,sep="\t")
##############################################################################################################################################################################

### ============== Apply some filters per family ================================================================================
#############################=1. Remove monomorphic markers across the data set ===================#########################################################

# Import Genotypes (A/C/T/G) for all RILs and Ras. Ras HH were removed. There are 178842 SNPs in the file.
RILs_ras_complete<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/KeepNucleotides/NAM_MNS_July2016_DP5_GQ30_mis80_5775ind.recode_transformedHETE_toNA_HH_hmp_RasHHremoved.txt",header=T,row.names=1)
#function to obtain only polymorphic markers
#Remove monomorphic markers and those segregating only between one homozyogus class and hete.
mono<-function(dat){
  GENOTYPES<-table(dat)
  # if there are more than 1 genotype after removing Hete, then it is polymorphic
  if (length(names(GENOTYPES)) == 1){monomorphic<-'yes'}
  if (length(names(GENOTYPES)) >1){
    #remove Hete genotypes
    if(length(which(names(GENOTYPES) == "HH")) >0){
      GENOTYPES_noHH<-GENOTYPES[-c(which(names(GENOTYPES) == "HH"))]
    }else{GENOTYPES_noHH<-GENOTYPES}
    
    # if there are still two genotypes then marker is polymorphic
    if (length(names(GENOTYPES_noHH)) == 1){monomorphic<-'yes'}else{monomorphic<-'no'}
  }
  
  return(monomorphic)
}

mono_count<-apply(RILs_ras_complete,1,mono)
print(paste("There were ",length(which(mono_count == "yes")), " SNPs removed because were monomorphic",sep="")) # 4 markers were removed

All_NAM_poly<-RILs_ras_complete[-c(which(mono_count == "yes")),] #37783 SNPs, 6080 samples

# write table for only polymorphic markers
write.table(All_NAM_poly,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_RasHHrm_Polymorphic.txt",quote=F,row.names=T,col.names=T,sep="\t")

# ============= Identify which SNPs on the per family basis , are closer than 100 bp. Make a list of them and remove them from the big file ==
# Import genotypes with only polymorphic markers. If the marker was seg only for hete and one homo it was removed as monomorphic (37,783  SNPs  for 6080 samples)
All_NAM_poly<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_RasHHrm_Polymorphic.txt",header=T,row.names=1)
# Import families info
Fam_info<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/smartPCA/80mis_allRIL_Parents/10missing/data/Pop_structure_propDiss_JULY2016.txt",header=T)

#Function to remove monomorphic SNPs if they are trully segregating ony for one variant
mono_heteCounts<-function(dat){
  GENOTYPES<-table(dat)
  # if there are more than 1 genotype after removing Hete, then it is polymorphic
  if (length(names(GENOTYPES)) <= 1){monomorphic<-'yes'}else{monomorphic<-'no'}
  return(monomorphic)
}


# Create a variable to retain the SNPs closer than 100bp in each family
SNPclose_allFam<-NULL

for (f in 1:dim(Fam_info)[1]){
  print(f)
  if (Fam_info[f,3] != "HR648_HR649_HR650" & Fam_info[f,3] !=  "HR632_HR651"){
    #Select all the individuals in a family
    Fam_genotypes<-All_NAM_poly [,grep(Fam_info[f,3], colnames(All_NAM_poly))]
  }
  
  if (Fam_info[f,3] == "HR648_HR649_HR650"){
    #Select all the individuals in a family
    Fam_genotypes<-All_NAM_poly [,grep("HR648|HR649|HR650", colnames(All_NAM_poly))]
  }
  
  if (Fam_info[f,3] == "HR632_HR651"){
    #Select all the individuals in a family
    Fam_genotypes<-All_NAM_poly [,grep("HR632|HR651", colnames(All_NAM_poly))]
  }
  

  #Remove non-segregating markers within a family
  mono_obs<-apply(Fam_genotypes, 1, mono_heteCounts)
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

print(length(unique(SNPclose_allFam)))
# Remove SNPs that are 100 bp close in any given family #497 SNPs total
All_NAM_poly_rm<-All_NAM_poly[!(row.names(All_NAM_poly) %in% SNPclose_allFam),]

write.table(All_NAM_poly_rm,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC.txt", quote=F,row.names=T,col.names=T,sep="\t")

##############################################################################################################################

### ======== = Convert to NA SNPs in a family that have excess of heterozygosity >0.20 =================
HETE<-function(dat){
  missing<-length(dat[is.na(dat)])
  AlleleAB_freq<-length(which(dat == 1))/(length(dat)-missing)
  return(AlleleAB_freq)
}

# Import file with SNPs closer than 100 bp removed
All_NAM_poly_rm <-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC.txt",header=T,row.names=1)

# Import Family names
Family<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/List_Fam_afterQC.txt",header=F)


for (n in 1:dim(Family)[1]){
  print(n)
  if (Family[n,1] != "HR632_HR651" & Family[n,1] != "HR648_HR649_HR650"){
    # grab all individuals in the family
    Position_indiv<-grep(Family[n,1], colnames(All_NAM_poly_rm))
    Genotypes_family<-All_NAM_poly_rm[, Position_indiv]
    # Estimate percentage heterozygosity at each SNP
    Hete_fam<-apply(Genotypes_family, 1, HETE)
    if(length(which(Hete_fam >0.2)) >0){
      All_NAM_poly_rm[which(Hete_fam >0.2), Position_indiv] <-NA
    } 			  		
  }
  
  if(Family[n,1] == "HR632_HR651"){
    # grab all individuals in the family
    Position_indiv<-grep("HR632|HR651", colnames(All_NAM_poly_rm))
    Genotypes_family<-All_NAM_poly_rm[, Position_indiv]
    # Estimate percentage heterozygosity at each SNP
    Hete_fam<-apply(Genotypes_family, 1, HETE)
    if(length(which(Hete_fam >0.2)) >0){
      All_NAM_poly_rm[which(Hete_fam >0.2), Position_indiv] <-NA
    } 	
  }
  
  if(Family[n,1] =="HR648_HR649_HR650"){
    # grab all individuals in the family
    Position_indiv<-grep("HR648|HR649|HR650", colnames(All_NAM_poly_rm))
    Genotypes_family<-All_NAM_poly_rm[, Position_indiv]
    # Estimate percentage heterozygosity at each SNP
    Hete_fam<-apply(Genotypes_family, 1, HETE)
    if(length(which(Hete_fam >0.2)) >0){
      All_NAM_poly_rm[which(Hete_fam >0.2), Position_indiv] <-NA
    } 	
  }
}


write.table(All_NAM_poly_rm,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH.txt", quote=F,row.names=T,col.names=T,sep="\t")

###############################################################################################################################################################
######################################### CONVERT TO NA HETES AND MINOR ALLELE IF SNP IS DISTORTED WITHIN A FAMILY ############################################
###############################################################################################################################################################
SEG_DIST<-function(dat){
  GENOTYPES<-table(dat)
  
  # if there is only one genotype, then do nothing
  if (dim(GENOTYPES)[1] == 1){dat<-dat}
  
  #if there are two genotypes: 
  if (dim(GENOTYPES)[1] == 2){
    # is any HH? convert all HH to NA
    if (length(which(names(GENOTYPES) == "HH")) == 1){dat[which(dat == "HH")] <-NA}
    # if neither is HH estimated seg distortion between the two alleles
    if (length(which(names(GENOTYPES) == "HH")) == 0){
      a<-GENOTYPES[1]
      b<-GENOTYPES[2]
       #calculate p-value
      total<-a+b
 
        DiffA<-((a/b - 1)^2)/1
        DiffB<-((b/a) -1)^2/1
        
        SumChiSq<-DiffA+DiffB
        p_value<-pchisq(SumChiSq,1, lower.tail=F)
        
        # If significant distorition then convert minor allele to NA
        SIGNIFICANCE<-0.01
        if(p_value <=0.01){
          #find minor allele
          if(a>b){mino<-names(b)}else{mino<-names(a)}
          dat[which(dat == mino)]<-NA
        }else{dat <-dat}
    }
  }
  # if there are three genotypes
  if (dim(GENOTYPES)[1] == 3){
    # Calculate chi2
    Hete_position<-which(names(GENOTYPES) == "HH")
    NoHetePosition<-which(names(GENOTYPES) != "HH")
    
    Allele1<-GENOTYPES[NoHetePosition[1]]*2 + ((GENOTYPES[Hete_position]))
    Allele2<-GENOTYPES[NoHetePosition[2]]*2 + ((GENOTYPES[Hete_position]))
    #calculate p-value
    total<-Allele1+Allele2
    
    DiffA<-((Allele1/Allele2 - 1)^2)/1
    DiffB<-((Allele2/Allele1) -1)^2/1
    
    SumChiSq<-DiffA+DiffB
    p_value<-pchisq(SumChiSq,1, lower.tail=F)
    
    # If significant distorition then convert minor allele to NA, as well as all Heterozygous calls
    SIGNIFICANCE<-0.01
    if(p_value <=0.01){
      #find minor allele
      if(Allele1>Allele2){mino<-names(Allele2)}else{mino<-names(Allele1)}
      dat[which(dat == mino)]<-NA
      dat[which(dat == "HH")]<-NA
    }else{dat <-dat}
  }
return(dat)
}

# 31076 SNPs for 6080 lines
All_NAM_poly_rm <-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH.txt",header=T,row.names=1)

# Import Family names
Family<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Fileter80missing/List_Fam_afterQC.txt",header=F)

# for each family, convert to NA sites where SNP is distorted or segregating only for one homo and a hete states.

for (n in 1:dim(Family)[1]){
  print(n)
  rm(SegDist_fam_t)
  if (Family[n,1] != "HR632_HR651" & Family[n,1] != "HR648_HR649_HR650"){
    # grab all individuals in the family
    Position_indiv<-grep(Family[n,1], colnames(All_NAM_poly_rm))
    Genotypes_family<-All_NAM_poly_rm[, Position_indiv]
    # Estimate segregation distortion at each SNP
    SegDist_fam<-apply(Genotypes_family, 1, SEG_DIST)
    SegDist_fam_t<-t(SegDist_fam)
  }
  
  if(Family[n,1] == "HR632_HR651"){
    # grab all individuals in the family
    Position_indiv<-grep("HR632|HR651", colnames(All_NAM_poly_rm))
    Genotypes_family<-All_NAM_poly_rm[, Position_indiv]
    # Estimate segregation distortion at each SNP
    SegDist_fam<-apply(Genotypes_family, 1, SEG_DIST)
    SegDist_fam_t<-t(SegDist_fam) 	
  }
  
  if(Family[n,1] =="HR648_HR649_HR650"){
    # grab all individuals in the family
    Position_indiv<-grep("HR648|HR649|HR650", colnames(All_NAM_poly_rm))
    Genotypes_family<-All_NAM_poly_rm[, Position_indiv]
    # Estimate segregation distortion at each SNP
    SegDist_fam<-apply(Genotypes_family, 1, SEG_DIST)
    SegDist_fam_t<-t(SegDist_fam)	
  }
  
  if(n ==1){
    All_fam_noDistortion<-SegDist_fam_t
  }else{
    if( identical(as.character(row.names(All_fam_noDistortion)), as.character(row.names(SegDist_fam_t))) == TRUE) {
      All_fam_noDistortion<-cbind(All_fam_noDistortion,SegDist_fam_t)
    }else{stop(paste("Error! Family ", as.character(Family[n,1])),"has SNPs in different order or missing; cannot be combined",sep="")}
  }
}

# Add Rasmusson

All_NAM_poly_NoSegDist<-cbind(as.character(All_NAM_poly_rm[,1]),All_fam_noDistortion)
colnames(All_NAM_poly_NoSegDist)[1]<-"Ras_consensus"

write.table(All_NAM_poly_NoSegDist,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_noSegDist.txt", quote=F,row.names=T,col.names=T,sep="\t")

###############################################################################################################################################################
#######################################################################
# Remove Samples missing more than 90% and those with heterozygosity >20%
#import data set with 31076SNPs  6080samples
All_NAM_poly_NoSegDist<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_noSegDist.txt", header=T,row.names=1)

MISSING<-function(dat){
  missing<-length(dat[is.na(dat)])/(length(dat))
  return(missing)
}

Missing_inSamples<-apply(All_NAM_poly_NoSegDist,2,MISSING)

print(paste("There are ",length(which(Missing_inSamples > 0.90) ) ," samples missing more than 90% of data",sep=""))
if (length(which(Missing_inSamples > 0.90)) >0){
  All_NAM_poly_rm_mis<-All_NAM_poly_NoSegDist [,-c(which(Missing_inSamples > 0.9))]
}else{All_NAM_poly_rm_mis<-All_NAM_poly_NoSegDist }

# Remove samples with hete >20%
HETE<-function(dat){
  missing<-length(dat[is.na(dat)])
  AlleleAB_freq<-length(which(dat == 'HH'))/(length(dat)-missing)
  return(AlleleAB_freq)
}

Hete_inSamples<-apply(All_NAM_poly_rm_mis,2,HETE)

print(paste("There are ",length(which(Hete_inSamples > 0.20) ) ," samples that have more than 20% hetes",sep=""))

if (length(which(Hete_inSamples > 0.20)) >0){
  All_NAM_poly_rm_mis_hete<-All_NAM_poly_rm_mis [,-c(which(Hete_inSamples > 0.2))]
}else{All_NAM_poly_rm_mis_hete<-All_NAM_poly_rm_mis }


# One more time remove SNPs monomorphic and samples missing more than 90%
mono_heteCounts<-function(dat){
  GENOTYPES<-table(dat)
  # if there are more than 1 genotype after removing Hete, then it is polymorphic
  if (length(names(GENOTYPES)) <= 1){monomorphic<-'yes'}else{monomorphic<-'no'}
  return(monomorphic)
}

SNPmono1<-apply(All_NAM_poly_rm_mis_hete,1,mono_heteCounts)

print(paste("There are ",length(which(SNPmono1 == "yes") ) ," SNPs that are monomorphic",sep=""))
if (length(which(SNPmono1 == "yes")) >0){
  All_NAM_poly_rm_mis_hete_mono<-All_NAM_poly_rm_mis_hete [-c(which(SNPmono1 == "yes")),]
}else{All_NAM_poly_rm_mis_hete_mono<-All_NAM_poly_rm_mis_hete }


# Remove samples missing more than 90% of SNPs
Missing_inSamples<-apply(All_NAM_poly_rm_mis_hete_mono,2,MISSING)

print(paste("There are ",length(which(Missing_inSamples >= 0.90) ) ," samples missing more than >= 90% of data",sep=""))

if (length(which(Missing_inSamples >= 0.90)) >0){
  All_NAM_poly_rm_mis_hete_mono_mis<-All_NAM_poly_rm_mis_hete_mono [,-c(which(Missing_inSamples >= 0.9))]
}else{All_NAM_poly_rm_mis_hete_mono_mis<-All_NAM_poly_rm_mis_hete_mono }


# Output RILs file to be use for imputations
write.table(All_NAM_poly_rm_mis_hete_mono_mis,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_noSegDist_rmIndmisHete_monoSNP.txt",quote=F,row.names=T,col.names=T,sep="\t")


######## ========== GET GBS PARENTAL GENOTYPES FOR THE SNPS TO BE IMPUTED, CALLs respect to Rasmusson ===================================
# Import GBS parents
GBS<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/MNS_July2016_Parents_DP5_GQ30.recode_transformedHETE_toNA_HH_hmp.txt",header=T)

# Import Rasmusson consensus 
Rasmusson<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/ConsensusRas_gbs_exome_gbsnames_inPriHmp.txt",header=F)

# Import RILs for imputation
All_NAM_poly_rm_mis_hete_mono_mis<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_noSegDist_rmIndmisHete_monoSNP.txt",header=T,row.names=1)


# Get only SNPs included in the file for imputation.Remove Parent PI328632 as it has the wrong genotypes. We will use the exome capture data for this parent instead.
GBS_im<-GBS[(GBS[,1] %in% row.names(All_NAM_poly_rm_mis_hete_mono_mis)),-c(which(colnames(GBS) == "PI328632"))]
Rasmusson_im<-Rasmusson[(Rasmusson[,1] %in% row.names(All_NAM_poly_rm_mis_hete_mono_mis)),]

#Sort SNPs based on RILs
GBS_im_or<-GBS_im[match(row.names(All_NAM_poly_rm_mis_hete_mono_mis), GBS_im[,1]),]
Rasmusson_im_or<-Rasmusson_im[match(row.names(All_NAM_poly_rm_mis_hete_mono_mis), Rasmusson_im[,1]),]

# Make sure All SNPs match with RILs
if ( identical(as.character(GBS_im_or[,1]), as.character(Rasmusson_im_or[,1])) == FALSE)stop("Rasmusson and GBS parents have no matching SNPs")

if ( identical(as.character(GBS_im_or[,1]), as.character(row.names(All_NAM_poly_rm_mis_hete_mono_mis))) == FALSE)stop("Rasmusson and GBS parents SNPs do not match RILs")

# Combine Parents and convert calls to be respect to Rass

PARENTS<-cbind(as.data.frame(Rasmusson_im_or[,2]), GBS_im_or[,-c(1:4)])
names(PARENTS)[1]<-"Ras_consensus"
row.names(PARENTS)<-GBS_im_or[,1]


#Combine Parents with RILs (remove the ras from the RILs)
PARENTS_RILSforimputation<-cbind(PARENTS,All_NAM_poly_rm_mis_hete_mono_mis[,-c(1)])

# Output RILs with Parents file to be use for imputations
write.table(PARENTS_RILSforimputation,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_noSegDist_rmIndmisHete_monoSNP_wPARENTS.txt",quote=F,row.names=T,col.names=T,sep="\t")

