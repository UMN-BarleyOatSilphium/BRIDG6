library(ape)

DATA<-read.table("/home/smithkp/agonzale/RemakeNAMgetDist/Final_filt_LDKNNI_noMis_noNAHHras_rasBased_MAF5_noLD.txt",header=T,row.names=1)
dim(DATA)
DATA[1:10,1:10]
SNP_gwas<-row.names(DATA)

write.table(SNP_gwas,"/home/smithkp/agonzale/RemakeNAMgetDist/SNP_gwas.txt",quote=F,row.names=F,col.names=F,sep="\t")
# Determine the gentic distance between each donor parent and Rasmusson
SNP_gwas<-read.table("/home/smithkp/agonzale/RemakeNAMgetDist/SNP_gwas.txt",header=FALSE)
QCdata<-read.table("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Input_forPriyanka/Input_allNuc_noRasBased/NAM_MNS_July2016_DP5_GQ30_mis80_6336ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_NOrasBased_Polymorphic_rm100SNPclose_NUC_naExcessHH_noSegDist_rmIndmisHete_monoSNP_wPARENTS.txt",header=T,row.names=1)

SNPsAfterQCforGWASwParents<-subset(QCdata,(row.names(QCdata) %in% SNP_gwas[,1]))

# import list of parents
Parents<-read.table("/home/smithkp/agonzale/RemakeNAMgetDist/Parents.txt")

SNPsAfterQCforGWASonlyParents<-SNPsAfterQCforGWASwParents[,(colnames(SNPsAfterQCforGWASwParents) %in% Parents[,1])]

# Import family Information
FamInfo<-read.csv("/home/smithkp/agonzale/RemakeNAMgetDist/Family_summary.csv")

################################################################################################
###### ----- FUNCTIONS -----####################################################################

# Convert to numeric: AA = 1, CC=2,TT=3,GG=4, NA=NA, all hete = 5
ConvertNum<-function(dat){
  dat[which(dat == "AA")]<-1
  dat[which(dat == "CC")]<-2
  dat[which(dat == "TT")]<-3
  dat[which(dat == "GG")]<-4
  dat[which(dat == "HH")]<- NA
  return(dat)
}

################################################################################################
# Genetic distance between each Donor parent and Rasmusson
#################################################################################################
DifferencesBtwParents<-matrix(NA, ncol=4,nrow= dim(SNPsAfterQCforGWASonlyParents)[2])
colnames(DifferencesBtwParents)<-c("ParentName","SNPcompared","SNPdiff","PropDiff")

for (d in 2:dim(SNPsAfterQCforGWASonlyParents)[2]){
  Parents_pair<-SNPsAfterQCforGWASonlyParents[,c(1,d)]
  Parents_pair<-as.data.frame(Parents_pair)
  
  DifferencesBtwParents[(d-1),1]<-names(Parents_pair)[2]
  
  #Remove SNPs that are missing in both or either parents
  BOTH_NA<-which(is.na(Parents_pair[,1]) | is.na(Parents_pair[,2]))
  if (length(BOTH_NA) >0){
    Parents_present<-Parents_pair[-c(BOTH_NA),]}else{Parents_present<-Parents_pair}
  
  # Find which SNPs are HETE in either parent or both parents
  BOTH_HH<-which(Parents_present[,1] == 'HH' | Parents_present[,2] == 'HH')
  if(length(BOTH_HH) >0){
    Parents_noBothHete<-Parents_present [-c(BOTH_HH),]}else{Parents_noBothHete<-Parents_present }
  
  
  GoodSNPs<-row.names(Parents_noBothHete)  # Number of SNPs compared
  
  DifferencesBtwParents[(d-1),2]<-length(GoodSNPs)
  
  ###=========convert to numeric to use correlation function ============
  dataConverted <- apply(Parents_noBothHete, 2, ConvertNum)
  DIFFERENCES_count<-length(which(dataConverted[,1] != dataConverted[,2]))
  
  DifferencesBtwParents[(d-1),3]<-DIFFERENCES_count
  
  DIFFERENCES<-dist.gene(t(Parents_noBothHete), method="percentage", pairwise.deletion = T, variance = FALSE)
  # amount differences based on sites compared
  DifferencesBtwParents[(d-1),4]<- round(DIFFERENCES, digits=3)
  
}

## add parent PI328632 from exome capture
HR655_pi<-read.table("/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter/Rasmusson_filtered/Consensus/Exome_PI328632.txt",header=F,row.names=1)
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

HR655parent<-subset(DATA,(row.names(DATA) %in% SNP_gwas[,1]))

sharedRas<-intersect(row.names(SNPsAfterQCforGWASonlyParents),row.names(HR655parent))
SNPwHR655<-SNPsAfterQCforGWASonlyParents[(row.names(SNPsAfterQCforGWASonlyParents) %in% sharedRas),]
HR655parent_shared<-subset(HR655parent,(row.names(HR655parent) %in% sharedRas))

# match order in Ras
HR655parent_shared<-cbind(row.names(HR655parent_shared),as.character(HR655parent_shared[,1]))
HR655parent_shared_or<-HR655parent_shared[match(row.names(SNPwHR655),as.character(HR655parent_shared[,1])),]
if (identical(as.character(HR655parent_shared_or[,1]),as.character(row.names(SNPwHR655)) ) == FALSE)stop("SNP order don't match")

PairParents<-cbind(HR655parent_shared_or,as.character(SNPwHR655[,1]))
colnames(PairParents)<-c("SNP","PI328032","Ras_consensus")
row.names(PairParents)<-PairParents[,1]
PairParents<-PairParents[,-c(1)]
Parents_pair<-PairParents

#### Estimate the percentage differentation from Rasmusson
DifferencesBtwParents[dim(DifferencesBtwParents)[1],1]<-"PI328632"
#Remove SNPs that are missing in both or either parents
BOTH_NA<-which(is.na(Parents_pair[,1]) | is.na(Parents_pair[,2]))
if (length(BOTH_NA) >0){
  Parents_present<-Parents_pair[-c(BOTH_NA),]}else{Parents_present<-Parents_pair}

# Find which SNPs are HETE in either parent or both parents
BOTH_HH<-which(Parents_present[,1] == 'HH' | Parents_present[,2] == 'HH')
if(length(BOTH_HH) >0){
  Parents_noBothHete<-Parents_present [-c(BOTH_HH),]}else{Parents_noBothHete<-Parents_present }


GoodSNPs<-row.names(Parents_noBothHete)  # Number of SNPs compared

DifferencesBtwParents[dim(DifferencesBtwParents)[1],2]<-length(GoodSNPs)

###=========convert to numeric to use correlation function ============
dataConverted <- apply(Parents_noBothHete, 2, ConvertNum)
DIFFERENCES_count<-length(which(dataConverted[,1] != dataConverted[,2]))

DifferencesBtwParents[dim(DifferencesBtwParents)[1],3]<-DIFFERENCES_count

DIFFERENCES<-dist.gene(t(Parents_noBothHete), method="percentage", pairwise.deletion = T, variance = FALSE)
# amount differences based on sites compared
DifferencesBtwParents[dim(DifferencesBtwParents)[1],4]<- round(DIFFERENCES, digits=3)

write.table(DifferencesBtwParents,"/home/smithkp/agonzale/RemakeNAMgetDist/OUTPUT_DPvsRas_genDist.txt",quote=F,row.names=F,col.names=T,sep="\t")

################################################################################################
############### -------- Count the number of segregating SNPs in each family ----###############
################################################################################################
DATA<-SNPsAfterQCforGWASwParents
# Import family Information
FamInfo<-read.csv("/Users/Mia/Downloads/Family_summary.csv")

PolyTest<-function(dat){
  if (dim(table(dat))<=1){test<-"mono"}else{test<-"poly"}
  return(test)
}


PolyFamResults <-matrix(NA,ncol=4 ,nrow = dim(FamInfo)[1] )
colnames(PolyFamResults) <-c("NAM_name","FamilyNo","No.RILsGWAS","NoSNPpoly")
for (x in 1:dim(FamInfo)[1]){
  #DiffRas <-DifferencesBtwParents[(which(DifferencesBtwParents[,1] == as.character(FamInfo[x,2]))),4]
  #PolyFamResults[x,4]<-DiffRas
  PolyFamResults[x,1]<- as.character(FamInfo[x,2])
  PolyFamResults[x,2]<- as.character(FamInfo[x,3])
  
  Fam<-FamInfo[x,10]
  #get genotypes for all members of this family
  FamGen<- DATA[,(grep(Fam, colnames(DATA)))]
  
  PolyFamResults[x,3]<-dim(FamGen)[2]
  
  # Assess which markers are monomorphic or polymorphic
  Assessing<-apply(FamGen,1,PolyTest)
  if (length(which(Assessing == "poly"))>0){
    CountPoly <-length(which(Assessing == "poly"))
  }else{CountPoly <- 0 }
  
  PolyFamResults[x,4]<-CountPoly
}

# Add the information for the families that were combined
FamCombined= c("HR632S|HR651S","HR648S|HR649S|HR650S")

PolyFamResults2 <-matrix(NA,ncol=4 ,nrow =length(FamCombined) )
colnames(PolyFamResults2) <-c("NAM_name","FamilyNo","No.RILsGWAS","NoSNPpoly")
)

for (x in 1:length(FamCombined)){
  PolyFamResults2[x,1]<-NA
  PolyFamResults2[x,2]<- as.character(FamCombined[x])
  
  Fam<-as.character(FamCombined[x])
  #get genotypes for all members of this family
  FamGen<- DATA[,(grep(Fam, colnames(DATA)))]
  
  PolyFamResults2[x,3]<-dim(FamGen)[2]
  
  # Assess which markers are monomorphic or polymorphic
  Assessing<-apply(FamGen,1,PolyTest)
  if (length(which(Assessing == "poly"))>0){
    CountPoly <-length(which(Assessing == "poly"))
  }else{CountPoly <- 0 }
  
  PolyFamResults2[x,4]<-CountPoly
}

AllResults<-rbind(PolyFamResults,PolyFamResults2)

write.table(AllResults,"/home/smithkp/agonzale/RemakeNAMgetDist/FamilySNPsegregating.txt",quote=F,row.names=F,col.names=T,sep="\t")

