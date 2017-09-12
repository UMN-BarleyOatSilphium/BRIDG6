# Author: Ana Poets
# Description: QTL statistics
# 1. number of families segregating where a significant marker segregates
# 2. QTL physical positions
# 3. Get largest and lowest allele effect at QTL 2.1
# 4. Table 2.2 : Genes within a QTL

########################################################################################################################################################################################################
rm (list =ls())
GenomeWide_QTLassigned<-read.table("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/QTLassignation_ANA/GenomeWide_QTLassigned_5e+06.xls",header=T)
# Get threshold for significant level
THR =  -log10(0.05 / ( nrow(GenomeWide_QTLassigned) * (1-0.05))) 

### 1. Determine how many families segregate for each of the significant variants
# Table of significant SNPs and their allele effects per family
SIG<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned$pval >= THR),c(((dim(GenomeWide_QTLassigned)[2])-87):dim(GenomeWide_QTLassigned)[2])]

# count how many families have an effect for the allele, if there's an effect estimated it means that the variant was segregating in that family.
SEG_Fam<-NULL
for(s in 1:dim(SIG)[1]){  
  SEG_Fam<-c(SEG_Fam,length(which(SIG[s,] != 0)))
}

# count how many families segreagate at each QTL
# Table of significant SNPs and their allele effects per family
SIGqtl<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned$pval >= THR),c(2,((dim(GenomeWide_QTLassigned)[2])-87):dim(GenomeWide_QTLassigned)[2])]

# List of the number of QTL assigned
QTL_all_sig<-names(table(SIGqtl[,1]))

#minimum number of families segregating for each QTL
QTLsegFamMin<-matrix(NA,ncol=2,nrow=length(QTL_all_sig))
for (q in 1:length(QTL_all_sig)){
  EFFECTS<-SIGqtl[(grep(QTL_all_sig[q],SIGqtl[,1])),]
  SEG_Fam_qtl<-NULL
  for(s in 1:dim(EFFECTS)[1]){  
    SEG_Fam_qtl<-c(SEG_Fam_qtl,length(which(EFFECTS[s,-1] != 0)))
  }
  QTLsegFamMin[q,1]<-QTL_all_sig[q]
  QTLsegFamMin[q,2]<-min(SEG_Fam_qtl)
}

# Report in manuscript
summary(as.numeric(as.matrix(QTLsegFamMin[,2])))

### 2. Interval (Mbp) representing each QTL
SIGqtl_size<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned$pval >= THR),c(1:2)]

INTERVAL<-matrix(NA,ncol=2,nrow=length(QTL_all_sig))

for (q in 1:length(QTL_all_sig)){
  SNPinQTL<-SIGqtl_size[(grep(QTL_all_sig[q],SIGqtl_size[,2])),]
  start_pos<-SNPinQTL[1,1]
  end_pos<-SNPinQTL[dim(SNPinQTL)[1],1]
  interval<-end_pos-start_pos
  INTERVAL[q,1]<-QTL_all_sig[q]
  INTERVAL[q,2]<-interval
}


#Replace QTL that had only one SNP for 5Mbp
INTERVAL[which(INTERVAL[,2]==0),2]<-c("5000000")
summary(as.numeric(as.matrix(INTERVAL[,2])))

# max and min allele effect at each significant marker
# Table of significant SNPs and their allele effects per family
SIGqtl<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned$pval >= THR),c(((dim(GenomeWide_QTLassigned)[2])-87):dim(GenomeWide_QTLassigned)[2])]

MaxMin_alleleEff<-matrix(NA, ncol=3,nrow=dim(SIGqtl)[1])
MaxMin_alleleEff[,1]<-row.names(SIGqtl)

for (q in 1:(dim(SIGqtl)[1])){
  MaxMin_alleleEff[q,2]<-min(SIGqtl[q,])
  MaxMin_alleleEff[q,3]<-max(SIGqtl[q,])
}

plot(MaxMin_alleleEff[,2],col="white", ylim=c(as.numeric(min(MaxMin_alleleEff[,2]))-5,as.numeric(max(MaxMin_alleleEff[,3]))))
lines(MaxMin_alleleEff[,2],col="blue")
lines(MaxMin_alleleEff[,3],col="red")


#########################################################################################################################################################################################################
### 3. Get largest and lowest allele effect at QTL 2.1
SIGqtl<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned$pval >= THR),c(2,((dim(GenomeWide_QTLassigned)[2])-87):dim(GenomeWide_QTLassigned)[2])]
qtl2_1<-SIGqtl[which(SIGqtl$QTL_assig == "2_1"),-1]
# number families with positive effects ast QTL 2.1
QTLmax_min_fam<-matrix(NA,ncol=2,nrow=88)
row.names(QTLmax_min_fam)<-colnames(qtl2_1)
for (f in 1:88){
  QTLmax_min_fam[f,1]<-max(qtl2_1[,f])
  QTLmax_min_fam[f,2]<-min(qtl2_1[,f])
}
# allele effect for a family is whichever (positive or negative) value is the largest absolute value. If both positive and negative effects are the same magnitude, chosee the positive one

Allele_effect<-function(dat){
  abs_val<-max(abs(dat[1]),abs(dat[2]))
  if(abs(dat[1]) == abs_val & abs(dat[2]) != abs_val){alleleEff<-dat[1]}
  if(abs(dat[2]) == abs_val & abs(dat[1]) != abs_val){alleleEff<-dat[2]}
  if(abs(dat[1]) == abs_val & abs(dat[2]) == abs_val){alleleEff<-dat[1]}
  return(alleleEff)
} 

QTLeffect<-as.data.frame(apply(QTLmax_min_fam,1,Allele_effect))
QTLmax_min_fam_effect<-cbind(QTLmax_min_fam,QTLeffect)
names(QTLmax_min_fam_effect)[3]<-"AllelesEffectQTL"

#families with largest positive effects at qtl2.1
length(which(QTLmax_min_fam_effect[,3] > 0)) #43

#families largest negative at qtl2.1
length(which(QTLmax_min_fam_effect[,3] < 0)) #45

###############################################################################################################################################################################
### 4. Table 2.2 Genes close or in QTL

# Import the position of known genes in barley
Annotations<-read.csv("~/Dropbox/GITHUB/BRIDG6/Datasets/Annotations/FT_annotations_select.csv",sep=",")

#Remove Genes in unknown chromosome as the positions of these SNPs are random to the chromosome.

Annotations_chr<-Annotations[-grep("UN", Annotations$Gene_Ch),]

## Get genome-wide cumulative positions
# cumulative positions until the end of each chromosome to be added to each chr. 
CHROM_length<-c(0,558535432,1326610456,2026321570,2673381728,3343411888,3926792401,4584016401)

Annotations_chr$Cumul_posStart<-rep(NA,dim(Annotations_chr)[1])
Annotations_chr$Cumul_posEnd<-rep(NA,dim(Annotations_chr)[1])
Annotations_chr$Cumul_posMid<-rep(NA,dim(Annotations_chr)[1])

for (c in 1:7){
  #cumul start
  CumPos_start<-Annotations_chr[grep(paste(c,"H",sep=""), Annotations_chr$Gene_Ch), 5] + CHROM_length[c]
  Annotations_chr[grep(paste(c,"H",sep=""), Annotations_chr$Gene_Ch),"Cumul_posStart"]<-CumPos_start
  #cumul end
  CumPos_end<-Annotations_chr[grep(paste(c,"H",sep=""), Annotations_chr$Gene_Ch), 6] + CHROM_length[c]
  Annotations_chr[grep(paste(c,"H",sep=""), Annotations_chr$Gene_Ch),"Cumul_posEnd"]<-CumPos_end
  #cumul mid point
  CumPos_mid<-Annotations_chr[grep(paste(c,"H",sep=""), Annotations_chr$Gene_Ch), 7] + CHROM_length[c]
  Annotations_chr[grep(paste(c,"H",sep=""), Annotations_chr$Gene_Ch),"Cumul_posMid"]<-CumPos_mid
  
}

# Summary annotations
Annotations_genes<-Annotations_chr[!is.na(Annotations_chr$Cumul_posStart),c(1,2,4,16,17,18)]

# Sumary QTL start and end positons
QTLsummary<-matrix(NA,ncol=2, nrow=length(QTL_all_sig))
row.names(QTLsummary)<-QTL_all_sig
colnames(QTLsummary)<-c("start","End")
for (z in 1:length(QTL_all_sig)){
    QTL_snps<-GenomeWide_QTLassigned[grep(QTL_all_sig[z],GenomeWide_QTLassigned$QTL_assig),1]
    if(length(QTL_snps) >0){
    QTLsummary[z,1]<-QTL_snps[1]
    
    QTLsummary[z,2]<-QTL_snps[length(QTL_snps)]
    }else{QTLsummary[z,1] <-QTL_snps[1] ; QTLsummary[z,2] <-QTL_snps[1] }
}


# Find the genes in the QTL and 1-5 Mbp from either side of a QTL 
DISTANCE<-5000000

#Table SNP near genes
Table_SNPnearest<-matrix(NA,nrow=1,ncol=10)
colnames(Table_SNPnearest)<-c("QTL","SNP","SNP -log(p)","QTL -log(p)","Nearest candidate gene","Distance to SNP","Chromosome","Minimum allele Effect","Maximum allele effect","Segregating Families")

for (k in 1:dim(QTLsummary)[1]){
  GENE<-Annotations_genes[which(Annotations_genes$Cumul_posMid >= (QTLsummary[k,"start"]- DISTANCE) &Annotations_genes$Cumul_posMid <= (QTLsummary[k,"End"] + DISTANCE)),]
  
  #find the SNP closest to each of the found genes
  if(dim(GENE)[1]>0){
    for (g in 1:dim(GENE)[1]){
      print (g)
      SNPnearest<-(GenomeWide_QTLassigned[which(abs(GenomeWide_QTLassigned[,1]-GENE[g,6])==min(abs(GenomeWide_QTLassigned[,1]-GENE[g,6]))),])
      # QTLname
      QTLpresent<-(GenomeWide_QTLassigned[which(abs(GenomeWide_QTLassigned[,1]-GENE[1,6])==min(abs(GenomeWide_QTLassigned[,1]-GENE[1,6]))),2])
      #Maximum lod score for any SNP in the QTL
      QTLinfo<-GenomeWide_QTLassigned[grep(QTLpresent, GenomeWide_QTLassigned$QTL_assig),]
      MAXLOD<-max(QTLinfo$pval)
      #LOD score for nearest SNP 
      SNP_LOD<-SNPnearest$pval
      #how close is the SNP to the mid point of the gene
      DISTANCE<-abs(GenomeWide_QTLassigned[,1]-GENE[g,6])[which(abs(GenomeWide_QTLassigned[,1]-GENE[g,6])==min(abs(GenomeWide_QTLassigned[,1]-GENE[g,6])))]
      #Abbreviation nearest gene
      NEAREST_GENE<-GENE[g,2]
      #chromosome. Keep the information from the first SNP in the QTL
      CHR<-GENE[1,3]
      # number of SNPs in QTL region
      SNP_region<-dim(QTLinfo)[1]
      # minumul allele effect for the SNP
      MIN_alleleEff<-min(SNPnearest[,c((dim(SNPnearest)[2] - 87):dim(SNPnearest)[2])])
      MAX_alleleEff<-max(SNPnearest[,c((dim(SNPnearest)[2] - 87):dim(SNPnearest)[2])])
      # number of families segregating for the allele
      NumFAM_SEG<-length(which(SNPnearest[,c((dim(SNPnearest)[2] - 87):dim(SNPnearest)[2])] != 0))
      
      #construct table
     Table_part<-data.frame(QTLpresent,row.names(SNPnearest), SNP_LOD,MAXLOD,NEAREST_GENE,round(DISTANCE/1000000,digits=3),CHR,MIN_alleleEff, MAX_alleleEff, NumFAM_SEG )
     names(Table_part)<-c("QTL","SNP","SNP -log(p)","QTL -log(p)","Nearest candidate gene","Distance to SNP","Chromosome","Minimum allele Effect","Maximum allele effect","Segregating Families")
     
     Table_SNPnearest<-rbind(Table_SNPnearest,Table_part)
    }
  }
  #print (paste( row.names(QTLsummary)[k], ":" ,as.character(GENE)))
}

# Remove first row of NA
Table_SNPnearest<-Table_SNPnearest[-1,]
write.table(Table_SNPnearest,"/Users/agonzale/Dropbox/SmithLab/NAM/write/Bridge_edited_files/Tables/Table2_2_Ana_qtlFloweringTime.xls",quote=F,row.names=F,col.names=T,sep="\t")
