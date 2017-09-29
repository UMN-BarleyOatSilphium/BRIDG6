# Author: Ana Poets
# Description: identify QTL at different Mbp distance and QTL statistics
####################################################################################################
rm(list=ls())

library(stringr)
# Determine length (bp) before and after a significant SNP to be defined as one QTL
LENGTH<-5000000

# Get output after GWAS

Polytest<-read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/my_gwas_PolyTest_out_80miss_BLUE.csv",header=T,row.names=1)

# SNP names
SNPnames<-read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/my_gwas_SNPs_out_80miss_BLUE.csv")

# Order of family names
Families<-read.csv("/Users/agonzale/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Input/phenos_LDKNNI.csv",sep=",")

# Import parent names
Parents_famNam<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/Pop_structure_propDiss_JULY2016_parentLessMissingData_grep.txt",header=T)

# Sort by family_grep
Parents_famNam_or<-Parents_famNam[order(Parents_famNam$Family_grep),]
# Add Parent name to each individual
for (i in 1:dim(Parents_famNam_or)[1]){
  Families[grep(Parents_famNam_or[i,"Family_grep"],Families$line_name),1]<-as.character(Parents_famNam_or[i,"Parent_lessMissingData"])
}


#select unique entries for total 88 families
Family_ID_uniq<-Families[!duplicated(Families$family),]

#sort by family entry order
Family_ID_uniq_or<-Family_ID_uniq[order(Family_ID_uniq[,"family"]),]

# Add parents names to the allele effect columns
names(Polytest)[c(((dim(Polytest)[2])-88 +1):dim(Polytest)[2]) ] <-Family_ID_uniq_or[,1]
# Add SNP names to the GWAS output

row.names(Polytest)<-sub("X","",SNPnames[,2])

# Fix physical positions genome-wide
Positions_phy_mbp <-as.data.frame(as.character(row.names(Polytest)))
Physical_positions <-  apply(Positions_phy_mbp,1, function(x) strsplit(as.character(x), "_")[[1]][2])
Physical_positions <-as.data.frame(Physical_positions)
row.names(Physical_positions)<-Positions_phy_mbp[,1]

for (i in 1:7){
  CHR1<-grep(paste(i,"H1", sep=""),row.names(Physical_positions))
  CHR2<-grep(paste(i,"H2", sep=""),row.names(Physical_positions))
  assign(paste("CHR",i,"_1",sep=""), CHR1)
  assign(paste("CHR",i,"_2",sep=""), CHR2)
}


List_chr_parts<-c("CHR1_1","CHR1_2","CHR2_1", "CHR2_2", "CHR3_1", "CHR3_2", "CHR4_1", "CHR4_2", "CHR5_1", "CHR5_2", "CHR6_1", "CHR6_2", "CHR7_1", "CHR7_2")

# Number of bp to add to each part to have a cumulative physical position
ADD_CHR<-c(0,312837513,558535432,952068106,1326610456,1720921089,2026321570,2381382776,2673381728,3054247210,3343411888,3638233958,3926792401,4252589917)

#Calculate joint physical positions for both parts of a chromosome, then calculate Cumulative positons genome wide
Positions_original<-Physical_positions[,1]

New_order_phy<-NULL
for (C in 1:14){
  CHR_pos<-as.numeric(as.character(Physical_positions[get(List_chr_parts[C]),1]))
  NEW_chr_pos<-CHR_pos + ADD_CHR[C]
  New_order_phy <-c(New_order_phy ,NEW_chr_pos)
}  
Positions_phy_mbp_cumul<-cbind(Positions_phy_mbp, New_order_phy)

# Add cumulative positions and a columns of NA for QTL assignation
Polytest_pos<-cbind(Positions_phy_mbp_cumul[,2],rep(NA,dim(Polytest)[1]),Polytest)
names(Polytest_pos)[1]<-"Cumul_Pos"
names(Polytest_pos)[2]<-"QTL_assig"

### Assign QTL

# Get threshold for significant level
THR =  -log10(0.05 / ( nrow(Polytest) * (1-0.05))) 

# Work one chromosome at the time

for (n in 1:7){
  print(n)
  # the first QTL is set to 2_1, use to digit numbers so 2.1 and 2.10 are differentiated
  x<-1
  Polytest_pos_chr<-Polytest_pos[grep(paste(n,"H",sep=""),row.names(Polytest_pos)),]
  
  # Get a vector of positions of p-values that are significant
  Pos_sig_pval<-Polytest_pos_chr[which(Polytest_pos_chr$pval >= THR),"Cumul_Pos"]
  
  if(length(Pos_sig_pval)>0){
    # Identify the position of the first significant marker and the position of +/- LENGTH
    Bottom<-Pos_sig_pval[1] - LENGTH
    Top<-Pos_sig_pval[1] + LENGTH
    
    START_QTL<-Pos_sig_pval[length(Pos_sig_pval)]
    while(!is.na(Bottom <= (Pos_sig_pval[length(Pos_sig_pval)]))){
      ### If the positions between Bottom and Top are NA for QTL assignation, then assign a new QTL. Otherwise, assign whichever QTL already exist for the first SNP.
      # SNPs in QTL
      SNP_qtl<-which(Polytest_pos_chr$Cumul_Pos >= Bottom & Polytest_pos_chr$Cumul_Pos <= Top)
      if (is.na(Polytest_pos_chr[SNP_qtl[1],"QTL_assig"])){QTLnumber<-(paste(n,"_",str_pad(x, 2, pad = "0"),sep="")) ; x<-(x+1)}else{QTL_assig<-Polytest_pos_chr[SNP_qtl[1],"QTL_assig"]}
      Polytest_pos_chr[SNP_qtl,"QTL_assig"]<-QTLnumber
      # move to the next significant p-value outside the last QTL
      Pos_sig_pval_new<-Pos_sig_pval[which(Pos_sig_pval >Top)[1]]
      
      Bottom<-Pos_sig_pval_new -LENGTH
      Top<-Pos_sig_pval_new + LENGTH
    }
  }else{Polytest_pos_chr<-Polytest_pos_chr}
  assign(paste("CHROM_",n,sep=""),Polytest_pos_chr)
  
}

GenomeWide_QTLassigned<-rbind(CHROM_1, CHROM_2,CHROM_3,CHROM_4,CHROM_5, CHROM_6, CHROM_7)

#write.table(GenomeWide_QTLassigned,paste("/Users/agonzale/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/QTLassignation_ANA/GenomeWide_QTLassigned_",LENGTH,".xls",sep=""),quote=F,row.names=T,col.names=T,sep="\t")

#########======== Table 1: average p-value for QTL and the total number of SNP in that QTL ============#####################

# List of the number of QTL assigned
QTL_all<-names(table(GenomeWide_QTLassigned[,2]))

QTLdescriptionTable<-matrix(NA,ncol=6,nrow=length(QTL_all))
colnames(QTLdescriptionTable)<-c("QTL","Maximum -log(p)","Number of significant SNP","Total number of SNP in QTL region (5Mbp)","Max. Families segregating","Min. Families Segregating")
for (i in 1:length(QTL_all)){
  p_values<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned$QTL_assig == QTL_all[i]),"pval"]
  QTLdescriptionTable[i,1]<-QTL_all[i]
  QTLdescriptionTable[i,2]<-round(max(p_values),2)
  QTLdescriptionTable[i,3]<-length(which(p_values >= THR))
  QTLdescriptionTable[i,4]<-length(p_values)
  
  # Get the max and minimum number of families for which a Significant SNP segregates in a QTL
  data<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned$QTL_assig == QTL_all[i]),]
  # get only significant SNPs in the QTL
  data_sig<-data[which(data$pval >= THR),(dim(data)[2]-87):dim(data)[2]]
  Segregation<-apply(data_sig,1, function(x) length(which(x >0)))
  QTLdescriptionTable[i,5]<-max(Segregation)
  QTLdescriptionTable[i,6]<-min(Segregation)
  
}

#write.table(QTLdescriptionTable,"/Users/agonzale/Documents/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/QTLassignation_ANA/Summary_QTL_max.xls",quote=F,row.names=F,col.names=T,sep="\t")


####### ====== visualization of the QTL location ==========================================================================
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

# Plot results
pdf("/Users/agonzale/Documents/SmithLab/NAM/write/Bridgs_edited_files_v2/ExtraFigs/GWAS_Genes.pdf",height=7, width=11)
par(mar=c(5,4,4,4))
#plot(GenomeWide_QTLassigned$Cumul_Pos/1000000, GenomeWide_QTLassigned$pval, cex=0.6, ylab="-log(p-value)",xlab="Physical Position (Mbp)",xaxt="n")
plot(GenomeWide_QTLassigned$Cumul_Pos/1000000, GenomeWide_QTLassigned$pval, cex=0.6, ylab="-log(p-value)",xlab="Physical Position (Mbp)")

abline(h=THR, col="red")

# Color each QTL identified with a different color
UNIQ_QTL<-unique(GenomeWide_QTLassigned[which(!is.na(GenomeWide_QTLassigned[,2])),2])
# create nQTL colors
palette(rainbow(length(UNIQ_QTL)))
#COLORS<-(palette(gray(seq(0,.9,len = 25)))) 
COLORS<-rep(c("blue","red"), 12)
# paint each QTL with a different color

for (p in 1:length(UNIQ_QTL)){
  QTLpoints<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned[,2] == UNIQ_QTL[p]),]
  points(QTLpoints$Cumul_Pos/1000000, QTLpoints$pval, cex=0.6, col=COLORS[p])
}

# Set chromosome boundaries
abline(v=c(558535432,1326610456,2026321570,2673381728,3343411888,3926792401)/1000000, col="gray", lty=2)
legend("topright",col=c("gray","black"),pch=c(NA,1),lty=c(2,NA),legend=c("Chromosome boundary","SNP"))
dev.off()
print(paste("Total number of QTL:",length(UNIQ_QTL), sep=" "))

# add the position of known genes associated with flowering time
axis(side=1, at=c(Annotations_genes$Cumul_posMid/1000000), labels = FALSE)
text(x=c(Annotations_genes$Cumul_posMid/1000000),  par("usr")[3], 
     labels = Annotations_genes$Abrevation, srt = 90, xpd = T, cex=0.5)

## qtlselect are the genes idetified in Table 2
axis(side=1, at=c(qtlselect/1000000), labels = FALSE)
text(x=c(qtlselect/1000000),  par("usr")[3],     labels = Table_SNPnearest[,5], srt = 90, xpd = T, cex=0.5)
