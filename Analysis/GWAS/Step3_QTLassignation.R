# Author: Ana Poets
# Desktop: identify QTL at different Mbp distance
####################################################################################################
rm(list=ls())

# Determine length (bp) before and after a significant SNP to be defined as one QTL
LENGTH<-3000000

# Get output after GWAS

Polytest<-read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/my_gwas_PolyTest_out_80miss_BLUE.csv",header=T,row.names=1)

# SNP names
SNPnames<-read.csv("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/my_gwas_SNPs_out_80miss_BLUE.csv")

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
  # the first QTL is set to 1
  x<-1
  Polytest_pos_chr<-Polytest_pos[grep(paste(n,"H",sep=""),row.names(Polytest_pos)),]
  
  # Get a vector of positions of p-values that are significant
  Pos_sig_pval<-Polytest_pos_chr[which(Polytest_pos_chr$pval >= THR),"Cumul_Pos"]

  if(length(Pos_sig_pval)>0){
  # Identify the position of the first significant marker and the position of +/- LENGTH
  Bottom<-Pos_sig_pval[1] - LENGTH
  Top<-Pos_sig_pval[1] + LENGTH

  while(!is.na(Bottom <= (Pos_sig_pval[length(Pos_sig_pval)]))){
    ### If the positions between Bottom and Top are NA for QTL assignation, then assign a new QTL. Otherwise, assign whichever QTL already exist for the first SNP.
    # SNPs in QTL
    SNP_qtl<-which(Polytest_pos_chr$Cumul_Pos >= Bottom & Polytest_pos_chr$Cumul_Pos <= Top)
    if (is.na(Polytest_pos_chr[SNP_qtl[1],"QTL_assig"])){QTLnumber<-paste(n,".",x,sep="") ; x<-(x+1)}else{QTL_assig<-Polytest_pos_chr[SNP_qtl[1],"QTL_assig"]}
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
  
write.table(GenomeWide_QTLassigned,"~/Desktop/GenomeWide_QTLassigned.xls",quote=F,row.names=T,col.names=T,sep="\t")
  

# Plot results

plot(GenomeWide_QTLassigned$Cumul_Pos/1000000, GenomeWide_QTLassigned$pval, cex=0.6, ylab="-lob(p-value)",xlab="Physical Position (Mbp)")
abline(h=THR, col="red")

# Color each QTL identified with a different color
UNIQ_QTL<-unique(GenomeWide_QTLassigned[which(!is.na(GenomeWide_QTLassigned[,2])),2])
# create nQTL colors
palette(rainbow(length(UNIQ_QTL)))
COLORS<-(palette(gray(seq(0,.9,len = 25)))) 

# paint each QTL with a different color

for (p in 1:length(UNIQ_QTL)){
  QTLpoints<-GenomeWide_QTLassigned[which(GenomeWide_QTLassigned[,2] == UNIQ_QTL[p]),]
  points(QTLpoints$Cumul_Pos/1000000, QTLpoints$pval, cex=0.6, col=COLORS[p])
}

# Set chromosome boundaries
abline(v=c(558535432,1326610456,2026321570,2673381728,3343411888,3926792401)/1000000, col="gray", lty=2)

print(paste("Total number of QTL:",length(UNIQ_QTL), sep=" "))
  
  