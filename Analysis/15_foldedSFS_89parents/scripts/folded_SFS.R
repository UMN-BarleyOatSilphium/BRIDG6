# Generate a folded site frequency spectrum for 9k and exon-capture data.

# Read the data files

GBS <- read.delim(file="/Users/mingqinshao/Downloads/maf0.003_meanDP5_psedumolecular_parents_MNS_July2016_production.recode.maf",header = T, sep="\t")
GBS$account_minor <- GBS$sample_NB*GBS$MAF

head(GBS)


# Bin them up in 2.5% frequency bins
breaks <- seq(0, 44, by=1)
GBS.sfs <- 100*table(cut(as.integer(GBS$account_minor), breaks=breaks, labels = FALSE))/nrow(GBS)
nrow(GBS.sfs)
GBS.sfs
# Make plot

pdf(file="~/Downloads/Folded_SFS_maf0.003_meanDP5__psedumolecular_parents_MNS_July2016_production_88.pdf", width=18, height=12)
par(mar=c(5,5.5,4,2))
barplot(t(GBS.sfs),
        beside=TRUE,
        col=c("gray"),
        cex.lab=1.4, 
        cex.axis=1.2,
        xlab="Minor Allele Count",
        ylab="Percentage SNPs",
        cex.main=1.4,
        ylim=c(0,5),
        main="Folded SFS using GBS data from 89 parents")
#barplot(t(exonSNP.sfs),
#        beside=TRUE,
#        col=c("blue"),
 #       xlab="Minor Allele Frequency",
 #       ylab="Proportion",
 #       cex.lab=1.2, 
 #       cex.axis=1.2,
  #      ylim=c(0,0.5),
  #      cex.main=1.4,
#        main="Folded SFS in exon captured SNP")

dev.off()

