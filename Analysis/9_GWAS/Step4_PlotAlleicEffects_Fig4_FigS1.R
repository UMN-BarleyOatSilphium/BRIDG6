# Author: Alex Ollhoff
# Co-author:Ana Poets
# Description: Plot alellic effects distribution for each QTL (Figures 4 and S1)
###############################################################################################################################################################################
rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape2)
# Look at allele effects
Fam_QTL_vis <- read.table("~/Dropbox/SmithLab/NAM/Analysis/WholeNAM_80mis_6060ind/Imputed_Output/Filtered_maf_mis_LD/GWAS/forGWAS/LDKNNI/Output/QTLassignation_ANA/GenomeWide_QTLassigned_5e+06.xls", header = T)

# Get threshold for significant level
THR =  -log10(0.05 / ( nrow(Fam_QTL_vis) * (1-0.05))) 

# Filter for SNPs and effects 
Few <- Fam_QTL_vis [,c(1,2,7, (dim(Fam_QTL_vis)[2]-87) : (dim(Fam_QTL_vis)[2]))]

# List of QTL names
QTLnames<-names(table(Fam_QTL_vis$QTL_assig))

for (q in 1:length(QTLnames)){
  Effects<-filter(Few,QTL_assig==QTLnames[q] & pval>=THR)
  assign(paste("Effects_",QTLnames[q],sep=""),Effects)
}

# Get largest allele effect in each family 

for (p in 1:length(QTLnames)){
  QTL_eff = get(paste("Effects_",QTLnames[p],sep="")) %>%
    summarise_each(funs(.[which.max(abs(.))]), -c(1:3))
  assign(paste("QTL_",QTLnames[p],sep=""),QTL_eff)
}

# Melt into long format
for (y in 1:length(QTLnames)){
  melted_qtl <- melt(get(paste("QTL_",QTLnames[y],sep="")))
  assign(paste("melted_",QTLnames[y],sep=""),melted_qtl)
}

# Get allele effects for all QTL
QTL_effects<-as.data.frame(names(Few)[-c(1:3)])
for (z in 1:length(QTLnames)){
  data<-get(paste("melted_",QTLnames[z],sep=""))
  QTL_effects <- cbind(QTL_effects, data$value)
  names(QTL_effects)[dim(QTL_effects)[2]]<-QTLnames[z]
}
names(QTL_effects)[1]<-"family"


# Rename families
parent_heading <- read.csv("~/Dropbox/GITHUB/BRIDG6/Datasets/Parents/NAM_parent_heading_nochecks_noduplicateparents.csv", header = T, na.strings = "NA", stringsAsFactors = F)

# Get family number and subpop
Parent_subpop <- select(parent_heading, c(new_fam, Pop_location))

# Rename the columns
colnames(Parent_subpop)[1] <- "family"

QTL_effects$family <- c(1:88)

# Get subpopulations
Effect_subpop <- full_join(Parent_subpop, QTL_effects, by = "family")

# melt
melted_eff <- melt(Effect_subpop, id.vars = c("family", "Pop_location"))

# Set values of zero to NA
melted_eff[melted_eff == 0] <- NA

pdf("~/Desktop/FigureS1.pdf",width=7,height=11)
ggplot(melted_eff, aes(value)) +
  labs(x = "Allele Effect", y = "Frequency") +
  facet_grid(variable ~ ., switch = "y", scales = "free") +
  scale_fill_manual(name = "BRIDG6 Diverse\nParent Subpopulation", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red", "grey55")) +
  theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks = element_blank(), panel.border = element_blank(), text = element_text(size = 8),
        legend.text = element_text(size = 14), strip.placement = "outside") +
  scale_x_continuous(expand = c(0,0)) +
  geom_histogram(binwidth = .2, aes(fill = factor(Pop_location)), na.rm = T)
dev.off()
# Get four top QTL
# Set values of zero to NA
melted_eff_4<-melted_eff[which(melted_eff$variable == "2_1" |melted_eff$variable=="7_2"|melted_eff$variable=="7_1"|melted_eff$variable=="3_1"),]
pdf("~/Desktop/Figure4.pdf",width=7,height=11)
ggplot(melted_eff_4, aes(value)) +
  labs(x = "Allele Effect", y = "Frequency") +
  facet_grid(variable ~ ., switch = "y", scales = "free") +
  scale_fill_manual(name = "BRIDG6 Diverse\nParent Subpopulation", values = c("goldenrod1", "violet", "olivedrab3", "turquoise1", "red", "grey55")) +
  theme(strip.background = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks = element_blank(), panel.border = element_blank(), text = element_text(size = 10),
        legend.text = element_text(size = 14) ) +
  scale_x_continuous(expand = c(0,0)) +
  geom_histogram(binwidth = .2, aes(fill = factor(Pop_location)), na.rm = T)
dev.off()