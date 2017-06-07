# Author: Alex Ollhoff
# Description: Map across all families using environment BLUPs and genotypes filtered by family with up to 80% missingness
# Date: March 16, 2017
#################################################################################################################################
# pbs file 
# vim Run_Step2_Map_80miss_fillin_GxE.sh

#!/bin/bash
#PBS    -l      walltime=20:00:00,mem=62gb,nodes=1:ppn=4     
#PBS    -o      /home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Output/mapping_output
#PBS    -e      /home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Output/mapping_error_output
#PBS    -N      NAM_mapping_80miss_MAF05
#PBS    -M      agonzale@umn.edu
#PBS    -m      abe 
#PBS    -r      n
#PBS    -q mesabi

cd /home/smithkp/agonzale/Projects/NAM/Scripts
module load R/3.3.1 
R CMD BATCH --slave Step2_Map_80miss_fillin_GxE.R

##########################################################################

# vim NAM_mapping_GxE.R

# Nested Association Mapping in the 6 row barley NSGC NAM Population using BLUEs and genotypes filtered by family with up to 80 missingness
# Alex Ollhoff
# April 24, 2017

library(NAM)

# Load genotypic data
# ref allele = 2 = major allele
# 0 = minor allele, this allows minor alleles to have different effects if stratification is provided
gen = read.csv("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Input/genos_GxE.csv")[,-1]
gen = data.frame(gen)

# Load phenotypic and family data
# phe1 = pheno in col 3, phe2 = fam in col 2
y = read.csv("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Input/phenos_GxE.csv", header = TRUE)
y = data.frame(y)
# get chrom
chr = data.frame(table(gsub("._.+$", "",colnames(gen))))[,2]

# GWAS on BLUEs
my_gwas = gwas2(y$DAP_BLUPs,gen,y$family,chr)

# export GWAS results
write.csv(my_gwas$MAP,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Output/my_gwas_MAP_out_80miss_BLUE.csv")

write.csv(my_gwas$Method,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Output/my_gwas_Method_out_80miss_BLUE.csv")

write.csv(my_gwas$PolyTest,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Output/my_gwas_PolyTest_out_80miss_BLUE.csv")

write.csv(my_gwas$SNPs,"/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Output/my_gwas_SNPs_out_80miss_BLUE.csv")


# Plot using FDR correction
pdf("/home/smithkp/agonzale/Projects/NAM/Analysis/Imputation/Imputed_May/GWAS_imp_pri/Output/PLOT_my_gwas_FDR_80miss_BLUE.pdf",width=7,height=5)
plot(my_gwas, FDR = 0.05, gtz = TRUE)
dev.off()


