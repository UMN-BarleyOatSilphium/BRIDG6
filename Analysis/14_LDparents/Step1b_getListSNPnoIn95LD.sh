#!/bin/bash

module load plink

cd /home/smithkp/agonzale/Projects/NAM/Analysis/LD_purge_rm100Close/AfterImp_maf5



WindowSize=10 #in variant count
stepSize=5 #variant count
r_squareThreshold=0.95

for i in {1..7}
do
plink --file CHR_${i} --indep-pairwise $WindowSize $stepSize $r_squareThreshold --out chr_${i}_pruned-snp-list
done

#move output files to another directory

mv *_pruned-snp-list* /home/smithkp/agonzale/Projects/NAM/Analysis/LD_purge_rm100Close/AfterImp_maf5/Output/

cd /home/smithkp/agonzale/Projects/NAM/Analysis/LD_purge_rm100Close/AfterImp_maf5/Output
cat *.prune.in >/home/smithkp/agonzale/Projects/NAM/Analysis/LD_purge_rm100Close/AfterImp_maf5/List_All_SNPs_noLD.txt