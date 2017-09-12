#!/bin/bash

#PBS -l walltime=90:00:00,mem=62gb,nodes=1:ppn=4
#PBS -N ConvertHete
#PBS -M agonzale@umn.edu
#PBS -m abe
#PBS -r n
#PBS -q small

# Convert all Homozygous that were segregating to NA (These are the sites that TASSEL converted to Homo even though there is some evidence of heterozygosity)
cd /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF

for i in `cat List_Families.txt`
do

./ConvertHomoToHete.py ${i}

done

mv *AlleleCount.txt /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/Processed_vcf/HomoHetetoNA_vcf/
mv *transformedHETE_toNA.vcf /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/Processed_vcf/HomoHetetoNA_vcf/


# Convert vcf files to HAPMAP
cd /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/Processed_vcf/HomoHetetoNA_vcf/

for j in `cat /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/List_Families.txt`
do
/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/Pipeline/GBarleyS/Pipeline_Scripts/VCF_processor_AP.py -i ${j}_transformedHETE_toNA.vcf -o ${j} -f

done

mv *_hmp.txt /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/Processed_vcf/HapMap_files/


# Convert vcf files to HAPMAP using genotype calls based on nucleotydes instead of A, B, AB. Use "HH" for heterozygote sites 
cd /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/Processed_vcf/HomoHetetoNA_vcf/

for j in `cat /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/List_Families.txt`
do
/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/Pipeline/GBarleyS/Pipeline_Scripts/Step1optional_VCF_processor_nucleotide_AP.py -i ${j}_transformedHETE_toNA.vcf -o ${j}_nucl -f

done

mv *_hmp.txt /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_filteredVCF/Processed_vcf/HapMap_files/HAPMAP_nucleotides/