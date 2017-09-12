#!/bin/bash
#PBS -l walltime=10:00:00,mem=62gb,nodes=1:ppn=4
#PBS -N FILTER_allNAM
#PBS -M agonzale@umn.edu
#PBS -m abe
#PBS -r n
#PBS -q small

VCF=/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Original_vcfs/MNS_July2016_production.vcf
OUTPUT=/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Filter

/home/smithkp/agonzale/vcftools-master/bin/vcftools --vcf $VCF --recode-INFO-all  --max-alleles 2 --min-alleles 2 --minDP 5 --minGQ 30 --max-missing 0.5 --maf 0.0003 --recode --remove-filtered-all --remove-indels  --keep /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Original_vcfs/List_all_NAM6_noDub.txt  --out $OUTPUT/MNS_July2016_NAM6ROW_DP5_GQ30_mis50
/home/smithkp/agonzale/vcftools-master/bin/vcftools --vcf $VCF --recode-INFO-all  --max-alleles 2 --min-alleles 2 --minDP 5 --minGQ 30 --max-missing 0.90 --maf 0.0003 --recode --remove-filtered-all --remove-indels  --keep /home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Original_vcfs/List_all_NAM6_noDub.txt  --out $OUTPUT/MNS_July2016_NAM6ROW_DP5_GQ30_mis10
cd $OUTPUT

# for filter with 50% missing data across NAM 6 row.
#/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Script/Step1a_ConvertHomoToHete.py MNS_July2016_NAM6ROW_DP5_GQ30_mis50.recode
/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Script/Step1optional_VCF_processor_nucleotide_AP.py -i MNS_July2016_NAM6ROW_DP5_GQ30_mis50.recode_transformedHETE_toNA.vcf -o MNS_July2016_NAM6ROW_DP5_GQ30_mis50.recode_transformedHETE_toNA_HH -f


#For filter with 10% missing data
/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Script/Step1a_ConvertHomoToHete.py MNS_July2016_NAM6ROW_DP5_GQ30_mis10.recode
/home/smithkp/agonzale/Shared/GBS_Projects/NAM_GBS_2_6row/TASSEL_072016/Script/Step1optional_VCF_processor_nucleotide_AP.py -i MNS_July2016_NAM6ROW_DP5_GQ30_mis10.recode_transformedHETE_toNA.vcf -o MNS_July2016_NAM6ROW_DP5_GQ30_mis10.recode_transformedHETE_toNA_HH -f