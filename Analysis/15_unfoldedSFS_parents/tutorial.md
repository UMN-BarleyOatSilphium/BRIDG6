# Objective: to make the folded SFS for 89 parents with gBS data

## Steps:

### 1 Extract the 89 parents SNPs from the raw vcf file:

 Ana sent me the [88 parents information](https://github.com/UMN-BarleyOatSilphium/BRIDG6/blob/master/Datasets/Parents/Parent_families_88.txt). I added the Russmusson to the list and use it to extract the parents SNPs from the raw data.

```
echo "module load vcftools_ML && vcftools --vcf /panfs/roc/scratch/llei/DRUM/BRIDGE6/MNS_July2016_production.vcf --keep /panfs/roc/scratch/llei/DRUM/BRIDGE6/89_parents.list --recode --recode-INFO-all --out /panfs/roc/scratch/llei/DRUM/BRIDGE6/parents_MNS_July2016_production.vcf"|qsub -l mem=22gb,nodes=1:ppn=7,walltime=24:00:00 -m abe -M llei@umn.edu -N gbs-filtering -q lab
```
### 2 reformat the vcf file:

```
sed -i 's/1H1/chr1H_part1/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/1H2/chr1H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/2H1/chr2H_part1/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/2H2/chr2H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/3H1/chr3H_part1/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/3H2/chr3H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/4H1/chr4H_part1/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/4H2/chr4H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/5H1/chr5H_part1/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/5H2/chr5H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/5H2/chr5H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/6H1/chr6H_part1/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/6H2/chr6H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/7H1/chr7H_part1/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/7H2/chr7H_part2/g' parents_MNS_July2016_production.vcf.recode.vcf
sed -i 's/UN/chrUn/g' parents_MNS_July2016_production.vcf.recode.vcf
```

### 3 Convert the parts to spedomolecular format:

```
python /panfs/roc/groups/9/morrellp/llei/Inversion/GATK_VariantRecalibrator/Convert_Parts_To_Pseudomolecules.py parents_MNS_July2016_production.vcf.recode.vcf >psedumolecular_parents_MNS_July2016_production.vcf.recode.vcf

sed -i 's/Schr/chr/g' psedumolecular_parents_MNS_July2016_production.vcf.recode.vcf

```

### 4 Keep SNPs Alex used for GWAS

```
vcftools --vcf psedumolecular_parents_MNS_July2016_production.vcf.recode.vcf --snps geno_GXE.sample --recode --recode-INFO-all --out filtered_psedumolecular_parents_MNS_July2016_production

```

### 5 Filtering based on the Alex's criteria for the progeny data:

```
vcftools --vcf filtered_psedumolecular_parents_MNS_July2016_production.recode.vcf --maf 0.003 --min-meanDP 5 --max-missing 0.1  --recode --recode-INFO-all --out maf0.003_meanDP5_psedumolecular_parents_MNS_July2016_production

```

### 6 Calculate the MAF and the count:

```
module load python3
python ~/Inversion/species_range_limits/script/VCF_MAF.py maf0.003_meanDP5_psedumolecular_parents_MNS_July2016_production.recode.vcf >maf0.003_meanDP5_psedumolecular_parents_MNS_July2016_production.recode.maf

```
### 7 Plot the folded SFS with `folded_SFS.R`
