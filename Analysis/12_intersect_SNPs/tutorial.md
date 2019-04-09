# Aim: to find the intersected SNPs among different datasets.
## Datasets:
- 9K iSelect chip SNPs
- SNPs from GBS data of BRIDG6 progines using for GWAS
- SNPs from exome capture data for BRIDG6 parents

## Convert 9k SNPs full to partial;
	#Just use python2:
```
python script/Pseudomolecules_to_Parts.py --vcf data/sorted_all_9k_masked_90idt.vcf >results/sorted_all_9k_masked_90idt_parts.vcf
```

## extract the SNPs_ids:
9k SNPS
```
grep -v "#" sorted_all_9k_masked_90idt_parts.vcf| cut -f 1,2 >9k.SNPid
```
Exome SNPs
```
grep -v "#" Barley_NAM_Regionlast26_All.raw.snps.indels.sorted.vcf| cut -f 1,2 >Barley_NAM_Regionlast26_All.SNPid
grep -v "#" Barley_GATK_Calling_NAM_genomewide.raw.snps.indels.vcf|cut -f 1,2 >Barley_GATK_Calling_NAM_genomewide.SNPid

cat Barley_NAM_Regionlast26_All.SNPid Barley_GATK_Calling_NAM_genomewide.SNPid|sort -k1,1 -k2,2n |unique >unique_NAM_parents.SNPsid
```
GBS

Extract the loci from the vcf files based on the SNP_ID
```
echo "module load vcftools_ML && vcftools --vcf /panfs/roc/scratch/llei/DRUM/BRIDGE6/MNS_July2016_production.vcf --snps /panfs/roc/groups/9/morrellp/llei/Alex_NAM_PROJECT/overlaps_datasets/results/genos_GxE.SNPid --recode --recode-INFO-all --out GWAS_MNS_July2016_production.vcf"|qsub -l mem=22gb,nodes=1:ppn=7,walltime=24:00:00 -m abe -M llei@umn.edu -N extract_GWAS_SNPs -q lab
	788538.nokomis0015.msi.umn.edu
grep -v "#" GWAS_MNS_July2016_production.vcf|cut -f 1,2 >GWAS_MNS_July2016_production.SNPid
```

## Run perl script to extract the intersect SNPs by two datasets:
```
perl .././script/intersect.pl 9k.SNPid unique_NAM_parents.SNPsid|grep "YES" >intersect_9k_NAM_parents

wc -l intersect_9k_NAM_parents 
6647 intersect_9k_NAM_parents


perl .././script/intersect.pl sorted_MNS_July2016_production_GWAS.SNPid unique_NAM_parents.SNPsid|grep "YES" >intersect_GWAS_NAM_parents

wc -l intersect_GWAS_NAM_parents
4068 intersect_GWAS_NAM_parents

perl .././script/intersect.pl 9k.SNPid sorted_MNS_July2016_production_GWAS.SNPid >intersect_9k_GWAS_NAM_parents

wc -l intersect_9k_GWAS_NAM_parents
45
```

## Run perl script tp extract the intersect SNPS by three SNPs

```
perl .././script/three_intersect.pl intersect_9k_GWAS unique_NAM_parents.SNPsid|grep "YES"|wc -l
40
```

## Then do a a little bit math for the numbers and use omnigraffle to draw the Venn diagram for the figure.

