#!/usr/bin/env python

#by Li Lei, 2017/06/02
#   A script to calculate the alt allele frequency in a VCF file
#   This is useful for calculate the minor allele frequency from vcf file
# adapted from Tomo Kono https://github.com/TomJKono/Misc_Utils/blob/master/VCF_To_Htable.py and https://github.com/TomJKono/Misc_Utils/blob/master/VCF_MAF.py
#usage: python3 VCF_MAF.py yourSNP.vcf >yourSNP.maf
import sys

#   A function to calculate the minor allele frequency
def MAF(x):
    #   get the set of alleles in the list
    genotypes = set(x)
    #   start counting up the frequencies
    freqs = []
    alleles = {}#create a dictionary for the alleles
    for g in genotypes:
        freqs.append(x.count(g)/float(len(x)))
        alleles[g] = x.count(g)
    return max(alleles, key=alleles.get), min(alleles, key=alleles.get), min(freqs)


#   Start iterating through the file
with open(sys.argv[1], 'r') as f:
    for line in f:
        #   ignore header lines
        if line.startswith('##'):
            continue
        #   This defines how many samples in the VCF
        elif line.startswith('#CHROM'):
             print ('Chrom\tPos\tsample_NB\tMajor\tMinor\tMAF')
        else:
            tmp = line.strip().split('\t')
            #   Parse out the relevant information
            chromosome = tmp[0]
            bp_pos = tmp[1]
            ref_allele = tmp[3]
            alt_alleles = tmp[4]
            
            genotypes = tmp[9:] 

            #INFO = tmp[7].split(';')
            #AN = INFO[2].split('=')
            #chr_nb = int(AN[1])/2
            #print (chromosome,bp_pos,ref_allele,alt_alleles,chr_nb)

            format = tmp[8].split(':')
            sample_genotypes = [x.split(':') for x in genotypes]
            #   check if AD is not in the format field
            #   if not, then skip it
            if 'AD' not in format:
                #notes = 'Missing Genotype Call'
                notes = 'Missing Genotype Call'
                maf = "NA"
                print (maf)
            else:
                notes = ''
                #   Which column is the AD?
                #AD_pos = format.index('AD')
                #   For each sample...
                g_column = []
                for g in genotypes:
                	#   In the genotype string, the first element (separated by :) is the actual genotype call
                    call = g.split(':')[0]
                	#   These are diploid calls, and we are assuming they are unphased
                	#   the are listed in the form allele1/allele2
                	#   with 0 = ref, 1 = alt1, 2 = alt2, and so on...
                    alleles = call.split('/')
                    individual_call = '' #define a dictionary for all of the alleles 
                    for x in alleles:
                        if x == '.': #ignor the missing data
                            continue
                    			#g_column.append( 'N')
                        		#individual_call += 'N'
                        else:
                            c = int(x)
                        		#   if it's 0, we just tack on the reference state
                            if c == 0:
                                g_column.append(ref_allele)
                            			#individual_call += ref_allele
                            else:
                            		#   Otherwise, we use it to alternate alleles
                                g_column.append(alt_alleles)
                            			#individual_call += alt_alleles
                	#   Then append the individual call to the column for the genotype matrix
                	#g_column.append(individual_call)
                	#print (len(g_column))
            #   Then, append that column to the genotype matrix
            #   If there is no variation in genotype calls (that is, all lines have the same genotype)
            #   then we don't care about it
            		#print (g_column)
                    unique_calls = set(g_column)
                    if len(unique_calls) <= 1:
                        major = "-"
                        minor = "-"
                        maf = "NA"
                    else:
                        major, minor, maf = MAF(g_column)
 

            #print ('\t'.join([chromosome, bp_pos, ref_allele, alt_alleles, chr_nb, maf]))
            print ('\t'.join([chromosome, bp_pos, str(int(len(g_column)/2)), major, minor, str(maf)]))
            #print (g_column)
            #print (leng_column)
