#!/usr/bin/env python

# This code was modified from Jeff Neyhart's VCF processor script.
# The main differences are that this script only converts from VCF file to a hapmap and
# it maintains the chromosome name as part of the SNP name.
# NOTE: this version convert genotypes using nucleotide states instead of A,B,AB. Heterozyote sites are marked as 'HH'

# Import modules
import subprocess # To spawn subprocesses
import sys # To take arguments
import re # To use regular expressions
import argparse # To get the arguments

#####
# Define functions
#####

def basic_hapmap(VCF):
	print "Writing the non-reformatted hapmap file using rrBLUP encoding."
	# Create filename
        filename = str(args.outfile) + '_hmp.txt'
        # Open handle for writing file
        handle = open(filename, 'w')
	
	# Lists for handling chromosome name
	chrom_l = [] # Empty list of chromosomes
	chrom_index = [] # Index of the chromosome positions

	# Start reading through the vcf file
       	for line in VCF:

               	if line.startswith('##'):
                       	continue
                # Look for the #CHROM line as this is the line that contains sample information
       	        elif line.startswith('#CHROM'):
               	        # Split the line on tabs
                       	tmp = line.strip().split('\t')
                        # Fine the format field
       	                format_field = tmp.index('FORMAT')
               	        # Get the samples out of the list
                        # Add 1 to the index, so "FORMAT" isn't included
                        samples = tmp[format_field + 1:]
			# First the header line with the information
               	        handle.write('rs\tallele\tchrom\tpos\t' + '\t'.join(samples) + '\n')

                # Now we have the sample names; let's move on to the genotype data
       	        else:
               	        # Create a new list to store the data
                       	toprint = []

                        tmp = line.strip().split('\t')

       	                # Assigning variable
               	        chrom = tmp[0]
                       	position = tmp[1]
                        ref_allele = tmp[3]
       	                alt_allele = tmp[4]

	#		# Handling chromosome names
	#		if chrom in chrom_l:
	#			pass
	#		else:
	#			chrom_l.append(chrom)
	#		# Assign the index of the chromosome within the unique list of chromosomes
	#		## as the name of that chromosome
	#		chrom_name = str(chrom_l.index(chrom) + 1)
			chrom_name = str(chrom)
			
               	        # The genotype data
                       	genotypes = tmp[9:]

                        # Create variable for the output file
       	                # Create the alleles variable
               	        alleles = ref_allele + '/' + alt_allele
       	                # The position variable was already created
               	        # Create a name for the SNP
                       	snp_id = chrom_name + '_' + position

                        # Append to the list each snp_id, alleles, etc
       	                toprint.append(snp_id)
               	        toprint.append(alleles)
                       	toprint.append(chrom_name)
                        toprint.append(position)

       	                for g in genotypes:
               	                # The genotype string is separated by :
                       	        # The first element of the genotype string is the genotype call
                               	call = g.split(':')[0]
                                # Genotypes are listed as allele1/allele2
       	                        # Assume the genotypes are unphased
               	                # 0 = ref, 1 = alt 1
                       	        # 0/0 = homo ref, 0/1 = het, 1/1 = homo alt
                                
				# Encode genotypes in rrBLUP format
				individual_call =''

	                        if call == '0/0': # If the call is 0/0, declare as 1
  	                        	individual_call += ref_allele + ref_allele
                	       	elif call == '0/1': # If the call is 0/1, declare as 0
                        	        individual_call += 'HH'
                                elif call == '1/0': # If the call is 1/0, declare as 0
                                        individual_call += 'HH'
                                elif call == '1/1': # If the call is 1/1, declare as -1
                                      	individual_call += alt_allele + alt_allele
	                        else:
        	                        individual_call += 'NA' # If it isn't any of the above, it its missing
               	                # Append the individual calls to the genotype matrix row
                       	        toprint.append(individual_call)

                        # Print the organized list
       	                handle.write('\t'.join(toprint) + '\n')
	
	print "File was written as " + filename
	# Close the handle
	handle.close()
##### End of function #####




#####
# Define the arguments
#####

# Description
DESC = """A Python program to convert a VCF file to a hapmap file in rrBLUP format {-1, 0, 1}. 
This tool is part of the GBarleyS pipeline. the -r flag. In this case, a reformatted hapmap file will be provided
in which the markers have been renamed to reflect their chromosome position. Note
that nothing is written to stdout; instead, files with the outfile name are written."""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input VCF file
parser.add_argument('-i',
                '--vcf_in',
                metavar = 'VCFIN',
                help = 'Input VCF file',
                required = True)
# Output file name
parser.add_argument('-o',
                '--outfile',
                metavar = 'OUTFILE',
                help = 'Output file basename (i.e. no extension)',
                required = True)
parser.add_argument('-f',
		'--hapmap',
		action = 'store_true',
		help = 'Boolean flag for whether a hapmap file should be exported',
		required = False)

# Parse the arguments
args = parser.parse_args()


#####
# Execute program functions
#####

# Print a statement if no flags are thrown
if args.hapmap:
	with open(args.vcf_in, 'r') as VCF:
		basic_hapmap(VCF)