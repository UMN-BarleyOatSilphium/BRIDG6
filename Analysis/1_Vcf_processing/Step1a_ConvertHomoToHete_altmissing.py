#!/usr/bin/env python

import sys

if len(sys.argv) == 1:
	print "Error: Input file was not determined. Usage: ConvertHomoToHete.py File_baseName"
	sys.exit(0)
	
#Create filename
inFile= str(sys.argv[1]) + '.vcf'
filename_out = str(sys.argv[1]) + '_transformedHETE_toNA.vcf'
filename_count = str(sys.argv[1]) + '_AlleleCount.txt'
#Open handle for writing file
INFILE=open(inFile , "r")
OUTPUT=open(filename_out, "a")
OUTPUT2=open(filename_count,"a")

#INFILE=open("HR202_DP5.vcf","r")
#OUTPUT=open("HR202_DP5_transformedHETE_toNA.vcf","a")
#OUTPUT2=open("AlleleCount_noRasAll_NA.txt","a")

CHANGES=0
AlleleCounts=[]

INFILE.seek(0)
head = [next(INFILE) for x in xrange(11)]
for info in head:
	newLine=info.strip("\n")
	print >>OUTPUT, newLine

for line in INFILE:
	GENOTYPES=[]
	LINE1=line.strip("\n")
	LINE=LINE1.split("\t")
	GENOTYPES.append("\t".join(LINE[0:9]))
	VALUES=LINE[9:]
	for item in VALUES:
		#print item
		if (("0/1" in item) or ("./." in item) or ("1/0" in item)):
			GENOTYPES.append(item)
			continue		
		pieces=item.split(":")
		#print pieces
		alleles=pieces[1].split(",")
		#print alleles
		# For some reason sometimes the number of alternative allele are not present, so assume it is 0
		if not (len(alleles) > 1):
			GENOTYPES.append(":".join(pieces))
		if (len(alleles) >1):
			if not (("0" in alleles[0]) or ("0" in alleles[1])):
				#print alleles
				#Change to:			
				pieces[0]="./."
				GENOTYPES.append(":".join(pieces))
				CHANGES=CHANGES+1
				AlleleCounts.append(alleles)
			else:
				GENOTYPES.append(":".join(pieces))
	print >>OUTPUT, ("\t".join(GENOTYPES))

print >>OUTPUT2, AlleleCounts
OUTPUT.close()
OUTPUT2.close()
INFILE.close()




