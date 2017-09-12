# Author: Ana Poets
# Description: Use a list of all the markers found in exome and gbs data combined. With header. This is the output of ConsensusRas_exome_gbs.R
# GBS SNPs are in the top of the file. Choose the GBS data all times, only  when GBS is NA then look at the call in the exome
# Format of input file: 
# 
# Exome[, which(names(Exome) == "Rusmusson")]
#chr7H_part1_45	GG
#chr7H_part1_63	CC
#chr7H_part1_94	CC
#chr7H_part1_112	GG
#chr7H_part1_131	AA
#chr7H_part1_145	GG
#chr7H_part1_168	CC
#chr7H_part1_185	AA
#chr7H_part1_192	GG


INPUT=open("Rasmusson_gbs_then_exome.txt",'r')
#INPUT.seek(0)

OUTPUT=open("ConsensusRas_gbs_exome.txt","a")



#Read file header
INPUT.readline()

# Create a dictionary for each SNP present in both exome and gbs data
years_dict ={}
for line in INPUT:
    item=line.strip('\n').split('\t')
    if item[0] in years_dict:
            years_dict[item[0]].append(item[1])
    else:
            years_dict[item[0]] = [item[1]]


# Convert dict to a list so we can investigate those SNPs with multiple hits
#lists
temp = []
dictList = []

for key,value in years_dict.iteritems():
	#print (key + value)
    temp= [key,value]
    dictList.append(temp)
    
CountDouble=0
CountTriple=0
SNPconsensus=[]
for SNP in dictList:
    # if there is only one value
    if len(SNP[1]) == 1:
        SNPconsensus.append((SNP[0],SNP[1][0]))
    # if there are two values for the same SNP
    if len(SNP[1]) ==2:
            CountDouble=CountDouble+1
            if SNP[1][0] == 'NA' and SNP[1][1] == 'NA':
                SNPconsensus.append((SNP[0],"NA"))
            else:
                if SNP[1][0] == SNP[1][1]:
                        SNPconsensus.append((SNP[0],SNP[1][0]))
                else:
                        if SNP[1][0] == "NA":
                            SNPconsensus.append((SNP[0],SNP[1][1]))
                        if SNP[1][1] == "NA":
                            SNPconsensus.append((SNP[0],SNP[1][0]))
                        #If the calls are different but are not NA, use the call from GBS
                        #if SNP[1][0] != "NA" and SNP[1][1] != "NA" and SNP[1][0] != SNP[1][1]:
                        if SNP[1][0] != "NA" and SNP[1][1] != "NA" and SNP[1][0] != SNP[1][1]:
                            SNPconsensus.append((SNP[0],SNP[1][0]))
    #If there are thres values for the Same SNP
    if len(SNP[1])==3:
        CountTriple=CountTriple+1
        NApresent=SNP[1].count('NA')
        # evaluate how many sites are NA. If two sites are NA, choose the genotype in the third site as the call. If one site is NA evaluate if the other two sites
        # are the same, if so choose either one, otherwise set to NA.
        if NApresent == 3:
            SNPconsensus.append((SNP[0],SNP[1][0]))
        if NApresent == 2:
            matches=(list(SNP[1]))
            for i in matches:
                if i !='NA':
                    SNPcall=SNP[0],i
                    SNPconsensus.append(SNPcall)
        if NApresent == 1:
            have_value =(list(SNP[1]))
            Genotypes=[]
            for j in have_value:
                if j != 'NA':
                    Genotypes.append(j)
            if Genotypes[0] == Genotypes[1]:
                SNPconsensus.append((SNP[0],Genotypes[1]))
            else:
                SNPconsensus.append((SNP[0],"NA"))
        if NApresent ==0:
            SNPconsensus.append((SNP[0],SNP[1][0]))
            #if SNP[1][0] == SNP[1][1] and SNP[1][1] == SNP[1][2]:
            #   SNPconsensus.append((SNP[0],SNP[1][0]))
            #else:
            #    SNPconsensus.append((SNP[0],"NA"))

for i in SNPconsensus:
    print >> OUTPUT, i[0],i[1]

INPUT.close()
OUTPUT.close()








