# Author: Ana M Poets
# Description: Calculate the frequency of donor parent at each family, and X^2 for segregation distortion. Expectation 1:1 frequency. Df = 2 -1 = 1 (two classes allowed AA and BB)
#  Use SNP windows to estimate frequency of Rasmusson alleles in each chromosome
################################################################################################################################################

rm(list=ls())

library (gtools)


# Get list for all families after QC, then find out which ones didn't have a good genetic map
FAMILY <-read.table("~/Documents/SmithLab/NAM/Analysis/HAPMAP_sharedSNP_nuc/Pri_hmp_Ras_gbsExome/50missing/List_Fam_afterQC.txt")

#Get the list of families that have complete clean genetic maps.
#FAMILY<-read.table("~/Documents/SmithLab/NAM/Analysis/GeneticMap/MSTmap/Plots/Gen_vs_Phy_pos_consecutive/Families_more50SNPs/List_fam_more50SNPs.txt")

# Import list of cross reference names between parents and families (only for those that passed QC ,no HR619)
CrosRef<-read.table("~/Documents/SmithLab/NAM/Analysis/HAPMAP_sharedSNP_nuc/Pri_hmp_Ras_gbsExome/50missing/List_crossRef_par_fam_QC.txt",header=F)


# Function to calculate frequency of "BB" allele
	COUNT_alleleBB<-function(dat){
		BBallele<-length(which(dat == "BB")) 
		if (BBallele  == 0 ){freq <- 0}
		if (BBallele >0){freq <- BBallele/(length(which(!is.na(dat))))}
		return(freq)
	}
	
# Get a the hapmap used to generate the genetic maps.
for (i in 1:(dim(FAMILY)[1])){
	HAPMAP<-read.table(paste("~/Documents/SmithLab/NAM/Analysis/HAPMAP_sharedSNP_nuc/Pri_hmp_Ras_gbsExome/50missing/HAPMAP_forNA_DCO_QC_uniqInd_rmCloseSNP2/", 	FAMILY[i,1], "_Hapmap_NA_DCO_QC_unique_rmCloseSNP2.txt", sep=""), header=T,row.names=1)
	
	# Follow the frequency of allele BB (donor parent). Since this hapmap has been conditionated in having both parents present and segregating, we don't have to follow parent AA
		# Remove parents from the file
		RIL<-HAPMAP[,-c(1,2)]
		
		FREQ<-(as.matrix(apply(RIL,1, COUNT_alleleBB)))
		FREQ<-cbind(row.names(FREQ), FREQ[,1])
		colnames(FREQ)<-c("SNP", as.character(FAMILY[i,1]))
		FREQ<-as.data.frame(FREQ)
		assign(paste("BB_FREQ_",i, sep=""), FREQ)
	
	}
	
# Combine the results for all families
# Make a vector comma separated with all the BB_FREQ_ variables, then to make a table of significant SNPs at each family
LIST_VECTORS<-list(NULL)
Variables_pos<-(grep("BB_FREQ_",ls()))
for (v in 1:length(Variables_pos)){
	POS<-as.numeric(Variables_pos[v])
	LIST_VECTORS[[v]] <-get(ls()[POS])
}

# The example in stackoverflow is: Reduce(function(x, y) merge(x, y, all=TRUE), list(df1, df2, df3))

MY_TABLE<-Reduce(function(x, y) merge(x, y, all=TRUE), LIST_VECTORS)

# Sort markers alphanumerically 

SNP_or_MY_TABLE<-mixedsort(as.character(MY_TABLE[,1]))

MY_TABLE_or<-MY_TABLE[match(SNP_or_MY_TABLE, MY_TABLE[,1]),]

#write.table(MY_TABLE_or, "~/Documents/SmithLab/NAM/Analysis/SNPdensity_&_segDist/88Fam_genmap_SNPdensityandFreq.txt", quote=F,row.names=F,col.names=T,sep="\t")

#===2. Bin chromosomes and estimate mean allele frequency of donor parent ===================================================================
# =========== Calculate the proportion of Rasmusson allele (1- donor allele freq ) in windows of ... SNPs =====================================
# Import table of frequency of donor allele
DATA<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Analysis/SNPdensity_and_segDist/88Fam_genmap_SNPdensityandFreq.txt", header=T, row.names=1 )

# Separate by chromosomes

CHR1<-DATA[grep("1H", row.names(DATA)),]
CHR2<-DATA[grep("2H", row.names(DATA)),]
CHR3<-DATA[grep("3H", row.names(DATA)),]
CHR4<-DATA[grep("4H", row.names(DATA)),]
CHR5<-DATA[grep("5H", row.names(DATA)),]
CHR6<-DATA[grep("6H", row.names(DATA)),]
CHR7<-DATA[grep("7H", row.names(DATA)),]
CHR_UN<-DATA[grep("UN", row.names(DATA)),]


# Find the positions of SNPs at the begining of 100 snp segments
RESULTS_mean<-matrix(NA, ncol=dim(DATA)[2])
colnames(RESULTS_mean)<-colnames(DATA)

CHRnames<-c("CHR1", "CHR2", "CHR3", "CHR4", "CHR5", "CHR6", "CHR7") # ignore the UNknown chrosmosome since we cannot group it by physical proximity


# Function to calculate the mean, setting it to NA if there are not values in that segment
GetMean<-function(dat){
	if (length(which(is.na(dat))) == length(dat)){mean_value<-NA}
	mean_value<-mean(dat[which(!is.na(dat))])
	return(mean_value)
}

for (n in 1:length(CHRnames)){
	print (n)
	# Devide each chromosome in 200 parts and get the average frequency of donor allele for that segment
	CHR_X<-get(CHRnames[n])
	# number of SNPs per segment to get 200 segments
	SNPsPerSegment<-round(dim(CHR_X)[1]/200, digits=0)
	
	# Total posible segments
	dim(CHR_X)[1]/SNPsPerSegment
	TotalSegments<-round(dim(CHR_X)[1]/SNPsPerSegment)
	
	EndPosition<-seq(TotalSegments , dim(CHR_X)[1], by= SNPsPerSegment)
	
	# if the last position is less than the SNPsperSegment from the end of the chromosome, change the position to be the end of the chromosome
	if ((dim(CHR_X)[1]  - EndPosition[length(EndPosition)]) <SNPsPerSegment) {EndPosition[length(EndPosition)]<-(dim(CHR_X)[1])}
	
	
	start<-1
	for (s in 1: length(EndPosition)){
		#print (s)
		start<-start
		end<-EndPosition[s]

			SEGMENT<-CHR_X[start:end,]
			#MEAN<-as.data.frame(t(apply(SEGMENT, 2, summary)[4,]))
			MEAN<-as.data.frame(t(apply(SEGMENT , 2, GetMean)))
			row.names(MEAN)<-paste(CHRnames[n], "_", s, sep="")
			# add values to our Results table
			if (s==1){RESULTS_mean[s,]<-MEAN}else{RESULTS_mean<-rbind(RESULTS_mean, MEAN)}
			start<-(end+1)	
	}
	
}

#
row.names(RESULTS_mean)[1]<-"CHR1_1"

write.table(RESULTS_mean, "~/Documents/SmithLab/NAM/Analysis/SNPdensity_&_segDist/88Fam_genmap_segDist_binned.txt", quote=F,row.names=T,col.names=T,sep="\t")

