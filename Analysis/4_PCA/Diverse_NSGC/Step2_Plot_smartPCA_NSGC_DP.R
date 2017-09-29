# Author: Ana M. Poets
# Description: Plot smartPCA results for NAM with and without parents
#################################################################################################

rm(list=ls())

#Calculate percentage of variance explained
		EVE<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Analysis/smartPCA/NSGC_diverseParents/output/NSGC_DP.eval")
		pc1<-round((EVE[1,]/sum(EVE)*100),2) #for PC1
		pc2<-round((EVE[2,]/sum(EVE)*100),2)#for PC2
		pc3<-round((EVE[3,]/sum(EVE)*100),2) #for PC3
		pc4<-round((EVE[4,]/sum(EVE)*100),2) #for PC4
		pc5<-round((EVE[5,]/sum(EVE)*100),2) #for PC5
		
		
		TOTAL<-sum(EVE)
		EXPLAIN<-(EVE[,1]/TOTAL)*100
		plot(EXPLAIN)
		DATA<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Analysis/smartPCA/NSGC_diverseParents/output/NSGC_DP.pca.evec",header=F,row.names=1)
		
		dim(DATA)
		
    Donor<-read.table("/Users/agonzale/Documents/SmithLab/NAM/Data/Alex/Data_SFS/NAM_donor_genotype.hmp.txt",header=T,row.names=1)
    
    plot(DATA[,1],DATA[,2],ylab=paste("PC2 (",pc2,"%)",sep=""),xlab=paste("PC1 (",pc1,"%)",sep=""))
    
    #select donor parents
    DP_pca<-DATA[(row.names(DATA) %in% colnames(Donor)),]
    points(DP_pca[,1],DP_pca[,2],col="red")
		