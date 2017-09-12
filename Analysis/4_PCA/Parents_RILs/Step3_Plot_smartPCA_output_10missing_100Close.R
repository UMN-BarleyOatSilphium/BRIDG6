# Author: Ana M. Poets
# Description: Plot smartPCA results for NAM with and without parents
#################################################################################################

rm(list=ls())


# Import family information, with RILs and parents to know which color to paint them
family_info<-read.csv("~/Dropbox/SmithLab/NAM/Data/Alex/Family_summary.csv",header=T)

# Separate family in populations
Admix<-subset(family_info, family_info$Pop_assignment == "Admixed")
coastal_m<-subset(family_info, family_info$Pop_assignment == "Coastal Mediterranean")
Asian<-subset(family_info, family_info$Pop_assignment == "Asian")
central_e<-subset(family_info, family_info$Pop_assignment == "Central European")
east_af<-subset(family_info, family_info$Pop_assignment == "East African")
unknown<-subset(family_info, family_info$Pop_assignment == "Unknown")

#get identifier for individuals in a family
List_indiv_admix<-c(as.character(Admix$Family), as.character(Admix$NAM_name))
List_indiv_coastal_m <-c(as.character(coastal_m$Family), as.character(coastal_m$NAM_name))
List_indiv_Asian <-c(as.character(Asian$Family), as.character(Asian$NAM_name))
List_indiv_central_e <-c(as.character(central_e$Family), as.character(central_e$NAM_name))
List_indiv_east_af <-c(as.character(east_af$Family), as.character(east_af$NAM_name))
List_indiv_unknown <-c(as.character(unknown$Family), as.character(unknown$NAM_name))

#Calculate percentage of variance explained
ANALYSIS<-c("NAM_all","NAM_RILs","NAM_RILs")
MAIN<-c("NAM RIL and Parents","NAM RIL","NAM RIL")
DIR<-"10missing_100Close"
for (i in 1:length(ANALYSIS)){

		EVE<-read.table(paste("~/Dropbox/SmithLab/NAM/Analysis/smartPCA/fromNAM80_afterQC_RIL_parents/",DIR,"/output/", ANALYSIS[i],".eval",sep=""))
		pc1<-round((EVE[1,]/sum(EVE)*100),2) #for PC1
		pc2<-round((EVE[2,]/sum(EVE)*100),2)#for PC2
		pc3<-round((EVE[3,]/sum(EVE)*100),2) #for PC3
		pc4<-round((EVE[4,]/sum(EVE)*100),2) #for PC4
		pc5<-round((EVE[5,]/sum(EVE)*100),2) #for PC5
		
		
		TOTAL<-sum(EVE)
		EXPLAIN<-(EVE[,1]/TOTAL)*100
		plot(EXPLAIN)
		DATA<-read.table(paste("~/Dropbox/SmithLab/NAM/Analysis/smartPCA/fromNAM80_afterQC_RIL_parents/",DIR,"/output/", ANALYSIS[i],".pca.evec",sep=""),header=F,row.names=1)
		
		dim(DATA)
		#DATA<-DATA[grep("PI|CIho|Ras",row.names(DATA)),]
		admix_indiv<-NULL
		for (p in 1:length(List_indiv_admix)){
			admix_indiv <-c(admix_indiv ,row.names(DATA)[grep(List_indiv_admix[p], row.names(DATA))])
		}
		
		coastalm_indiv<-NULL
		for (p in 1:length(List_indiv_coastal_m)){
			coastalm_indiv <-c(coastalm_indiv ,row.names(DATA)[grep(List_indiv_coastal_m[p], row.names(DATA))])
		}
		
		asian_indiv<-NULL
		for (p in 1:length(List_indiv_Asian)){
			asian_indiv <-c(asian_indiv ,row.names(DATA)[grep(List_indiv_Asian[p], row.names(DATA))])
		}
		
		centrale_indiv<-NULL
		for (p in 1:length(List_indiv_central_e)){
			centrale_indiv <-c(centrale_indiv ,row.names(DATA)[grep(List_indiv_central_e[p], row.names(DATA))])
		}
		
		eastaf_indiv<-NULL
		for (p in 1:length(List_indiv_east_af)){
			eastaf_indiv <-c(eastaf_indiv ,row.names(DATA)[grep(List_indiv_east_af[p], row.names(DATA))])
		}
		unknown_indiv<-NULL
		for (p in 1:length(List_indiv_unknown)){
			unknown_indiv <-c(unknown_indiv ,row.names(DATA)[grep(List_indiv_unknown[p], row.names(DATA))])
		}
		
		#Select the first 5 PC
		DATA<-DATA[,1:5]
		
		pdf(paste("~/Dropbox/SmithLab/NAM/Analysis/smartPCA/fromNAM80_afterQC_RIL_parents/",DIR,"/Plots/",ANALYSIS[i],"PC1vsPC2.pdf",sep=""),width=7,height=5)
		plot(DATA[,1],DATA[,2], ylab=paste("PC2 (",pc2,"%)",sep=""),xlab=paste("PC1 (",pc1,"%)",sep=""),col="white", cex=0.6, main= paste("PC1 vs PC2: ", MAIN[i], sep=""))
		points(DATA[admix_indiv,1],DATA[admix_indiv,2], col="orange", cex=0.6)
		points(DATA[coastalm_indiv,1],DATA[coastalm_indiv,2], col="cyan", cex=0.6)
		points(DATA[asian_indiv,1],DATA[asian_indiv,2], col="magenta", cex=0.6)
		points(DATA[centrale_indiv,1],DATA[centrale_indiv,2], col="green", cex=0.6)
		points(DATA[eastaf_indiv,1],DATA[eastaf_indiv,2], col="pink", cex=0.6)
		points(DATA[unknown_indiv,1],DATA[unknown_indiv,2], col="gray", cex=0.6)
		
		#change point type for parents
		points(DATA[admix_indiv[grep("PI|CIho", admix_indiv)],1],DATA[admix_indiv[grep("PI|CIho", admix_indiv)],2], col="black", cex=0.6, pch=24, bg="orange")
		points(DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],1],DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],2], col="black", cex=0.6, pch=24, bg="cyan")
		points(DATA[asian_indiv[grep("PI|CIho", asian_indiv)],1],DATA[asian_indiv[grep("PI|CIho", asian_indiv)],2], col="black", cex=0.6, pch=24, bg="magenta")
		points(DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],1],DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],2], col="black", cex=0.6, pch=24, bg="green")
		points(DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],1],DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],2], col="black", cex=0.6, pch=24, bg="pink")
		points(DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],1],DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],2], col="black", cex=0.6, pch=24, bg="gray")
		points(DATA[grep("Ras",row.names(DATA)),1],DATA[grep("Ras",row.names(DATA)),2], col="black", cex=0.6,pch=19)
		
		if (i == 1){
		points(DATA[1,1],DATA[1,2], col="black", cex=0.8,pch=19)
		legend("bottomright",pch=c(1,1,1,1,1,1,2,19),col=c("orange","magenta","green","cyan","pink","gray","black","black"), c("Admixed","Asian","Central European","Coastal Mediterranean","East Africa","Unknown","Donor Parent","Rasmusson"), cex=0.6)
		}
		
		if (i ==2){
			legend("topright",pch=c(1,1,1,1,1,1,2,19),col=c("orange","magenta","green","cyan","pink","gray","black","black"), c("Admixed","Asian","Central European","Coastal Mediterranean","East Africa","Unknown","Donor Parent","Rasmusson"), cex=0.6)
		}
		dev.off()
		
		pdf(paste("~/Dropbox/SmithLab/NAM/Analysis/smartPCA/fromNAM80_afterQC_RIL_parents/",DIR,"/Plots/",ANALYSIS[i],"PC1vsPC3.pdf",sep=""),width=7,height=5)
		plot(DATA[,1],DATA[,3], ylab=paste("PC3 (",pc3,"%)",sep=""),xlab=paste("PC1 (",pc1,"%)",sep=""),col="white", cex=0.6, main= paste("PC1 vs PC3: ", MAIN[i], sep=""))
		points(DATA[admix_indiv,1],DATA[admix_indiv,3], col="orange", cex=0.6)
		points(DATA[coastalm_indiv,1],DATA[coastalm_indiv,3], col="cyan", cex=0.6)
		points(DATA[asian_indiv,1],DATA[asian_indiv,3], col="magenta", cex=0.6)
		points(DATA[centrale_indiv,1],DATA[centrale_indiv,3], col="green", cex=0.6)
		points(DATA[eastaf_indiv,1],DATA[eastaf_indiv,3], col="pink", cex=0.6)
		points(DATA[unknown_indiv,1],DATA[unknown_indiv,3], col="gray", cex=0.6)
		points(DATA[grep("Ras",row.names(DATA)),1],DATA[grep("Ras",row.names(DATA)),3], col="black", cex=0.6,pch=19)

		
		#change point type for parents
		points(DATA[admix_indiv[grep("PI|CIho", admix_indiv)],1],DATA[admix_indiv[grep("PI|CIho", admix_indiv)],3], col="black", cex=0.6, pch=24, bg="orange")
		points(DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],1],DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],3], col="black", cex=0.6, pch=24, bg="cyan")
		points(DATA[asian_indiv[grep("PI|CIho", asian_indiv)],1],DATA[asian_indiv[grep("PI|CIho", asian_indiv)],3], col="black", cex=0.6, pch=24, bg="magenta")
		points(DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],1],DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],3], col="black", cex=0.6, pch=24, bg="green")
		points(DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],1],DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],3], col="black", cex=0.6, pch=24, bg="pink")
		points(DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],1],DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],3], col="black", cex=0.6, pch=24, bg="gray")
		if (i == 1){
		points(DATA[1,1],DATA[1,3], col="black", cex=0.8,pch=19)
		legend("bottomright",pch=c(1,1,1,1,1,1,2,19),col=c("orange","magenta","green","cyan","pink","gray","black","black"), c("Admixed","Asian","Central European","Coastal Mediterranean","East Africa","Unknown","Donor Parent","Rasmusson"), cex=0.6)
		}
		
		if (i ==2){
			legend("topleft",pch=c(1,1,1,1,1,1,2,19),col=c("orange","magenta","green","cyan","pink","gray","black","black"), c("Admixed","Asian","Central European","Coastal Mediterranean","East Africa","Unknown","Donor Parent","Rasmusson"), cex=0.6)
		}
		dev.off()
		
		pdf(paste("~/Dropbox/SmithLab/NAM/Analysis/smartPCA/fromNAM80_afterQC_RIL_parents/",DIR,"/Plots/",ANALYSIS[i],"PC1vsPC4.pdf",sep=""),width=7,height=5)
		plot(DATA[,1],DATA[,4], ylab=paste("PC4 (",pc4,"%)",sep=""),xlab=paste("PC1 (",pc1,"%)",sep=""),col="white", cex=0.6, main= paste("PC1 vs PC4: ", MAIN[i], sep=""))
		points(DATA[admix_indiv,1],DATA[admix_indiv,4], col="orange", cex=0.6)
		points(DATA[coastalm_indiv,1],DATA[coastalm_indiv,4], col="cyan", cex=0.6)
		points(DATA[asian_indiv,1],DATA[asian_indiv,4], col="magenta", cex=0.6)
		points(DATA[centrale_indiv,1],DATA[centrale_indiv,4], col="green", cex=0.6)
		points(DATA[eastaf_indiv,1],DATA[eastaf_indiv,4], col="pink", cex=0.6)
		points(DATA[unknown_indiv,1],DATA[unknown_indiv,4], col="gray", cex=0.6)
		points(DATA[grep("Ras",row.names(DATA)),1],DATA[grep("Ras",row.names(DATA)),3], col="black", cex=0.6,pch=19)

		#change point type for parents
		points(DATA[admix_indiv[grep("PI|CIho", admix_indiv)],1],DATA[admix_indiv[grep("PI|CIho", admix_indiv)],4], col="black", cex=0.6, pch=24, bg="orange")
		points(DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],1],DATA[coastalm_indiv[grep("PI|CIho", coastalm_indiv)],4], col="black", cex=0.6, pch=24, bg="cyan")
		points(DATA[asian_indiv[grep("PI|CIho", asian_indiv)],1],DATA[asian_indiv[grep("PI|CIho", asian_indiv)],4], col="black", cex=0.6, pch=24, bg="magenta")
		points(DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],1],DATA[centrale_indiv[grep("PI|CIho", centrale_indiv)],4], col="black", cex=0.6, pch=24, bg="green")
		points(DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],1],DATA[eastaf_indiv[grep("PI|CIho", eastaf_indiv)],4], col="black", cex=0.6, pch=24, bg="pink")
		points(DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],1],DATA[unknown_indiv[grep("PI|CIho", unknown_indiv)],4], col="black", cex=0.6, pch=24, bg="gray")
		
		if (i == 1){
		points(DATA[1,1],DATA[1,4], col="black", cex=0.8,pch=19)
		legend("bottomright",pch=c(1,1,1,1,1,1,2,19),col=c("orange","magenta","green","cyan","pink","gray","black","black"), c("Admixed","Asian","Central European","Coastal Mediterranean","East Africa","Unknown","Donor Parent","Rasmusson"), cex=0.6)
		}
		
		if (i ==2){
			legend("topleft",pch=c(1,1,1,1,1,1,2,19),col=c("orange","magenta","green","cyan","pink","gray","black","black"), c("Admixed","Asian","Central European","Coastal Mediterranean","East Africa","Unknown","Donor Parent","Rasmusson"), cex=0.6)
		}
		
		dev.off()
		
		

}
