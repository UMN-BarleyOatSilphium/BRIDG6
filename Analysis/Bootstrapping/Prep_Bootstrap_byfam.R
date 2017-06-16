# Bootstrap sampling mapping panels for the NAM

# load package
require(NAM)

# import files for processing
y <- read.csv("~/Documents/PhD/NAM/BLUPs/fam_fixed/fam_BLUEs_nob_no19_famcorrect.csv", header=T)
gen <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Genotypes/NAM_MNS_July2016_DP5_GQ30_mis80_6060ind.recode_transformedHETE_toNA_HH_hmp_withRAShomo_rasBased_rm100Close_naExcessHH_NAonlyHeteHomo_naMAF_imputed_naMAFagain.csv", header = T)
LD <- read.table("~/Documents/PhD/NAM/NAM_mapping/Genotypes/List_SNPtoKeep_noLD80.txt")

# Add X so list matches column names in gen
LD_X <- sub("^", "X", as.character(LD$V1))

# See what I'm working with
class(y)
class(gen)
dim(y)
dim(gen)

# Filter for markers not in LD, this is going to give an error message but it's okay, it's just the markers with unknown location which are added in the next step
gen_noLD <- select(gen, one_of(LD_X))

# Add Unk markers 
#gen_Unk <- select(gen, contains("UN"))

# Combine markers not in LD and Unk
#gen_noLD_Unk <- cbind(gen_noLD, gen_Unk)

# Add line names back
gen_noLD$X <- gen$X

# Remove Ras from the first row
gen_noRas <- gen_noLD[-1,]

# Filter for rows present in both genotypes and phenotypes
gen_1 <- gen_noRas[match(y$line_name, gen_noRas$X, nomatch = NA, incomparables = F),]

# Remove rows where all genos are NA in phenos
gen_2 <- gen_1[rowSums(is.na(gen_1)) != ncol(gen_1),]

# filter for rows present in both in phenos
y_1 <- y[match(gen_2$X, y$line_name, nomatch = NA, incomparables = F),]
y_2 <- y_1[rowSums(is.na(y_1)) != ncol(y_1),]

# View(gen_2)
# View(y_2)

# Remove line names
gen_naked <- select(gen_2, -(X))

# Marker QC + imputation
adjusted_genotypes = snpQC( gen=gen_naked, MAF=0.05, impute=FALSE)
rownames(adjusted_genotypes) = y_2$line_name
dim(adjusted_genotypes)

# 36,212 markers

# getting chr
chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
class(chr)

# Make file for MSI
write.csv(adjusted_genotypes, "~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/genos_MAF05.csv")
write.csv(y_2, "~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/phenos_MAF05.csv")
### Prep loops for mapping below, then copy and edit in map_bootstrap script ###

# write.csv(adjusted_genotypes, "~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/genos_80miss_byfam.csv")
# write.csv(y_2, "~/Documents/PhD/NAM/NAM_mapping/80miss_byfam/phenos_80miss_byfam.csv")


# Import files for testing loop
y_1 <- read.csv("~/Documents/PhD/NAM/NAM_mapping/Bootstrapping/phenos_80miss_forbootstrap.csv", header=T)

y_2 <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/phenos_MAF05.csv", header=T)[,-1]
gen_2 <- read.csv("~/Documents/PhD/NAM/NAM_mapping/80miss_MAF05/genos_MAF05.csv", header = T) 
#gen_test <- as.data.frame(adjusted_genotypes)
#y_test <- as.data.frame(y_2)

# Store a list of families
SAMPLE_NAMES<-c(1:88)

# Empty matrix to store results
RESULTS_SIGNIF<-NULL
SNPs_significant<-NULL
# for b bootstrappings do:
for (i in 1:10){
  # Select 20 individuals randomly from the population
  Sample200 <- SAMPLE_NAMES[ sample(1:length(SAMPLE_NAMES), 2, replace = T)]
  
  # Filter for phenotypes present in the sample
  Sampled_phen<-y_2[(y_2[,3] %in% Sample200),]
  
  # Filter for genotypes present in the sample
  Sampled_gen<-gen_2[(gen_2[,1] %in% Sampled_phen$line_name),]
 
  # Strip off first column 
  gen_naked <- Sampled_gen[,-1]
  
  # QC for MAF of 1/250 (0.00047)
  adjusted_genotypes <- snpQC(gen = gen_naked, MAF = 0.003, impute = FALSE)
  rownames(adjusted_genotypes) <- Sampled_phen$line_name
  head(Sampled_phen)
  dim(adjusted_genotypes)
  
  # Get chrom 
  chr = data.frame(table(gsub("._.+$", "",colnames(adjusted_genotypes))))[,2]
  class(chr) # should be integer
  
  # Need to write to 10 different csv... 
  my_gwas = gwas2(Sampled_phen$BLUE,adjusted_genotypes,Sampled_phen$family,chr)
  number_of_markers = nrow(adjusted_genotypes)

    #= Saving plots ====
  pdf(paste("~/Desktop/GWAS_",i,".pdf",sep=""),width=7,height=5)
  plot( my_gwas, FDR = 0.05)
  dev.off()
  
  #= Assign output to a variable and safe also to a file
    RESULTS<-my_gwas$PolyTest
    if (length(which(RESULTS$pval >(log(number_of_markers*0.05)))) >0){
      Significant_Results<-RESULTS[which(RESULTS$pval >(log(number_of_markers*0.05))),]
      SNPs_significant<-c(SNPs_significant,row.names(Significant_Results))
      RESULTS_SIGNIF<-c(RESULTS_SIGNIF,length(which(RESULTS$pval > (log(number_of_markers*0.05))))) # change significance threshold
      
    
    }
    assign(paste("NAM_",i,sep=""),my_gwas$PolyTest)
  write.table(my_gwas$PolyTest, paste("~/Desktop/NAM_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  write.table(my_gwas$SNPs, paste("~/Desktop/NAM_",i,".txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  
   }
write.table(as.data.frame(RESULTS_SIGNIF), "~/Desktop/Sig_count_round.txt",quote=F,row.names=F,col.names=F,sep="\t")

# The number of times that each significant SNP was found significant 
write.table(as.data.frame(table(SNPs_significant)), "~/Desktop/Count_timesSigSNP.txt",quote=F,row.names=F,col.names=F,sep="\t")

# 1888 markers


log(6000*0.05)


### Concatonate results of single family analyses ###

for ( i in 1:(dim(Family_list)[1])){
  if(i == 19)next
  Fam_QTL <- read.table(paste("~/Documents/PhD/NAM/NAM_mapping/By_fam/QTL_", Family_list[i,1], ".txt", sep = ""), header = T)
  Fam_QTL_lod<-(Fam_QTL[,c(2,7)])
  names(Fam_QTL_lod)<-c("SNP",paste("LOD_HR",Family_list[i,1],sep=""))
  assign(paste("NAM_", Family_list[i,1], sep = ""), Fam_QTL_lod)
}

LIST_VECTORS<-list(NULL)
Variables_pos<-(grep("NAM_",ls()))
for (v in 1:length(Variables_pos)){
  POS<-as.numeric(Variables_pos[v])
  LIST_VECTORS[[v]] <-get(ls()[POS])
}
MY_TABLE<-Reduce(function(...) merge(..., all=TRUE), LIST_VECTORS)
write_csv(MY_TABLE, "~/Documents/PhD/NAM/NAM_mapping/By_fam/Table_1_try.csv")
