#!/usr/bin/env R

#This file takes an input .vcf file and produces two output files. Both are tab-delimited tables. The first is a table of locus positions that are entirely invariant across all samples. The second is a table of locus positions that are only variant in one sample. Files are saved to the directory of the input .vcf file, and can then be used with vcftools to remove the unwanted loci.
#A .vcf file produced by whichever variant-calling software may contain sites where all individuals have the same genotype. They will still be considered variant, however, if any allele in that genotype is non-reference. SNP-based Bayesian softwares also use only variant sites, so you may want to remove completely homogeneous ones ahead of time for that purposes as well. Singleton loci can also be good to ID ahead of time.


##################################
######## Initial Checks ##########
##################################


#Read in the file name
input.file <- commandArgs(trailingOnly=T)

#Is only one .vcf file being used as input?
if(length(input.file)>1){
  print("Try only running one .vcf file at a time. Have some mercy on this puny script.")
  stop()
}

#Is the file a .vcf file?
if(substr(input.file,
          start=nchar(input.file)-2,
          stop=nchar(input.file)) != "vcf"){
  print("This script requires a .vcf file with the correct file extension (.vcf) as input.")
  stop()
}

#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("vcfR","dplyr","magrittr")
if(sum(is.installed(packages)) < length(packages)){
  print("Something is missing, check that dplyr, vcfR, and magrittr are all installed.")
  stop()
} else {
  print("Loading packages.")
}
#This is ugly as hell, but using lapply() to load a vector of packages makes it impossible to suppress all messages
#Ugly on this end, neater when run through bash
suppressMessages(library("vcfR"))
suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))

##################################
######## Data Preparation ########
##################################


print(paste0("Input file = ",input.file))

#Read in the .vcf file
input <- vcfR::read.vcfR(input.file,verbose=F)

#Create an object with the locus positions and genotypes of the loci, also keeping the FORMAT field
geno <- cbind(input@fix[,1:2],input@gt)

#Pare down genotype fields to only the genotype value
for (i in 1:nrow(geno)){
  
  #the first three characters in the string should be e.g. "0/0"
  
  geno[i,4:ncol(geno)] <- substr(start=1,stop=3,geno[i,4:ncol(geno)])
  
}

#Replace all / with | for consistency with the remaining genotype values
geno <- gsub("/","|",geno)

#Replace missing data with NA
geno[geno == ".|."] <- NA
geno[geno == ".:."] <- NA

#Change all "1|0" values to "0|1" so that a site with all hets is considered invariant
geno[geno == "1|0"] <- "0|1"

##################################
## Identifying Invariant Sites ###
##################################


#Produce indices of the sites with the number of unique genotypes and number of missing individuals
for (i in 1:nrow(geno)){
  
  if(i == 1){
    ident<-c()
    empty<-c()
  }
  
  ident <- c(ident,geno[i,4:ncol(geno)] %>% as.vector %>% unique %>% length)
  empty <- c(empty,geno[i,4:ncol(geno)] %>% is.na %>% sum)
  
}

#Select idents without missing data
badloci.full <- geno[ident==1 & empty == 0,]

#Select rows with missing data
missing <- geno[empty>0,]

#Loop over rows to remove empties, then flag by whether length(unique()) == 1
for (i in 1:nrow(missing)){
  
  if(i==1){
    missing.ident<-c()
  }
  
  genos <- missing[i,4:ncol(missing)] %>% as.vector
  genos <- genos[!is.na(genos)]
  missing.ident <- c(missing.ident,length(unique(genos))==1)
  
}

#Use this index to isolate rows that have missing data AND homogenous genotypes
badloci.missing <- missing[missing.ident,]

#Concatenate full list of invariant loci, retaining only locus ID information
badloci.ident <- rbind(badloci.full,badloci.missing)
badloci.ident <- badloci.ident[,1:2]

#Save this as a table to the directory
out.file <- paste0(input.file,".badloci.ident")
write.table(badloci.ident,file=out.file,sep="\t",row.names=F,quote=F)


##################################
## Identifying Singleton Sites ###
##################################


#Remove the badloci.ident sites from the geno object
#To use the subset() function you need to have a single vector of values to check
#Combine the CHROM and POS columns into a unique ID column for the geno and badloci.ident objects
geno <- cbind(geno,paste0(geno[,1],".",geno[,2]))
badloci.ident <- cbind(badloci.ident,paste0(badloci.ident[,1],".",badloci.ident[,2]))

#Remove the values in the geno ID column that also appear in the badloci.ident ID columm
#Stored as new object; the varian sites fromthe geno object
geno.var <- subset(geno,!(geno[,ncol(geno)] %in% badloci.ident[,ncol(badloci.ident)]))

#And remove the temporary ID column
geno <- geno[,1:(ncol(geno)-1)]
geno.var <- geno.var[,1:(ncol(geno.var)-1)]
badloci.ident <- badloci.ident[,1:(ncol(badloci.ident)-1)]

#Loop over the rows of geno.var and flag rows where sum(xtabs(~genos) > 1) == 1
#Condition is that the the number of genotypes with >1 representative is only 1

for (i in 1:nrow(geno.var)){
  
  if(i==1){
    singletons<-c()
  }
  
  genos <- geno.var[i,]
  genos <- genos[!is.na(genos)]
  
  check <- sum(xtabs(~genos)>1)==1
  singletons <- c(singletons,check)
  
}

#Pare down geno.var dataset to only parsimony uninformative sites
#Retaining only tbhe CHROM and POS values
badloci.singleton <- geno.var[singletons,][,1:2]

#Save the results
out.file <- paste0(input.file,".badloci.singleton")
write.table(badloci.singleton,file=out.file,sep="\t",row.names=F,quote=F)

#Print some final messages
print(paste0("The total number of invariant sites is ",nrow(badloci.ident),". Locus positions saved to ",paste0(input.file,".badloci.ident")))
print(paste0("The total number of singleton sites is ",nrow(badloci.singleton),". Locus positions saved to ",paste0(input.file,".badloci.singleton")))
