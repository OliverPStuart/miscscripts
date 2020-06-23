#!/usr/bin/env R


#This script takes a .vcf file of ONLY SNPs, ideally from FreeBayes, and gives you a structure file formatted for infocalc.
#The format is something like this:
#loc1 loc2  loc3  loc4
#ind  pop.num pop.name  region  country loc1.geno loc2.geno loc3.geno loc4.geno


##################################
######## Initial Checks ##########
##################################


#Read in the file name
inputs <- commandArgs(trailingOnly=T)

#Is there an input?
if(length(inputs) != 2){
  print("This script needs exactly two inputs; a vcf file and a table indicating the population of each sample. Whichever other file you specify as input will be interpreted as the table, regardless of its extension.")
  stop()
}

#Is there a .vcf file being used as input?
if(identical(grep(".vcf",inputs),integer(0)) == T){
  print("One of the inputs needs to be a vcf file; check the file extension of your vcf.")
  stop()
}

#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("vcfR")
if(sum(is.installed(packages)) < length(packages)){
  print("Something is missing, check that the vcfR is installed.")
  stop()
} else {
  print("Loading packages.")
}

suppressMessages(library("vcfR"))


##################################
######## Data Preparation ########
##################################


vcf.file <- inputs[grep(".vcf",inputs)] ; print(paste0("Input .vcf file = ",vcf.file))
table.file <- inputs[inputs != vcf.file] ; print(paste0("Input table file = ",table.file))

#Read in .vcf as character lines
vcf <- readLines(vcf.file)

#Remove header, selecting only the rows beginning from the header, removing meta-information
vcf <- vcf[grep("CHROM",vcf):length(vcf)]

#Split apart columns of the data table by the tab separator
vcf <- strsplit(vcf,"\t")

#Unlist the split data starting from the second list element, make a dataframe
geno <- data.frame(matrix(unlist(vcf[2:length(vcf)]),nrow=length(vcf)-1, byrow=T),
                   stringsAsFactors=FALSE)

#Make headers for the dataframe
colnames(geno) <- unlist(vcf[1])

#Remove extraneous columns to match what vcfR gives you
geno <- geno[,c(1,2,9:ncol(geno))]

#Find sites with full information, i.e. those that have more than just GT in the FORMAT field
complete <- geno[grep("DP",geno[,3]),]

#Pare down genotype fields to only the genotype value
for (i in 1:nrow(complete)){
  
  #The first three characters in the string should be e.g. "0/0"
  complete[i,4:ncol(complete)] <- substr(start=1,
                                         stop=3,
                                         complete[i,4:ncol(complete)])
  
}

#Replace missing data with NA
complete <- gsub(".:.","-9|-9",as.matrix(complete))

#Replace all / with | for consistency with the remaining genotype values
complete <- gsub("/","|",complete)

#Replace rows with full information in original genotype object with new, pared down rows
geno[grep("DP",geno[,3]),] <- complete

#Select only genotypes
geno <- geno[,4:ncol(geno)]

#Replace NAs with characters appropriate for structure
geno[is.na(geno)==T|geno=="."] <- "-9|-9"
geno[is.na(geno)==T|geno==".|."] <- "-9|-9"
geno[is.na(geno)==T|geno=="./."] <- "-9|-9"


##################################
### Writing the structure file ###
##################################


#Make the output file connection
output.file <- paste0(vcf.file,".struc")

#Get the number of samples
n <- length(vcf[[1]])-9
#The 9 refers to the first nine fields of a standard .vcf file
#CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT

#The first line is a list of all loci
loci <- unlist(vcf)[seq(from=1,
                        to=length(unlist(vcf))-length(vcf[[1]])+1,
                        by=length(unlist(vcf[1])))][-1]
cat(paste(loci,collapse=" "),sep="\n",file=output.file)

#Now, time to make the prepending information
samples <- colnames(geno)
species <- read.table(table.file,header=T)

#Check that sample names match up
if(sum(samples %in% species$id) != length(samples)){
  
  print("There are some samples in the vcf file unaccounted for in your table file; check your inputs and run again. Also check that the column header of the species map table is 'id'.")
  unlink(output.file)
  stop()
  
}

#Add a new index column to the species table
index <- seq(from=1,to=length(unique(species$species)),by=1)
index <- data.frame(index=index,species=unique(species$species))
species <- merge(species,index,by="species")

#Write the prepend
species$prepend <- paste(species$id,species$index,1)
samples <- data.frame(id=samples,order=seq(from=1,to=length(samples)))

#Order the table by the order of the samples in the vcf file
final <- merge(species,samples,by="id")
final <- final[order(final$order),]$prepend

#Writing in the loci, one individual at a time, adding the prepend (which is now correctly ordered) to the beginning
for(i in 1:ncol(geno)){
  
  genos <- unlist(strsplit(geno[,i],"|",fixed=T))
  geno1 <- c(final[i],genos[seq(from=1,to=length(genos)-1,by=2)])
  geno2 <- c(final[i],genos[seq(from=2,to=length(genos),by=2)])
  
  cat(paste(geno1,collapse=" "),sep="\n",file=output.file,append=T)
  cat(paste(geno2,collapse=" "),sep="\n",file=output.file,append=T)
  
}

#And we're done!