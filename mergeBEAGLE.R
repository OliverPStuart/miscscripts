#!/usr/bin/env R


#This script will take two beagle files, unzipped, and merge them
#Data importation is done strictly as characters so there shouldn't be any loss of information by implicit number conversation.
#The script will check that your files have identical alleles at their intersecting markers.
  #If they don't, it will stop and tell you.
#Unfortunately, due to the nature of the beast, the number of loci in the final BEAGLE file == the number in the smaller input file.


##################################
######## Initial Checks ##########
##################################


#Read in the file name
input.files <- commandArgs(trailingOnly=T)[1:2]
outfile <- commandArgs(trailingOnly=T)[3]

#Are there two beagle files being used as input?
if(length(input.files[grep("beagle",input.files)]) != 2){
  print("This script needs two .beagle files in it.")
  stop()
}


##################################
######## Data Preparation ########
##################################

file1 <- read.delim(input.files[1],check.names=F,colClasses="character")
samples1 <- (ncol(file1)-3)/3 ; loci1 <- nrow(file1)
file2 <- read.delim(input.files[2],check.names=F,colClasses="character")
samples2 <- (ncol(file2)-3)/3 ; loci2 <- nrow(file2)
print(paste0("Your files have ",samples1," samples with ",loci1, " loci and ",samples2," samples with ",loci2," loci respectively."))

comb <- merge(file1,file2,by="marker")


##################################
##### Checking Allele Columns ####
##################################


x1 <- comb$allele1.x
x2 <- comb$allele2.x
y1 <- comb$allele1.y
y2 <- comb$allele2.y

if((sum(x1 == y1) / length(x1)) != 1 & (sum(x2 == y2) / length(x2)) != 1){
  
  print("Your input files do not reference the same alleles at their intersecting loci. Script stopping.")
  stop()
  
}

comb <- comb[,-which(names(comb) %in% c("allele1.y","allele2.y"))]
names(comb)[2:3] <- c("allele1","allele2")
names(comb) <- gsub("\\.[0-9]","",names(comb))


##################################
########## File Writing ##########
##################################


write.table(comb,paste0(outfile,".beagle"),col.names=T,row.names=F,quote=F,sep="\t")

