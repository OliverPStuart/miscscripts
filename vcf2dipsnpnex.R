#!/usr/bin/env R


#This script takes a .vcf file of ONLY SNPs, ideally from FreeBayes, and gives you a nexus file in diploid format:
#Sample1 012221?1010
#This is the format required for SNAPP.
#Reference and alternate homozygotes are coded as 0 and 2, heterozygotes as 1, missing as ?


##################################
######## Initial Checks ##########
##################################


#Read in the file name
input.file <- commandArgs(trailingOnly=T)
output.file <- paste0(input.file,".nex")

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
packages <- c("vcfR")
if(sum(is.installed(packages)) < length(packages)){
  print("Something is missing, check that the required packages are all installed",
        paste0(packages,collapse=" ,"),".")
  stop()
} else {
  print("Loading packages.")
}
#This is ugly as hell, but using lapply() to load a vector of packages makes it impossible to suppress all messages
#Ugly on this end, neater when run through bash
suppressMessages(library("vcfR"))


##################################
######## Data Preparation ########
##################################


print(paste0("Input file = ",input.file))

#Read in .vcf as character lines
vcf <- readLines(input.file)

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

#Pare down genotype fields to only the genotype value
for (i in 1:nrow(geno)){
  
  #The first three characters in the string should be e.g. "0/0"
  geno[i,4:ncol(geno)] <- substr(start=1,stop=3,geno[i,4:ncol(geno)])
  
}

#Replace missing data with NA
geno[geno == ".|." | geno == "./."] <- NA
geno[geno == ".:."] <- NA

#Replace all / with | for consistency with the remaining genotype values
geno <- gsub("/","|",as.matrix(geno))

#One final call to remove weird missing data
geno[geno != "0|0" & geno != "0|1" & geno != "1|0" & geno != "1|1"] <- NA

#Select only genotypes
geno <- geno[,4:ncol(geno)]

#Replace NAs with characters appropriate for NEXUS
geno[is.na(geno)==T|geno=="."] <- "?|?"


##################################
##### Writing the nexus file #####
##################################


#Print the header information to the output file
#Nexus header
cat("#NEXUS",sep="\n",file=output.file)
#New line
cat("\n",file=output.file,append=T)
#Begin data block
cat("BEGIN DATA;",sep="\n",file=output.file,append=T)
#Add dimension information, how many taxa and characters
cat(paste0("\tDIMENSIONS NTAX=",ncol(geno)," NCHAR=",nrow(geno),";"),sep="\n",file=output.file,append=T)
#Add format information, always SNP
cat(paste0("\tFORMAT DATATYPE=SNP MISSING=? GAP=- ;"),sep="\n",file=output.file,append=T)
#Begin matrix block
cat("MATRIX",sep="\n",file=output.file,append=T)

#For each column of the genotype object
for (sample in 1:ncol(geno)){
  
  #Generate a new object of the ith individual's genotypes
  ind <- geno[,sample]

  #Make empty vectors to store genotypes
  hap <- c()
  
  #Then at every site
  for(site in 1:length(ind)){
    
    
    if(ind[site] == "1|1"){hap <- c(hap,"2")}
    else if(ind[site] == "0|0"){hap <- c(hap,"0")}
    else if(ind[site] == "0|1" | ind[site] == "1|0"){hap <- c(hap,"1")}
    else if(ind[site] == "?|?"){hap <- c(hap,"?")}
    
    #Then, once we've done that for all sites
    if(site == length(ind)){
      
      #Collapse the thing into a single character vector
      
      #And append them to the nexus file
      cat(paste0(colnames(as.data.frame(geno))[sample],
                 "\t",paste0(hap,collapse="")),sep="\n",append=T,file=output.file)
      
    }
    
  }
  
}

#Add the final semi-colon to indicate the matrix block is finished
cat(";",sep="\n",append=T,file=output.file)

#Add necessary final line
cat("END;",sep="\n",append=T,file=output.file)

#Print some final messages
print(paste0("File has been saved in current directory as ",output.file))