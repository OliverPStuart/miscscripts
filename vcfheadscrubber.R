#!/usr/bin/env R

#This script is a great way to reduce the file-size of your .vcf file after filtering.
#If you're working with a non-model organism, your .vcf might have variants called from reads mapped to a large collection of contigs.
#After filtering, your .vcf file may have many contigs in the header that do not represent any variants.
#This file takes all the variants in your .vcf file, checks which contigs they come from, then removes all others from the header.
#It assumes a contig naming convention of locus.xxxxx, where xxxxx represents a number.


##################################
######## Initial Checks ##########
##################################


#Read in the file name, make an output file name.
input.file <- commandArgs(trailingOnly=T)
output.file <- paste0("SCRUBBED.",input.file)

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
packages <- c("stringr","magrittr")
if(sum(is.installed(packages)) < length(packages)){
  print("Something is missing, check that stringr, and magrittr are all installed.")
  stop()
} else {
  print("Loading packages.")
}
#This is ugly as hell, but using lapply() to load a vector of packages makes it impossible to suppress all messages
#Ugly on this end, neater when run through bash
suppressMessages(library("stringr"))
suppressMessages(library("magrittr"))


##################################
######## Data Preparation ########
##################################


#Read in vcf
vcf <- readLines(input.file)

#Get list of contigs, comma appended to avoid substring matches
  #(Is that even a concern for the %in% operator if it checks matches of entire strings?)
vcf.all.contigs <- vcf[grep("contig",vcf)]
vcf.all.contigs <- str_extract(vcf.all.contigs,"locus.([0-9]+),")

#Get list of the contigs with loci in this vcf, again comma appended
vcf.loc.contigs <- paste0(str_extract(vcf[-grep("#",vcf)],"locus.([0-9]+)"),",")

#Select rows of vcf header where the contig is found in the list of variants
vcf.header.red <- vcf[grep("contig",vcf)][vcf.all.contigs %in% vcf.loc.contigs]

#Make up new header

#Section before contig list
vcf.head.1 <- vcf[1:grep("##contig",vcf)[1] - 1]

#Section after contig list
vcf.head.3 <- vcf[(tail(grep("##contig",vcf),n=1) + 1):tail(grep("#",vcf),n=1)]

#Paste all together
vcf.header.full <- c(vcf.head.1,vcf.header.red,vcf.head.3)

#Paste full header onto the .vcf lines without #, namely the variant representation
vcf.full <- c(vcf.header.full,vcf[-grep("#",vcf)])

#Write the .vcf file
writeLines(vcf.full,output.file)
