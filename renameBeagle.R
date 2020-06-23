#!/bin/env R

# Author:  Oliver Stuart
# Date:    31/03/2020
# Project: Locust extinction

#This script takes a beagle and the list of files used to create it and renames the samples in the beagle file after the input files.
#This is useful since ANGSD output renames all individuals to "Ind1", "Ind2", etc.


##################################
######## Initial Checks ##########
##################################


#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("optparse")
if(sum(is.installed(packages)) < length(packages)){
  print(paste0("Something is missing, check that all required packages (",paste(packages,collapse=", "),") are installed."))
  stop()
}
suppressMessages(library("optparse"))

#Bring in the inputs
option_list = list(
  make_option(c("-b", "--beagle"), type="character", default=NULL, 
              help="beagle file to be modified", metavar="character"),
  make_option(c("-l", "--list"), type="character",  default=NULL,
              help="the list of bam files used to create the beagle file", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$beagle) | is.null(opt$list)){
  print_help(opt_parser)
  stop("Provide a beagle file and its corresponding input list.", call.=FALSE)
}


##################################
######### Preparing data #########
##################################


#Read in beagle file
input.beagle <- read.delim(opt$beagle,
                           check.names=F,
                           header=T,
                           colClasses="character")
#Get associated list
input.list <- opt$list

#Prepare sample name data.frame for selecting spretus data in beagle file
sample.names <- readLines(input.list)

#Below are custom string edits designed specifically for this project
#Since my bam files have a few different suffixes denoting different data densities
#I need several different modifications
#sample.names <- gsub(".o6.bam","",sample.names)
#sample.names <- gsub(".o6.sub.bam","",sample.names)

#Repeat each three times
sample.names <- rep(sample.names,each=3)

#Replace the column names
colnames(input.beagle)[4:ncol(input.beagle)] <- sample.names


##################################
######### Write outputs ##########
##################################


#Write the beagle file
write.table(input.beagle,gsub(".beagle.gz",".renamed.beagle",opt$beagle),col.names=T,row.names=F,quote=F,sep="\t")

#Zip it
system(paste0("gzip ",gsub(".beagle.gz",".renamed.beagle",opt$beagle)))