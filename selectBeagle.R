#!/bin/env R

# Author:  Oliver Stuart
# Date:    15/06/2020
# Project: Locust extinction

#This script subsets a beagle file with a list of individuals
#Beagle file has columns marker\tallele1\tallele2\tIND01\tIND01\tIND01\t... etc.
#List is one name per line


##################################
######## Initial Checks ##########
##################################


#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("optparse","dplyr","magrittr")
if(sum(is.installed(packages)) < length(packages)){
  print(paste0("Something is missing, check that all required packages (",paste(packages,collapse=", "),") are installed."))
  stop()
}
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))

#Bring in the inputs
option_list = list(
  make_option(c("-b", "--beagle"), type="character", default=NULL, 
              help="beagle file to be thinned", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="outfile name", metavar="character"),
  make_option(c("-i", "--inds"), type="character", default=NULL, 
              help="list of individuals", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$beagle) | is.null(opt$inds)){
  print_help(opt_parser)
  stop("Provide a beagle file and a list of individuals.", call.=FALSE)
}


##################################
############ le script ###########
##################################


#Define files
file <- opt$beagle
outfile<- paste0(opt$out,".beagle")
inds <- opt$inds

#Get names
names <- c("marker","allele1","allele2",readLines(inds))

#Read in first line
header <- unlist(strsplit(readLines(file)[1],"\t"))

#For both, add some boundiong characters to avoid partial matches at beginning/end string
#e.g. for if Ind1 also matches Ind10, Ind11, etc.
names <- paste0(":",names,":")
header <- paste0(":",header,":")

#Read in whole beagle file (bottleneck)
beagle <- read.delim(file,
                     header=T,
                     check.names=F,
                     stringsAsFactors = F)
beagle <- beagle[, header %in% names ]
names <- gsub(":","",names)
colnames(beagle) <- c(names[1:3],rep(names[4:length(names)],each=3))

#Write output
write.table(beagle,outfile,quote=F,row.names=F,col.names=T,sep="\t")
system(paste0("gzip ",outfile))