#!/bin/env R

# Author:  Oliver Stuart
# Date:    15/06/2020
# Project: Locust extinction

#This script subsets a beaglr file with a list of sites
#Assumes an ANGSD created beagle with format of "marker" column as CHR_POS 

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
  make_option(c("-s", "--sites"), type="character", default="out", 
              help="list of sites", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$beagle) | is.null(opt$sites)){
  print_help(opt_parser)
  stop("Provide a beagle file and a list of sites.", call.=FALSE)
}


##################################
############ le script ###########
##################################


#Define files
file <- opt$beagle
outfile<- paste0(opt$out,".beagle")
sites <- opt$sites

#Get sites
s <- read.table(sites,header=F,
                stringsAsFactors = F)
s <- paste0(s$V1,"_",s$V2)

#Read in whole beagle file (bottleneck)
beagle <- read.delim(file,
                     header=T,
                     check.names=F,
                     stringsAsFactors = F)

#Thin
beagle <- beagle[beagle$marker %in% s,]

#Write and gzip the beagle
write.table(beagle,outfile,quote=F,row.names=F,col.names=T,sep="\t")
system(paste0("gzip ",outfile))