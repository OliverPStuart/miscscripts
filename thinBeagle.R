#!/bin/env R

# Author:  Oliver Stuart
# Date:    15/06/2020
# Project: Locust extinction

#This script subsets a beagle file from a list of sites
#Assumes an ANGSD created beagle with format of "marker" column as CHR_POS 
#Sites is standard CHR\tPOS format, no header

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
              help="outfile name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$beagle)){
  print_help(opt_parser)
  stop("Provide a beagle file.", call.=FALSE)
}


##################################
############ le script ###########
##################################


#Define files
file <- opt$beagle
outfile<- paste0(opt$out,".beagle")

#Get marker columns
loci <- read.delim(file,
                   header=T,
                   check.names=F,
                   stringsAsFactors = F)[,1]

#For later
total <- length(loci)

#Split string and coerce into dataframe
loci <- unlist(strsplit(loci,"_"))
loci <- data.frame(CHR = loci[seq(from=1,to=length(loci)-1,by=2)],
                   POS = loci[seq(from=2,to=length(loci),by=2)])

#Group and sample randomly
loci <- loci %>% group_by(CHR) %>% sample_n(1) %>% as.data.frame

#For later
thin <- nrow(loci)

#Reformat to original marker columns
loci <- paste0(loci$CHR,"_",loci$POS)

#Read in whole beagle file (bottleneck)
beagle <- read.delim(file,
                     header=T,
                     check.names=F,
                     stringsAsFactors = F)

#Thin
beagle <- beagle[beagle$marker %in% loci,]

#Print some messages
print(paste0(total," loci thinned to ",thin," loci, writing beagle now"))

#Write and gzip the beagle
write.table(beagle,outfile,quote=F,row.names=F,col.names=T,sep="\t")
system(paste0("gzip ",outfile))