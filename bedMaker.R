#!/usr/bin/env R


# Author: Oliver Stuart
# Date: 3/4/2020

#This script takes a list of contigs and their lengths and spits out a bed file targeting reads on those contigs


##################################
######## Initial Checks ##########
##################################


#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("optparse","dplyr")
if(sum(is.installed(packages)) < length(packages)){
  print(paste0("Something is missing, check that all required packages (",paste(packages,collapse=", "),") are installed."))
  stop()
}
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))

#Bring in the inputs
option_list = list(
  make_option(c("-c", "--contig"), type="character", default=NULL, 
              help="tab separated table of contigs\tlengths", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$contig)){
  print_help(opt_parser)
  stop("Provide a contig file.", call.=FALSE)
}


##################################
######### Preparing data #########
##################################


#Read table
d <- read.delim(opt$contig,header=T,
                colClasses="character")
colnames(d) <- c("contig","lengths")
w <- data.frame(contig=d$contig,
                start=rep(0,times=nrow(d)),
                stop=as.numeric(d$length)-1)


##################################
########## Writing data ##########
##################################


#Write the bed file
write.table(w,file=paste0(opt$contig,".bed"),col.names=F,row.names=F,quote=F,sep="\t")
