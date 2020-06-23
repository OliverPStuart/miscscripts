#!/usr/bin/env R

#This script takes the fasta file output by vsearch and gives you a subset fasta of clusters with a max specified depth.
#Inputs are the file, the depth required, and the length of contig required. This assumes that sequences are on a single line.


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
  make_option(c("-f", "--fasta"), type="character", default=NULL, help="fasta of vsearch cluster contigs, one line per sequence", metavar="character"),
  make_option(c("-d", "--depth"), type="integer",  default=NULL, help="minimum depth", metavar="integer"),
  make_option(c("-l", "--length"), type="integer", default=NULL, help="mininum length", metavar="integer"),
  make_option(c("-m", "--max"), type="integer", default=NULL, help="maximum depth", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="out.fa", 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$fasta) | is.null(opt$depth) | is.null(length)){
  print_help(opt_parser)
  stop("Provide a single-line fasta file, a depth, and a length.", call.=FALSE)
}


##################################
######### Preparing data #########
##################################


#Read in file
fasta <- readLines(opt$fasta)

#Length of object
n <- length(fasta)

#Make a manipulable data.frame object
d <- data.frame(name=fasta[seq(from=1,to=n-1,by=2)],
                seq=fasta[seq(from=2,to=n,by=2)],stringsAsFactors=F)

#Make variables for depth and length of contigs
d$depth <- as.numeric(gsub(".*=","",d$name)) ; d$len <- nchar(d$seq)

#Pare down data.frame object
if(is.null(opt$max)){
  d <- d[d$depth > opt$depth-1 & d$len > opt$length-1,]
} else {
  d <- d[d$depth > opt$depth-1 & d$len > opt$length-1 & d$depth < opt$max+1,]
}


##################################
########## Write Ouput ###########
##################################


#Write file
writeLines(paste(rbind(d$name,d$seq)),opt$out)
