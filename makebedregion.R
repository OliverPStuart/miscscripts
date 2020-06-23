#!/usr/bin/env R


#This file takes in a table of variant information and gives you a .bed file of regions around those variants at a given distance.

#Three inputs:
	#Contig information, tabular (name\tlength)
	#SNP positions on contigs, tabular (name\tpos)
	#Value for padding, how many bases on either side of all SNPs?


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
  make_option(c("-c", "--contigs"), type="character", default=NULL, 
              help="contig information in forms of name\tlength", metavar="character"),
  make_option(c("-s", "--sitelist"), type="character", default=NULL,
              help="locus file name, must have contig and position", metavar="character"),
  make_option(c("-p", "--padding"), type="integer", default=NULL,
              help="integer value for the padding around each SNP", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$contigs) | is.null(opt$sitelist) | is.null(opt$padding)){
  print_help(opt_parser)
  stop("Provide files for contig information, site information, and width of region.", call.=FALSE)
}

#Brief warning
print("Check that your inputs are exactly as specified in this script's header.")


##################################
######### Preparing data #########
##################################


#Read in contig and site information
contigs <- read.table(opt$contigs,header=F) ; colnames(contigs) <- c("CHR","LEN")
sites <- read.table(opt$sitelist,header=F) ; colnames(sites) <- c("CHR","POS")
contigs <- merge(contigs,sites,by="CHR")

#Lower bound
contigs$START <- contigs$POS-opt$padding
#Floor (boundary stopd at position 1)
contigs$START[contigs$START < 1] <- 1

#Upper bound (boundary stops at contig's length)
contigs$STOP <- contigs$POS+opt$padding
#Ceiling
contigs$STOP[contigs$STOP > contigs$LEN] <- contigs$LEN[contigs$STOP > contigs$LEN]


##################################
########## Write Ouput ###########
##################################


write.table(contigs[,c(1,4,5)],opt$out,quote=F,sep='\t',
            col.names=F,row.names=F)
