#!/usr/bin/env R

# Author:  Oliver Stuart
# Date:    19/06/2020
# Project: Locust extinction

#This script renames a .nwk file with a list of names
#The pipeline this fits into makes distance matrices from .beagle files
#And then fastme is used to make a NJ tree from the matrix, and then RAxML is used to paint bootstraps onto the tree
#The list of names needs to be one name per line, same order as in the original .beagle file


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
  make_option(c("-n", "--nwk"), type="character", default=NULL, 
              help="newick tree to be renamed", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", 
              help="outfile name", metavar="character"),
  make_option(c("-i", "--inds"), type="character", default=NULL, 
              help="list of individuals", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$nwk) | is.null(opt$inds)){
  print_help(opt_parser)
  stop("Provide a beagle file and a list of individuals.", call.=FALSE)
}


##################################
############ le script ###########
##################################


#Define files
file <- opt$nwk
outfile<- paste0(opt$out,".nwk")
inds <- opt$inds

#Get names and vector ordered vector for matching ordered names in nwk file
names <- readLines(inds)
names <- data.frame(id=names,
                    ord=0:(length(names)-1))

#Get species to make labels with complete information
species <- read.table("/Volumes/Alter/grasshoppers/lowcov/analysis/variants/speciesmap.txt",
                      header=T)
species <- merge(names,species,by="id")
species$lab <- paste0(species$id,"_",species$species)
species$ind <- paste0("Ind_",species$ord)
species <- species[order(species$ord),c(5,4)]
names <- species

#The following is to control for overlapping names
#e.g. searching for a hypothetical "ind_1" would also return "ind_15", "ind_113", etc.
names[,1] <- paste0(as.character(names[,1]),":")
names[,2] <- paste0(as.character(names[,2]),":")

#Read in tree
tree <- readLines(file)

#For all rows, search the original name (Ind##) and replace with the new label
for(i in 1:nrow(names)){

  tree <- gsub(names[i,1],names[i,2],tree)
  
}

#Write outfile
writeLines(tree,outfile)
