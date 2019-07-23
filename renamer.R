#!/usr/bin/env R

#This file renames your tree sequence names.
#It takes a tab-delimited text file (first column is current names, second column is new names) and a tree file in Newick.
  #IN THAT ORDER.

inputs <- commandArgs(trailingOnly=T)

names <- as.data.frame(read.delim(inputs[1],header=F))

#The following is to control for overlapping names
#e.g. searching for a hypothetical "ind_1" would also return "ind_15", "ind_113", etc.
names[,1] <- paste0(as.character(names[,1]),":")
names[,2] <- paste0(as.character(names[,2]),":")
tree <- readLines(inputs[2])

for(i in 1:nrow(names)){

  tree <- gsub(names[i,1],names[i,2],tree)

}

writeLines(tree,inputs[2])