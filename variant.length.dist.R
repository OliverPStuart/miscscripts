#!/usr/bin/env R

#This script takes a list of called variants and spits out a plot of the distribution of their lengths

#The input needs to be taken directly from a .vcf file, first using awk to remove headers ('!/#/')
#Then using cut to get the alternate allele genotype(s) field

#Error to be resolved: account for the fact that indels may also be wrt the reference
#Probably unimportant for now, gross information content is roughly the same
#These variants will be undetected provided that they're length invariant across samples

input <- commandArgs(trailingOnly=T)

#Check the number of inputs, if > 1, stop
if(length(input)>1){
  print("Make sure you only have one input file.")
  stop()
}

file <- input
library(magrittr)

#Read in the file
dist <- as.character(read.delim(file,header=F)[,1])

#Split apart the genotypes
dist <- strsplit(dist,",")

#Remove indels
#Identify list elements where there are genotypes with different lengths
dist.indels <- lapply(FUN=nchar,dist) %>% lapply(FUN=unique) %>% lapply(FUN=length)
#Select variants that are not these
dist <- dist[unlist(dist.indels)==1]

#Loop over the list and count the number of characters in each genotype
for(i in 1:length(dist)){
  
  if(i==1){lengths<-c()}
  
  lengths <- c(lengths,nchar(dist[[i]][1]))
  
}

#Open plotting device, a .tiff
tiff(file=paste0(file,"poly.length.dist.tiff"),width=1000,height=1000,units="px")

#Make a histogram
hist(lengths,breaks=max(lengths),
     main="frequency of polymorphism lengths\nweak filter (bin=1 extends up to ~1000)")

#Close the plotting device
dev.off()
