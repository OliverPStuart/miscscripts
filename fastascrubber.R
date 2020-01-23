#!/usr/bin/env R


#This script is a great way to reduce the file-size of your .fa file after filtering.
#If you're working with a non-model organism, your .vcf might have variants called from reads mapped to a large collection of contigs.
#After filtering, your .fa file may have many contigs on which there are no variants in the .vcf.
#This file takes all the variants in your .fa file, checks which contigs they come from, then removes all others from the header.
#It assumes a contig naming convention of locus.xxxxx, where xxxxx represents a number.


library(stringr)
library(magrittr)

input.vcf <- commandArgs(trailingOnly=T)[grep(".vcf",commandArgs(trailingOnly=T))]
input.fa <- commandArgs(trailingOnly=T)[grep(".fa",commandArgs(trailingOnly=T))]

fa <- readLines(input.fa)
vcf <- readLines(input.vcf)

vcf.contigs <- paste0(">",str_extract(vcf[-grep("#",vcf)],"locus.([0-9]+)"))

contig.names <- c(1:length(fa))[fa %in% vcf.contigs]
sequences <- contig.names + 1

fa.red <- fa[c(rbind(contig.names,sequences))]

writeLines(fa.red,paste0("SCRUBBED.",input.fa))
