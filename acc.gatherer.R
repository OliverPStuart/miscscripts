#!/usr/bin/env R

#This is a basic script which takes in the start and end numbers and letter code of a series of GenBank accessions.
#The output is a text file with all the accessions from start to end on separate rows.
#This can be submitted to GenBank to download a custom set of records.

print("Arguments have to be in the right order: firstnumber secondnumber lettercode outputname")
print("Otherwise the output will be nonsensical.")

inputs <- commandArgs(trailingOnly=T)

seqs <- paste0(inputs[3],seq(from=inputs[1],to=inputs[2]))

write.table(seqs, file=inputs[4],sep="\t",row.names=F,col.names=F,quote=F)
