#!/usr/bin/env R


#This script will take in a few inputs:
	#A reference file in fasta format, one line per sequence
	#A list of sites from that fasta file with contig, position, ref and alt allele fields
	#A replacement specification: do you want to replace the positions with the alternate or a random third character?
	#And an output filename

#This file was motivated as a way to generate alternative reference files to map to for the Melanoplus project
#I was finding reference bias in my variant calls, and in the past, one way to get around bias in the mapping stage is to switch up the alleles at the variant positions so that reads with alternates don't get downwardly biased mapping scores

#This script is designed for RADseq data =<130 bp, and was originally used .bam files mapped with bowtie


##################################
######## Initial Checks ##########
##################################


#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("optparse")
if(sum(is.installed(packages)) < length(packages)){
  print(paste0("Something is missing, check that all required packages (",paste(packages,collapse=", "),") are installed."))
  stop()
} else {
  print("Loading packages.")
}
suppressMessages(library("optparse"))

#Bring in the inputs
option_list = list(
  make_option(c("-f", "--fastaref"), type="character", default=NULL, 
              help="reference file to modify in fasta format", metavar="character"),
  make_option(c("-s", "--sitelist"), type="character",  default=NULL,
              help="file with site information, must have contig, position, ref allele, and alt allele fields", metavar="character"),
  make_option(c("-r", "--replacement"), type="character", default="alt", 
              help="type of modification to make, takes values alt (replace with alternate allele) or other (replace with a random third allele)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$fastaref) | is.null(opt$sitelist)){
  print_help(opt_parser)
  stop("You have to provide a reference fasta and a list of sites.", call.=FALSE)
}

#Brief warning
print("This script assumes that your fasta file is separated into two lines per contig/locus. If the sequence of a contig is separated over several lines this script will not work. It was designed for a fasta file of RADloci, so <150 bases really.")


##################################
######### Preparing data #########
##################################


#Read file in of positions and alt/ref
loci <- read.table(opt$sitelist,header=F)

colnames(loci) <- c("CHR","POS","REF","ALT")
loci$REF <- as.character(loci$REF) ; loci$ALT <- as.character(loci$ALT) ; loci$CHR <- as.character(loci$CHR)

#Read reference
ref <- readLines(opt$fastaref)

#Padding to make character matching easier
loci$CHR <- paste0(">",loci$CHR,":")
ref[seq(from=1,to=length(ref)-1,by=2)] <- paste0(ref[seq(from=1,to=length(ref)-1,by=2)],":")

#Make structures for selecting random nucleotides
nuc <- c("A","G","T","C")
'%ni%' <- Negate('%in%')


##################################
############ The Loop ############
##################################


#For replacement with alternative allele
if(opt$replacement == "alt"){
  
  for(i in 1:nrow(loci)){
  
    if(i == 1){ref.mod <- ref}
    substr(ref.mod[grep(loci$CHR[1],ref.mod)+1],loci$POS[i],loci$POS[i]) <- loci$ALT[i]
    
  }
  
} else {
#For replacement with a randomly chosen third allele
  
  for(i in 1:nrow(loci)){
    
    if(i == 1){ref.mod <- ref}
    
    #Sample randomly one nucleotide not in either reference or alternate
    rep <- sample(nuc[nuc %ni% loci[i,3:4]],1)
    substr(ref.mod[(grep(loci$CHR[i],ref.mod)+1)],loci$POS[i],loci$POS[i]) <- rep
    
  }
  
}


##################################
########## Write Ouput ###########
##################################


ref.mod <- gsub(":","",ref.mod)
writeLines(ref.mod,opt$out)