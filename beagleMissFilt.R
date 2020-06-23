#!/bin/env R


# Author:  Oliver Stuart
# Date:    26/03/2020
# Project: Locust extinction

#This script takes a beagle file from the locust extinction project and filters out loci missing among spretus individuals.
#It can be modified to examine other individuals, but for this project we're only looking at the samples from M. spretus.
#The inputs are a beagle file and a double which specifies the maximum site-specific missingness for filtering.
#The beagle file and the list used to make it MUST be in the same directory. The naming conventions is:
  #names.for.file.beagle.gz
  #names.for.file.list


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
  make_option(c("-b", "--beagle"), type="character", default=NULL, 
              help="beagle file to be modified", metavar="character"),
  make_option(c("-m", "--miss"), type="integer",  default=NULL,
              help="a percentage value for missingness filtering", metavar="integer")
  ); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$beagle) | is.null(opt$miss)){
  print_help(opt_parser)
  stop("Provide a beagle file and a missingness threshold.", call.=FALSE)
}


##################################
######### Preparing data #########
##################################


#Read in beagle file
input.beagle <- read.delim(opt$beagle,
                           check.names=F,
                           header=T,
                           colClasses="character")
output.beagle <-  gsub("\\.beagle.gz","\\.filt\\.beagle",opt$beagle)
#Get associated list
input.list <- gsub("\\.beagle.gz","\\.list",opt$beagle)

#Prepare sample name data.frame for selecting spretus data in beagle file
sample.names <- readLines(input.list)
sample.names <- gsub(".o6.bam","",sample.names)
sample.names <- gsub(".o6.sub.bam","",sample.names)
sample.names <- strsplit(sample.names,split="/")
sample.names <- unlist(lapply(sample.names,`[[`,length(sample.names[[1]])))
species.map <- read.table("~/Desktop/grasshoppers/data/speciesmaps/speciesmap.txt",header=T,stringsAsFactors=F)
sample.names <- merge(data.frame(id=sample.names,ord=1:length(sample.names)),species.map,by="id")
sample.names <- sample.names[order(sample.names$ord),]


##################################
### Define missingmess function ##
##################################


indMissGet <- function(beagle){
  
  #make indices of the number of individuals and the number of sites
  n_ind <- (ncol(beagle) - 3) / 3
  n_site <- nrow(beagle)
  
  #for each individual
  for(i in 1:n_ind){
    
    #make an empty vector to store individual level missing data to begin with
    if(i == 1){ind_miss <- c()}
    
    #for each site
    for(j in 1:n_site){
      
      #make an empty vector to store site level missingness to begin with
      if(j == 1){site_miss <- c()}
      
      #if all entries for this particular individual == 0.333333
      #honey that's a missing data point
      if(sum(beagle[j,(1+i*3):(3+i*3)] == 0.333333) == 3){
        
        #add a 1 to the site missingness vector
        site_miss <- c(site_miss,1)
        
      } else {
        
        #if not, there's something there, yay, nice, cool
        site_miss <- c(site_miss,0)
        
      }
      
      #once finished looping over sites
      if(j == n_site){
        
        #take the mean of the missingness vector to get the proportion of missing sites
        ind_miss <- c(ind_miss,mean(site_miss))
        
      }
      
    }
    
    #once finished looping over individuals
    if(i == n_ind){
      
      #make a nice data.frame object to analyse
      return(data.frame(id = colnames(beagle)[seq(from=4,to=ncol(beagle)-2,by=3)],
                        miss = ind_miss,
                        stringsAsFactors=F))
      
    }
    
  }
  
}


siteMissGet <- function(beagle){
  
  #make indices of the number of individuals and the number of sites
  n_ind <- (ncol(beagle) - 3) / 3
  n_site <- nrow(beagle)
  
  #for each site
  for(i in 1:n_site){
    
    #make an empty vector to store individual level missing data to begin with
    if(i == 1){site_miss <- c()}
    
    #for each site
    for(j in 1:n_ind){
      
      #make an empty vector to store site level missingness to begin with
      if(j == 1){ind_miss <- c()}
      
      #if all entries for this particular individual == 0.333333
      #honey that's a missing data point
      if(sum(beagle[i,(1+j*3):(3+j*3)] == 0.333333) == 3){
        
        #add a 1 to the site missingness vector
        ind_miss <- c(ind_miss,1)
        
      } else {
        
        #if not, there's something there, yay, nice, cool
        ind_miss <- c(ind_miss,0)
        
      }
      
      #once finished looping over sites
      if(j == n_ind){
        
        #take the mean of the missingness vector to get the proportion of missing sites
        site_miss <- c(site_miss,mean(ind_miss))
        
      }
      
    }
    
    #once finished looping over individuals
    if(i == n_site){
      
      #make a nice data.frame object to analyse
      return(data.frame(site = as.character(beagle$marker),
                        miss = site_miss,
                        stringsAsFactors=F))
      
    }
    
  }
  
}


##################################
###### Beagle manipulation #######
##################################


#Get index of spretus files in data.frame
#Get the sample names in the beagle file based on their position
inds <- paste0("Ind",grep("spretus",sample.names$species) - 1)
#Subset the beagle file based on these names
spretus <- cbind(input.beagle[,c(1:3)],input.beagle[,colnames(input.beagle) %in% inds])

#Calculate the site-specific missingness for the spretus samples in the beagle file
spr.site.miss <- siteMissGet(spretus)

#Select the sites below the missingness threshold
good.loci <- spr.site.miss[spr.site.miss$miss < (opt$miss/100),]

#Filter out any site with all missing data in spretus
beagle.filt <- input.beagle[input.beagle$marker %in% good.loci$site,]


##################################
######### Write outputs ##########
##################################


#Write the beagle file
write.table(beagle.filt,output.beagle,col.names=T,row.names=F,quote=F,sep="\t")

#Write the report
rep <- paste0(output.beagle,".report")
cat("input file properties",sep="\n",file=rep)
cat(paste0("name = ",opt$beagle),sep="\n",file=rep,append=T)
cat(paste0("n_inds = ",(ncol(input.beagle)/3)-1),sep="\n",file=rep,append=T)
cat(paste0("n_loci = ",nrow(input.beagle)),sep="\n",file=rep,append=T)
cat(paste0("n_loci after filtering = ",nrow(beagle.filt)),sep="\n",file=rep,append=T)
cat(paste0("mean site missingness in spretus before filtering = ",mean(spr.site.miss$miss)),sep="\n",file=rep,append=T)
cat(paste0("mean site missingness in spretus before filtering = ",mean(good.loci$miss)),sep="\n",file=rep,append=T)

