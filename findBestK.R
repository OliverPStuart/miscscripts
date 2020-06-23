#!/bin/env R

# Author: Oliver Stuart
# Date: 27/03/2020

#This script runs NGSadmix on a beagle file and uses the Evanno method to find the best K values
#Inputs are:
  #a beagle file
  #a minor allele frequency threshold
  #a maximum K to analyse
  # and the number of replicate analyses
#The call to ngsadmix in the first loop should be modified for wherever ngsadmix is in the system
#All outputs will be saved in a subdiretory named for the minor allele frequency and the beagle file input


start.time <- Sys.time()


##################################
######## Initial Checks ##########
##################################


#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("optparse","stringr","dplyr","ggplot2")
if(sum(is.installed(packages)) < length(packages)){
  print(paste0("Something is missing, check that all required packages (",paste(packages,collapse=", "),") are installed."))
  stop()
}
suppressMessages(library("optparse"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))

#Bring in the inputs
option_list = list(
  make_option(c("-b", "--beagle"), type="character", default=NULL, 
              help="beagle file to be analysed", metavar="character"),
  
  make_option(c("-m", "--maf"), type="integer",  default=5,
              help="minor allele frequency threshold", metavar="integer"),
  
  make_option(c("-r", "--Nrep"), type="integer",  default=10,
              help="the number of replicate analyses", metavar="integer"),
  
  make_option(c("-k", "--maxK"), type="integer",  default=10,
              help="the maximum value of K to use in the analyse", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Do the inputs actually exist?
if (is.null(opt$beagle)){
  print_help(opt_parser)
  stop("Provide a beagle file.", call.=FALSE)
}

if(opt$Nrep < 10){
  
  print("It is strongly recommended to run the analysis for at least 10 reps. Otherwise the program may crash.")
  
}


##################################
######## Running analyses ########
##################################


#Get variables
beagle <- opt$beagle
maf <- opt$maf
maxK <- opt$maxK
Nrep <- opt$Nrep

#Make output directory
out.dir <- paste0(getwd(),"/",tail(unlist(strsplit(beagle,split="/")),n=1),".",maf/100,".bestK")
dir.create(path=out.dir)

#Load the beagle file just to get some stats from it
Nind <- (ncol(read.delim(beagle)) - 3) / 3

#Perform ngsadmix analyses
for (j in 1:maxK) {
  
  for (i in 1:Nrep) {
    
    system(paste0("/Users/oliverstuart/ngsadmix/ngsadmix -likes ", beagle," -minMaf ", maf/100," -K ", 
                 j, " -o ",out.dir, "/allindsK", j, "rep", i, " -minInd ", ceiling(Nind*0.7), 
                 sep = ""))
    
  }
  
}


##################################
## Making datasets for analysis ##
##################################


#Make a list of all the log files
logs <- list.files(out.dir)[grep("\\.log",list.files(out.dir))]

#Extract from them the information for a data.frame with columns K, rep, like
Ks <- as.numeric(gsub("K","",str_extract(logs,"K([0-9]+)")))
Ns <- as.numeric(gsub("rep","",str_extract(logs,"rep([0-9]+)")))
options(digits=10)
for(i in 1:length(logs)){
  if(i == 1){likes <- c()}
  
  likes <- as.numeric(c(likes,
             str_extract(readLines(paste0(out.dir,"/",logs[i]))[9],
                                    "-([0-9]+)\\.([0-9]+)")))
  
}

#Very inelegant reshuffling of a dataframe please don't look Sasha
likes <- data.frame(K = Ks, rep = Ns, likes = likes)
likes <- likes[order(likes$K,likes$rep),]

#Briefly make a plot of likelihoods values for each K
likelihoodplot <- ggplot(likes,aes(x=as.factor(K),y=likes)) + geom_boxplot() +
  theme_linedraw() + labs(x="K",y="log Likelihoods")

#Reshuffle the data.frame some more to put it into usable format
likes <- data.frame(t(as.data.frame(split(likes$likes,likes$K))),
                    row.names = paste0("K",1:maxK),stringsAsFactors=F)
colnames(likes) <- paste0("Rep",1:Nrep)




##################################
###### Analyse the results #######
##################################


#A loop which iterates over columns and calculates L"(K) for all K values (except 1 and maxK)

#For all reps
for(i in 1:Nrep){
  
  #For every K values
  for(k in 1:maxK){
    
    #At the first instance make an empty vector to contain the new L"(K) values for that rep
    if(k == 1){newcol <- c()}
    
    #If K is either 1 or the maximum
    if(k == 1 | k == maxK){
      
      #The whole column is NA
      newcol <- c(newcol,NA)
      
      #Other wise (K = 2:(maxK-1))
    } else {
      
      #The change in K is the absolute value of the next K minus twice the sum of this plus the previous K
      LK <- abs(likes[k+1,i]-2*likes[k,i]+likes[k-1,i])
      #And plonk it onto the data.frame
      newcol <- c(newcol,LK)
      
    }
    
    #At the last K for this rep
    if(k == maxK){
      
      #Put the vector onto the original likes data.frame
      likes <- cbind(likes,newcol)
      
    }
    
  }
  
  #Now, in the last instance
  if(i == Nrep){
    
    #Rename the columns to make it pretty
    colnames(likes)[(1+Nrep):ncol(likes)] <- paste0("L''(K) Rep",1:Nrep)
  
    #And since dplyr doesn't want to behave we're using another loop
    for(k in 1:maxK){
      
      #In the first instance
      if(k == 1){
        
        #Make some placeholder vectors
        mean <- c()
        std <- c()
        deltaK <- c()
        
      }
      
      #Extract and unlist the L''(K) values
      values <- unlist(likes[k,(1+Nrep):ncol(likes)])
      
      #Take their means and standard deviations and add them to the respective
      mean <- c(mean,mean(values))
      std <- c(std,sd(values))
      
      #Calculate delta K following Evannos et al.
      deltaK <- mean / std
      
      #At the last instances
      if(k == maxK){
        
        #Create new columns for the new values
        likes$Mean <- mean
        likes$SD <- std
        likes$deltaK <- deltaK
        
      }
      
    }
    
  }
  
}


##################################
####### Write the outputs ########
##################################


#Write the likelihood output table
write.table(likes,file=paste0(out.dir,"/likelihood.table"),
            quote=F,sep="\t",
            row.names=T,col.names=T)

#Write a small output report
logfile <- paste0(out.dir,"/out.log")
cat("Analysis log",sep="\n",file=logfile)
cat("#########################",sep="\n",file=logfile,append=T)
cat(paste0("File analysed = ",tail(unlist(strsplit(beagle,split="/")),n=1)),sep="\n",file=logfile,append=T)
cat(paste0("Time taken = ",round(Sys.time() - start.time,digits=3)," seconds"),sep="\n",file=logfile,append=T)
cat("#########################",sep="\n",file=logfile,append=T)
cat("Best K = ",which.max(likes$deltaK),sep="\n",file=logfile,append=T)

#Make a figure of deltaK to print out
plotting <- data.frame(K = 1:maxK,
                       deltaK = likes$deltaK,
                       mean = likes$Mean,
                       std = likes$SD)

deltaKplot <- ggplot(plotting,aes(x=K,y=deltaK)) + geom_point(na.rm=TRUE) + geom_line(na.rm=TRUE) + 
  theme_linedraw() + theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(1:(Nrep-1)),limits=c(1,(maxK-1))) +
  scale_y_continuous(limits=c(0,
                              plotting$deltaK[which.max(plotting$deltaK)]+plotting$deltaK[which.max(plotting$deltaK)]*0.02),
                     expand=c(0,0))


if(is.element("patchwork", installed.packages()[,1])){
  
  suppressMessages(library(patchwork))
  outplot <- likelihoodplot | deltaKplot
  
  png(paste0(out.dir,"/summary.plot.png"),res=300,height=8,width=8,units="in")
  plot(outplot)
  dev.off()
  
} else {
  
  png(paste0(out.dir,"/likelihood.plot.png"),res=300,height=8,width=4,units="in")
  plot(likelihoodplot)
  dev.off
  
  png(paste0(out.dir,"/deltaK.plot.png"),res=300,height=8,width=4,units="in")
  plot(deltaKplot)
  dev.off()
  
}
