#!/bin/env R

#This script has functions to estimate the individual and site level proportions of missing data from a .beagle file.
#A site is considered missing if all three genotype likelihoods are calculated as 0.333333.
#This is equivalent to no information being present for that site as all possibilities are equally likely.

#Test beagle file
#b <- read.delim("/Users/oliverstuart/Desktop/grasshoppers/data/beaglefiles/subsample.sansprbru.sangred.beagle.gz")[1:5,1:15]

indMissGet <- function(beagle){
  
  #make indices of the number of individuals and the number of sites
  n_ind <- (ncol(beagle) - 3) / 3
  n_site <- nrow(beagle)
  
  #for each individual
  for(i in 1:n_ind){
    
    #make an empty vector to store individual level missing data to begin with
    if(i == 1){ind_miss <- c()}
    
    print(paste0("processing sample ",i," of ",n_ind," total"))
    
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
    #as well as a vector for printing progress
    if(i == 1){
      site_miss <- c()
      prog <- seq(from=0,to=n_site,by=1000)[-1]
      }
    
    if(i %in% prog){
    print(paste0("reached site ",i," of ",n_site," total"))
    }
    
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

missTableGet <- function(beagle){
  
  n_ind <- (ncol(beagle) - 3) / 3
  n_site <- nrow(beagle)
  
  for(i in 1:n_ind){
    
    print(paste0("processing sample ",i," of ",n_ind," total"))
    
    #make an empty data.frame to store individual level missing data to begin with
    if(i == 1){ miss <- data.frame(CHROM_POS=beagle$marker,
                                   REF=beagle$allele1,
                                   ALT=beagle$allele2)}
    
    #for each site
    for(j in 1:n_site){
      
      #make an empty vector to store site level missingness to begin with
      if(j == 1){site_miss <- c()}
      
      #if all entries for this particular individual == 0.333333
      #honey that's a missing data point
      if(sum(beagle[j,(1+i*3):(3+i*3)] == 0.333333) == 3){
        
        #add a 1 to the site missingness vector
        site_miss <- c(site_miss,0)
        
      } else {
        
        #if not, there's something there, yay, nice, cool
        site_miss <- c(site_miss,1)
        
      }
      
      if(j == n_site){
        
        miss <- cbind(miss,site_miss)
        
      }
      
    }
    
    #once finished looping over individuals
    if(i == n_ind){
      
      #make a nice data.frame object to analyse
      colnames(miss)[4:ncol(miss)] <- colnames(beagle)[seq(from=4,to=ncol(beagle)-2,by=3)]
      prop <- c()
      for(k in 1:nrow(miss)){
        
        prop <- c(prop,1-(sum(miss[k,4:ncol(miss)])/n_ind))
        
      }
      
      miss$PROP <- prop
      return(miss)
      
    }
    
  }
  
}


