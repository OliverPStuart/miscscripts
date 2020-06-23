###coupon collecting simulation
###used to estimate required depth of sequencing for biallelic sites in polyploids

library(ggplot2)
library(RColorBrewer)

#number of simulations to perform
n_sim <- 100000

#minor frequency
maf <- 1/6
#in the case of getting required depth for polyploid sequencing, this is 1/ploidy
#the idea being that the minimum number of copies of a minor allele at a site in a polyploid is 1
#so this sets the probability of sampling that allele at that site

#a vector of minimum depths of the minor allele to investigate
mins <- c(5,10,15,20,30)

for(m in 1:length(mins)){
  
  if(m == 1){depths <- c()}
  
  for(i in 1:n_sim){
    
    #make empty vector of waiting times 
    if(i == 1){times <- c()}
    
    #make empty draws vector
    draws <- c()
    
    #iterate until the number of minor alleles in draws reaches minimum frequency
    while(sum(draws == 1) < mins[m]){
      
      draws <- c(draws,rbinom(n=1,size=1,prob=maf))
      
    }
    
    times <- c(times,length(draws))
  
  }
  
  depths <- c(depths,times)
  
  if(m == length(mins)){
    d <- data.frame(depth=depths,
                    req=c(rep(mins,each=n_sim)))
    
    p <- ggplot(d,aes(x=depth,group=as.factor(req),fill=as.factor(req))) + 
      geom_density(alpha=0.6) + 
      scale_fill_brewer(palette="PRGn",name="minimum\nrequired\nminor allele\ndepth") +
      theme_linedraw() +
      theme(axis.ticks=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank()) +
      scale_x_continuous(limits=c(0,max(d$depth)),expand=c(0,0)) +
      scale_y_continuous(limits=c(0,NA),expand=c(0,0)) +
      labs(x="expected required coverage")
    plot(p)
    
  }
  
}

