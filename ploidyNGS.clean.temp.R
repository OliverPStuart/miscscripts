###cleaning up the output of ngsploidy
#considering only biallelic sites

setwd("/Volumes/Alter/BGILHISIRNA19112019/omesandindexes/transcriptome/alignments/alignedtofull/ploidyngs/")

library(dplyr)
library(magrittr)
library(ggplot2)

d <- read.table("M2_depth100.tab",header=F,
                stringsAsFactors = F)
d$V5 <- paste0(d$V1,"_",d$V2)
#shockingly, using dplyr to summarise then filter is faster than a loop
biallelic <- d %>% group_by(V5) %>% filter(V3 == "First" | V3 == "Second") %>%
  summarise(sum(V4)) %>% filter(`sum(V4)` == 100) %>% as.data.frame

print(paste0(nrow(biallelic)," of ",nrow(d)," (",
             round(nrow(biallelic)/nrow(d),digits=4)*100,
             "%) total sites analysed were biallelic"))

d <- d[d$V5 %in% biallelic$V5,]
d <- d[d$V3 == "Second" | d$V3 == "First",]
d <- data.frame(POS=d$V5,ALL=d$V3,PROP=d$V4)

ggplot(d[d$ALL == "First",],aes(x=PROP)) + 
  geom_histogram(fill="grey",colour="black") +
  scale_x_continuous(expand=c(0,0),limits=c(49,100)) +
  geom_vline(xintercept=66.66,colour="red",lty="dashed") +
  geom_vline(xintercept=83.33,colour="red",lty="dashed")

