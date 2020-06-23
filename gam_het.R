gam_het <- function(x) {
  
  gam <- as.data.frame(t(combn(x,m=2)),stringsAsFactors=F)
  return(sum(gam$V1 != gam$V2)/nrow(gam))
  
}
