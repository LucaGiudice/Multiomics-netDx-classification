sim_SPER <- function(input_m) {
  #install.packages("Hmisc")
  require(Hmisc)
  input_m=as.matrix(input_m)
  SPER_m <- rcorr(input_m,type=c("spearman"))
  SPER_m=abs(SPER_m$r)
  return(SPER_m)
}
