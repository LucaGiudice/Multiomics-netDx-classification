sim_HOE <- function(input_m) {
  #install.packages("Hmisc")
  require(Hmisc)
  input_m=as.matrix(input_m)
  HOE_m <- hoeffd(input_m)
  HOE_m=abs(HOE_m$D)
  return(HOE_m)
}
