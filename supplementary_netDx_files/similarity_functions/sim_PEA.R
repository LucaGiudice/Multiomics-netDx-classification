#Function that takes a matrix: genes (measures) x patients
#Returns a simmetric similarity matrix of the weighted jaccard between all the patients
#It has been made to work with smo_pred2.R

sim_PEA <- function(input_m) {
  #install.packages("HiClimR")
  require(HiClimR)
  PEA_m <- abs(fastCor(input_m, upperTri = FALSE))
  #PEA_m[is.nan(PEA_m)] <- 0
  return(PEA_m)
}
  

  
