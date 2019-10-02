#Function that takes a matrix: genes (measures) x patients
#Returns a simmetric similarity matrix of the weighted jaccard between all the patients
#It has been made to work with smo_pred2.R

sim_WJ <- function(input_m) {
  require(data.table)
  require(matrixStats)
  
  samples <- 1:ncol(input_m)
  comb <- CJ(samples, samples)
  comb[, i := .I]
  comb <- melt(comb, 'i')
  setorder(comb, value)
  v2 <- paste0("V", 1:2)
  comb[, variable2 := v2 , keyby = i]
  comb2 <- dcast(comb, i ~ variable2, value.var = 'value')
  combUnique <- unique(comb2, by = c('V1', 'V2'))
  
  XX <- apply(combUnique[, -'i'], 1, function(x) {
    x2 <- rowRanges(input_m, cols = x)
    s <- colSums2(x2)
    s[1] / s[2]
  })
  
  set(combUnique, j = 'xx', value = XX)
  rez2 <- merge(comb2, combUnique[, -'i'], by = c('V1', 'V2'), all.x = T)
  setorder(rez2, i)
  rez2 <- array(rez2$xx, dim = rep(ncol(input_m), 2))
  rownames(rez2) <- colnames(input_m)
  colnames(rez2) <- colnames(input_m)
  return(rez2)
}
  

  