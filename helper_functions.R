#Takes marker matrix as 0,1,2 and returns Van Raden GRM
m_to_g <- function(M) {
  pVec <- colMeans(M)/2 
  
  P <- matrix(rep(2 * (pVec - .5), nrow(M)), 
              ncol = ncol(M), byrow = T)
  
  Z <- (allM - 1) - as.matrix(P)
  
  scaleFactor <- sum(2 * (pVec * (1-pVec)))
  
  G <- (Z %*% t(Z)) / scaleFactor
}

#Takes prcomp object and returns number of PCs needed to accumulate some proportion of variance
get_PC_n <- function(PCs, propVar) {
  PCs_var <- PCs$sdev^2 / sum(PCs$sdev^2)
  
  PCs_cumSum <- 0
  PCs_N <- 0
  while (PCs_cumSum < propVar) {
    PCs_N <- PCs_N + 1
    PCs_cumSum <- PCs_cumSum + PCs_var[PCs_N]
  }
  
  return(PCs_N)
}