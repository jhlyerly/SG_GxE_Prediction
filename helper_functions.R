m_to_g <- function(M) {
  pVec <- colMeans(M)/2 
  
  P <- matrix(rep(2 * (pVec - .5), nrow(M)), 
              ncol = ncol(M), byrow = T)
  
  Z <- (allM - 1) - as.matrix(P)
  
  scaleFactor <- sum(2 * (pVec * (1-pVec)))
  
  G <- (Z %*% t(Z)) / scaleFactor
}