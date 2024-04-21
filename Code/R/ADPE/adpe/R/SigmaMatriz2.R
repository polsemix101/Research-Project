#' Sigma Matrix 2
#'
#' This is an internal function
#' Computes the covariance matrix Sigma2
#'
#'

SigmaMatriz2 = function(m,q){
  k = factorial(m)

  Q.1 = Q.1(q)
  Q.2 = Q.2(q)

  M = Q.1 + Q.2
  Dq = diag(q)


  SigmaMatriz2 = Dq - (2*m - 1) * q %*% t(q) + M + t(M)

  return(SigmaMatriz2)
}
