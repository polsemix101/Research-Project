#' Function that computes the covariance matrix for a time series
#'
#' @param TS time series
#' @param emb embedding dimension
#' @param ent type of entropy: S (Shannon), T (Tsallis, \eqn{beta} is required), R (Rényi, \eqn{beta} is required), or F (Fisher)
#' @param beta parameter for the Tsallis and Rényi entropies
#'


SigmaqM2 <- function(TS, emb, ent, beta){
  # Compute Sigma matrix
  S <- SigmaM2(TS, emb)

  # Compute vector of probabilities
  qvec <- OPprob(TS, emb)

  # Compute matrix of partial derivatives
  J <- PDmatrixM(q = qvec, ent, beta)
  
  return(J %*% S %*% t(J))
}



