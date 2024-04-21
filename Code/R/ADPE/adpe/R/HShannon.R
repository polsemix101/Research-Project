#' Asymptotic mean of the Shannon Permutation Entropy
#'
#' This function computes the asymptotic uncorrected mean of the Shannon Permutation Entropy
#' of a vector of proportions
#'
#' @param p a vector of proportions sorted according to our notation
#' @keywords asymptotic distribution, permutation entropy, Shannon entropy
#'
#' @examples
#' X <- rnorm(1000) # a time series
#' q1 <- pdc::codebook(X, m=3, t=1) # computes the ordinal patterns
#' q <- c(q1[6],q1[5],q1[4],q1[2],q1[3],q1[1]) # sorts the patterns to match our notation
#' mu <- HShannon(q) # computes the asymptotic uncorrected mean
#' nu2 <- sigma.n.q(n, m=3, q) # computes the asymptotic variance
#'
#' sigma2 <- v.n.q(n, m, q) # computes the asymptotic variance under the Multinomail law

HShannon <- function(p){

  prob <- p[p > 0]
  N <- length(p)
  H <- -sum(prob * log(prob)) / log(N)

  return(H)
}
