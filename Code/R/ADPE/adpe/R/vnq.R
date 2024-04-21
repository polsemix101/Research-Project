#' Asymptotic variance of the Multinomial model
#' 
#' This function computes the asymptotic variance of the Multinomial model
#' with proportions q. It is provided as a comparison with the asymptotic
#' variance of the Ordinal Patterns (that embeds the serial correlation)
#' 
#' @param n The number of possible bins (patterns)
#' @param m The number of observed patterns
#' @param q The observed proportions (in no specific order)
#' 
#' #' X <- rnorm(1000) # a time series
#' q1 <- pdc::codebook(X, m=3, t=1) # computes the ordinal patterns
#' v.n.q() # COMPLETAR

v.n.q = function(n, m, q){
  k = factorial(m)
  Term2 = 0
  for (j in 1:(k-1)){
    for (i in (j+1):k){
      Term2 = Term2 + q[i]*q[j]*(1+log(q[i]))*(1+log(q[j]))
    }
  }
  v = (1/n) * sum(q*(1-q)*(log(q)+1)^2) - (2/n)*Term2
  return(v)
  
}
