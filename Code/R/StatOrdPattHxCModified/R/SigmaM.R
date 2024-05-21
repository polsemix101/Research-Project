#' Function that computes the covariance matrix of ordinal patterns
#' from a time series, given embedding dimension
#'
#' @param TS time series
#' @param emb embedding dimension
#' @returns A covariance matrix

path = "./StatOrdPattHxCModified/R/"
files = c("Bandt-PompeModified.R")
for(i in files){
  source(paste(path,i,sep=""))  
}

SigmaM <- function(TS, emb){
  # Find OP probabilities
  q = formationPatternM(TS, emb,output=0)

  # Find sum of Q matrices
  k <- factorial(emb)
  Q_lag <- matrix(0, nrow = k, ncol = k)

  for (l in 1:(emb-1)){
    Qaux <- QmatrixM(formationPatternM(TS,emb,l,output=1))
    Q_lag <- Q_lag + Qaux + t(Qaux)
  }
  return(diag(q) - (2 * emb - 1) * q %*% t(q) + Q_lag)
}



