#' Function that computes the asymptotic variance of ordinal patterns' entropies
#'
#' @param TS time series
#' @param emb embedding dimension
#' @param ent type of entropy: S (Shannon), T (Tsallis, \eqn{beta} is required), R (Rényi, \eqn{beta} is required), or F (Fisher)
#' @param beta the parameter for the Tsallis and Rényi entropies
#' @returns The asymptotic variance
#' 
#' @export
#' 
#' @examples
#' # Time series
#' set.seed(1234567890, kind="Mersenne-Twister")
#' x <- rnorm(1000) # white noise
#' y <- mov.av(x, order=51) # smoothed with moving averages 
#' # Asymptotic variances of their Shannon entropies
#' sigma2q(x, emb=4, ent="S")
#' sigma2q(y, emb=4, ent="S")
#' 

path = "./StatOrdPattHxCModified/R/"
files = c("Sigmaq.R","Bandt-PompeModified.R")
for(i in files){
  source(paste(path,i,sep=""))  
}

sigma2qM <- function(TS, emb, ent, beta){

  # Find the number of OP
  n <- length(TS) - emb + 1
  # Compute the covariance matrix
  S <- SigmaqM(TS, emb, ent, beta)
  # Compute the
  switch(ent,

         # Shannon
         "S" = {
           sig <- sum(S) / (log(factorial(emb)))^2
         },

         # Tsallis
         "T" = {
           sig <- sum(S) / (1 - (factorial(emb))^(beta - 1))^2
         },

         # Rényi
         "R" = {
           sig <- sum(S)
         },

         # Fisher
         "F" = {
           sig <- 16 * sum(S)
         }
  )

  return(sig)
}


