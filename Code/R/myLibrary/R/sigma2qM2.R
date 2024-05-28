#' Function that computes the asymptotic variance of ordinal patterns' entropies
#'
#' @param TS time series
#' @param emb embedding dimension
#' @param ent type of entropy: S (Shannon), T (Tsallis, \eqn{beta} is required), R (Rényi, \eqn{beta} is required), or F (Fisher)
#' @param beta the parameter for the Tsallis and Rényi entropies
#' @returns The asymptotic variance
#' 
#' @export


sigma2qM2 <- function(TS, emb, ent, beta){

  # Find the number of OP
  n <- length(TS) - emb + 1

  # Compute the covariance matrix
  S <- SigmaqM2(TS, emb, ent, beta)

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


