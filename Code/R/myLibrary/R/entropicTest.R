#' Function that applies of the hypothesis test for two time series,
#' a given embedding dimension, and
#' a given entropy ent = S (Shannon), T (Tsallis, beta is required),
#' R (Rényi, beta is required), or F (Fisher information measure)

#' @param TSx a time series
#' @param TSy another time series
#' @param emb the embedding dimension
#' @param ent the type of entropy: S (Shannon), T (Tsallis, \eqn{\beta} is required), R (Rényi, \eqn{\beta} is required), or F (Fisher)
#' @param beta the beta parameter, only required for Tsallis and Rényi entropies
#' @param pvalue (optional, default 0.05) the \eqn{p}-value
#'
#' @export


entropicTest <- function(TSx, TSy, emb, ent, beta, pvalue=0.05){

  pv <- pval(TSx, TSy, emb, ent, beta)

  if(pv < pvalue){
    print(paste0("The null hypothesis is rejected with p-value = ", pv, "."))
  }
  else{
    print(paste0("The null hypothesis is not rejected with p-value = ", pv, "."))
  }
}





