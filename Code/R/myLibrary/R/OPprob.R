#' Function that finds the Ordinal Patterns probabilities in a time series for a given embedding dimension
#'
#' @usage OPprob(TS, emb)
#' @param TS time series of length \eqn{n-m+1}
#' @param emb embedding dimension \eqn{m}
#' @returns a sequence of \eqn{n} patterns
#'
#' @import tibble
#' @import dplyr
#' @import prodlim
#' @importFrom stats pnorm
#'
#' @export



OPprob <- function(TS, emb){

  op <- tibble(
    OP = factor(OPseq(TS, emb), levels = 1:factorial(emb))
  )

  fr <- op %>% count(OP, .drop = FALSE)
  return(fr$n / sum(fr$n))
}


