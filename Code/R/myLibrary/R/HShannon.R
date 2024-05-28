#' Function that computes the normalized Shannon entropy (using natural logarithms) of a probability function
#'
#' Given a probability function \eqn{\bm p=(p_1,p_2,\dots,p_k)},
#' its normalized Shannon entropy is
#' \deqn{S^{(S)}[\bm p] = -\frac{1}{\ln k}\sum_{i=1}^k p_i \ln p_i.}
#' We adopt the convention \eqn{0\ln 0 = 0}.
#'
#' @usage HShannon(p)
#' @param p a probability function
#' @returns A value in \eqn{[0,1]}
#'
#' @export


HShannon <- function(p){

  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){

    prob <- p[p > 0]
    N <- length(p)
    H <- -sum(prob * log(prob)) / log(N)

    return(H)
  } else {
    print("ERROR: Not a valid probability function")
  }
}

