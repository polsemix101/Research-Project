#' Statistical Complexity of a probability function
#'
#' Given a probability function \eqn{\bm p}, its statistical complexity
#' is the product of its entropy \eqn{H(\bm p)} and its distance to an
#' equilibrium distribution \eqn{\bm u}.
#'
#' The normalized Jensen-Shannon distance to the uniform distribution
#' \eqn{\bm u=(1/k, 1/k,\dots,1/k)}, where \eqn{k\geq 2} is the size of \eqn{p},
#' is \eqn{Q=Q'/\max{Q'}}, where
#' \deqn{Q'(\bm{p}, \bm{u}) = \sum_{\ell=1}^{k} \Big[\big(p(\ell)-\frac{1}{k}\big) \ln\big(k p(\ell)\big)\Big],}
#' and
#' \deqn{\max Q' = -2\Big[\frac{k!+1}{k!} \ln(k!+1) - 2 \ln(2k!) + \ln k!\Big].}
#'
#' @export
#'
#' @param p a probability function
#' @usage StatComplexity(p)
#' @returns a value in \eqn{[0,1]}


StatComplexity <- function(p){
  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){

    k <- length(p)
    kfactorial <- factorial(k)

    maxQprime <- -2 * (
      (kfactorial+1)/kfactorial * log(kfactorial+1) -
        2 * log(2*kfactorial) +
        log(kfactorial)
    )

    Qprime <- sum(
      (p-1/k)*log(k*p)
    )

    HS <- HShannon(p)

    return(Qprime/maxQprime*HS)

  } else {
    print("ERROR: Not a valid probability function")
  }
}
