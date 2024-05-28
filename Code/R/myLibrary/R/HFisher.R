#' Function that computes the normalized Fisher information measure of a probability function
#'
#' Given a probability function \eqn{\bm p=(p_1,p_2,\dots,p_k)},
#' its Fisher information measure is
#' \deqn{H^{(F)}[\bm p] = 4  \sum_{i=1}^{k-1} (\sqrt{p_{i+1}} - \sqrt{p_i})^2.}
#'
#' @usage HFisher(p)
#' @param p a probability function
#' @returns a value in \eqn{[0,1]}
#'
#' @export

HFisher <- function(p){
  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){
    p1 <- p[-1]
    p2 <- p[-length(p)]

    H <- 4 * sum((sqrt(p1) - sqrt(p2))^2)

    return(H)} else {
      print("ERROR: Not a valid probability function")
    }
}


