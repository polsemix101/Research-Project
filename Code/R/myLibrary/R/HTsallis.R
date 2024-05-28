#' Function that computes normalized Tsallis entropy of index beta of a probability function
#'
#' Given a probability function \eqn{\bm p=(p_1,p_2,\dots,p_k)}, its
#' Tsallis entropy with index \eqn{\beta \in \mathbbm{R}\setminus\{1\}}$ is
#' \deqn{
#' S^{(T_{\beta})}[\bm p] = \sum_{i=1}^k \frac{p_i - p_i^\beta}{\beta-1}.}
#'
#' @param p a probability function
#' @param beta the parameter beta for the Tsallis entropy; default 1.5
#' @usage HTsallis(p, beta)
#' @returns a value in \eqn{[0,1]}
#'
#' @export
#'

HTsallis <- function(p, beta=1.5){
  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){
  N <- length(p)
  H <- sum(p - p^beta) / (1 - N^(1-beta))

  return(H)} else {
    print("ERROR: Not a valid probability function")
  }
}

