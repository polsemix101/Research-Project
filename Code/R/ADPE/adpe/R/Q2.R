#' Second-lag transition matrix
#'
#' This function computes the second-lag transition matrix from
#' the observed marginal proportions q
#'
#' @param q the observed marginal proportions

Q.2 = function(q){
  Q.2.vec = c(0.05*q[1],0.05*q[1],0.15*q[1],0.15*q[1],0.30*q[1],0.30*q[1],
              0.15*q[2],0.15*q[2],0.20*q[2],0.20*q[2],0.15*q[2],0.15*q[2],
              0.05*q[3],0.05*q[3],0.15*q[3],0.15*q[3],0.30*q[3],0.30*q[3],
              0.30*q[4],0.30*q[4],0.15*q[4],0.15*q[4],0.05*q[1],0.05*q[4],
              0.15*q[5],0.15*q[5],0.20*q[5],0.20*q[5],0.15*q[5],0.15*q[5],
              0.30*q[6],0.30*q[6],0.15*q[6],0.15*q[6],0.05*q[6],0.05*q[6])
  return(matrix(Q.2.vec,ncol = 6))
}
