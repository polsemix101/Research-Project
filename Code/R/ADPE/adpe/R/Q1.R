#' First-lag transition matrix
#'
#' This function computes the first-lag transition matrix from
#' the observed marginal proportions q
#'
#' @param q the observed marginal proportions

Q.1 = function(q){
  #solo m = 3
  Q.1.vec = c(0.25*q[1], 0.25*q[1],0,0.5*q[1],0,0,
              0,0,0.25*q[2],0,0.25*q[2],0.5*q[2],
              0.25*q[3],0.5*q[3],0,0.25*q[3],0,0,
              0,0,0.25*q[4],0,0.5*q[4],0.25*q[4],
              0.5*q[5],0.25*q[5],0,0.25*q[5],0,0,
              0,0,0.5*q[6],0,0.25*q[6],0.25*q[6])

  return(matrix(Q.1.vec,ncol = 6))
}
