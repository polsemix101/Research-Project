#' Sigma Matrix 2
#'
#' This is an internal function
#' Computes the covariance matrix Sigma
#'
#'


Sigmaq = function(m, q){
  #Ecuacion 13

  k = factorial(m)
  S = rep(0,k^2)

  Sigmaq = matrix(S,ncol = k)

  SigmaMatriz = SigmaMatriz2(m,q)
  #Calculo de Sigmaq Ec. 13

  for (i in seq(1,k)){
    for (j in seq(1,k)){
      if(j == i){
        Sigmaq[i,j] <- (log(q[i]) + 1)^2 * SigmaMatriz[i,i]
      }else{
        Sigmaq[i,j] <- (log(q[i]) + 1) * (log(q[j]) + 1) * SigmaMatriz[i,j]
      }
    }
  }


  return(Sigmaq)
}

