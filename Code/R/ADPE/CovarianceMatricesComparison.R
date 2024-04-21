#Comparaci?n de la matriz de covarianza asint?tica, asumiendo que 
# los patrones ordinales son independientes (que no lo son) y asumiendo que no. 


#m = 3 (embedding), k = 6
# Calculamos Q^1 y Q^2
#Matriz de probabilidades que relacionan a los patrones
# para l = 1 o sea un desplazamiento
Q.1 = function(q){
  #solo m = 3
  Q.1.vec = c(0.25*q[1], 0.25*q[1],0,0.5*q[1],0,0,
              0,0,0.25*q[2],0,0.25*q[2],0.5*q[2],
              0.25*q[3],0.5*q[3],0,0.25*q[3],0,0,0,
              0,0.25*q[4],0,0.5*q[4],0.25*q[4],
              0.5*q[5],0.25*q[5],0,0.25*q[5],0,0,0,
              0,0.5*q[6],0,0.25*q[6],0.25*q[6])
  
  return(matrix(Q.1.vec,ncol = 6))
}
##############################
#Lo mismo que la anterior pero para l=2 o sea dos desplazamientos

Q.2 = function(q){
  Q.2.vec = c(0.05*q[1],0.05*q[1],0.15*q[1],0.15*q[1],0.30*q[1],0.30*q[1],
              0.15*q[2],0.15*q[2],0.20*q[2],0.20*q[2],0.15*q[2],0.15*q[2],
              0.05*q[3],0.05*q[3],0.15*q[3],0.15*q[3],0.30*q[3],0.30*q[3],
              0.30*q[4],0.30*q[4],0.15*q[4],0.15*q[4],0.05*q[1],0.05*q[4],
              0.15*q[5],0.15*q[5],0.20*q[5],0.20*q[5],0.15*q[5],0.15*q[5],
              0.30*q[6],0.30*q[6],0.15*q[6],0.15*q[6],0.05*q[6],0.05*q[6])
  return(matrix(Q.2.vec,ncol = 6))
}

#Covarianza real
CovarianceMatricesEq8 = function(m, q, k, epsilon){
  
  D.q = q * diag(k)
  qxqt = q %*% t(q)
  
  #necesito Q^1 y Q^2
  
  
  Q.1 = Q.1(q)
  
  Q.2 = Q.2(q)
  
  M = Q.1 + Q.2
  #Ecuaci?n 7
  Sigma = D.q -(2*m-1)*qxqt + D.q %*% M + t(M) %*% D.q
  
  return(Sigma)
  
  
  
}
########################################
#Ecuaci?n 28 Suponiendo la independencia
SEq28 = function(m, q, k, epsilon){
  D.q = q * diag(k)
  qxqt = q %*% t(q)
  S = D.q - qxqt
  return(S)
}
#Comparacion
ComparaMatrices = function (m,q,k,epsilon){
  Sigma = CovarianceMatricesEq8(m, q, k, epsilon)
  S = SEq28(m,q,k,epsilon)

  Compara = as.vector(S - Sigma)

  return( sqrt(sum(Compara*Compara)))

}
# 1- Caso1: probabilidad de alcanzar el patron pi_i: Equiprobable q_i = 1/k

m = 3
k = factorial(m)
q1 = rep(1/k, k)

difer1 = ComparaMatrices(m,q1,k,epsilon = 0)

# 2- Caso 2: probabilidad de alcanzar el patron pi_i: P_2, q_i = 1/k, i=1,...k-2
#q_{k-1} = 1/k  + epsilon, q_k = 1/k - epsilon, 0<= epsilon <= 1/k

m = 3
k = factorial(m)
epsilon = 1/(2*k)
q2 = rep(1/k, k)
q2[k-1] <- q2[k-1] + epsilon
q2[k] <- q2[k] - epsilon

difer2 = ComparaMatrices(m,q2,k,epsilon)

###############################################

# 3- Caso 3probabilidad de alcanzar el patr?n pi_i: P_H, q_i = 1/k -  epsilon, i=1,...k/2
#q_i = 1/k  + epsilon, i = k/2 +1...k , 0<= epsilon <= 1/k

m = 3
k = factorial(m)
epsilon = 1/(2*k)
q3 = rep(1/k, k)
for (i in seq(1,k/2))
  q3[i] <- q3[i] - epsilon

for (i in seq(k/2 + 1,k))
  q3[i] <- q3[i] + epsilon

difer3 = ComparaMatrices(m,q3,k,epsilon)

########################################
# 4- probabilidad de alcanzar el patron lineal pi_i: P_L, q_i = i/sum(1+...+k),
#i=1,...k

m = 3
k = factorial(m)
q4 = rep(0, k)
for (i in seq(1,k))
  q4[i] = (2*i)/(k*(k+1))
epsilon = 0 

difer4 = ComparaMatrices(m,q4,k,epsilon)

###############################################
