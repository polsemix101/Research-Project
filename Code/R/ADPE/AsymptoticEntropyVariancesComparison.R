#Este programa mide la diferencia entre las varianzas asintoticas 
# de la entropia de Shannon.
#sigma.n.q es la varianza asintotica considerando la realidad 
#existente de secuencias de patrones no independientes
#sigma.n.p es la varianza asintotica asumiendo independencia 
#en la secuencia de patrones, lo cual es una simplificacion

#n es la cantidad de elementos de la serie
#m es el tamano de los patrones, este programa vale solamente para m = 3
#q es la probabilidad de aparicion de un patron
# k es m!, la cantidad de permutaciones
# tao es el embedding time delay, en este caso tao = 1
#Sigmaq es la matriz de covarianzas asint?tica del vector h(q_sombrero)
#######################################
#Matriz de probabilidades que relacionan a los patrones
# para l = 1 o sea un desplazamiento
#Q.1 y Q.2 son las matrices de probabilidades de que aparezca un patrón, 
#conociendo el anterior.

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
#################################

SigmaMatriz = function(m,q){
  #Calculo de la matriz Sigma Ec.14 y 15
  k = factorial(m)
  S = rep(0,k^2)
  SigmaMatriz = matrix(S,ncol = k)
  Q.1 = Q.1(q)
  Q.2 = Q.2(q)
  
  M = Q.1 + Q.2
  
  for (i in seq(1,k)){
    for (j in seq(1,k)){
      if(j == i){
        SigmaMatriz[i,j] <- q[i] - (2*m - 1) * q[i]^2 + 2 * M[i,i]
      }else{
        SigmaMatriz[i,j] <- -(2*m - 1) * q[i] *q[j] + (M[i,j] + M[j,i])
      }
    }
  }
  return(SigmaMatriz)
}

#################################
#calculoSigmamatriz con la formula de la Eq. 8
SigmaMatriz2 = function(m,q){
  k = factorial(m)
  
  Q.1 = Q.1(q)
  Q.2 = Q.2(q)
  
  M = Q.1 + Q.2
  Dq = diag(q)
  

  SigmaMatriz2 = Dq - (2*m - 1) * q %*% t(q) + M + t(M)
 
  return(SigmaMatriz2)
}

##################################################
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



#######################################
#sigma.n.q es la varianza asintótica de la entropía de qsombrero
#Calculo de sigma.n.q usando la primera parte de la formula 19
sigma.n.q = function(n, m, q){
  S = Sigmaq(m,q)
  k = factorial(m)
 
  Term2 = 0
  for(i in 1:(k-1)){
    for (j in (i+1):k){
      Term2 = Term2 + S[i,j]
    }
  } 
  sigma.n.q = (1/n)* sum(diag(S)) + (2/n) * Term2 
  return(sigma.n.q)
}

#Calculo de sigma.n.q usando la segunda parte de la formula 19
sigma.n.q2 = function(n, m, q){
  k = factorial(m)
  Q.1 = Q.1(q)
  Q.2 = Q.2(q)
  M = Q.1 + Q.2
  Term1 = sum((log(q)+1)^2*(q-(2*m -1 )*q^2 + 2*diag(M)))
  
  
  Term2 = 0
  for (i in 1:(k-1)){
    for (j in (i+1): k){
      Term2 = Term2 + (log(q[i]) + 1) * (log(q[j]) + 1) * ((2*m - 1)* q[i]*q[j]-
                                            (M[i,j]+ M[j,i]))
    }
  }
 
  sigma.n.q2 = (1/n) * Term1 - (2/n)* Term2
  return(sigma.n.q2)
}


#calculo de v.n.q ecuacion 31

v.n.q = function(n, m, q){
  k = factorial(m)
  Term2 = 0
  for (j in 1:(k-1)){
    for (i in (j+1):k){
      Term2 = Term2 + q[i]*q[j]*(1+log(q[i]))*(1+log(q[j]))
    }
  }
  v = (1/n) * sum(q*(1-q)*(log(q)+1)^2) - (2/n)*Term2
  return(v)
  
}

### Ahora la diferencia sigma.n.q - v.n.q
Diference = function(n, m, q){
  
  k = factorial(m)
  Q.1 = Q.1(q)
  Q.2 = Q.2(q)
  M = Q.1 + Q.2
  Term1 =  sum((log(q) + 1)^2 * ((2*m - 2) * q^2 -2 * diag(M)))
  Term2 = 0
  for (j in 1:(k-1)){
    for (i in (j+1):k){
      Term2 = Term2 + (1+log(q[i]))*(1+log(q[j])) * ((2 *m -2) * q[i]*q[j] - 
                                                       (M[i,j]+ M[j,i]))
    }
  }
  
  Diference = - Term1/n - 2 * Term2/n
    
    return(Diference)
  
  
}

##########################
#####
#Ejemplos de uso
# 1- probabilidad de alcanzar el patr?n pi_i: Equiprobable q_i = 1/k
#Calculo Sigma cuando los patrones no son independientes 

m = 3
k = factorial(m)
q1 = rep(1/k, k)
epsilon = 0


sigma11 = sigma.n.q(100,m,q1)
s21 = sigma.n.q2(100,m,q1)
#ok, me dan diferentes pero muy cercanas a cero

nu.multinomial1 = v.n.q(100,m,q1)


##################################################

# 2- probabilidad de alcanzar el patr?n pi_i: P_2, q_i = 1/k, i=1,...k-2
#q_{k-1} = 1/k  + epsilon, q_k = 1/k - epsilon, 0<= epsilon <= 1/k

m = 3
k = factorial(m)
epsilon = 1/(2*k)
q2 = rep(1/k, k)
q2[k-1] <- q2[k-1] + epsilon
q2[k] <- q2[k] - epsilon


sigma12 = sigma.n.q(1000,m,q2)
s22 = sigma.n.q2(1000,m,q2)
#ok, me dan igual

nu.multinomial2 = v.n.q(1000,m,q2)

#deberia  ser sigma12 - nu.multinomial2 = diferencia2
diferencia2 = Diference(1000, m, q2)

sigma12 - nu.multinomial2

# 3- probabilidad de alcanzar el patr?n pi_i: P_H, q_i = 1/k -  epsilon, i=1,...k/2
#q_i = 1/k  + epsilon, i = k/2 +1...k , 0<= epsilon <= 1/k

m = 3
k = factorial(m)
epsilon = 1/(2*k)
q3 = rep(1/k, k)
for (i in seq(1,k/2))
  q3[i] <- q3[i] - epsilon

for (i in seq(k/2 + 1,k))
  q3[i] <- q3[i] + epsilon

sigma13 = sigma.n.q(1000,m,q3)
s23 = sigma.n.q2(1000,m,q3)
#ok, me dan igual

nu.multinomial3 = v.n.q(1000,m,q3)
#deberia  ser sigma13 - nu.multinomial3 = diferencia3
diferencia3 = Diference(1000, m, q3)

sigma13 - nu.multinomial3

########################################
# 4- probabilidad de alcanzar el patr?n lineal pi_i: P_L, q_i = i/sum(1+...+k),
#i=1,...k

m = 3
k = factorial(m)
q4 = rep(0, k)
for (i in seq(1,k))
  q4[i] = (2*i)/(k*(k+1))
epsilon = 0 

sigma14 = sigma.n.q(100,m,q4)
s24 = sigma.n.q2(100,m,q4)
#ok, me dan igual

nu.multinomial4 = v.n.q(100,m,q4)
diferencia4 = Diference(100, m, q4)

sigma14 - nu.multinomial4

###############################################

#La diferencia segun la Ecuacion 32 s = nu + diferencia
# AR: Es la misma que Diference, ¿por qué la hacemos de nuevo?

diferencia = function(n, m, q){
  k = factorial(m)
  Q.1 = Q.1(q)
  Q.2 = Q.2(q)
  M = Q.1 + Q.2
  
  Term2 = 0
  for (j in 1:(k-1)){
    for (i in (j+1):k){
      Term2 = Term2 + (1+log(q[i]))*(1+log(q[j])) * ((2*m-2)* q[i]*q[j] - 
                                                       (M[i,j]+ M[j,i]))
      }
  }

  diferencia = - (1/n) * sum((log(q)+1)^2*((2*m-2)*q^2-2*diag(M))) - (2/n)*Term2
  return(diferencia)
  
}

Diference(100,m,q4) + nu.multinomial4
diferencia(100,m,q4) + nu.multinomial4
sigma14
