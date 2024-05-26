########################################################################################################
# Author: Eduarda Chagas
# Date : Jun 18, 2020
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################
# Packages and sources ---------------------------------------------------------------------------------
if(!require(gtools)){
  install.packages("gtools")
  require(gtools)
} 
library(dplyr)

# Auxiliar Function -------------------------------------------------------------------------------------

define.symbols <- function(D){
  d = c(1:D)
  symbol = matrix(unlist(permutations(n = D, r = D, v = d)), nrow = factorial(D), ncol = D, byrow = FALSE)
  # symbol = matrix(unlist(permutations(n = D)), nrow = factorial(D), ncol = D, byrow = FALSE)
  symbol = symbol - 1
  symbol
}

FP <- function(n, dimension, delay){
  dyn.load("FormationPatterns.so")
  p <- .Call("FormationPatterns", n, dimension, delay)
  p = t(p) + 1
  return(p)
}

percentual.equalities <- function(patterns){
  n.patterns = dim(patterns)[1]
  n.duplicated = 0
  for(i in 1:n.patterns){
    if(length(which(duplicated(patterns[i,]) == TRUE))){
      n.duplicated = n.duplicated + 1
      print(length(which(duplicated(patterns[i,]) == TRUE)))
    }
  }
  return(n.duplicated/n.patterns)
}


formationPattern <- function(series, D, tau, option){
  
  i = 1
  n = length(series)
  p_patterns = elements = matrix(nrow = n, ncol = D)
  index = c(0:(D-1))
  
  for(s in seq(1, length(series)-(D-1)*tau, by = 1)){
    # the indices for the subsequence
    ind = seq(s, s+(D-1)*tau, by = tau)
    elements[i,] = series[ind]
    p_patterns[i,] = index[order(elements[i,])]
    i = i + 1
  }
  
  if(option == 0){
    p_patterns = na.omit(p_patterns)
    return(p_patterns[1:(i-1),])
  }else if(option == 1){
    elements = na.omit(elements)
    return(elements[1:(i-1),])    
  }
}

#Only used for magnus' version
identicalValues <- function(elements){
  ranked = rank(elements)
  unique = unique(ranked)
  len = length(ranked)
  weight = 1
  perm=1
  if (length(unique)!=len){
    permList = c()
    indexList = c()
    for (i in 1:length(unique)){
      index = which(ranked == unique[i])
      n = length(index)
      perm = perm * factorial(n)
      perms = permutations(n,n,v=(unique[i]-(n-1)/2):(unique[i]+(n-1)/2))
      permList = append(permList,list(perms))
      indexList = append(indexList,list(index))
    }
    weight = 1/perm
    finalPatterns = matrix(rep(0,len),nrow=1,ncol=len)
    for (i in 1:length(permList)){
      nPerms = nrow(permList[[i]])
      temp = indexList[i][[1]]
      nrowBeforeExtension = nrow(finalPatterns)
      patterns = finalPatterns
      if (nPerms!=1){
        for(o in 1:(nPerms-1)){
          finalPatterns = rbind(finalPatterns,patterns)
        }
      }
      for (j in 1:length(temp)){
        value = permList[[i]][,j]
        for (p in 1:nPerms){
          for(k in 1:nrowBeforeExtension){
            rowIndex = k+(p-1)*nrowBeforeExtension
            finalPatterns[rowIndex,temp[j]]=value[p]
          }
        }
      }
    }
    return(cbind(rep(weight,perm),finalPatterns))
  } else {
    return(c(weight,ranked))
  }
  
}

#Only used for magnus' version
ordinal_pattern <-function(weightAndPatterns){
  d = ncol(weightAndPatterns)-1
  possiblePatterns = permutations(d,d,v=1:d)
  weight = weightAndPatterns[,1]
  df = as.data.frame(weightAndPatterns[,-c(1)])
  weights = c()
  for (i in 1:factorial(d)){
    row = possiblePatterns[i,]
    index = 1:(nrow(df))
    for(j in 1:d){
      index = intersect(index,which(df[,j]==row[j]))
    }
    sum = sum(weight[index])
    weights = append(weights,sum)
  }
  return(weights)
}

QmatrixM <- function(seq,lag=1){
  d = ncol(seq)-1
  n = nrow(seq)
  i=1
  Q = matrix(0,nrow=factorial(d),ncol=factorial(d))
  possiblePatterns = permutations(d,d,v=1:d)
  inp = c() #index of new pattern
  while(i<=n){
    np = round(seq[i,1]**(-1)) #number of possible patterns at pattern i
    inp = append(inp,i)
    i=i+np
  }
  
  for(i in 1:(length(inp)-lag)){
    w1=seq[inp[i],1] #weight of each pattern in first pattern index
    n_p1 = round(w1**(-1)) #number of possible patterns starting from i
    p1 = seq[i:(i+n_p1-1),(2:(d+1))] #possible patterns 1
    
    p1index = c()
    row = p1
    for(j in 1:n_p1){
      if(n_p1>1){
        row = p1[j,]
      }
      index = 1:factorial(d)
      for(k in 1:d){
        index = intersect(index,which(row[k]==possiblePatterns[,k]))
      }
      p1index = append(p1index,index)
    }
    
    
    j = inp[i+lag] #find possible pattern 2
  
    w2=seq[j,1]
    n_p2 = round(w2**(-1))
    p2 = seq[j:(j+n_p2-1),(2:(d+1))]
    p2index = c()
    row = p2
    for(j in 1:n_p2){
      if(n_p2>1){
        row = p2[j,]
      }
      index = 1:factorial(d)
      for(k in 1:d){
        index = intersect(index,which(row[k]==possiblePatterns[,k]))
      }
      p2index = append(p2index,index)
    }
    
    
    wQ = w1*w2 #weight of each transition in Q - they will all always be the same value
    for(i1 in 1:length(p1index)){
      for(i2 in 1:length(p2index)){
        Q[p1index[i1],p2index[i2]] = Q[p1index[i1],p2index[i2]]+wQ
      }
    }
  }
  return(Q/sum(Q))
}


#Only used for magnus' version, splits ties evenly among possible patterns
formationPatternM <- function(series, D, tau=1,output=0){ #Output 0 gives the ordinal pattern distribution, output 1 gives the sequences of patterns
  i=1
  n = length(series)
  p_patterns = matrix(nrow = 0, ncol = D+1)
  elements = matrix(nrow = 0, ncol = D)
  index = c(0:(D-1))
  
  for(s in seq(1, length(series)-(D-1)*tau, by = 1)){
    # the indices for the subsequence
    ind = seq(s, s+(D-1)*tau, by = tau)
    elements = rbind(elements,series[ind])
    p_patterns = rbind(p_patterns,identicalValues(elements[i,]))
    i=i+1
  }
  p_patterns = na.omit(p_patterns)
  if(output==1){
    return(p_patterns)
  } else {
    temp = ordinal_pattern(p_patterns)
    return(temp/sum(temp))
  }
}


# Bandt-Pompe function ---------------------------------------------------------------------------------

bandt.pompe.optimizated <- function(series, dimension, delay){
  dyn.load("BandtPompe.so")
  elements = formationPattern(series, dimension, delay, 1)
  element.size = dim(elements)[1]
  probability <- .Call("BandtPompe", elements, dimension, element.size)
  return(probability)
}

bandt.pompe <- function(serie, dimension, delay){  
  fat = factorial(dimension)
  probability = rep(0, fat)
  symbols = define.symbols(dimension)
  p_patterns = formationPattern(serie, dimension, delay, 0)
  n_symbols = dim(p_patterns)[1]
  for(j in 1:n_symbols){
    for(i in 1:fat){
      if(all(p_patterns[j,] == symbols[i,])){ 
        probability[i] = probability[i] + 1
        break
      }
    }
  }
  return(probability/n_symbols)
}


# Shannon Entropy function -------------------------------------------------------------------------------
shannon.entropy <- function(prob){
  entropy = prob * log(prob)
  entropy[is.nan(entropy)] = 0
  return(-sum(entropy))
}

shannon.entropy.normalized <- function(prob){
  entropy = (shannon.entropy(prob)/log(length(prob)))
  entropy[is.nan(entropy)] = 0
  return(entropy)
}

# Jensen Divergence function ------------------------------------------------------------------------------
jensen.divergence <- function(prob){
  cc = rep(1/length(prob),length(prob))
  s_p = shannon.entropy(prob)
  s_q = shannon.entropy(cc)
  s_pq = shannon.entropy((prob + cc)/2)
  divergence = sum(s_pq - (s_p/2) - (s_q/2))
  return(divergence)
}

# Statistical Complexity function ---------------------------------------------------------------------------

constant <- function(prob){
  k = (0.5)/length(prob)
  a1 = (0.5 + k) * log(0.5 + k)
  a2 = (length(prob) - 1) * k * log(k)
  a3 = (1 - 0.5) * log(length(prob))
  b = -1/(a1 + a2 + a3)
  return(b)
}

Ccomplexity<-function(prob){
  cc = jensen.divergence(prob) * constant(prob) * shannon.entropy.normalized(prob)
  return(cc)
}

# Trozos functions ---------------------------------------------------------------------------

cotas <- function(dimension, option = 0){
  
  if(option == 0){
    cx = readingMPR(dimension,1)
    cy = readingMPR(dimension,2)
    
  }else{
    cx = readingMPR(dimension,3)
    cy = readingMPR(dimension,4)
    
  }
  cotas.data = data.frame("cx" = cx, "cy" = cy)
  return(cotas.data)
}


readingMPR<-function(dimension, option=0){
  if(dimension == 3){ 
    continua = "../../Data/trozos/continuaN6.txt"
    trozo = "../../Data/trozos/trozosN6.txt"
  }
  if(dimension == 4){ 
    continua = "../../Data/trozos/continuaN24.txt"
    trozo = "../../Data/trozos/trozosN24.txt"
  }
  if(dimension == 5){ 
    continua = "../../Data/trozos/continuaN120.txt"
    trozo = "../../Data/trozos/trozosN120.txt"
  }
  if(dimension == 6){ 
    continua = "../../Data/trozos/continuaN720.txt"
    trozo = "../../Data/trozos/trozosN720.txt"
  }
  if(dimension == 36){ 
    continua = "../../Data/trozos/continuaN36.txt"
    trozo = "../../Data/trozos/trozosN36.txt"
  }
  if(dimension == 576){ 
    continua = "../../Data/trozos/continuaN576.txt"
    trozo = "../../Data/trozos/trozosN576.txt"
  }
  if(dimension == 14400){ 
    continua = "../../Data/trozos/continuaN14400.txt"
    trozo = "../../Data/trozos/trozosN14400.txt"
  }
  if(dimension == 518400){ 
    continua = "../../Data/trozos/continuaN518400.txt"
    trozo = "../../Data/trozos/trozosN518400.txt"
  }
  curva1x = read.table(continua, stringsAsFactors=FALSE, fileEncoding="latin1")[,1]
  if(option==1) return(curva1x)
  curva1y = read.table(continua, stringsAsFactors=FALSE, fileEncoding="latin1")[,2]
  if(option==2) return(curva1y)
  curva2x = read.table(trozo, stringsAsFactors=FALSE, fileEncoding="latin1")[,1]
  if(option==3) return(curva2x)
  curva2y = read.table(trozo, stringsAsFactors=FALSE, fileEncoding="latin1")[,2]
  if(option==4) return(curva2y)
}

histogram <- function(series, D, tau = 1, max_y, title, x_label){
    fat = factorial(D)
    p.patterns = formationPattern(series, D, tau, 0)
    n.symbols = dim(p.patterns)[1]
    symbol = define.symbols(D)
    index.rep = array(0, n.symbols)
    
    for(i in 1:n.symbols){
      for(j in 1:fat){
        if(all(p.patterns[i,] == symbol[j, ])){
          index.rep[i]=j
          break
        }
      }
    }
    
    index.rep = index.rep[1:n.symbols]
    index.rep = data.frame(i = index.rep)
    
    p = ggplot(index.rep) +
      geom_histogram(aes(x = i, y = ..density..), binwidth = 1, fill = "grey", color = "black") + 
      ggtitle(title) + 
      ylim(0, max_y) +
      labs(x=x_label, y="") +
      theme_few(base_size = 12, base_family = "sans") +  
      theme(plot.title = element_text(hjust=0.5)) 
  return(p)
}

hist_prop <- function(series, D, tau = 1){
  fat = factorial(D)
  p.patterns = formationPattern(series, D, tau, 0)
  n.symbols = dim(p.patterns)[1]
  symbol = define.symbols(D)
  index.rep = array(0, n.symbols)
  
  for(i in 1:n.symbols){
    for(j in 1:fat){
      if(all(p.patterns[i,] == symbol[j, ])){
        index.rep[i]=j
        break
      }
    }
  }
  
  index.rep = index.rep[1:n.symbols]
  index.rep = data.frame(i = index.rep)
  return(index.rep)
}

histD3 <- function(series){

  p.patterns = formationPattern(series, D = 3, tau = 1, 0)
  n.symbols = dim(p.patterns)[1]
  symbol = matrix(c(0,1,2,
                    0,2,1,
                    1,0,2,
                    1,2,0,
                    2,0,1,
                    2,1,0), ncol = 3, byrow = TRUE)
  index.rep = array(0, n.symbols)

  for(i in 1:n.symbols){
    for(j in 1:6){
      if(all(p.patterns[i,] == symbol[j, ])){
        index.rep[i]=j
        break
      }
    }
  }

  index.rep = index.rep[1:n.symbols]
  index.rep = data.frame(i = index.rep)
  return(index.rep)
}