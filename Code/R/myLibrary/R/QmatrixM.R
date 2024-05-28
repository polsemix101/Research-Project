#' @export

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