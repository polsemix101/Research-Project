#' @export



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