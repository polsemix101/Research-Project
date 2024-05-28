#' @export


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