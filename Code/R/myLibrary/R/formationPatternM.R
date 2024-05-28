#' @export

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