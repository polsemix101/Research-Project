#' @export

  
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