# Shannon Entropy --------------------------------------------------------------
shannon.entropy <- function(prob){
  entropy = prob * log2(prob)
  entropy[is.nan(entropy)] = 0
  return(-sum(entropy))
}

shannon.entropy.ln <- function(prob){
  entropy = prob * log(prob)
  entropy[is.nan(entropy)] = 0
  return(-sum(entropy))
}

shannon.entropy.log10 <- function(prob){
  entropy = prob * log(prob,10)
  entropy[is.nan(entropy)] = 0
  return(-sum(entropy))
}
shannon.entropy.normalized<- function(prob){
  entropy = (shannon.entropy(prob)/log2(length(prob)))
  entropy[is.nan(entropy)] = 0
  return(entropy)
}

# Jensen Divergence function ---------------------------------------------------
jensen.divergence <- function(prob){
  cc = rep(1/length(prob),length(prob))
  s_p = shannon.entropy(prob)
  s_q = shannon.entropy(cc)
  s_pq = shannon.entropy((prob + cc)/2)
  divergence = sum(s_pq - (s_p/2) - (s_q/2))
  return(divergence)
}
# Tsallis entropy------------------------------------------------
tsallis.entropy <- function(prob,q){
  entropy = prob - prob^q
  entropy[is.nan(entropy)] = 0
  return((1/(q-1))*sum(entropy))
}


# Statistical Complexity function ----------------------------------------------
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