.logLH0_c <- function(a, a1, a2, x, n){
  # Fixed parameters
  a1 <- a1
  b1 <- a2
  # Free parameters
  p_het <- a[1]
  pi_err <- a[2]
  dbb <- choose(n, x) * beta(x+a1, n-x+b1)/beta(a1,b1)
  db <-  dbinom(x, n, pi_err)
  p_tot <- p_het * dbb + (1 - p_het) * db
  l1 <- -sum(log(p_tot))
  return(l1)
}
# Alternate model
.logLH1_c <- function(a, x, n){
  # Free parameters
  a1 <- a[1]
  b1 <- a[2]
  p_het <- a[3]
  pi_err <- a[4]
  dbb <- choose(n, x) * beta(x+a1, n-x+b1)/beta(a1,b1)
  db <-  dbinom(x, n, pi_err)
  p_tot <- p_het * dbb + (1 - p_het) * db
  l1 <- -sum(log(p_tot))
  return(l1)
}