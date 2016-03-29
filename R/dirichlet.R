rdirichlet <- function(n, shape){
  K <- length(shape)
  out <- matrix(NA_real_, K, n)
  for(i in 1:n){
    gamma_draws <- rgamma(K, shape, scale = 1)
    out[,i] <- gamma_draws / sum(gamma_draws)
  }
  return(t(out[-K,]))
}


