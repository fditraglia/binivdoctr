#' Simulation DGP for binary example with independent instrument
#'
#' @param a0 Up error rate
#' @param a1 Down error rate
#' @param b ATE
#' @param c Intercept
#' @param n Sample size
#' @param d Controls non-compliance
#' @param rho Error correlation
#'
#' @return Samples from model
#' @export
#'
#' @examples
binDGP <- function(a0, a1, b = 1, c = 0, n = 1000, d = 0.15, rho = 0.5){
  n_treat <- ceiling(n/2)
  n_control <- n - n_treat
  z <- c(rep(0, n_control), rep(1, n_treat)) # offer of treatment
  errors <- MASS::mvrnorm(n, mu = c(0, 0),
                          Sigma = matrix(c(1, rho, rho, 1), 2, 2, byrow = TRUE))
  g0 <- qnorm(d)
  g1 <- qnorm(1 - d) - qnorm(d)
  Tstar <- as.numeric(g0 + g1 * z + errors[,2] > 0) #select into treatment
  y <- c + b * Tstar + errors[,1]
  #mis-classification
  Tobs <- (1 - Tstar) * rbinom(n, 1, a0) + Tstar * rbinom(n, 1, 1 - a1)
  return(data.frame(Tobs, y, z, Tstar))
}
