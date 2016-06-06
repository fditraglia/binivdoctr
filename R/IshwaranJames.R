#' Gibbs sampler for Dirichlet Process Mixture of Normals
#'
#' @param X Vector of input data.
#' @param nDraws Number of draws from the Gibbs sampler (exluding burn-in).
#' @param burnIn Number of burn-in draws.
#' @param maxClust Maximum number of mixture components.
#' @param A Prior parameter.
#' @param xi_sq Prior parameter.
#' @param eta1 Prior parameter.
#' @param eta2 Prior parameter.
#' @param nu1 Prior parameter.
#' @param nu2 Prior parameter.
#'
#' @return List of MCMC draws for mu (means of mixture components), sigma_sq
#' (variances of mixture components), p (mixing probabilities), alpha, theta,
#' and K (assignment of observations to clusters).
#' @export
#'
#' @examples
drawDPM <- function(X, nDraws = 5000, burnIn = 1000, maxClust = 50,
                    A = 1000, xi_sq = NULL, eta1 = 2, eta2 = 2,
                    nu1 = 2, nu2 = 2){

  # Default to "Empirical Bayes" for xi_sq
  if(is.null(xi_sq)) xi_sq <- (4 * sd(X))^2

  nObs <- length(X)

  # object to store the draws
  draws <- list(mu = matrix(NA_real_, nrow = maxClust, ncol = nDraws),
                sigma_sq = matrix(NA_real_, nrow = maxClust, ncol = nDraws),
                K = matrix(NA_integer_, nrow = nObs, ncol = nDraws),
                p = matrix(NA_real_, nrow = maxClust, ncol = nDraws),
                alpha = rep(NA_real_, nDraws),
                theta = rep(NA_real_, nDraws))

  # initialize the sampler
  alpha <- eta1 / eta2 # prior mean
  theta <- mean(X)
  K <- factor(sample(1:maxClust, nObs, replace = TRUE), 1:maxClust)
  clusters <- tapply(X, K, c)
  empty <- sapply(clusters, is.null)
  mu <- tapply(X, K, mean)
  sigma_sq <- tapply(X, K, var)
  # Handle clusters with only one observation
  sigma_sq[!empty & is.na(sigma_sq)] <- var(X)
  obsPerClust <- tapply(X, K, length)
  obsPerClust[is.na(obsPerClust)] <- 0
  p <- (1 / nObs) * obsPerClust

  # loop over draws of the sampler
  for(i in 1:(burnIn + nDraws)){

    # update mu
    s_sq <- rep(NA_real_, maxClust)
    m <- rep(NA_real_, maxClust)
    s_sq[!empty] <- 1 / (obsPerClust[!empty] / sigma_sq[!empty] + 1 / xi_sq)
    m[!empty] <- s_sq[!empty] * (obsPerClust[!empty] * (theta / xi_sq) +
      sapply(clusters[!empty], sum) / sigma_sq[!empty])
    m[empty] <- theta
    s_sq[empty] <- xi_sq
    mu <- rnorm(maxClust, mean = m, sd = sqrt(s_sq))

    # update sigma_sq
    rho <- rep(NA_real_, maxClust)
    rho[!empty] <- mapply(function(y, z) sum((y - z)^2 / 2),
                          clusters[!empty], mu[!empty])
    rho[empty] <- 0
    sigma_sq <- 1 / rgamma(maxClust, shape = nu1 + obsPerClust / 2,
                           rate = nu2 + rho)

    # update cluster membership
    pi <- matrix(rep(p / sqrt(sigma_sq), times = nObs) *
                 exp(-0.5 * (rep(X, each = maxClust) - rep(mu, times = nObs))^2 /
                 rep(sigma_sq, times = nObs)), nrow = maxClust, ncol = nObs,
                 byrow = FALSE)
    K <- factor(apply(pi, 2, function(x) sample(1:maxClust, size = 1, prob = x)),
                1:maxClust)
    obsPerClust <- tapply(X, K, length)
    obsPerClust[is.na(obsPerClust)] <- 0
    clusters <- tapply(X, K, c)
    empty <- sapply(clusters, is.null)

    # update p
    V <- c(rbeta(maxClust - 1, shape1 = 1 + obsPerClust[-maxClust],
               shape2 = alpha + cumsum(obsPerClust)[-1]), 1)
    log_remain <- c(0, cumsum(log(1 - V[-maxClust])))
    p <- exp(log(V) + log_remain)


    # update alpha
    alpha <- rgamma(1, shape = maxClust + eta1 - 1,
                    rate = eta2 - log_remain[(maxClust - 1)])

    # update theta
    s_theta_sq <- 1 / (maxClust / xi_sq + 1 / A)
    m_theta <- (s_theta_sq / xi_sq) * sum(mu)
    theta <- rnorm(1, mean = m_theta, sd = sqrt(s_theta_sq))

    # store the draws
    if(i > burnIn){
      drawIndex <- i - burnIn
      draws$mu[,drawIndex] <- mu
      draws$sigma_sq[,drawIndex] <- sigma_sq
      draws$K[,drawIndex] <- K
      draws$p[,drawIndex] <- p
      draws$alpha[drawIndex] <- alpha
      draws$theta[drawIndex] <- theta
    }
  }
  return(draws)
}

#' Process results of drawDPM to drop empty clusters
#'
#' @param draws Output of drawDPM
#'
#' @return List of lists with values for non-empty clusters
#' @export
#'
#' @examples
onlyActive <- function(draws){
  nDraws <- with(draws, ncol(K))
  active <- with(draws, lapply(1:nDraws, function(col) unique(K[,col])))
  pActive <- with(draws, lapply(1:nDraws,
                              function(i) p[active[[i]],i] / sum(p[active[[i]],i])))
  muActive <- with(draws, lapply(1:nDraws, function(i) mu[active[[i]],i]))
  sigma_sqActive <- with(draws, lapply(1:nDraws,
                                       function(i) sigma_sq[active[[i]],i]))
  out <- list(mu = muActive, sigma_sq = sigma_sqActive, K = draws$K,
              p = pActive, alpha = draws$alpha, theta = draws$theta)
  return(out)
}

#' Generate reduced form draws for particular cell
#'
#' @param y_cell Observations for one of the four cells.
#'
#' @return Reduced form draws
#' @export
#'
#' @examples
cellOutput <- function(y_cell){
  s_cell <- sd(y_cell)
  draws <- drawDPM(y_cell / s_cell, A = 10000, nu1 = 2, nu2 = 2)
  draws$mu <- draws$mu * s_cell
  draws$sigma_sq <- draws$sigma_sq * s_cell^2
  draws<- onlyActive(draws)
  means <- mapply(function(p, m) sum(p * m), p = draws$p, m = draws$mu)
  vars <-  mapply(function(p, m, s_sq) sum(p * s_sq) + sum(p * m^2) - sum(p * m)^2,
                       p = draws$p, m = draws$mu, s_sq = draws$sigma_sq)
  out <- list(mu = draws$mu, p = draws$p, sigma_sq = draws$sigma_sq, cellMean = means,
              cellVar = vars)
  return(out)
}

alphaBoundsIndep <- function(p0, p1, draw00, draw01, draw10, draw11){

}


drawObsIJ <- function(y_name, T_name, z_name, data, nDraws = 1000){

  # Do tons of stuff...

  # Pack things up and return
  obsDraws <- data.frame(n = rep(n, nDraws),
                         p = rep(p, nDraws),
                         q = rep(q, nDraws),
                         p0 = rep(p0, nDraws),
                         p1 = rep(p1, nDraws),
                         a0upper = a0upper, #VECTOR!
                         a1upper = a1upper, #VECTOR!
                         p00 = rep(p00, nDraws),
                         p01 = rep(p01, nDraws),
                         p10 = rep(p10, nDraws),
                         p11 = rep(p11, nDraws),
                         yb00 = yb00,
                         yb01 = yb01,
                         yb10 = yb10,
                         yb11 = yb11,
                         yt00 = yt00,
                         yt01 = yt01,
                         yt10 = yt10,
                         yt11 = yt11,
                         s2_00 = s2_00, #VECTOR!
                         s2_01 = s2_01, #VECTOR!
                         s2_10 = s2_10, #VECTOR!
                         s2_11 = s2_11, #VECTOR!
                         s_zT_upper = rep(s_zT_upper, nDraws),
                         N1 = N1,
                         N2 = rep(N2, nDraws),
                         N3 = N3,
                         N4 = rep(N4, nDraws),
                         beta_iv = beta_iv)
  return(obsDraws)
}

