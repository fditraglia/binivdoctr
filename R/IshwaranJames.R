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

