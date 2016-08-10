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

#' Generate reduced form draws for particular (Tobs, Z) cell
#'
#' @param y_cell Observations for one of the four cells.
#' @param nDraws Number of MCMC draws
#' @param burnIn Number of burn-in MCMC draws
#'
#' @return Reduced form draws after removing "inactive" clusters
#' @export
#'
#' @examples
cellOutput <- function(y_cell, nDraws = 5000, burnIn = 1000){
  s_cell <- sd(y_cell)
  draws <- drawDPM(y_cell / s_cell, nDraws = nDraws, burnIn = burnIn)
  draws$mu <- draws$mu * s_cell
  draws$sigma_sq <- draws$sigma_sq * s_cell^2
  draws<- onlyActive(draws)
  means <- mapply(function(p, m) sum(p * m), p = draws$p, m = draws$mu)
  vars <-  mapply(function(p, m, s_sq) sum(p * s_sq) + sum(p * m^2) - sum(p * m)^2,
                       p = draws$p, m = draws$mu, s_sq = draws$sigma_sq)
  out <- list(mu = draws$mu, p = draws$p, sigma_sq = draws$sigma_sq,
              cellMean = means, cellVar = vars)
  return(out)
}

#' Calculate upper bounds for a0 and a1 based on independence assumption
#'
#' @param p0 P(Tobs=1|Z=0)
#' @param p1 P(Tobs=1|Z=1)
#' @param draw00 List with mu, p, sigma_sq from DPM for (Tobs=0,Z=0) cell
#' @param draw01 List with mu, p, sigma_sq from DPM for (Tobs=0,Z=1) cell
#' @param draw10 List with mu, p, sigma_sq from DPM for (Tobs=1,Z=0) cell
#' @param draw11 List with mu, p, sigma_sq from DPM for (Tobs=1,Z=1) cell
#'
#' @return List with upper bounds for a0 and a1
#' @export
#'
#' @examples
alphaBoundsIndep <- function(p0, p1, draw00, draw01, draw10, draw11){

  F00 <- function(x){
    sum(draw00$p * pnorm(x, mean = draw00$mu, sd = sqrt(draw00$sigma_sq)))
  }

  F10 <- function(x){
    sum(draw10$p * pnorm(x, mean = draw10$mu, sd = sqrt(draw10$sigma_sq)))
  }

  F01 <- function(x){
    sum(draw01$p * pnorm(x, mean = draw01$mu, sd = sqrt(draw01$sigma_sq)))
  }

  F11 <- function(x){
    sum(draw11$p * pnorm(x, mean = draw11$mu, sd = sqrt(draw11$sigma_sq)))
  }

  F0 <- function(x){
    (1 - p0) * F00(x) + p0 * F10(x)
  }

  F1 <- function(x){
    (1 - p1) * F01(x) + p1 * F11(x)
  }

  B1_a0_k0 <- function(x){
    p0 * F10(x) / F0(x)
  }

  B2_a0_k0 <- function(x){
    p0 * (1 - F10(x)) / (1 - F0(x))
  }

  B1_a0_k1 <- function(x){
    p1 * F11(x) / F1(x)
  }

  B2_a0_k1 <- function(x){
    p1 * (1 - F11(x)) / (1 - F1(x))
  }

  # Only look at regions where there is data
  l_k0 <- max(min(y00), min(y10))
  u_k0 <- min(max(y00), max(y10))
  l_k1 <- max(min(y01), min(y11))
  u_k1 <- min(max(y01), max(y11))

  # Upper bound for a0
  a0upper <- min(optimize(B1_a0_k0, lower = l_k0, upper = u_k0)$obj,
                 optimize(B2_a0_k0, lower = l_k0, upper = u_k0)$obj,
                 optimize(B1_a0_k1, lower = l_k1, upper = u_k1)$obj,
                 optimize(B2_a0_k1, lower = l_k1, upper = u_k1)$obj)

  # Upper bounds for  a1
  a1upper <- min(optimize(function(x) 1 - B1_a0_k0(x), lower = l_k0, upper = u_k0)$obj,
                 optimize(function(x) 1 - B2_a0_k0(x), lower = l_k0, upper = u_k0)$obj,
                 optimize(function(x) 1 - B1_a0_k1(x), lower = l_k1, upper = u_k1)$obj,
                 optimize(function(x) 1 - B2_a0_k1(x), lower = l_k1, upper = u_k1)$obj)

  upperBounds <- list(a0 = a0upper, a1 = a1upper)
  return(upperBounds)
}


#' Draw reduced form parameters from Ishwaran & James DPM
#'
#' @param y_name Name of outcome
#' @param T_name Name of treatment
#' @param z_name Name of instrument
#' @param data Dataframe
#' @param nDraws Number of MCMC draws
#' @param burnIn Number of burn-in draws
#'
#' @return Matrix of reduced form draws
#' @export
#'
#' @examples
drawObsIJ <- function(y_name, T_name, z_name, data,
                      nDraws = 1000, burnIn = 1000){

  # Extract data for named columns
  y <- get(y_name, data)
  Tobs <- get(T_name, data)
  z <- get(z_name, data)

  # Outcomes in the four (Tobs, Z) cells
  y00 <- y[Tobs == 0 & z == 0]
  y01 <- y[Tobs == 0 & z == 1]
  y10 <- y[Tobs == 1 & z == 0]
  y11 <- y[Tobs == 1 & z == 1]

  # Quantities on which we condition (no sampling uncertainty!)
  n <- length(Tobs)
  p <- mean(Tobs)
  q <- mean(z)
  p0 <- mean(Tobs[z == 0])
  p1 <- mean(Tobs[z == 1])
  p00 <- mean(Tobs == 0 & z == 0)
  p01 <- mean(Tobs == 0 & z == 1)
  p10 <- mean(Tobs == 1 & z == 0)
  p11 <- mean(Tobs == 1 & z == 1)

  # Fit DPM to each (Tobs, Z) cell
  draws00 <- cellOutput(y00, nDraws = nDraws, burnIn = burnIn)
  draws01 <- cellOutput(y01, nDraws = nDraws, burnIn = burnIn)
  draws10 <- cellOutput(y10, nDraws = nDraws, burnIn = burnIn)
  draws11 <- cellOutput(y11, nDraws = nDraws, burnIn = burnIn)

  # Compute upper bounds for a0, a1
  alpha0upper <- rep(NA_real_, nDraws)
  alpha1upper <- rep(NA_real_, nDraws)
  for(i in 1:nDraws){
    temp00 <- lapply(draws00, function(x) x[[i]])
    temp01 <- lapply(draws01, function(x) x[[i]])
    temp10 <- lapply(draws10, function(x) x[[i]])
    temp11 <- lapply(draws11, function(x) x[[i]])
    tempBounds <- alphaBoundsIndep(p0, p1, temp00, temp01, temp10, temp11)
    alpha0upper[i] <- tempBounds$a0
    alpha1upper[i] <- tempBounds$a1
  }

  # Unpack Cell means and variances from DPM draws
  yb00 <- draws00$cellMean
  yb01 <- draws01$cellMean
  yb10 <- draws10$cellMean
  yb11 <- draws11$cellMean

  s2_00 <- draws00$cellVar
  s2_01 <- draws01$cellVar
  s2_10 <- draws10$cellVar
  s2_11 <- draws11$cellVar

  # Implied y-tildes
  yt00 <- (1 - p0) * yb00
  yt01 <- (1 - p1) * yb01
  yt10 <- p0 * yb10
  yt11 <- p1 * yb11

  # Calculate Wald estimates at posterior draws from DPM
  beta_iv <- ((yt01 + yt11) - (yt00 + yt10)) / (p1 - p0)

  # No controls for the time being!
  s_zT_upper <- 1 / cov(z, Tobs) # Used to compute beta, Doesn't vary across draws
  C2 <- C4 <- 0 # Don't vary across draws even if there are covariates
  C1 <- C3 <- rep(0, nDraws) # Vary across draws when there are covariates


  # Pack things up and return
  obsDraws <- data.frame(n = rep(n, nDraws),
                         p = rep(p, nDraws),
                         q = rep(q, nDraws),
                         p0 = rep(p0, nDraws),
                         p1 = rep(p1, nDraws),
                         a0upper = alpha0upper, #VECTOR!
                         a1upper = alpha1upper, #VECTOR!
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
                         C1 = C1,
                         C2 = rep(C2, nDraws),
                         C3 = C3,
                         C4 = rep(C4, nDraws),
                         beta_iv = beta_iv)
  return(obsDraws)
}

