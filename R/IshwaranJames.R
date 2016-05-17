# Gibbs Sampler for Ishwaran and James

# Some fake data
set.seed(2837)
nObs <- 100
x1 <- -1 + rnorm(nObs)
x2 <- 1 + rnorm(nObs)
y <- sample(0:1, nObs, replace = TRUE, prob = c(0.65, 0.35))
X <- y * x1 + (1 - y) * x2


# Prior hyperparameters
A <- 1000
xi_sq <- (4 * sd(X))^2
nu1 <- nu2 <- 2
eta1 <- eta2 <- 2

# Maximum number of clusters
maxClust <- 50


# initialize randomly
set.seed(2617)
K <- factor(sample(1:50, nObs, replace = TRUE), 1:maxClust) #Keep track of empties!
clusters <- tapply(X, K, c)
empty <- sapply(clusters, is.null)
mu <- tapply(X, K, mean)
sigma_sq <- tapply(X, K, var)
# Handle clusters with only one observation
sigma_sq[!empty & is.na(sigma_sq)] <- var(X)
obsPerClust <- tapply(X, K, length)
obsPerClust[is.na(obsPerClust)] <- 0
p <- (1 / nObs) * obsPerClust
theta <- mean(X)
alpha <- rgamma(1, shape = eta1, rate = eta2)

# object to store the draws
nDraws <- 5000
draws <- list(mu = matrix(NA_real_, nrow = maxClust, ncol = nDraws), 
              sigma_sq = matrix(NA_real_, nrow = maxClust, ncol = nDraws), 
              K = matrix(NA_integer_, nrow = nObs, ncol = nDraws), 
              p = matrix(NA_real_, nrow = maxClust, ncol = nDraws), 
              alpha = rep(NA_real_, nDraws), 
              theta = rep(NA_real_, nDraws))


# The Gibbs Sampler!
for(i in 1:nDraws){
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
  V <- rbeta(maxClust - 1, shape1 = 1 + obsPerClust[-maxClust], 
             shape2 = alpha + cumsum(obsPerClust)[-1])
  log_remain <- cumsum(log(1 - V))
  p <- rep(NA_real_, maxClust)
  p[-maxClust] <- c(V[1], exp(log(V)[-1] + log_remain[-(maxClust - 1)]))
  p[maxClust] <- 1 - sum(p[-maxClust])
  
  # update alpha
  alpha <- rgamma(1, shape = maxClust + eta1 - 1, 
                  rate = eta2 - log_remain[(maxClust - 1)])
  
  # update theta
  s_theta_sq <- 1 / (maxClust / xi_sq + 1 / A)
  m_theta <- (s_theta_sq / xi_sq) * sum(mu)
  theta <- rnorm(1, mean = m_theta, sd = sqrt(s_theta_sq))
  
  # store the draws
  draws$mu[,i] <- mu
  draws$sigma_sq[,i] <- sigma_sq
  draws$K[,i] <- K
  draws$p[,i] <- p
  draws$alpha[i] <- alpha
  draws$theta[i] <- theta
  
}


plot(draws$alpha, type = 'l')
plot(draws$theta, type = 'l')
plot(apply(draws$K, 2, function(x) length(unique(x))), type = 'l')

table(draws$K[,5000])
round(draws$p[,5000], 4)

plot_pdf <- function(n){
  pdf <- function(x, n){
    sum(draws$p[,n] * dnorm(x, mean = draws$mu[,n], sd = sqrt(draws$sigma_sq[,n])))
  }
  
  phi <- Vectorize(function(x) pdf(x, n), "x")
  curve(phi(x), from = -5, to = 5)
}

plot_pdf(sample(1000:5000, 1))

plot_cdf <- function(n){
CDF <- function(x, n){
  sum(draws$p[,n] * pnorm(x, mean = draws$mu[,n], sd = sqrt(draws$sigma_sq[,n])))
}

Phi <- Vectorize(function(x) CDF(x, n), "x")
curve(Phi(x), from = -5, to = 5, n = 1000)
  
}

plot_cdf(sample(1000:5000, 1))

# evaluate at posterior mean
mu_pm <- apply(draws$mu[,-c(1:1000)], 1, mean)
sigma_pm <- sqrt(apply(draws$sigma_sq[,-c(1:1000)], 1, mean))
p_pm <- apply(draws$p[,-c(1:1000)], 1, mean)

f_pm <- function(x){
  sum(p_pm * dnorm(x, mean = mu_pm, sd = sigma_pm))
}

F_pm<- function(x){
  sum(p_pm * pnorm(x, mean = mu_pm, sd = sigma_pm))
}

phi_pm <- Vectorize(f_pm, "x")
Phi_pm<- Vectorize(F_pm, "x")

curve(Phi_pm(x), from = -5, to = 5, n = 1000)
curve(phi_pm(x), from = -5, to = 5, n = 1000)









