# Gibbs Sampler for Ishwaran and James

# Some fake data
set.seed(2837)
n <- 100
x1 <- -1 + 0.5 * rnorm(n)
x2 <- 1 + 1.5 * rnorm(n)
y <- sample(0:1, n, replace = TRUE, prob = c(0.65, 0.35))
X <- y * x1 + (1 - y) * x2


# Prior hyperparameters
A <- 1000
xi_sq <- (4 * sd(X))^2
nu1 <- nu2 <- 2
eta1 <- eta2 <- 2

# Maximum number of clusters
M <- 50


# initialize randomly with only 10 nonempty clusters
set.seed(2617)
K <- factor(sample(c(1, 5, 7, 9, 13, 25, 36, 40), n, replace = TRUE), 1:M) #Keep track of empties!
clusters <- tapply(X, K, c)
empty <- sapply(clusters, is.null)
n_empty <- sum(empty)
mu <- tapply(X, K, mean)
sigma_sq <- tapply(X, K, var)
n_k <- tapply(X, K, length)
n_k[is.na(n_k)] <- 0
p <- (1 / n) * n_k


# Another possible way?
#clusters <- lapply(1:M, function(x) data[which(K == x)])
#mu <- sapply(clusters, mean)
#sigma_sq <- sapply(clusters, var)
#p <- (1 / n) * sapply(clusters, length)

# update mu
s_sq <- rep(NA_real_, M)
s_sq[!empty] <- 1 / (n_k[!empty] / sigma_sq[!empty] + 1 / xi_sq)
#m[!empty] <- 

m[empty] <- theta
s_sq[empty] <- xi_sq 
rnorm(M, mean = theta, sd = sqrt(xi_sq))



# update sigma_sq

# update cluster membership
pi <- matrix(rep(p / sqrt(sigma_sq), times = n) * 
             exp( -0.5 * (rep(X, each = M) - rep(mu, times = n))^2 / 
             rep(sigma_sq, times = n)), nrow = M, ncol = n, byrow = FALSE)
K <- factor(apply(pi, 2, function(x) sample(1:M, size = 1, prob = x)), 1:M)
n_k <- tapply(X, K, length)
n_k[is.na(n_k)] <- 0

# update p
V <- rbeta(M - 1, shape1 = 1 + n_k[-M], shape2 = cumsum(n_k)[-1])
log_remain <- cumsum(log(1 - V))
p <- c(V[1], exp(log(V)[-1] + log_remain[-(M - 1)]))

# update alpha
alpha <- rgamma(1, shape = M + eta1 - 1, rate = eta2 - log_remain[(M - 1)])

# Update theta
s_theta_sq <- 1 / (M / xi_sq + 1 / A)
m_theta <- (s_theta_sq / xi_sq) * sum(mu)
theta <- rnorm(1, mean = m_theta, sd = sqrt(s_theta_sq))






















