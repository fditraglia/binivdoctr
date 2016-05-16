# Gibbs Sampler for Ishwaran and James

# Some fake data
set.seed(2837)
n <- 500
x1 <- -1 + 0.5 * rnorm(n)
x2 <- 1 + 1.5 * rnorm(n)
y <- sample(0:1, n, replace = TRUE, prob = c(0.65, 0.35))
data <- y * x1 + (1 - y) * x2


# Prior hyperparameters
A <- 1000
xi_sq <- (4 * sd(data))^2
nu1 <- nu2 <- 2
eta1 <- eta2 <- 2

# Maximum number of clusters
M <- 50


# initialize randomly with only 10 nonempty clusters
set.seed(2617)
K <- factor(sample(1:10, n, replace = TRUE), 1:M) #Keep track of empties!
clusters <- tapply(data, K, c)
empty <- sapply(clusters, is.null)
clusters[!empty]
mu <- tapply(data, K, mean)
sigma_sq <- tapply(data, K, var)
n_k <- tapply(data, K, length)
p <- (1 / n) * n_k


# Another possible way?
#clusters <- lapply(1:M, function(x) data[which(K == x)])
#mu <- sapply(clusters, mean)
#sigma_sq <- sapply(clusters, var)
#p <- (1 / n) * sapply(clusters, length)

# update mu
s_sq <- 1 / (n_k / sigma_sq + 1 / xi_sq)


# update sigma_sq

# update cluster membership

# update p

# update alpha

# Update theta
