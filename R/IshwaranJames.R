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

# initialize randomly
set.seed(2617)
K <- sample(1:M, n, replace = TRUE)

# update mu

# update sigma_sq

# update cluster membership

# update p

# update alpha

# Update theta
