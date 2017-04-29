setwd('C:/Users/mmm/Documents/GITHUB_FOLDERS/binivdoctr/R')

source('binCDFs.R')
source('dirichlet.R')
source('identifiedSet.R')
source('IshwaranJames.R')
source('plotFunctions.R')
source('SamplingUncertainty.R')
source('simDGPs.R')
source('summary.R')

library('nleqslv')
library('rootSolve')

# context("Example Datasets")

# test_that("Test an interior min/max for dz_tilde ", {
# 
#   
#   
# })

L <- 100 # Size of grid

obs <- list()

obs$p <- 0.3878092
obs$q <- 0.4749366
obs$p0 <- 0.2438337
obs$p1 <- 0.5469806
obs$p00 <- 0.3970353
obs$p01 <- 0.2151555
obs$p10 <- 0.1280281
obs$p11 <- 0.2597811
obs$yb00 <- 0.9371405
obs$yb01 <- -1.063885
obs$yb10 <- 1.010211
obs$yb11 <- -0.3593873
obs$yt00 <- 0.7086341
obs$yt01 <- -0.4819606
obs$yt10 <- 0.2463235
obs$yt11 <- -0.1965779
obs$s2_00 <- 0.6952148
obs$s2_01 <- 1.233553
obs$s2_10 <- 0.1491803
obs$s2_11 <- 0.06342888
obs$C1 <- 0
obs$C2 <- 0
obs$C3 <- 0
obs$C4 <- 0
dxs <-  -10 # -2.561388

a0upper <- get_alpha_bounds(obs)[1]
a1upper <- get_alpha_bounds(obs)[2]

a0_grid <- seq(0,a0upper,length.out = L)
a1_grid <- seq(0,a1upper,length.out = L)


dzs_matrix <- outer(a0_grid, a1_grid, function(a0, a1)
                                      get_dz_tilde(dxs,a0,a1,obs))

coordinates_argmin <- which(dzs_matrix==min(dzs_matrix),arr.ind=TRUE)
coordinates_argmax <- which(dzs_matrix==max(dzs_matrix),arr.ind=TRUE)


initial_guess <- expand.grid(seq(0,a0upper,length.out = 5),
                             seq(0,a1upper,length.out = 5))

result <- searchZeros(as.matrix(initial_guess),gradient_fn)
result$x

gradient_fn(coordinates_argmin)
gradient_fn(coordinates_argmax)

# result$x



obs <- list()
# obs$n <- length(Tobs)
obs$p <- p_vec[1:200]
obs$q <- q_vec[1:200]
obs$p0 <- p_0_vec[1:200]
obs$p1 <- p_1_vec[1:200]
obs$p00 <- p_00_vec[1:200]
obs$p01 <- p_01_vec[1:200]
obs$p10 <- p_10_vec[1:200]
obs$p11 <- p_11_vec[1:200]
obs$yb00 <- y_00_vec[1:200]
obs$yb01 <- y_01_vec[1:200]
obs$yb10 <- y_10_vec[1:200]
obs$yb11 <- y_11_vec[1:200]
obs$yt00 <- obs$yb00*(1-obs$p0)
obs$yt01 <- obs$yb01*(1-obs$p1)
obs$yt10 <- obs$yb10*(obs$p0)
obs$yt11 <- obs$yb11*(obs$p1)
obs$s2_00 <- s_00_vec[1:200]
obs$s2_01 <- s_01_vec[1:200]
obs$s2_10 <- s_10_vec[1:200]
obs$s2_11 <- s_11_vec[1:200]

obs$C1 <- rep(0,200)
obs$C2 <- rep(0,200)
obs$C3 <- rep(0,200)
obs$C4 <- rep(0,200)

dxs <- dxs_vec[1:200]

##

# -1.7038968 -1.7038968

