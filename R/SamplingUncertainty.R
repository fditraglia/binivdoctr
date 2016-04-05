getMomentsForSampler <- function(y_name, T_name, z_name, controls = NULL, data){
  # Extract data for named columns
  y <- get(y_name, data)
  Tobs <- get(T_name, data)
  z <- get(z_name, data)
  n <- length(y)

  ybar <- c(mean(y[Tobs == 0 & z == 0]),
            mean(y[Tobs == 0 & z == 1]),
            mean(y[Tobs == 1 & z == 0]),
            mean(y[Tobs == 1 & z == 1]))
  ynames <- c("yb00", "yb01", "yb10", "yb11")
  names(ybar) <- ynames

  s2y <- c(var(y[Tobs == 0 & z == 0]),
           var(y[Tobs == 0 & z == 1]),
           var(y[Tobs == 1 & z == 0]),
           var(y[Tobs == 1 & z == 1]))

  ny <- c(sum(Tobs == 0 & z == 0),
          sum(Tobs == 0 & z == 1),
          sum(Tobs == 1 & z == 0),
          sum(Tobs == 1 & z == 1))

  Syy <- diag(n * s2y / ny)
  rownames(Syy) <- colnames(Syy) <- ynames
  out_var <- Syy / n
  out_mean <- ybar

  if(!is.null(controls)){

    xiRegFormula <- y ~ Tobs + z + Tobs:z
    omega <- lm(xiRegFormula)$residuals
    Q <- matrix(c(1, 0, 0, 0,
                  1, 0, 1, 0,
                  1, 1, 0, 0,
                  1, 1, 1, 1), 4, 4, byrow = TRUE)
    A <- model.matrix(xiRegFormula)
    Sigma_AA_inv <- solve(crossprod(A) / n)

    second_stage <- reformulate(c(T_name, controls), response = y_name)
    first_stage <- reformulate(c(z_name, controls))
    bReg <- AER::ivreg(second_stage, first_stage, data)
    gamma_iv <- coefficients(bReg)[-c(1,2)]
    rho <- residuals(bReg)
    R <- model.matrix(first_stage, data)
    W <- model.matrix(second_stage, data)
    Sigma_RW_inv <- solve(crossprod(R, W) / n)

    Xi_rr <- t(R) %*% diag(rho^2) %*%  R / (n - 1)
    Xi_rw <- t(R) %*% diag(rho * omega) %*% A / (n - 1)

    H_yy <- Syy
    H_bb <- Sigma_RW_inv %*% Xi_rr %*% t(Sigma_RW_inv)
    H_by <- Sigma_RW_inv %*% Xi_rw %*% t(Q %*% Sigma_AA_inv)
    H <- rbind(cbind(H_bb, H_by),
               cbind(t(H_by), H_yy))
    S <- H[-c(1,2), -c(1,2)]

    out_mean <- c(gamma_iv, ybar)
    out_var <- S / n
  }
  return(list(mean = out_mean, var = out_var))
}


drawObs <- function(y_name, T_name, z_name, controls = NULL, data,
                    nDraws = 1000){

  # Draw reduced form parameters
  moments <- getMomentsForSampler(y_name, T_name, z_name, controls, data)
  RFdraws <- MASS::mvrnorm(nDraws, moments$mean, moments$var)

  gamma_iv <- RFdraws[,-which(colnames(RFdraws) %in%
                                c("yb00", "yb01", "yb10", "yb11"))]
  yb00 <- RFdraws[,colnames(RFdraws) == "yb00"]
  yb01 <- RFdraws[,colnames(RFdraws) == "yb01"]
  yb10 <- RFdraws[,colnames(RFdraws) == "yb10"]
  yb11 <- RFdraws[,colnames(RFdraws) == "yb11"]

  # Extract data for named columns
  y <- get(y_name, data)
  Tobs <- get(T_name, data)
  z <- get(z_name, data)

  # Quantities for which we do not propagate sampling uncertainty
  n <- length(Tobs)
  p <- mean(Tobs)
  q <- mean(z)
  p0 <- mean(Tobs[z == 0])
  p1 <- mean(Tobs[z == 1])
  p00 <- mean(Tobs == 0 & z == 0)
  p01 <- mean(Tobs == 0 & z == 1)
  p10 <- mean(Tobs == 1 & z == 0)
  p11 <- mean(Tobs == 1 & z == 1)
  s2_00 <- var(y[Tobs == 0 & z == 0])
  s2_01 <- var(y[Tobs == 0 & z == 1])
  s2_10 <- var(y[Tobs == 1 & z == 0])
  s2_11 <- var(y[Tobs == 1 & z == 1])

  # These are now vectors since yb00 etc are vectors
  yt00 <- (1 - p0) * yb00
  yt01 <- (1 - p1) * yb01
  yt10 <- p0 * yb10
  yt11 <- p1 * yb11

  if(!is.null(controls)){
    # No intercept since we "project it out" by working with Cov matrix below
    x <- model.matrix(reformulate(controls, intercept = FALSE), data)
    Sigma <- rbind(cbind(cov(z, Tobs), cov(z, x)),
                       cbind(cov(x, Tobs), cov(x)))
    Sigma_inv <- solve(Sigma)
    s_zT_upper <- Sigma_inv[1,1] # We'll need this to back out the implied beta
    s_xT_upper <- matrix(Sigma_inv[-1,1], ncol(x), 1)

    # These don't vary across draws
    N2 <- drop(cov(z, x) %*% s_xT_upper)
    N4 <- drop(cov(Tobs, x) %*% s_xT_upper)

    # These will be vectors: they vary across draws
    N1 <- drop(cov(z, x) %*% t(gamma_iv) / var(z))
    N3 <- drop(cov(Tobs, x) %*% t(gamma_iv))

  }else{
    N2 <- N4 <- 0 # Don't vary across draws
    N1 <- N3 <- rep(0, nDraws) # Varies across draws
  }

  # Pack things up and return somehow...
  obsDraws <- data.frame(n = rep(n, nDraws),
                         p = rep(p, nDraws),
                         q = rep(q, nDraws),
                         p0 = rep(p0, nDraws),
                         p1 = rep(p1, nDraws),
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
                         s2_00 = rep(s2_00, nDraws),
                         s2_01 = rep(s2_01, nDraws),
                         s2_10 = rep(s2_10, nDraws),
                         s2_11 = rep(s2_11, nDraws),
                         s_zT_upper = rep(s_zT_upper, nDraws),
                         N1 = N1,
                         N2 = rep(N2, nDraws),
                         N3 = N3,
                         N4 = rep(N4, nDraws))
  return(obsDraws)
}


draw_dz_tilde <- function(y_name, T_name, z_name, controls = NULL, data,
                          dTstar_tilde_range,
                          drawAlphas = function() rdirichlet(1, c(1, 1, 5)),
                          nRF = 10, nIS = 1, maxIter = nIS * 100){

  RFdraws <- drawObs(y_name , T_name, z_name, controls, data , nDraws  = nRF)
  for(i in 1:nrow(RFdraws)){
    obs <- as.list(RFdraws[i,])

    drawCounter <- IScounter <- 0

    while(IScounter < nIS){

      dTstar_tilde <- runif(1, dTstar_tilde_range[1], dTstar_tilde_range[2])
      alphas <- drawAlphas()
      a0 <- alphas[1]
      a1 <- alphas[2]
      dz_tilde <- get_dz_tilde_check(dTstar_tilde, a0, a1, obs)

      drawCounter <- drawCounter + 1
      if(!is.na(dz_tilde))

    }

  }
}


















