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
  ynames <- c("y00", "y01", "y10", "y11")
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
