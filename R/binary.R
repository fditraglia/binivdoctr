#' Calculate summary statistics of observables
#'
#' @description This function calculates summary statistics for later use in the
#' partial identification exercise for a binary, mis-measured, endogenous
#' regressor based on a binary instrument that may itsef be endogenous.
#' @param y_name Name of the outcome variable.
#' @param T_name Name of the potentially endogenous and mis-measured binary
#' treatment variable.
#' @param z_name Name of the potentially endogenous binary instrument.
#' @param controls Vector of names of exogenous control regressors.
#' @param data Name of dataframe whose columns contain the data for the outcome,
#' regressor, instrument, and controls.
#'
#' @return A list of summary statistics: probabilities, conditional means and
#' variances, etc., in notation corresponding to the paper.
#' @export
#' @examples
#' AngristControls <- c("age", "svy", "sex2", "phone", "hsvisit",
#'                      "d1995", "djamundi", "dmonth1", "dmonth2",
#'                      "dmonth3", "dmonth4", "dmonth5", "dmonth6", "dmonth7",
#'                      "dmonth8", "dmonth9", "dmonth10", "dmonth11",
#'                      "strata1", "strata2", "strata3", "strata4", "strata5")
#' AngristY <- "totalRepeats"
#' AngristT <- "usesch"
#' AngristZ <- "vouch0"
#' getObs(AngristY, AngristT, AngristZ, AngristControls, angrist)
#'
getObs <- function(y_name, T_name, z_name, controls = NULL, data){
  # Extract data for named columns
  y <- get(y_name, data)
  Tobs <- get(T_name, data)
  z <- get(z_name, data)

  obs <- list()
  obs$p <- mean(Tobs)
  obs$q <- mean(z)
  obs$p0 <- mean(Tobs[z == 0])
  obs$p1 <- mean(Tobs[z == 1])
  obs$p00 <- mean(Tobs == 0 & z == 0)
  obs$p01 <- mean(Tobs == 0 & z == 1)
  obs$p10 <- mean(Tobs == 1 & z == 0)
  obs$p11 <- mean(Tobs == 1 & z == 1)
  obs$yb00 <- mean(y[(Tobs == 0) & (z == 0)])
  obs$yb01 <- mean(y[(Tobs == 0) & (z == 1)])
  obs$yb10 <- mean(y[(Tobs == 1) & (z == 0)])
  obs$yb11 <- mean(y[(Tobs == 1) & (z == 1)])
  obs$yt00 <- with(obs, (1 - p0) * yb00)
  obs$yt01 <- with(obs, (1 - p1) * yb01)
  obs$yt10 <- with(obs, p0 * yb10)
  obs$yt11 <- with(obs, p1 * yb11)
  obs$s2_00 <- var(y[Tobs == 0 & z == 0])
  obs$s2_01 <- var(y[Tobs == 0 & z == 1])
  obs$s2_10 <- var(y[Tobs == 1 & z == 0])
  obs$s2_11 <- var(y[Tobs == 1 & z == 1])

  if(!is.null(controls)){
    second_stage <- reformulate(c(T_name, controls), response = y_name)
    first_stage <- reformulate(c(z_name, controls))
    gamma_iv <- coefficients(AER::ivreg(second_stage, first_stage, data))[-c(1,2)]
    # No intercept since we "project it out" by working with Cov matrix below
    x <- model.matrix(reformulate(controls, intercept = FALSE), data)
    Sigma <- rbind(cbind(cov(z, Tobs), cov(z, x)),
                       cbind(cov(x, Tobs), cov(x)))
    Sigma_inv <- solve(Sigma)
    obs$s_zT_upper <- Sigma_inv[1,1]
    s_xT_upper <- matrix(Sigma_inv[-1,1], ncol(x), 1)
    obs$N1 <- drop(cov(z, x) %*% gamma_iv / var(z))
    obs$N2 <- drop(cov(z, x) %*% s_xT_upper)
    obs$N3 <- drop(cov(Tobs, x) %*% gamma_iv)
    obs$N4 <- drop(cov(Tobs, x) %*% s_xT_upper)
  }else{
    # Case without covariates!
    obs$N1 <- obs$N2 <- obs$N3 <- obs$N4 <- 0
  }
  return(obs)
}


