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
  obs$n <- length(Tobs)
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
    obs$s_zT_upper <- Sigma_inv[1,1] # need this to back out implied beta
    s_xT_upper <- matrix(Sigma_inv[-1,1], ncol(x), 1)
    obs$C1 <- drop(cov(z, x) %*% gamma_iv / var(z))
    obs$C2 <- drop(cov(z, x) %*% s_xT_upper)
    obs$C3 <- drop(cov(Tobs, x) %*% gamma_iv)
    obs$C4 <- drop(var(z) * cov(Tobs, x) %*% s_xT_upper)
  }else{
    # Case without covariates!
    obs$s_zT_upper <- 1 / (cov(z, Tobs))
    obs$C1 <- obs$C2 <- obs$C3 <- obs$C4 <- 0
  }
  return(obs)
}


#' Calculate endogeneity of binary regressor
#'
#' @param dz_tilde Endogeneity of binary instrument.
#' @param a0 Probability of observing T = 1 when true treatment status is 0.
#' @param a1 Probability of observing T = 0 when true treatment status is 0.
#' @param obs List of obvservables calculated using get Obs.
#'
#' @return Endogneity of binary regressor implied by specified values of
#' instrument invalidity and mis-classification probabilities.
#' @export
#'
get_dTstar_tilde <- function(dz_tilde, a0, a1, obs){

  h1 <- with(obs, (1 - q) * yt00 + q * yt01)
  h2 <- with(obs, (1 - q) * yt10 + q * yt11)
  h <- with(obs, (h1 - a1 * (h1 + h2)) / (1 - p - a1))

  g <- with(obs, (yt01 - yt00) - a1 * ((yt01 - yt00) + (yt11 - yt10)))
  D <- with(obs, ((1 - a0) * yt10 - a0 * yt00)/(p0 - a0) -
              ((1 - a0) * yt11 - a0 * yt01)/(p1 - a0))

  Ff <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
  S <- with(obs, (p1 - p0) / (1 - a0 - a1))
  B <- with(obs, (g - (p0 - p1) * h) / (1 - a0 - a1) -
              ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (1 - a0 - a1)))
  S_tilde <- with(obs, S / (S * Ff * C4 - C2 + 1))
  B_tilde <- with(obs, (S * Ff * C3 + B - C1) / (S * Ff * C4 - C2 + 1))

  return(with(obs, -B_tilde / S_tilde + (1 / S_tilde) * dz_tilde))
}

#' Calculate endogeneity of binary instrument
#'
#' @param dTstar_tilde Endogeneity of binary regressor.
#' @param a0 Probability of observing T = 1 when true treatment status is 0.
#' @param a1 Probability of observing T = 0 when true treatment status is 1.
#' @param obs List of observables calculated by \code{getObs}.
#'
#' @return Value of instrument endogeneity implied by the data and specified values of regressor endogeneity and mis-classification probabilities.
#' @export
#'
get_dz_tilde <- function(dTstar_tilde, a0, a1, obs){

  h1 <- with(obs, (1 - q) * yt00 + q * yt01)
  h2 <- with(obs, (1 - q) * yt10 + q * yt11)
  h <- with(obs, (h1 - a1 * (h1 + h2)) / (1 - p - a1))

  g <- with(obs, (yt01 - yt00) - a1 * ((yt01 - yt00) + (yt11 - yt10)))
  D <- with(obs, ((1 - a0) * yt10 - a0 * yt00)/(p0 - a0) -
              ((1 - a0) * yt11 - a0 * yt01)/(p1 - a0))

  Ff <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
  S <- with(obs, (p1 - p0) / (1 - a0 - a1))
  B <- with(obs, (g - (p0 - p1) * h) / (1 - a0 - a1) -
              ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (1 - a0 - a1)))
  S_tilde <- with(obs, S / (S * Ff * C4 - C2 + 1))
  B_tilde <- with(obs, (S * Ff * C3 + B - C1) / (S * Ff * C4 - C2 + 1))

  return(with(obs, B_tilde + S_tilde * dTstar_tilde))
}


#' Calculate sharp bounds for a0 and a1
#'
#' @param obs List of observables calculated by \code{getObs}.
#'
#' @return Vector with two elements: the first is the upper bound for a0 and the
#' second is the upper bound for a1. (The lower bounds for each are both zero.)
#' @export
#' @examples
#' afghanControls <- c("headchild", "age",  "yrsvill",  "farsi",  "tajik",
#'                    "farmers",  "agehead",  "educhead",  "nhh",  "land",
#'                    "sheep", "distschool", "chagcharan")
#' afghanY <- "testscore"
#' afghanT <- "enrolled"
#' afghanZ <- "buildschool"
#' afghanObs <- getObs(afghanY, afghanT, afghanZ, afghanControls, afghan)
#' get_alpha_bounds(afghanObs)
get_alpha_bounds <- function(obs){

  A0 <- with(obs, p0 * (1 - p0) * (yb10 - yb00)^2 + (1 - p0) * s2_00 + p0 * s2_10)
  A1 <- with(obs, p1 * (1 - p1) * (yb11 - yb01)^2 + (1 - p1) * s2_01 + p1 * s2_11)
  Q0 <- with(obs, p0^2 * s2_10 - (1 - p0)^2 * s2_00)
  Q1 <- with(obs, p1^2 * s2_11 - (1 - p1)^2 * s2_01)

  quad_a0_0 <- with(obs, c(p0^2 * s2_10, -Q0 - A0, A0))
  quad_a0_1 <- with(obs, c(p1^2 * s2_11, -Q1 - A1, A1))
  roots_a0_0 <- Re(polyroot(quad_a0_0))
  roots_a0_1 <- Re(polyroot(quad_a0_1))
  a0_upper <- min(c(roots_a0_0, roots_a0_1))

  quad_a1_0 <- with(obs, c((1 - p0)^2 * s2_00, Q0 - A0, A0))
  quad_a1_1 <- with(obs, c((1 - p1)^2 * s2_01, Q1 - A1, A1))
  roots_a1_0 <- Re(polyroot(quad_a1_0))
  roots_a1_1 <- Re(polyroot(quad_a1_1))
  a1_upper <- min(c(roots_a1_0, roots_a1_1))

  return(c(a0upper = a0_upper, a1upper = a1_upper))

}
