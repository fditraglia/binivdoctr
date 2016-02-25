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


get_dTstar_tilde <- function(dz_tilde, a0, a1, obs){

  h1 <- with(obs, (1 - q) * yt00 + q * yt01)
  h2 <- with(obs, (1 - q) * yt10 + q * yt11)
  h <- with(obs, (h1 - a1 * (h1 + h2)) / (1 - p - a1))

  g <- with(obs, (yt01 - yt00) - a1 * ((yt01 - yt00) + (yt11 - yt10)))
  D <- with(obs, ((1 - a0) * yt10 - a0 * yt00)/(p0 - a0) -
              ((1 - a0) * yt11 - a0 * yt01)/(p1 - a0))

  K_first <- g / (p0 - p1)
  K_second <- ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (p0 - p1))
  K <- K_first - K_second - h

  F1 <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
  F2 <- with(obs, q * (1 - q) * F1)
  S <- with(obs, (1 - a0 - a1) / (p0 - p1))

  return(with(obs, K - F1 - S * N1 + (F2 + S * (N2 - 1)) * dz_tilde))
}

get_dz_tilde <- function(dTstar_tilde, a0, a1, obs){

  h1 <- with(obs, (1 - q) * yt00 + q * yt01)
  h2 <- with(obs, (1 - q) * yt10 + q * yt11)
  h <- with(obs, (h1 - a1 * (h1 + h2)) / (1 - p - a1))

  g <- with(obs, (yt01 - yt00) - a1 * ((yt01 - yt00) + (yt11 - yt10)))
  D <- with(obs, ((1 - a0) * yt10 - a0 * yt00)/(p0 - a0) -
              ((1 - a0) * yt11 - a0 * yt01)/(p1 - a0))

  K_first <- g / (p0 - p1)
  K_second <- ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (p0 - p1))
  K <- K_first - K_second - h

  F1 <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
  F2 <- with(obs, q * (1 - q) * F1)
  S <- with(obs, (1 - a0 - a1) / (p0 - p1))

  denom <- (F2 + S * (N2 - 1))
  return(with(obs, (F1 + S * N1 - K) / (denom) + (1 / (denom)) * dTstar_tilde))
}




get_m11_bounds <- function(a0, a1, obs){
  D <- with(obs, ((1 - a0) * yt10 - a0 * yt00)/(p0 - a0) -
              ((1 - a0) * yt11 - a0 * yt01)/(p1 - a0))
  m00 <- with(obs, ((1 - a1) * yt00 - a1 * yt10) / (1 - p0 - a1))
  m01 <- with(obs, ((1 - a1) * yt01 - a1 * yt11) / (1 - p1 - a1))
  ybStar10 <- with(obs, ((1 - a0) * yt10 - a0 * yt00) / (p0 - a0))
  ybStar11 <- with(obs, ((1 - a0) * yt11 - a0 * yt01) / (p1 - a0))
  V00 <- with(obs, a1 * (1 - a0) * (p0 - a0) * (1 - p0 - a1) /
    ((1 - p0)^2 * (1 - a0 - a1)^2))
  V10 <- with(obs, a1 * (1 - a0) * (p1 - a0) * (1 - p1 - a1) /
    ((1 - p1)^2 * (1 - a0 - a1)^2))
  V01 <- with(obs, a0 * (1 - a1) * (p0 - a0) * (1 - p0 - a1) /
    (p0^2 * (1 - a0 - a1)^2))
  V11 <- with(obs, a0 * (1 - a1) * (p1 - a0) * (1 - p1 - a1) /
    (p1^2 * (1 - a0 - a1)^2))

  # Restrictions for m10
  c00 <- -1 * (ybStar10(a0) * (ybStar10(a0) - 2 * m00(a1)) -
    s2_00 / V00(a0, a1))
  c10 <- -1 * (ybStar10(a0) * (ybStar10(a0) - 2 * m00(a1)) -
    s2_10 / V10(a0, a1))
  # Since these are the same quadratic up to a vertical shift only need
  # to look at the one with a *smaller* intercept
  c0 <- min(c00, c10)
  b0 <- -2 * m00(a1)
  disc0 <- b0^2 - 4 * 1 * c0
  if(disc0 > 0){
    roots_m10 <- sort(0.5 * (-b0 + c(-1, 1) * sqrt(disc0)))
  }else{
    roots_m10 <- c(NA, NA)
  }
  # Restrictions for m11
  c01 <- -1 * (ybStar11(a0) * (ybStar11(a0) - 2 * m01(a1)) -
    s2_01 / V01(a0, a1))
  c11 <- -1 * (ybStar11(a0) * (ybStar11(a0) - 2 * m01(a1)) -
    s2_11 / V11(a0, a1))
  # Since these are the same quadratic up to a vertical shift only need
  # to look at the one with a *smaller* intercept
  c1 <- min(c01, c11)
  b1 <- -2 * m00(a1)
  disc1 <- b1^2 - 4 * 1 * c1
  if(disc1 > 0){
    roots_m11 <- sort(0.5 * (-b1 + c(-1, 1) * sqrt(disc1)))
  }else{
    roots_m11 <- c(NA, NA)
  }
  roots <- list(m10 = roots_m10, m11 = roots_m11)


  if(all(is.na(unlist(roots)))){
    return(list(NA))
  }else if(all(is.na(roots$m10)) & all(!is.na(roots$m11))){
    return(list(roots$m11))
  }else if(all(!is.na(roots$m10)) & all(is.na(roots$m11))){
    return(list(roots$m10 - D))
  }else{
    set1 <- roots$m10 - D
    set2 <- roots$m11
    if(min(set1) <= min(set2)){
      L <- set1
      O <- set2
    }else{
      O <- set1
      L <- set2
    }
    if(min(L) <= min(O) & max(L) >= max(O)){
     return(list(L))
    }else if(min(L) <= min(O) & max(L) <= max(O) & min(O) <= max(L)){
      return(list(c(min(L), max(O))))
    }else{
      return(list(L, O))
    }
  }
}




