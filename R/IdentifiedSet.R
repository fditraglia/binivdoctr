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
    obs$N1 <- drop(cov(z, x) %*% gamma_iv / var(z))
    obs$N2 <- drop(cov(z, x) %*% s_xT_upper)
    obs$N3 <- drop(cov(Tobs, x) %*% gamma_iv)
    obs$N4 <- drop(cov(Tobs, x) %*% s_xT_upper)
  }else{
    # Case without covariates!
    obs$s_zT_upper <- 1 / (cov(z, Tobs))
    obs$N1 <- obs$N2 <- obs$N3 <- obs$N4 <- 0
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

  K_first <- with(obs, g / (p0 - p1))
  K_second <- with(obs, ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (p0 - p1)))
  K <- K_first - K_second - h

  F1 <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
  F2 <- with(obs, (n / (n-1)) * q * (1 - q) * F1) # var(z) in R divides by n - 1
  S <- with(obs, (1 - a0 - a1) / (p0 - p1))

  return(with(obs, K - F1 * N3 - S * N1 + (F2 * N4 + S * (N2 - 1)) * dz_tilde))
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

  K_first <- with(obs, g / (p0 - p1))
  K_second <- with(obs, ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (p0 - p1)))
  K <- K_first - K_second - h

  F1 <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
  F2 <- with(obs, (n / (n-1)) * q * (1 - q) * F1) # var(z) in R divides by n - 1
  S <- with(obs, (1 - a0 - a1) / (p0 - p1))

  denom <- with(obs, (F2 * N4 + S * (N2 - 1)))
  return(with(obs, (F1 * N3 + S * N1 - K) / (denom) + (1 / (denom)) *
                dTstar_tilde))
}


#' Bounds on m11 from variance conditions
#'
#' @param a0 Probability of observing T = 1 when true treatment status is 0.
#' @param a1 Probability of observing T = 0 when true treatment status is 1.
#' @param obs List of observables calculated using \code{getObs}.
#'
#' @return A list of bounds for m11 or NA if there are no bounds.
#' @export
#'
#' @examples
#' noCov <- getObs("totalRepeats", "usesch", "vouch0", data = angrist)
#' get_m11_bounds(0.2, 0.4, noCov)
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
  c00 <- with(obs, -1 * (ybStar10 * (ybStar10 - 2 * m00) - s2_00 / V00))
  c10 <- with(obs, -1 * (ybStar10 * (ybStar10 - 2 * m00) - s2_10 / V10))
  # Since these are the same quadratic up to a vertical shift only need
  # to look at the one with a *smaller* intercept
  c0 <- min(c00, c10)
  b0 <- -2 * m00
  disc0 <- b0^2 - 4 * 1 * c0
  if(disc0 > 0){
    roots_m10 <- sort(0.5 * (-b0 + c(-1, 1) * sqrt(disc0)))
  }else{
    roots_m10 <- c(NA, NA)
  }
  # Restrictions for m11
  c01 <- with(obs, -1 * (ybStar11 * (ybStar11 - 2 * m01) - s2_01 / V01))
  c11 <- with(obs, -1 * (ybStar11 * (ybStar11 - 2 * m01) - s2_11 / V11))
  # Since these are the same quadratic up to a vertical shift only need
  # to look at the one with a *smaller* intercept
  c1 <- min(c01, c11)
  b1 <- -2 * m00
  disc1 <- b1^2 - 4 * 1 * c1
  if(disc1 > 0){
    roots_m11 <- sort(0.5 * (-b1 + c(-1, 1) * sqrt(disc1)))
  }else{
    roots_m11 <- c(NA, NA)
  }
  roots <- list(m10 = roots_m10, m11 = roots_m11)
  # Case I: no bounds from var conditions.
  if(all(is.na(unlist(roots)))){
    return(list(NA))
  # Case II: bounds for m11 only from var conditions.
  }else if(all(is.na(roots$m10)) & all(!is.na(roots$m11))){
    return(list(roots$m11))
  # Case II: bounds for m10 only from var conditions
  }else if(all(!is.na(roots$m10)) & all(is.na(roots$m11))){
    return(list(roots$m10 - D)) # Implied bound for m11
  # Case IV: bounds for both m10 and m11 from var conditions
  }else{
    set1 <- roots$m10 - D # Implied bound for m11
    set2 <- roots$m11
    # Call set with smallest lower bound L (for lower) and the other O
    if(min(set1) <= min(set2)){
      L <- set1
      O <- set2
    }else{
      O <- set1
      L <- set2
    }
    # Case (a): O is a subset of L
    if(min(L) <= min(O) & max(L) >= max(O)){
     return(list(L))
    # Case (b): O is not a subset of L but O and L overlap
    }else if(min(L) <= min(O) & max(L) <= max(O) & min(O) <= max(L)){
      return(list(c(min(L), max(O))))
    # Case (c): O and L are disjoint
    }else{
      return(list(L, O))
    }
  }
}


#' If constraint is satisfied, calculate endogneity of binary regressor
#'
#' @param dz_tilde Endogeneity of binary instrument.
#' @param a0 Probability of observing T = 1 when true treatment status is 0.
#' @param a1 Probability of observing T = 0 when true treatment status is 1.
#' @param obs List of observables calculated by \code{getObs}.
#'
#' @return If the variance constraint is satistfied, return the value of regressor endogeneity implied by the data and specified values of instrument endogeneity and mis-classification probabilities. Otherwise return \code{NA}.
#' @export
#'
#' @examples
#' AngristControls <- c("age", "svy", "sex2", "phone", "hsvisit",
#'                      "d1995", "djamundi", "dmonth1", "dmonth2",
#'                      "dmonth3", "dmonth4", "dmonth5", "dmonth6", "dmonth7",
#'                      "dmonth8", "dmonth9", "dmonth10", "dmonth11",
#'                      "strata1", "strata2", "strata3", "strata4", "strata5")
#' AngristY <- "totalRepeats"
#' AngristT <- "usesch"
#' AngristZ <- "vouch0"
#' AngristObs <- getObs(AngristY, AngristT, AngristZ, AngristControls, angrist)
#' get_dTstar_tilde_check(0, 0.05, 0.05, AngristObs)
#' get_dTstar_tilde_check(0, 0.3, 0.4, AngristObs)
get_dTstar_tilde_check <- function(dz_tilde, a0, a1, obs){
  m11_bounds <- get_m11_bounds(a0, a1, obs)
  if(all(is.na(unlist(m11_bounds)))){ # No bounds
    return(get_dTstar_tilde(dz_tilde, a0, a1, obs))
  }else{ # Bounds
    g <- with(obs, (yt01 - yt00) - a1 * ((yt01 - yt00) + (yt11 - yt10)))
    D <- with(obs, ((1 - a0) * yt10 - a0 * yt00)/(p0 - a0) -
                ((1 - a0) * yt11 - a0 * yt01)/(p1 - a0))
    dz <- with(obs, N1 - (N2 - 1) * dz_tilde)
    m11 <- with(obs, (g - (1 - a0 - a1) * dz - (p0 - a0) * D) / (p0 - p1))
    constraint_fails <- sapply(m11_bounds, function(x)
      (m11 >= x[1] & m11 <= x[2]))
    if(any(constraint_fails)){
      return(NA)
    }else{
      return(get_dTstar_tilde(dz_tilde, a0, a1, obs))
    }
  }
}


#' If constraint is satisfied, calculate endogneity of binary instrument
#'
#' @param dTstar_tilde Endogeneity of binary regressor.
#' @param a0 Probability of observing T = 1 when true treatment status is 0.
#' @param a1 Probability of observing T = 0 when true treatment status is 1.
#' @param obs List of observables calculated by \code{getObs}.
#'
#' @return If the variance constraint is satistfied, return the value of instrument endogeneity implied by the data and specified values of instrument endogeneity and mis-classification probabilities. Otherwise return \code{NA}.
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
#' AngristObs <- getObs(AngristY, AngristT, AngristZ, AngristControls, angrist)
#' get_dz_tilde_check(-0.2, 0.1, 0.1, AngristObs)
#' get_dz_tilde_check(-0.1, 0.1, 0.2, AngristObs)
get_dz_tilde_check <- function(dTstar_tilde, a0, a1, obs){
  m11_bounds <- get_m11_bounds(a0, a1, obs)
  if(all(is.na(unlist(m11_bounds)))){ # No bounds
    return(get_dz_tilde(dTstar_tilde, a0, a1, obs))
  }else{ # Bounds
    h1 <- with(obs, (1 - q) * yt00 + q * yt01)
    h2 <- with(obs, (1 - q) * yt10 + q * yt11)
    h <- with(obs, (h1 - a1 * (h1 + h2)) / (1 - p - a1))
    D <- with(obs, ((1 - a0) * yt10 - a0 * yt00)/(p0 - a0) -
                ((1 - a0) * yt11 - a0 * yt01)/(p1 - a0))
    F1 <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
    F2 <- with(obs, (n / (n-1)) * q * (1 - q) * F1) # var(z) in R divides by n-1
    dTstar <- with(obs, dTstar_tilde + F1 * N3 -
      F2 * N4 * get_dz_tilde(dTstar_tilde, a0, a1, obs))
    m11 <- with(obs, ((p - a0) * (dTstar + h) - (1 - q) * (p0 - a0) * D) /
      ((1 - q) * (p0 - a0) + q * (p1 - a0)))
    constraint_fails <- sapply(m11_bounds, function(x)
      (m11 >= x[1] & m11 <= x[2]))
    if(any(constraint_fails)){
      return(NA)
    }else{
      return(get_dz_tilde(dTstar_tilde, a0, a1, obs))
    }
  }
}


