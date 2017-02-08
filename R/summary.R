#' Title
#'
#' @param y_name
#' @param T_name
#' @param z_name
#' @param data
#' @param controls
#' @param robust
#'
#' @return
#' @export
#'
#' @examples
get_summary_stats <- function(y_name, T_name, z_name, data, controls = NULL,
                          robust = FALSE) {
  first_stage <- reformulate(c(z_name, controls), response = NULL)
  second_stage <- reformulate(c(T_name, controls), response = y_name)
  OLS <- lm(second_stage, data)
  IV <-  AER::ivreg(second_stage, first_stage, data)

  b_OLS <- coef(OLS)[T_name]
  b_IV <- coef(IV)[T_name]

  if (robust) {
    se_OLS <- sqrt(diag(sandwich::vcovHC(OLS, type = 'HC0')))[T_name]
    se_IV <- sqrt(diag(sandwich::vcovHC(IV, type = 'HC0')))[T_name]
  } else {
    se_OLS <- sqrt(diag(vcov(OLS)))[T_name]
    se_IV <- sqrt(diag(vcov(IV)))[T_name]
  }

  obs <- getObs(y_name, T_name, z_name, controls, data)
  alpha_bounds <- get_alpha_bounds(obs)

  list(n = nrow(data),
       b_OLS = b_OLS,
       se_OLS = se_OLS,
       b_IV = b_IV,
       se_IV = se_IV,
       a0_upper = alpha_bounds[1],
       a1_upper = alpha_bounds[2])
}

#' Title
#'
#' @param draws
#' @param level
#'
#' @return
#' @export
#'
#' @examples
get_HPDI <- function(draws, level = 0.9){
  interval <- coda::HPDinterval(coda::as.mcmc(draws), level)
  lower <- interval[[1]]
  upper<- interval[[2]]
  return(data.frame(lower = lower, median = median(draws), upper = upper))
}

#' Title
#'
#' @param draws
#'
#' @return
#' @export
#'
#' @examples
summarize_dz_draws <- function(draws){

  a0_bounds <- t(sapply(draws$IS, function(x) range(x$a0)))
  a0_upper <- a0_bounds[,2]
  a1_bounds <- t(sapply(draws$IS, function(x) range(x$a1)))
  a1_upper <- a1_bounds[,2]

  b_bounds_analyt <- get_beta_bounds(draws$dTstar_tilde_range[1],
                              draws$dTstar_tilde_range[2],
                              a0_upper+0.15, a1_upper+0.15,draws$RF,
                              draws$evaluateInterior)

  b_bounds_sim <- t(sapply(draws$IS, function(x) range(x$beta)))
  
  b_lower <- pmin(b_bounds_analyt[,1],b_bounds_sim[,1],na.rm=TRUE)
  b_upper <- pmax(b_bounds_analyt[,2],b_bounds_sim[,2],na.rm=TRUE)
  
  b_lower <- b_bounds_analyt[,1]
  b_upper <- b_bounds_analyt[,2]
  
  #dz_bounds <- t(sapply(draws$IS, function(x) range(x$dz_tilde)))
  dz_bounds_analyt <- get_dz_tilde_bounds(draws$dTstar_tilde_range[1],
                                    draws$dTstar_tilde_range[2],
                                    a0_upper, a1_upper,
                                    draws$RF,
                                    draws$evaluateInterior)

  dz_bounds_sim <- t(sapply(draws$IS, function(x) range(x$dz_tilde)))
  
  dz_lower <- pmin(dz_bounds_analyt[,1],dz_bounds_sim[,1],na.rm=TRUE)
  dz_upper <- pmax(dz_bounds_analyt[,2],dz_bounds_sim[,2],na.rm=TRUE)

  IS <- do.call(rbind, draws$IS)
  b_draws <- IS$beta
  dz_draws <- IS$dz_tilde

  rbind(b_lower = get_HPDI(b_lower),
        b_upper = get_HPDI(b_upper),
        dz_lower = get_HPDI(dz_lower),
        dz_upper = get_HPDI(dz_upper),
        a0_upper = get_HPDI(a0_upper),
        a1_upper = get_HPDI(a1_upper),
        dz_bayes = get_HPDI(dz_draws),
        b_bayes = get_HPDI(b_draws))
}


#' Title
#'
#' @param draws
#'
#' @return
#' @export
#'
#' @examples
summarize_dTstar_draws <- function(draws){
  
  # b_bounds <- t(sapply(draws$IS, function(x) range(x$beta)))
  b_bounds <- t(sapply(draws$IS, function(x) range(x$beta)))
  b_lower <- b_bounds[,1]
  b_upper <- b_bounds[,2]
  
  dTstar_bounds <- t(sapply(draws$IS, function(x) range(x$dTstar_tilde)))
  dTstar_lower <- dTstar_bounds[,1]
  dTstar_upper <- dTstar_bounds[,2]

  a0_bounds <- t(sapply(draws$IS, function(x) range(x$a0)))
  a0_upper <- a0_bounds[,2]
  a1_bounds <- t(sapply(draws$IS, function(x) range(x$a1)))
  a1_upper <- a1_bounds[,2]

  IS <- do.call(rbind, draws$IS)
  b_draws <- IS$beta
  dTstar_draws <- IS$dTstar_tilde

  rbind(b_lower = get_HPDI(b_lower),
        b_upper = get_HPDI(b_upper),
        dTstar_lower = get_HPDI(dTstar_lower),
        dTstar_upper = get_HPDI(dTstar_upper),
        a0_upper = get_HPDI(a0_upper),
        a1_upper = get_HPDI(a1_upper),
        dTstar_bayes = get_HPDI(dTstar_draws),
        b_bayes = get_HPDI(b_draws))
}
