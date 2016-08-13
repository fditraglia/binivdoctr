HPDI <- function(draws, level = 0.9){
  interval <- coda::HPDinterval(coda::as.mcmc(draws), level)
  lower <- interval[[1]]
  upper<- interval[[2]]
  return(data.frame(median = median(draws), lower = lower, upper = upper))
}

#' Summarize posterior draws from draw_dz_tilde
#'
#' @param draws
#' @param digits
#' @return
#' @export
#'
#' @examples
summarize_dz_draws <- function(draws, digits = 2){
  b_bounds <- t(sapply(draws$IS, function(x) range(x$beta)))
  b_lower <- b_bounds[,1]
  b_upper <- b_bounds[,2]

  dz_bounds <- t(sapply(draws$IS, function(x) range(x$dz_tilde)))
  dz_lower <- dz_bounds[,1]
  dz_upper <- dz_bounds[,2]

  a0_bounds <- t(sapply(draws$IS, function(x) range(x$a0)))
  a0_upper <- a0_bounds[,2]
  a1_bounds <- t(sapply(draws$IS, function(x) range(x$a1)))
  a1_upper <- a1_bounds[,2]

  IS <- do.call(rbind, draws$IS)
  b_draws <- IS$beta
  dz_draws <- IS$dz_tilde

  out <- rbind(b_lower = HPDI(b_lower),
               b_upper = HPDI(b_upper),
               dz_lower = HPDI(dz_lower),
               dz_upper = HPDI(dz_upper),
               a0_upper = HPDI(a0_upper),
               a1_upper = HPDI(a1_upper),
               dz_bayes = HPDI(dz_draws),
               b_bayes = HPDI(b_draws))
  return(round(out, digits))
}


#' Summarize posterior draws from draw_dTstar_tilde
#'
#' @param draws
#' @param digits
#'
#' @return
#' @export
#'
#' @examples
summarize_dTstar_draws <- function(draws, digits = 2){
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

  out <- rbind(b_lower = HPDI(b_lower),
               b_upper = HPDI(b_upper),
               dTstar_lower = HPDI(dTstar_lower),
               dTstar_upper = HPDI(dTstar_upper),
               a0_upper = HPDI(a0_upper),
               a1_upper = HPDI(a1_upper),
               dTstar_bayes = HPDI(dTstar_draws),
               b_bayes = HPDI(b_draws))
  return(round(out, digits))
}
