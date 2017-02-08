# #### This file contains three functions:
# 1. make_I(): Initial tex line for each example
# 2. make_II_III(): Line for each prior with each example
# 3. makeExample(): Carries out numerical calculations and uses the above
#    functions as building blocks to generate whole section of table devoted
#    to example

# ------------------------ Generic Helper Functions
myformat <- function(x){
  x <- ifelse(is.na(x), NA, format(round(x, 2), n_digits = 2, nsmall = 2))
}
format_est <- function(est) {
  paste0('$', myformat(est), '$')
}
format_se <- function(se) {
  paste0('$(', myformat(se), ')$')
}
format_HPDI <- function(bounds) {
  paste0('$[', myformat(bounds[1]), ',', myformat(bounds[2]), ']$')
}

make_tex_row <- function(char_vec, shift = 0) {
  out <- paste0(char_vec, collapse = ' & ')
  if(identical(shift, 0)){
    out <- paste(out, '\\\\')
  } else {
    out <- paste(paste0(rep('&', shift), collapse = ''), out, '\\\\')
  }
  return(out)
}

# ------------------------ Functions Specific to Binary Case
make_I <- function(stats, example_name) {
  example_name <- paste(example_name, paste0('($n=', stats$n, '$)'))
  est <- with(stats, sapply(c(b_OLS, b_IV, a0_upper, a1_upper), format_est))
  est <- make_tex_row(c(example_name, est))
  se <- with(stats, sapply(c(se_OLS, se_IV), format_se))
  se <- make_tex_row(se, shift = 1)
  paste(est, se, sep = '\n')
}

make_II_III <- function(stats, prior_name) {

  if('dTstar_bayes' %in% rownames(stats)) {
    medians <- c(format_est(stats['dTstar_lower', 'median']),
                 format_est(stats['dTstar_upper', 'median']),
                 format_est(stats['b_lower', 'median']),
                 format_est(stats['b_upper', 'median']),
                 format_est(stats['dTstar_bayes', 'median']),
                 format_est(stats['b_bayes', 'median']))

    HPDIs <- c(format_HPDI(stats['dTstar_lower', 'lower'],
                           stats['dTstar_lower', 'upper']),
               format_HPDI(stats['dTstar_upper', 'lower'],
                           stats['dTstar_upper', 'upper']),
               format_HPDI(stats['b_lower', 'lower'],
                           stats['b_lower', 'upper']),
               format_HPDI(stats['b_upper', 'lower'],
                           stats['b_upper', 'upper']),
               format_HPDI(stats['dTstar_bayes', 'lower'],
                           stats['dTstar_bayes', 'upper']),
               format_HPDI(stats['b_bayes', 'lower'],
                           stats['b_bayes', 'upper']))
  } else {
    medians <- c(format_est(stats['dz_lower', 'median']),
                 format_est(stats['dz_upper', 'median']),
                 format_est(stats['b_lower', 'median']),
                 format_est(stats['b_upper', 'median']),
                 format_est(stats['dz_bayes', 'median']),
                 format_est(stats['b_bayes', 'median']))

    HPDIs <- c(format_HPDI(stats['dz_lower', 'lower'],
                           stats['dz_lower', 'upper']),
               format_HPDI(stats['dz_upper', 'lower'],
                           stats['dz_upper', 'upper']),
               format_HPDI(stats['b_lower', 'lower'],
                           stats['b_lower', 'upper']),
               format_HPDI(stats['b_upper', 'lower'],
                           stats['b_upper', 'upper']),
               format_HPDI(stats['dz_bayes', 'lower'],
                           stats['dz_bayes', 'upper']),
               format_HPDI(stats['b_bayes', 'lower'],
                           stats['b_bayes', 'upper']))
  }
  row1 <- paste('\\hspace{2em}', prior_name, make_tex_row(medians, shift = 5))
  row2 <- make_tex_row(HPDIs, shift = 5)
  paste(row1, row2, sep = '\n')
}

#' Generates table of parameter estimates given user restrictions and data
#'
#' @param y_name Character string with the column name of the dependent variable
#' @param T_name Character string with the column name of the endogenous regressor(s)
#' @param z_name Character string with the column name of the instrument(s)
#' @param data Data frame
#' @param controls Vector of character strings specifying the exogenous variables
#' @param robust Indicator for heteroskedasticity-robust standard errors
#' @param r_TstarU_restriction Matrix of desired mins and maxes for r_TstarU (must be same dimensions as k_restriction)
#' @param k_restriction Matrix of desired minx and maxes for kappa (must be same dimensions as r_TstarU_restriction)
#' @param n_draws Number of draws when generating frequentist-friendly draws of the covariance matrix
#' @param n_RF_draws Number of reduced-form draws
#' @param n_IS_draws Number of draws on the identified set
#' @param resample Indicator of whether or not to resample using magnification factor
#' @param Jeffreys Indicator of whether to draw covariance matrix using the Jeffreys prior to ensure positive definiteness
#' @param example_name Character string describing the example
#' @export
makeExample <- function(y_name, T_name, z_name, data, controls = NULL,
                        robust = FALSE,
                        dTstar_tilde_range = NULL,
                        dz_tilde_range = NULL,
                        a0_restriction = NULL,
                        a1_restriction = NULL,
                        n_draws = 5000, n_RF_draws = 1000,
                        n_IS_draws = 1000, resample = FALSE, Jeffreys = FALSE,
                        example_name) {

  # if (is.null(r_TstarU_restriction)) {
  #   r_TstarU_restriction <- matrix(c(-1, 1), nrow = 1)
  # }
  # if (is.null(k_restriction)) {
  #   k_restriction <- matrix(c(0, 1), nrow = 1)
  # }
  # if (nrow(r_TstarU_restriction) != nrow(k_restriction)) {
  #   if (is.null(r_TstarU_restriction)) {
  #     r_TstarU_restriction <- matrix(c(-1, 1), nrow = nrow(k_restriction), byrow = TRUE)
  #   } else if (is.null(k_restriction)) {
  #     k_restriction <- matrix(c(0, 1), nrow = nrow(r_TstarU_restriction), byrow = TRUE)
  #   } else {
  #     stop("Dimension mismatch between r_TstarU_restriction and k_restriction.
  #          Please make sure that if there are restrictions on kappa or r_TstarU
  #          that they are of the same dimension so that all examples are accounted
  #          for.")
  #   }
  #   }

  ## NOTE: INCLUDE AN OPTION FOR WHETHER THE USER WANTS TO SPECIFY A VALID
  ## INSTRUMENT OR NOT. IF THE ARGUMENT FOR BOUNDS ON DTSTAR_U. PROBABLY
  ## BETTER TO INCLUDE THIS INSIDE THE

  summary_stats <- get_summary_stats(y_name, T_name, z_name, data, controls,robust)
  obs <- getObs(y_name, T_name, z_name, controls, data)

  headline <- make_I(stats_I, example_name)
  nExamples <- nrow(dTstar_tilde_range)
  exampleTex <- NULL

  for (i in 1:nExamples) {

    draws <- draw_dTstar_tilde_valid(y_name, T_name, z_name, controls,
                                     data, a0bound, a1bound,
                                     nRF = n_RF_draws, nIS = n_IS_draws)

    bayes <-  summarize_dz_draws(draws)

    ## Introduce conditionals on what information to report.

    newRow <- make_II_III(stats, paste0("$(\\kappa, \\rho_{T^*u}) \\in (",
                                        k_restriction[i, 1], ",",
                                        k_restriction[i, 2],
                                        ") \\times [",
                                        r_TstarU_restriction[i, 1], ",",
                                        r_TstarU_restriction[i, 2], "]"))


    exampleTex <- paste(exampleTex, newRow, sep = "/n")
  }
  output <- paste(headline, exampleTex, sep = "/n")
  return(output)
    }
