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
format_HPDI <- function(lower, upper) {
  paste0('$[', myformat(lower), ',', myformat(upper), ']$')
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
#' @param dTstar_tilde_range Matrix of desired mins and maxes for dTstar_tilde_range (must be same dimensions as k_restriction)
#' @param dz_tilde_range Matrix of desired mins and maxes for dTstar_tilde_range (must be same dimensions as r_TstarU_restriction)
#' @param a0bound
#' @param a1bound
#' @param n_draws Number of draws when generating frequentist-friendly draws of the covariance matrix
#' @param n_RF_draws Number of reduced-form draws
#' @param n_IS_draws Number of draws on the identified set
#' @param resample Indicator of whether or not to resample using magnification factor
#' @param Jeffreys Indicator of whether to draw covariance matrix using the Jeffreys prior to ensure positive definiteness
#' @param example_name Character string describing the example
#' @export
makeExample <- function(y_name, T_name, z_name, data, controls = NULL,
                        dTstar_tilde_range = NULL,
                        dz_tilde_range = NULL,
                        a0bound= NULL, a1bound = NULL,
                        n_draws = 5000, n_RF_draws = 1000,
                        n_IS_draws = 1000, resample = FALSE, Jeffreys = FALSE,
                        example_name,
                        option,
                        evaluateInterior=FALSE) {

  # Baseline table object

  table_obj <- list(y_name = y_name,
                    T_name = T_name,
                    z_name = z_name,
                    data = data,
                    controls = controls,
                    n_draws = n_draws,
                    n_RF_draws = n_RF_draws,
                    n_IS_draws = n_IS_draws,
                    resample = resample,
                    Jeffreys = Jeffreys,
                    example_name = example_name,
                    evaluateInterior = evaluateInterior)

  # Obtain summary statistics

  summary_stats <- do.call(get_summary_stats, table_obj)
  headline <- make_I(summary_stats, example_name)

  nExamples <- nrow(dTstar_tilde_range)
  exampleTex <- NULL

  for (i in 1:nExamples) {

    ## Define appropriate labels

    #print(i)

    text_PequalPstar <- ""
    text_A0equalA1 <- ""
    text_a0bound <- ""
    text_a1bound <- ""
    option_default <- NULL

    if(length(option) >= i) {

      #print("here we go")
      #print(option[i])

      if(!is.na(option[i])) {

        if(option[i] == "PequalPstar") {

          text_PequalPstar <- c("$\\;p=p^*$")

        } else if (option[i] == "A0equalA1") {

          text_A0equalA1 <- c("$\\;\\alpha_0=\\alpha_1$")

        }

      } else {

        option_default <- option[i]

      }

    }

    if( length(a0bound) >= i) {

      if( a0bound[i] == 0 ) {
        text_a0bound <- "$\\;\\alpha_0 = 0$"
      } else if(a0bound[i] < 1 ){
        text_a0bound <- paste0("$\\;\\alpha_0 < $",a0bound[i])
      }
    }

    if( length(a1bound) >= i ) {

      if( a1bound[i] == 0 ) {
        text_a1bound <- "$\\;\\alpha_1 = 0$"
      } else if( a1bound[i] < 1 ) {
        text_a1bound <- paste0("$\\;\\alpha_1 < $",a1bound[i])
      }
    }

    ## Obtain summary statistics

    temp_table_obj <- c(table_obj,
                        list(dTstar_tilde_range = dTstar_tilde_range[i,],
                             a0bound = a0bound[i],
                             a1bound = a1bound[i],
                             option = option[i]))

    if(!is.na(dTstar_tilde_range[i,1])) {

      temp_summary <- summarize_dz_draws(do.call(draw_dz_tilde, temp_table_obj))

      newRow <- make_II_III(temp_summary, paste0("$\\delta_{T^*} \\in [",
                                          dTstar_tilde_range[i,1],
                                          ",",
                                          dTstar_tilde_range[i,2],
                                          "]$",
                                          text_a0bound,text_a1bound,
                                          text_PequalPstar,text_A0equalA1))

    } else {

      temp_summary <- summarize_dTstar_draws(do.call(draw_dTstar_tilde_valid,
                                                     temp_table_obj))

      newRow <- make_II_III(temp_summary, paste0("$\\delta_z=0$",
                                          text_a0bound,text_a1bound))
    }

    ## Introduce conditionals on what information to report.

    exampleTex <- paste(exampleTex, newRow, sep = "")
  }
  output <- paste(headline, exampleTex, sep = "")
  return(output)
}

# verifiyArguments <- function(y_name, T_name, z_name, data, controls = NULL,
#                              robust = FALSE,
#                              dTstar_tilde_range = NULL,
#                              dz_tilde_range = NULL,
#                              a0bound= NULL, a1bound = NULL,
#                              n_draws = 5000, n_RF_draws = 1000,
#                              n_IS_draws = 1000, resample = FALSE, Jeffreys = FALSE,
#                              example_name,
#                              option) {
#
#   if (is.null(dTstar_tilde_range) & is.null(dz_tilde_range)) {
#     stop("Execution was stopped because no argument was entered
#          for dTstar_tilde_range or dz_tilde_range.
#          Please specify a range.")
#   } else if (is.null(dTstar_tilde_range) & !is.null(dz_tilde_range))  {
#
#     if(dim(dTstar_tilde_range)[2] > 3 ){
#
#       stop("Execution was stopped because dTstar_tilde_range has 3 or more
#            columns. dTstar_tilde_range range should be a Mx2 matrix, where
#            M is the number of priors being evaluated.")
#     }
#
#     if( is.null(a0bound) ) {
#       a0bound <- matrix(c(1),nrow=dim(dTstar_tilde_range)[1])
#
#     } else if(dim(a0bound)[2] != 1) {
#       stop("Execution was stopped because a0bound is not a vector")
#
#     } else if(dim(a0bound)[1] != dim(dim(dTstar_tilde_range)[1])) {
#
#       stop("Execution was stopped because dTstar_tilde_range and a0 bound
#            do not have the same number of rows.")
#
#     }
#
#     if( is.null(a1bound) ) {
#       a1bound <- matrix(c(1),nrow=dim(dTstar_tilde_range)[1])
#
#     } else if(dim(a1bound)[2] != 1) {
#       stop("Execution was stopped because a1bound is not a vector")
#
#     } else if(dim(a1bound)[1] != dim(dim(dTstar_tilde_range)[1])) {
#
#       stop("Execution was stopped because dTstar_tilde_range and a1 bound
#            do not have the same number of rows.")
#
#     }
#
#
#     } else if (!is.null(dTstar_tilde_range) & is.null(dz_tilde_range))  {
#
#
#     } else {
#
#
#
#     }
#   is.null(dTstar_tilde_range) & is.null(dz_tilde_range)
#
#
#
#
#
#     }
#
# }

table_header_fn <- function(){

  table_header <- '\\begin{tabular}{lcccccccccc}
                 \\hline \\hline
  &\\multicolumn{4}{c}{(I) Summary Statistics}
  &\\multicolumn{4}{c}{(II) Frequentist-Friendly}
  &\\multicolumn{2}{c}{(III) Fully Bayesian} \\\\
  \\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-11}
  & OLS & IV & $\\bar{\\alpha}_0$ & $\\bar{\\alpha}_1$
  & $\\underline{\\delta}_{T^*/z}$ & $\\bar{\\delta}_{T^*/z}$
  & $\\underline{\\beta}$ & $\\bar{\\beta}$
  & $\\delta_{T^*/z}$ & $\\beta$ \\\\'

}

table_footer_fn <- function() {

  table_footer <- '\\hline
  \\end{tabular}'

}

makeTable <- function(file,...) {

  table_examples <- c()
  list_examples <- list(...)

  for( i in 1:length(list(...))) {

   table_examples <- c(table_examples,list_examples[[i]],"\\\\")

  }

  cat(table_header_fn(),'\\\\',
      table_examples,
      table_footer_fn(),
      sep = '\n',file=file)

}
