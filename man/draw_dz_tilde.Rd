% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SamplingUncertainty.R
\name{draw_dz_tilde}
\alias{draw_dz_tilde}
\title{Make draws for dz_tilde}
\usage{
draw_dz_tilde(y_name, T_name, z_name, controls = NULL, data,
  dTstar_tilde_range, a0bound = NULL, a1bound = NULL, nRF = 500,
  nIS = 500, option = NULL)
}
\arguments{
\item{PequalPstar}{}
}
\description{
Make draws for dz_tilde
}
\examples{
afghanControls <- c("headchild", "age",  "yrsvill",  "farsi",  "tajik",
                   "farmers",  "agehead",  "educhead",  "nhh",  "land",
                   "sheep", "distschool", "chagcharan")
foo <- draw_dz_tilde(y_name = "testscore", T_name = "enrolled",
                     z_name = "buildschool", controls = afghanControls,
                     data = afghan, dTstar_tilde_range = c(0, 1))
}

