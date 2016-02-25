# Function to calculate all the summary statistics of the data needed to use
# our method for the identified set *without* sampling uncertainty
getObs <- function(y_name, T_name, z_name, controls, data){
  # Extract data for named columns
  y <- get(y_name, data)
  Tobs <- get(T_name, data)
  z <- get(z_name, data)

  # Calculate all the observables
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

  second_stage <- reformulate(c(T_name, controls), response = y_name)
  first_stage <- reformulate(c(z_name, controls))
  obs$gamma_iv <- coefficients(ivreg(second_stage, first_stage, data))[-c(1,2)]
  # No intercept since we "project it out" by working with Cov matrix below
  x <- model.matrix(reformulate(controls, intercept = FALSE), data)
  obs$D_ETobs <- colMeans(x[Tobs == 1,]) - colMeans(x[Tobs == 0,])
  invert_me <- rbind(cbind(cov(z, Tobs), cov(z, x)),
                     cbind(cov(x, Tobs), cov(x)))
  obs$s_xT_upper <- matrix(solve(invert_me)[-1,1], ncol(x), 1)

  return(obs)
}

#AngristControls <- c("age", "svy", "sex2", "phone", "hsvisit",
#                     "d1995", "djamundi", "dmonth1", "dmonth2",
#                     "dmonth3", "dmonth4", "dmonth5", "dmonth6", "dmonth7",
#                     "dmonth8", "dmonth9", "dmonth10", "dmonth11",
#                     "strata1", "strata2", "strata3", "strata4", "strata5")
#AngristY <- "totalRepeats"
#AngristT <- "usesch"
#AngristZ <- "vouch0"
#
#getObs(AngristY, AngristT, AngristZ, AngristControls, angrist)
