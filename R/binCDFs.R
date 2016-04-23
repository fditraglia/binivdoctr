#' Make plot of bounds on a0 and a1 implied by observable CDFs
#'
#' @param inData Dataframe with columns named Tobs (binary), y, and z (binary)
#' @param inc Increment for tau
#' @param frac Location of dashed lines: quantiles of y distribution
#'
#' @return Make plot of bounds on a0 and a1
#' @export
#'
#' @examples
plotAlphaBounds <- function(inData, inc = 0.01, frac = 0.05){
  Tobs <- inData$Tobs
  z <- inData$z
  y <- inData$y
  p0 <- mean(Tobs[z == 0])
  p1 <- mean(Tobs[z == 1])
  y00 <- y[Tobs == 0 & z == 0]
  y01 <- y[Tobs == 0 & z == 1]
  y10 <- y[Tobs == 1 & z == 0]
  y11 <- y[Tobs == 1 & z == 1]
  F00 <- ecdf(y00)
  F01 <- ecdf(y01)
  F10 <- ecdf(y10)
  F11 <- ecdf(y11)

  F00tilde <- function(tau){
    (1 - p0) * F00(tau)
  }
  F01tilde <- function(tau){
    (1 - p1) * F01(tau)
  }
  F10tilde <- function(tau){
    p0 * F10(tau)
  }
  F11tilde <- function(tau){
    p1 * F11(tau)
  }

  # Bounds for a0
  B0_z0_a0 <- function(tau){
  F10tilde(tau) / (F10tilde(tau) + F00tilde(tau))
  }

  B1_z0_a0 <- function(tau){
  (1 - F10tilde(tau) - (1 - p0)) / (1 - F10tilde(tau) - F00tilde(tau))
  }

  B0_z1_a0 <- function(tau){
  F11tilde(tau) / (F11tilde(tau) + F01tilde(tau))
  }

  B1_z1_a0 <- function(tau){
  (1 - F11tilde(tau) - (1 - p1)) / (1 - F11tilde(tau) - F01tilde(tau))
  }

  # Bounds for a1
  B0_z0_a1 <- function(tau){
  F00tilde(tau) / (F00tilde(tau) + F10tilde(tau))
  }

  B1_z0_a1 <- function(tau){
  (1 - F00tilde(tau) - p0) / (1 - F00tilde(tau) - F10tilde(tau))
  }

  B0_z1_a1 <- function(tau){
  F01tilde(tau) / (F01tilde(tau) + F11tilde(tau))
  }

  B1_z1_a1 <- function(tau){
  (1 - F01tilde(tau) - p1) / (1 - F01tilde(tau) - F11tilde(tau))
  }

  tau_range_z0 <- c(max(min(y00), min(y10)), min(max(y00), max(y10)))
  tau_z0 <- seq(tau_range_z0[1], tau_range_z0[2], inc)
  tau_range_z1 <- c(max(min(y01), min(y11)), min(max(y01), max(y11)))
  tau_z1 <- seq(tau_range_z1[1], tau_range_z1[2], inc)

  plot_min <- min(tau_z0, tau_z1)
  plot_max <- max(tau_z0, tau_z1)

  # Set up side-by-side plot with legend beneath
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), oma = c(4, 1, 1, 1))

  # Make plot for a0
  a0_B0_z0 <- B0_z0_a0(tau_z0)
  a0_B1_z0 <- B1_z0_a0(tau_z0)
  a0_B0_z1 <- B0_z1_a0(tau_z1)
  a0_B1_z1 <- B1_z1_a0(tau_z1)
  a0_min <- min(c(a0_B0_z0, a0_B1_z0, a0_B0_z1, a0_B1_z1))
  a0_max <- max(c(a0_B0_z0, a0_B1_z0, a0_B0_z1, a0_B1_z1))

  plot(tau_z0, a0_B0_z0, type = 'l', xlab = expression(tau),
       ylab = expression(bar(alpha[0])), lwd = 2, ylim = c(a0_min, a0_max),
       xlim = c(plot_min, plot_max), col = "blue",
       main = substitute(paste(alpha[0], " Upper Bounds")))
  abline(v = quantile(y, frac), lwd = 2, lty = 2)
  abline(v = quantile(y, 1 - frac), lwd = 2, lty = 2)

  points(tau_z0, a0_B1_z0, type = 'l', col = "red", lwd = 2)

  points(tau_z1, a0_B0_z1, type = 'l', col = "green", lwd = 2)

  points(tau_z1, a0_B1_z1, type = 'l', col = "orange", lwd = 2)

  # Make plot for a1
  a1_B0_z0 <- B0_z0_a1(tau_z0)
  a1_B1_z0 <- B1_z0_a1(tau_z0)
  a1_B0_z1 <- B0_z1_a1(tau_z1)
  a1_B1_z1 <- B1_z1_a1(tau_z1)
  a1_min <- min(c(a1_B0_z0, a1_B1_z0, a1_B0_z1, a1_B1_z1))
  a1_max <- max(c(a1_B0_z0, a1_B1_z0, a1_B0_z1, a1_B1_z1))

  plot(tau_z0, a1_B0_z0, type = 'l', xlab = expression(tau),
       ylab = expression(bar(alpha[1])), lwd = 2, ylim = c(a1_min, a1_max),
       xlim = c(plot_min, plot_max), col = "blue",
       main = substitute(paste(alpha[1], " Upper Bounds")))
  abline(v = quantile(y, frac), lwd = 2, lty = 2)
  abline(v = quantile(y, 1 - frac), lwd = 2, lty = 2)

  points(tau_z0, a1_B1_z0, type = 'l', col = "red", lwd = 2)

  points(tau_z1, a1_B0_z1, type = 'l', col = "green", lwd = 2)

  points(tau_z1, a1_B1_z1, type = 'l', col = "orange", lwd = 2)

  # Add legend beneath side-by-side plots
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = c("B0, z0", "B1, z0", "B0, z1", "B1, z1"), xpd = TRUE,
         horiz = TRUE, inset = c(0,0), bty = "n",
         fill = c("blue", "red", "green", "orange"), cex = 1.3)

  par(op)
}



