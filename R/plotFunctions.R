#' Valid IV plot
#'
#' @param obs
#' @param a0_range
#' @param a1_range
#' @param theta
#' @param phi
#' @param forceZero
#'
#' @return
#' @export
#'
#' @examples
plot_valid_iv <- function(obs, alpha_grid = 0.01,
                          theta = 0, phi = 15,
                          forceZero = FALSE, xlab = "a0",
                          ylab = "a1",
                          zlab = "dTstar",
                          col = "dodgerblue2"){

  alpha_bounds <- get_alpha_bounds(obs)
  a0upper <- alpha_bounds[1]
  a1upper <- alpha_bounds[2]
  a0_range <- seq(0, a0upper, by = alpha_grid)
  a1_range <- seq(0, a1upper, by = alpha_grid)

  f <- function(a0, a1){
    get_dTstar_tilde(0, a0, a1, obs)
  }
  dTstar <- outer(a0_range, a1_range, f)
  dTstar_range <- range(dTstar, na.rm = TRUE)

  # Optionally, force the plot to include zero on z-axis
  if(forceZero){
    sameSign <- prod(sign(dTstar_range)) == 1
    if(sameSign){
      bothNegative <- sign(dTstar_range[1]) == -1
      if(bothNegative){
        dTstar_range <- c(dTstar_range[1], 0)
      }else{
        dTstar_range <- c(0, dTstar_range[2])
      }
    }
  }
  persp(a0_range, a1_range, dTstar, phi = phi, theta = theta,
        zlim = dTstar_range, ticktype = "detailed", xlab = xlab,
        ylab = ylab, zlab = zlab, col = col)
}



#' Fixed dTstar plot
#'
#' @param obs
#' @param dTstar
#' @param a0_range
#' @param a1_range
#' @param theta
#' @param phi
#' @param forceZero
#'
#' @return
#' @export
#'
#' @examples
plot_fixed_dTstar <- function(obs, dTstar, alpha_grid = 0.005, theta = 0,
                              phi = 15, forceZero = TRUE, xlab = "a0",
                              ylab = "a1", zlab = "dz",
                              col = "dodgerblue2"){

  alpha_bounds <- get_alpha_bounds(obs)
  a0upper <- alpha_bounds[1]
  a1upper <- alpha_bounds[2]
  a0_range <- seq(0, a0upper, by = alpha_grid)
  a1_range <- seq(0, a1upper, by = alpha_grid)

  f <- function(a0, a1){
      get_dz_tilde(dTstar, a0, a1, obs)
  }
  dz <- outer(a0_range, a1_range, f)
  dz_range <- range(dz, na.rm = TRUE)

  # Optionally, force the plot to include zero on z-axis
  if(forceZero){
    sameSign <- prod(sign(dz_range)) == 1
    if(sameSign){
      bothNegative <- sign(dz_range[1]) == -1
      if(bothNegative){
        dz_range <- c(dz_range[1], 0)
      }else{
        dz_range <- c(0, dz_range[2])
      }
    }
  }
  persp(a0_range, a1_range, dz, phi = phi, theta = theta, zlim = dz_range,
        ticktype = "detailed", xlab = xlab, ylab = ylab, zlab = zlab,
        col = col)
}



