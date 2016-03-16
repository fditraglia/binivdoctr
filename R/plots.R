plot_valid_iv <- function(obs, a0_range = seq(0, 0.2, 0.01),
                          a1_range = seq(0, 0.2, 0.01), theta = 0, phi = 15,
                          forceZero = FALSE){
  f <- function(a0, a1){
    get_dTstar_tilde_check(0, a0, a1, obs)
  }
  fv <- Vectorize(f)
  dTstar <- outer(a0_range, a1_range, fv)

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
        zlim = dTstar_range, ticktype = "detailed", xlab = "a0", ylab = "a1",
        zlab = "dTstar")
}



plot_fixed_dTstar <- function(obs, dTstar, a0_range = seq(0, 0.2, 0.01),
                              a1_range = seq(0, 0.2, 0.01), theta = 0,
                              phi = 15, forceZero = TRUE){
  f <- function(a0, a1){
    get_dz_tilde_check(dTstar, a0, a1, obs)
  }
  fv <- Vectorize(f)
  dz <- outer(a0_range, a1_range, fv)

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
        ticktype = "detailed", xlab = "a0", ylab = "a1", zlab = "dz")
}



