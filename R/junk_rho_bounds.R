library('nleqslv')

get_dz_tilde_bounds <- function(a0_upper, a1_upper,
                                dTstar_tilde_lower,dTstar_tilde,upper,
                                obs) {
  
  # a0_upper <- 0.5
  # a1_upper <- 0.5
  
  # Obtain bounds for a0 and a1 implied by observables, without impoising user beliefs.
  unrestricted <- get_alpha_bounds(obs)
  
  # Storing the more restrictive bounds for each parameter
  # (either user or unrestricted)
  a0_upper_bound <- ifelse(a0_upper >= unrestricted['a0upper'],
                           unrestricted['a0upper'],
                           a0_upper)
  a1_upper_bound <- ifelse(a1_upper >= unrestricted['a1upper'],
                           unrestricted['a1upper'],
                           a1_upper)
  
  # Finding dz_tilde bounds based on restricted bounds above
  # Candidate Set I - Corner for both a0 and a1
  set1a <- candidate1(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set1b <- candidate1(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)
  
  # Candidate Set II - Corner for a0 and interior for a1
  set2a <- candidate2(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set2b <- candidate2(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)
  
  # Candidate Set III - Interior for a0, corner for rho_TstarU
  set3a <- candidate3(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set3b <- candidate3(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)
  
  # Candidate Set IV - Interior for both a0 and a1
  set4a <- candidate4(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set4b <- candidate4(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)
  
  # Finally: overall max and min
  r_uz_max <- pmax(set1a$max_corner, set1b$max_corner,
                   set2a$max_edge,set2b$max_edge,
                   set3a$r_uz$max_edge,set3b$r_uz$max_edge na.rm = TRUE)
  
  r_uz_min <- pmin(set1$min_corner, set2$min_edge, set3$r_uz$min_edge, na.rm = TRUE)
  data.frame(min = r_uz_min, max = r_uz_max)
  
}


# Evaluates the corners given user bounds. Vectorized wrt a0 and a1 bounds.
candidate1 <- function(dTstar_tilde_bound, a0_upper, a1_upper, obs) {
  
  # dTstar_tilde_lower <- 1
  # dTstar_tilde_upper <- 2
  # dTstar_tilde_bound <- 1
  
  corner1 <- get_dz_tilde(dTstar_tilde_bound,0,0,obs)
  corner2 <- get_dz_tilde(dTstar_tilde_bound,a0_upper,0,obs)
  corner3 <- get_dz_tilde(dTstar_tilde_bound,0,a1_upper,obs)
  corner4 <- get_dz_tilde(dTstar_tilde_bound,a0_upper,a1_upper,obs)
  min_corner <- pmin(corner1, corner2, corner3, corner4, na.rm = TRUE)
  max_corner <- pmax(corner1, corner2, corner3, corner4, na.rm = TRUE)
  ans <- list(min_corner = min_corner,
              max_corner = max_corner)
  return(ans)
}

# Evaluates the edge where a0 is on the boundary. Vectorized wrt a0 and a1 bounds.
candidate2 <- function(dTstar_tilde_bound, a0,a1upper, obs) {
  
  a1upper <- get_alpha_bounds(obs)[2]
  a0 <- get_alpha_bounds(obs)[1]
  dTstar_tilde_bound <- dxs
  
  g0 <- with(obs,yt01-yt00)
  g1 <- with(obs,(yt01-yt00)-(yt11-yt10))
  h0 <- with(obs,(1-q)*yt00+q*yt01)
  h1 <- with(obs,((1-q)*yt00+q*yt01)+((1-q)*yt01+q*yt11))
  Delta01 <- with(obs,yt10)
  Delta11 <- with(obs,yt10+yt00)
  Delta02 <- with(obs,yt11)
  Delta12 <- with(obs,yt11+yt10)
  DeltaA0 <- with(obs,((Delta01-Delta11*a0)/(p-a0))-((Delta02-Delta12*a0)/(p-a0)))
  
  # Getting quadratic coefficients
  c <- g0-(1-a0)*g1-(p0-p1)*h1-((p0-a0)*(p1-a0)/(p-a0))*DeltaA0+
    (p1-p0)*dTstar_tilde_bound
  b <- 2*(1-p)*(-g0+2*(1-p)*(1-a0)*g1-(p1-p0)*dTsartilde_bound)+
    2*(1-p)*((p0-a0)*(p1-a0)/(p-a0))*DeltaA0 +
    (p0-p1)*(((1-p)*h1+h0)+(h0-(1-p)*h1))
  a <- (1-(p^2))*(g0-(1-a0)*g1+(p1-p0)*dTstar_tilde_bound)+
    (1-(p^2))*((p0-a0)*(p1-a0)/(p-a0))*DeltaA0-
    (p0-p1)*((1-p)-(1-a0)*h0-(1-a0)*(1-p)*h1)
  coefs <- cbind(a, b, c)
  
  # Getting roots of cubic
  all_roots <- apply(coefs, 1, polyroot)
  
  # Getting roots of cubic
  all_roots <- apply(coefs, 1, polyroot)
  real_roots <- Re(unlist(all_roots))
  real_roots <- real_roots[(real_roots >= 0) & (real_roots <= a0upper)]
  
  # Evaluating all combinations of roots and bounds for r_TstarU
  dz_tilde <- outer(real_roots, a1,
                    function(a0,a1) get_dz_tilde(dTstar_tilde_bound,a1,a0,obs))
  
  # Getting max and min if there are real roots
  min_edge <- ifelse(length(real_roots) > 0, min(dz_tilde), NA)
  max_edge <- ifelse(length(real_roots) > 0, max(dz_tilde), NA)
  ans <- list(dz_tilde = list(min_edge = min_edge, max_edge = max_edge),
              a0_roots = real_roots)
  return(ans)
  
}

# Evaluates the edge where a1 is on the boundary.
candidate3 <- function(dTstar_tilde_bound,a0upper,a1, obs) {
  
  a0upper <- get_alpha_bounds(obs)[1]
  a1 <- get_alpha_bounds(obs)[2]
  dTstar_tilde_bound <- 2
  
  g0 <- with(obs,yt01-yt00)
  g1 <- with(obs,(yt01-yt00)-(yt11-yt10))
  h0 <- with(obs,(1-q)*yt00+q*yt01)
  h1 <- with(obs,((1-q)*yt00+q*yt01)+((1-q)*yt01+q*yt11))
  Delta01 <- with(obs,yt10)
  Delta11 <- with(obs,yt10+yt00)
  Delta02 <- with(obs,yt11)
  Delta12 <- with(obs,yt11+yt10)
  p0 <- with(obs,p0)
  p1 <- with(obs,p1)
  p  <- with(obs,p)
  h_a1 <- (h0-h1*a1)/(1-p-a1)
  
  # Getting quadratic coefficients
  d <- (Delta02-Delta01)-(Delta12-Delta11)
  c <- (g0-g1*a1)-(p0-p1)*h_a1+Delta01*((p-a0)+p1)+
    Delta02*(-(p-a1)-p0)+
    2*Delta11*p1-2*p0*Delta12+(p1-p0)*dTstar_tilde_bound
  b <- -2*p*(g0-g1*a1)-(p0-p1)+2*p*(p0-p1)*h_a1+
    (p*p1*Delta11-(p+p1)*Delta01)-(p*p0*Delta12-(p+p0)*Delta02)-
    (1-a1)*(Delta01*(p+p1)+Delta11*(p-p0)-Delta02*(p+p0)-Delta12*(p-p0))+
    Delta01*(p1-p*p1-p)-Delta02*(p0-p*p0-p)-
    (p+p1)*(p-p0)*dTstar_tilde_bound
  a <- (p^2)*(g0-g1*a1)-(p^2)*(p0-p1)*h_a1-
    p*p1*Delta01+p*p0*Delta02-
    (1-a1)*(Delta01*(p1-p*p1-p)-Delta02*(p0-p*p0-p))+
    (p^2)*(p1-p0)*dTstar_tilde_bound
  
  coefs <- cbind(a, b, c, d)
  
  # Getting roots of cubic
  all_roots <- apply(coefs, 1, polyroot)
  real_roots <- Re(unlist(all_roots))
  real_roots <- real_roots[(real_roots >= 0) & (real_roots <= a0upper)]
  
  # Evaluating all combinations of roots and bounds for r_TstarU
  dz_tilde <- outer(real_roots, a1,
                    function(a0,a1) get_dz_tilde(dTstar_tilde_bound,a1,a0,obs))
  
  # Getting max and min if there are real roots
  min_edge <- ifelse(length(real_roots) > 0, min(dz_tilde), NA)
  max_edge <- ifelse(length(real_roots) > 0, max(dz_tilde), NA)
  ans <- list(dz_tilde = list(min_edge = min_edge, max_edge = max_edge),
              a0_roots = real_roots)
  return(ans)
  
}


candidate4 <- function(dTstar_tilde_bound,obs) {
  
  a0upper <- get_alpha_bounds(obs)[1]
  a1upper <- get_alpha_bounds(obs)[2]
  dTstar_tilde_bound <- dxs
  
  g0 <- with(obs,yt01-yt00)
  g1 <- with(obs,(yt01-yt00)+(yt11-yt10))
  h0 <- with(obs,(1-q)*yt00+q*yt01)
  h1 <- with(obs,((1-q)*yt00+q*yt01)+((1-q)*yt10+q*yt11))
  Delta01 <- with(obs,yt10)
  Delta11 <- with(obs,yt10+yt00)
  Delta02 <- with(obs,yt11)
  Delta12 <- with(obs,yt11+yt01)
  p0 <- with(obs,p0)
  p1 <- with(obs,p1)
  p  <- with(obs,p)
  
  gradient_fn <- function(vector_a) {
    
    #   vector_a <- a
    a0 <- vector_a[1]
    a1 <- vector_a[2]
    
    h <- (h0-h1*a1)/(1-p-a1)
    DeltaA0 <- with(obs,((Delta01-Delta11*a0)/(p0-a0))-
                      ((Delta02-Delta12*a0)/(p1-a0)))
    dDeltaA0Term_da0 <- (1/(p-a0)^2)*(p0-a0)*(p1-a0)*DeltaA0+
      (1/(p-a0))*((p1-a0)*(-Delta11)-(Delta01-a0*Delta11)-
                    (p0-a0)*(-Delta12)+(Delta02-a0*Delta12))
    dh_da1 <- (1/(1-p-a1)^2)*(h0-h1*a1)+(1/(1-p-a1))*(-h1)
    dg_da1 <- -g1
    
    dB_da0 <- (1/(1-a0-a1)^2)*(g0-g1*a1)-(1/(1-a0-a1)^2)*(p0-p1)*h-
      (1/(1-a0-a1)^2)*((p0-a0)*(p1-a0)/(p-a0))*DeltaA0-
      (1/(1-a0-a1))*dDeltaA0Term_da0
    
    dS_da0 <- (1/(1-a0-a1)^2)*(p1-p0)
    
    dB_da1 <- (1/(1-a0-a1)^2)*(g0-g1*a1)-(1/(1-a0-a1)^2)*(p0-p1)*h-
      (1/(1-a0-a1)^2)*((p0-a0)*(p1-a0)/(p-a0))*DeltaA0+
      (1/(1-a0-a1))*dg_da1-(1/(1-a0-a1))*(p0-p1)*dh_da1
    
    dS_da1 <- (1/(1-a0-a1)^2)*(p1-p0)
    
    dFUNC_da0 <- dB_da0 + dS_da0*dTstar_tilde_bound
    dFUNC_da1 <- dB_da1 + dS_da1*dTstar_tilde_bound
    
    # get_dz_tilde(dTstar_tilde_bound,a0,a1,obs)
    # B_manual + S_manual*dTstar_tilde_bound
    # B_manual <- (1/(1-a0-a1))*(g0-g1*a1)-(1/(1-a0-a1))*(p0-p1)*h_a1-
    #              (1/(1-a0-a1))*((p0-a0)*(p1-a0)/(p-a0))*DeltaA0
    # S_manual <- (1/(1-a0-a1))*(p1-p0)
    
    #     cond1_c <- g0-(1-a0)*g1-(p0-p1)*h1-((p0-a0)*(p1-a0)/(p-a0))*DeltaA0+
    #               (p1-p0)*dTstar_tilde_bound
    #     cond1_b <- 2*(1-p)*(-g0+2*(1-p)*(1-a0)*g1-(p1-p0)*dTstar_tilde_bound)+
    #                2*(1-p)*((p0-a0)*(p1-a0)/(p-a0))*DeltaA0 +
    #                (p0-p1)*(((1-p)*h1+h0)+(h0-(1-p)*h1))
    #     cond1_a <- (1-(p^2))*(g0-(1-a0)*g1+(p1-p0)*dTstar_tilde_bound)+
    #                (1-(p^2))*((p0-a0)*(p1-a0)/(p-a0))*DeltaA0-
    #                (p0-p1)*((1-p)-(1-a0)*h0-(1-a0)*(1-p)*h1)
    #
    #     cond2_d <- (Delta02-Delta01)-(Delta12-Delta11)
    #     cond2_c <- (g0-g1*a1)-(p0-p1)*h_a1+Delta01*((p-a0)+p1)+
    #                 Delta02*(-(p-a1)-p0)+
    #                 2*Delta11*p1-2*p0*Delta12+(p1-p0)*dTstar_tilde_bound
    #     cond2_b <- -2*p*(g0-g1*a1)-(p0-p1)+2*p*(p0-p1)*h_a1+
    #                 (p*p1*Delta11-(p+p1)*Delta01)-(p*p0*Delta12-(p+p0)*Delta02)-
    #                 (1-a1)*(Delta01*(p+p1)+Delta11*(p-p0)-Delta02*(p+p0)-Delta12*(p-p0))+
    #                 Delta01*(p1-p*p1-p)-Delta02*(p0-p*p0-p)-
    #                 (p+p1)*(p-p0)*dTstar_tilde_bound
    #     cond2_a <- (p^2)*(g0-g1*a1)-(p^2)*(p0-p1)*h_a1-
    #                 p*p1*Delta01+p*p0*Delta02-
    #                 (1-a1)*(Delta01*(p1-p*p1-p)-Delta02*(p0-p*p0-p))+
    #                 (p^2)*(p1-p0)*dTstar_tilde_bound
    #
    #     cond1_coefs <- cbind(cond1_a, cond1_b, cond1_c)
    #     cond2_coefs <- cbind(cond2_a, cond2_b, cond2_c, cond2_d)
    
    #     gradient <- c(cond1_coefs%*%c(1,a1,a1^2),cond2_coefs%*%c(1,a0,a0^2,a0^3))
    
    gradient <- cbind(dFUNC_da0,dFUNC_da1)
    
    return(gradient)
    
  }
  
  a <-  c(a0_grid[coordinates_argmin[1]],a1_grid[coordinates_argmin[2]])
  
  gradient_fn(a)
  gradient_fn(c(0,0))
  
  
  a0_grid <- seq(0,a0upper,length.out = L)
  a1_grid <- seq(0,a1upper,length.out = L)
  
  dzs_matrix <- matrix(NA,100,100)
  
  for(i in 1:100) {
    for(j in 1:100) {
      
      dzs_matrix[i,j] <- gradient_fn(c(a0_grid[i],a1_grid[j]))
      
    }
  }
  
  
  gradient_fn(a)
  
  initial_guess <- expand.grid(seq(0,a0upper,length.out = 10),
                               seq(0,a1upper,length.out = 10))
  
  nleqslv(c(a0upper,a1upper),gradient_fn)
  
  result <- searchZeros(as.matrix(initial_guess),gradient_fn)
  
  result <- searchZeros(t(as.matrix(a)),gradient_fn)
  result$x
  
  b <- c(a0_grid[coordinates_argmax[1]],a1_grid[coordinates_argmax[2]])
  gradient_fn(b)
  
  a0upper
  a1upper
  
  ## image2D(x=a0_grid,y=a1_grid,z=dzs_matrix,xlab="a0",ylab="a1")
  ## COMMAND TO FIND ROOTS
  
}
