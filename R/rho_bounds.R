
get_dz_tilde_bounds <- function(dTstar_tilde_lower,dTstar_tilde_upper,
                                a0_upper, a1_upper,
                                obs) {
  
  # Obtain bounds for a0 and a1 implied by observables,
  # without imposing user beliefs.
  
  # a0_upper <- 1
  # a1_upper <- 1
  # dTstar_tilde_lower <- dxs_vec[44]
  # dTstar_tilde_upper <- dxs_vec[44]
  
  vectorized_alpha_bounds <- function(obs) {
    get_alpha_bounds(as.list(obs))
  }
  
  unrestricted <- t(apply(as.data.frame(obs),1,vectorized_alpha_bounds))

  # Storing the more restrictive bounds for each parameter
  # (either user or unrestricted)
  a0_upper_bound <- ifelse(a0_upper >= unrestricted[,'a0upper'],
                           unrestricted[,'a0upper'],
                           a0_upper)
  a1_upper_bound <- ifelse(a1_upper >= unrestricted[,'a1upper'],
                           unrestricted[,'a1upper'],
                           a1_upper)

  # Finding dz_tilde bounds based on restricted bounds above
  # Candidate Set I - Corner for both a0 and a1
  print('dz_tilde bounds: Evaluating corners')
  set1a <- candidate1(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set1b <- candidate1(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)

  print('dz_tilde bounds: Evaluating edges')
  # Candidate Set II - Corner for a0 and interior for a1
  set2a <- candidate2(dTstar_tilde_lower, rep(0,length(a0_upper_bound)),
                      a0_upper_bound,a1_upper_bound, obs)
  set2b <- candidate2(dTstar_tilde_upper, a0_upper_bound,
                      a0_upper_bound,a1_upper_bound, obs)
  
  # Candidate Set III - Interior for a0, corner for a1
  set3a <- candidate3(dTstar_tilde_lower, rep(0,length(a1_upper_bound)),
                      a0_upper_bound,a1_upper_bound, obs)
  set3b <- candidate3(dTstar_tilde_upper, a1_upper_bound,
                      a0_upper_bound,a1_upper_bound, obs)
  
  print('dz_tilde bounds: Evaluation interior solution')
  # Candidate Set IV - Interior for both a0 and a1
  set4a <- candidate4(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set4b <- candidate4(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)
  
  # Finally: overall max and min
  dz_tilde_min <- pmin(set1a$min_corner, set1b$max_corner,
                       set2a$min_edge,set2b$max_edge,
                       set3a$min_edge,set3b$max_edge,
                       set4a$min_int,set4b$max_int,
                       na.rm = TRUE)
  
  dz_tilde_max <- pmax(set1a$max_corner, set1b$max_corner,
                   set2a$max_edge,set2b$max_edge,
                   set3a$max_edge,set3b$max_edge,
                   set4a$max_int,set4b$max_int,
                   na.rm = TRUE)

  data.frame(dz_tilde_min = dz_tilde_min, dz_tilde_max = dz_tilde_max)
  
}

# Evaluates the corners given user bounds. Vectorized wrt a0 and a1 bounds.
candidate1 <- function(dTstar_tilde_bound, a0_upper, a1_upper, obs) {
  
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

candidate2 <- function(dTstar_tilde_bound, a0,
                       a0_upper_bound,a1_upper_bound, obs) {
  
  # dTstar_tilde_bound <- dxs[44]
  obs$a0_upper_bound <- a0_upper_bound
  obs$a1_upper_bound <- a1_upper_bound
  obs$a0_edge <- a0

  # obs_df <- obs
  # obs <- obs_df
  
  solve_foc <- function(obs) {
    
    #obs <- as.list(as.data.frame(obs_df)[44,])
    obs <- as.list(obs)
    
    gradient_fn <- function(a1) {
      
      gradient_dztilde <- foc_dztilde_a1(dTstar_tilde_bound,
                                         with(obs,a0_edge),a1,obs)
      return(gradient_dztilde)
      
    }

    a1_roots <- uniroot.all(gradient_fn,interval=c(0,with(obs,a1_upper_bound)))

    valid_roots <- matrix(NA,length(a1_roots),2)
    valid_roots[,1] <- with(obs,a0_edge)
    valid_roots[,2] <- a1_roots
    
  
    if(length(valid_roots) != 0) {
      
      value <- apply(valid_roots,1,
                     function(a) get_dz_tilde(dTstar_tilde_bound,
                                              a[1],a[2],obs))

      min_int <- min(t(value), na.rm = TRUE)
      max_int <- max(t(value), na.rm = TRUE)
      
    } else {
      
      min_int <- NA
      max_int <- NA
      
    }
    
    return(c(min_int,max_int)) # rbind(obs$p,obs$q)
    
  }
  
  results <- t(apply(as.data.frame(obs),1,solve_foc))
  
  return(list(min_edge=as.vector(results[,1]),
              max_edge=as.vector(results[,2])))
  
}

candidate3 <- function(dTstar_tilde_bound, a1,
                       a0_upper_bound,a1_upper_bound, obs) {
  
  # dTstar_tilde_bound <- dxs[44]
  obs$a0_upper_bound <- a0_upper_bound
  obs$a1_upper_bound <- a1_upper_bound
  obs$a1_edge <- a1
  
  # obs_df <- obs
  # obs <- obs_df
  
  solve_foc <- function(obs) {
    
    #obs <- as.list(as.data.frame(obs_df)[44,])
    obs <- as.list(obs)
    
    gradient_fn <- function(a0) {
      
      gradient_dztilde <- foc_dztilde_a0(dTstar_tilde_bound,a0,
                                         with(obs,a1_edge),obs)
      return(gradient_dztilde)
      
    }
    
    a0_roots <- uniroot.all(gradient_fn,interval=c(0,with(obs,a0_upper_bound)))
    
    valid_roots <- matrix(NA,length(a0_roots),2)
    valid_roots[,1] <- a0_roots
    valid_roots[,2] <- with(obs,a1_edge)
    
    
    if(length(valid_roots) != 0) {
      
      value <- apply(valid_roots,1,
                     function(a) get_dz_tilde(dTstar_tilde_bound,
                                              a[1],a[2],obs))
      
      min_int <- min(t(value), na.rm = TRUE)
      max_int <- max(t(value), na.rm = TRUE)
      
    } else {
      
      min_int <- NA
      max_int <- NA
      
    }
    
    return(c(min_int,max_int)) # rbind(obs$p,obs$q)
    
  }
  
  results <- t(apply(as.data.frame(obs),1,solve_foc))
  
  return(list(min_edge=as.vector(results[,1]),
              max_edge=as.vector(results[,2])))
  
}

candidate4 <- function(dTstar_tilde_bound,
                       a0_upper_bound, a1_upper_bound, obs) {

  # Add bounds to thhe obs object.
  # dTstar_tilde_bound <- dxs_vec[44]
  obs$a0_upper_bound <- a0_upper_bound
  obs$a1_upper_bound <- a1_upper_bound
  
  # obs_df <- obs
  # obs <- obs_df
  
  solve_foc <- function(obs) {

    #obs <- as.list(as.data.frame(obs_df)[44,])
    obs <- as.list(obs)
        
    gradient_fn <- function(a) {
      
      a0 <- a[1]
      a1 <- a[2]
      
      gradient_dztilde <- cbind(foc_dztilde_a0(dTstar_tilde_bound,a0,a1,obs),
                                foc_dztilde_a1(dTstar_tilde_bound,a0,a1,obs))
      return(gradient_dztilde)
      
    }
    
    # gradient_fn(a)
    
    # Use a grid of initial guesses in the identified set (currently only one)
    initial_guess <- expand.grid(seq(with(obs,a0_upper_bound)*(1/2),
                                     with(obs,a0_upper_bound)*(1/2),
                                     length.out = 1),
                                 seq(with(obs,a1_upper_bound)*(1/2),
                                     with(obs,a1_upper_bound)*(1/2),
                                     length.out = 2))

    result <- searchZeros(as.matrix(initial_guess),gradient_fn)
    valid_roots <- result$x[(result$x[,1]>=0) &
                            (result$x[,1]<=with(obs,a0_upper_bound)) &
                            (result$x[,2]>=0) &
                            (result$x[,2]<=with(obs,a1_upper_bound)),]
    
    if(length(valid_roots) != 0) {
      
      if(length(valid_roots == 2)) {
        value <- apply(t(as.matrix(valid_roots)),1,
                       function(a) get_dz_tilde(dTstar_tilde_bound,
                                                a[1],a[2],obs))
      } else {
        value <- apply(as.matrix(valid_roots),1,
                       function(a) get_dz_tilde(dTstar_tilde_bound,
                                                a[1],a[2],obs))
      }
      

      min_int <- min(t(value), na.rm = TRUE)
      max_int <- max(t(value), na.rm = TRUE)
      
    } else {

      min_int <- NA
      max_int <- NA
            
    }
    
    return(c(min_int,max_int)) # rbind(obs$p,obs$q)
    
  }
  
  # ptm <- proc.time()
  results <- t(apply(as.data.frame(obs),1,solve_foc))
  # proc.time()-ptm
  
  return(list(min_int=as.vector(results[,1]),
              max_int=as.vector(results[,2])))
  
}

foc_dztilde_a0 <- function(dTstar_tilde,a0,a1,obs) {

  g1 <- with(obs,yt01-yt00)
  g2 <- with(obs,(yt01-yt00)+(yt11-yt10))
  h1 <- with(obs,(1-q)*yt00+q*yt01)
  h2 <- with(obs,((1-q)*yt00+q*yt01)+((1-q)*yt10+q*yt11))
  Delta01 <- with(obs,yt10)
  Delta11 <- with(obs,yt10+yt00)
  Delta02 <- with(obs,yt11)
  Delta12 <- with(obs,yt11+yt01)

  g <- g1-g2*a1
  h <- with(obs,(h1-h2*a1)/(1-p-a1))
  D <- with(obs,((Delta01-Delta11*a0)/(p0-a0))-
                      ((Delta02-Delta12*a0)/(p1-a0)))
  dD_Term_da0 <- with(obs,(1/(p-a0)^2)*(p0-a0)*(p1-a0)*D+
                          (1/(p-a0))*((p1-a0)*(-Delta11)-(Delta01-a0*Delta11)-
                          (p0-a0)*(-Delta12)+(Delta02-a0*Delta12)))

  dB_da0 <- with(obs,(1/(1-a0-a1)^2)*g-(1/(1-a0-a1)^2)*(p0-p1)*h-
                     (1/(1-a0-a1)^2)*((p0-a0)*(p1-a0)/(p-a0))*D-
                     (1/(1-a0-a1))*dD_Term_da0)

  dS_da0 <- with(obs,(1/(1-a0-a1)^2)*(p1-p0))

  if(with(obs,C1==0 & C2==0 & C3==0 & C4==0)) {
    foc_a0 <- with(obs,dB_da0+dS_da0*dTstar_tilde)
  } else {

    S <- with(obs, (p1 - p0) / (1 - a0 - a1))
    B <- with(obs, (g - (p0 - p1) * h) / (1 - a0 - a1) -
                ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (1 - a0 - a1)))
    Ff <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))

    dStilde_dS  <- with(obs,(1/(Ff*C4*S-C2+1)^2)*(-C2+1))
    dStilde_dFf <- with(obs,(1/(Ff*C4*S-C2+1)^2)*(-C4*(S^2)))
    dBtilde_dS  <- with(obs,(1/(Ff*C4*S-C2+1)^2)*((Ff*C4*S-C2+1)*(Ff*C3)-
                                         (S*Ff*C3+B-C1)*(Ff*C4)))
    dBtilde_dFf <- with(obs,(1/(Ff*C4*S-C2+1)^2)*((Ff*C4*S-C2+1)*(S*C3)-
                                         (S*Ff*C3+B-C1)*(S*C4)))
    dBtilde_dB <- with(obs,(1/(Ff*C4*S-C2+1)))
    dFf_da0     <- with(obs,(1/(p-a0)^2)*((1-a0-a1)/(1-p-a1))-
                   (1/(p-a0))*(1/(1-p-a1)))

    dStilde_da0 <- (dStilde_dS*dS_da0)+(dStilde_dFf*dFf_da0)
    dBtilde_da0 <- (dBtilde_dS*dS_da0)+(dBtilde_dFf*dFf_da0)+(dBtilde_dB*dB_da0)

    foc_a0 <- dBtilde_da0+dStilde_da0*dTstar_tilde

  }

  return(foc_a0)

}

foc_dztilde_a1 <- function(dTstar_tilde,a0,a1,obs) {

  g1 <- with(obs,yt01-yt00)
  g2 <- with(obs,(yt01-yt00)+(yt11-yt10))
  h1 <- with(obs,(1-q)*yt00+q*yt01)
  h2 <- with(obs,((1-q)*yt00+q*yt01)+((1-q)*yt10+q*yt11))
  Delta01 <- with(obs,yt10)
  Delta11 <- with(obs,yt10+yt00)
  Delta02 <- with(obs,yt11)
  Delta12 <- with(obs,yt11+yt01)

  g <- g1-g2*a1
  h <- with(obs,(h1-h2*a1)/(1-p-a1))
  D <- with(obs,((Delta01-Delta11*a0)/(p0-a0))-
              ((Delta02-Delta12*a0)/(p1-a0)))
  dD_Term_da0 <- with(obs,(1/(p-a0)^2)*(p0-a0)*(p1-a0)*D+
                          (1/(p-a0))*((p1-a0)*(-Delta11)-(Delta01-a0*Delta11)-
                          (p0-a0)*(-Delta12)+(Delta02-a0*Delta12)))
  dh_da1 <- with(obs,(1/(1-p-a1)^2)*(h1-h2*a1)+(1/(1-p-a1))*(-h2))
  dg_da1 <- -g2

  dB_da1 <- with(obs,(1/(1-a0-a1)^2)*g-(1/(1-a0-a1)^2)*(p0-p1)*h-
                     (1/(1-a0-a1)^2)*((p0-a0)*(p1-a0)/(p-a0))*D+
                     (1/(1-a0-a1))*dg_da1-(1/(1-a0-a1))*(p0-p1)*dh_da1)

  dS_da1 <- with(obs,(1/(1-a0-a1)^2)*(p1-p0))

  if(with(obs,C1==0 & C2==0 & C3==0 & C4==0)) {
    foc_a1 <- dB_da1+dS_da1*dTstar_tilde
  } else {
    
    S <- with(obs, (p1 - p0) / (1 - a0 - a1))
    B <- with(obs, (g - (p0 - p1) * h) / (1 - a0 - a1) -
                   ((p0 - a0) * (p1 - a0) * D) / ((p - a0) * (1 - a0 - a1)))
    Ff <- with(obs, (1 - a0 - a1) / ((p - a0) * (1 - p - a1)))
    
    dStilde_dS  <- with(obs,(1/(Ff*C4*S-C2+1)^2)*(-C2+1))
    dStilde_dFf <- with(obs,(1/(Ff*C4*S-C2+1)^2)*(-C4*(S^2)))
    dBtilde_dS  <- with(obs,(1/(Ff*C4*S-C2+1)^2)*((Ff*C4*S-C2+1)*(Ff*C3)-
                                                    (S*Ff*C3+B-C1)*(Ff*C4)))
    dBtilde_dFf <- with(obs,(1/(Ff*C4*S-C2+1)^2)*((Ff*C4*S-C2+1)*(S*C3)-
                                                    (S*Ff*C3+B-C1)*(S*C4)))
    dBtilde_dB <- with(obs,(1/(Ff*C4*S-C2+1)))
    dFf_da1     <- with(obs,(1/(1-p-a1)^2)*((1-a0-a1)/(p-a0))-
                            (1/(1-p-a1))*(1/(p-a0)))
    
    dStilde_da1 <- (dStilde_dS*dS_da1)+(dStilde_dFf*dFf_da1)
    dBtilde_da1 <- (dBtilde_dS*dS_da1)+(dBtilde_dFf*dFf_da1)+(dBtilde_dB*dB_da1)
    
    foc_a1 <- dBtilde_da1+dStilde_da1*dTstar_tilde
    
  }

  return(foc_a1)

}

# a <- c(a0_grid[coordinates_argmin[1]],a1_grid[coordinates_argmin[2]])
# b <- c(a0_grid[coordinates_argmax[1]],a1_grid[coordinates_argmax[2]])
# 
# foc_dztilde_a0(dxs,a[1],a[2],obs)[44]
# foc_dztilde_a0(dxs,b[1],b[2],obs)[44]
# foc_dztilde_a1(dxs,a[1],a[2],obs)[44]
# foc_dztilde_a1(dxs,b[1],b[2],obs)[44]
