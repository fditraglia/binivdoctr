library('plot3D')
library('numDeriv')
library('rootSolve')
library('nleqslv')

setwd('C:/Users/mmm/Documents/GITHUB_FOLDERS/binivdoctr/R')

source('binCDFs.R')
source('dirichlet.R')
source('identifiedSet.R')
source('IshwaranJames.R')
source('plotFunctions.R')
source('SamplingUncertainty.R')
source('simDGPs.R')
source('summary.R')
source('rho_bounds.R')

get_beta_bounds <- function(dTstar_tilde_lower,dTstar_tilde_upper,
                                a0_upper, a1_upper,
                                obs) {
  
  # Obtain bounds for a0 and a1 implied by observables,
  # without imposing user beliefs.
  
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
  print('Beta bounds: Evaluating corners')
  set1a <- candidate1(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set1b <- candidate1(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)
  
  print('Beta bounds: Evaluating edges')
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
  
  print('Beta bounds: Evaluation interior solution')
  # Candidate Set IV - Interior for both a0 and a1
  set4a <- candidate4(dTstar_tilde_lower, a0_upper_bound, a1_upper_bound, obs)
  set4b <- candidate4(dTstar_tilde_upper, a0_upper_bound, a1_upper_bound, obs)
  
  # Finally: overall max and min
  beta_min <- pmin(set1a$min_corner, set1b$max_corner,
                       set2a$min_edge,set2b$max_edge,
                       set3a$min_edge,set3b$max_edge,
                       set4a$min_int,set4b$max_int,
                       na.rm = TRUE)
  
  beta_max <- pmax(set1a$max_corner, set1b$max_corner,
                       set2a$max_edge,set2b$max_edge,
                       set3a$max_edge,set3b$max_edge,
                       set4a$max_int,set4b$max_int,
                       na.rm = TRUE)
  
  data.frame(beta_min = beta_min, beta_max = beta_max)
  
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
                       a0_upper_bound, a1_upper_bound,
                       beta_iv,obs) {
  
  # Add bounds to thhe obs object.
  # dTstar_tilde_bound <- dxs_vec[44]
  obs$a0_upper_bound <- a0_upper_bound
  obs$a1_upper_bound <- a1_upper_bound
  obs$beta_iv <- 
  
  # obs_df <- obs
  # obs <- obs_df
  
  solve_foc <- function(obs) {
    
    #obs <- as.list(as.data.frame(obs_df)[44,])
    obs <- as.list(obs)
    
    gradient_fn <- function(a) {
      
      a0 <- a[1]
      a1 <- a[2]
      
      foc_beta_a0 <- with(obs,(1-a0-a1)*
                        foc_dztilde_a0(dTstar_tilde_bound,a0,a1,obs)-
                        (beta_iv-
                        (1/(p1-p0))*get_dz_tilde(dTstar_tilde_bound,a0,a1,obs))
      
      foc_beta_a1 <- with(obs,(1-a0-a1)*
                        foc_dztilde_a1(dTstar_tilde_bound,a0,a1,obs)-
                        (beta_iv-
                        (1/(p1-p0))*get_dz_tilde(dTstar_tilde_bound,a0,a1,obs))
      
      gradient_beta <- cbind(foc_beta_a0,foc_beta_a1)

      return(gradient_beta)
      
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
                       function(a) get_beta(dTstar_tilde_bound,a[1],a[2],obs,
                                            beta_iv))
      } else {
        value <- apply(as.matrix(valid_roots),1,
                       function(a) get_beta(dTstar_tilde_bound,a[1],a[2],obs,
                                            beta_iv))
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

foc_beta_a0 <- function(dTstar_tilde,a0,a1,obs) {
  
  with((1-a0-a1)*(b_iv)
  
}