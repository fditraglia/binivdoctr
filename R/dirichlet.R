rdirichlet <- function(n, shape){
  K <- length(shape)
  gamma_draws <- matrix(rgamma(n * K, shape, scale = 1), K)
  out <- t(t(gamma_draws) / colSums(gamma_draws))
  return(out[-K,])
}

ddirichlet <- function(x, shape, log.ret = FALSE){
  K <- length(shape)
  stopifnot(is.numeric(x))
  if(is.vector(x)){
    if(length(x) == K){
      stopifnot(identical(sum(x), 1))
      x <- matrix(x, ncol = length(x))
    }else{
      stopifnot(length(x) == (K - 1))
      stopifnot(sum(x) < 1)
      x <- c(x, 1 - sum(x))
      x <- matrix(x, ncol = K)
    }
  }else if(is.matrix(x)){
    if(ncol(x) == K){
      stopifnot(identical(rowSums(x), 1))
    }else{
      stopifnot(ncol(x) == (K - 1))
      stopifnot(all(rowSums(x) < 1))
      x <- cbind(x, 1 - rowSums(x))
    }
  }else{
    stop("x must be a vector or matrix")
  }
  logD <- drop(log(x) %*% (shape - 1) + lgamma(sum(shape)) -  sum(lgamma(shape)))
  if(log.ret){
    return(logD)
  }else{
    return(exp(logD))
  }
}

image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

dirichletPlot <- function(scale, inc = 0.005){
  stopifnot(length(scale) == 3)
  a0 <- a1 <- seq(0.01, 0.99, inc)
  a <- expand.grid(a0 = a0, a1 = a1)
  aRestrict <- subset(a, a0 + a1 < 1)
  densRestrict <- ddirichlet(as.matrix(aRestrict), scale)
  aRestrict <- cbind(aRestrict, densRestrict)
  names(aRestrict) <- c("a0", "a1", "dens")
  aViolate <- subset(a, a0 + a1 >= 1)
  densViolate <- rep(NA, nrow(aViolate))
  aViolate <- cbind(aViolate, densViolate)
  names(aViolate) <- c("a0", "a1", "dens")
  a <- rbind(aRestrict, aViolate)
  a <- a[order(a$a1, a$a0),]
  dens<- matrix(a$dens, nrow = length(a0), ncol = length(a1), byrow = FALSE)

  layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths=c(4,1), heights=c(4,4))
  pal <- colorRampPalette(c("black", "red", "yellow"), space="rgb")
  mybreaks <- seq(min(dens, na.rm = TRUE), max(dens, na.rm = TRUE), length.out=100)
  image(a0, a1, dens, col = pal(length(mybreaks) - 1), breaks = mybreaks)
  image.scale(dens, col=pal(length(mybreaks)-1), breaks = mybreaks, horiz = FALSE,
              xlab = "", ylab = "")
}

dirichletInteractive <- function(){
  manipulate::manipulate(dirichletPlot(c(s1, s2, s3)),
             s1 = slider(1, 10, 1, step = 0.25),
             s2 = slider(1, 10, 1, step = 0.25),
             s3 = slider(1, 10, 2, step = 0.25))
}
