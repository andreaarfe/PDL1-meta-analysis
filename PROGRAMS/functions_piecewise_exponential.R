
# Function to compute the survival function associated of a piecewise exponential distribution
surv.piecewise.exp <- function(times,lambda,breaks){
  delta <- c(diff(breaks),Inf)
  surv <- vector(mode="numeric",length = length(times))
  for(i in seq_along(times)){
    tt <- pmin(delta,pmax(0,times[i]-breaks))
    surv[i] <- exp(-sum(tt*lambda))
  }
  return(surv)
}

# Function to sample from a piecewise exponential distribution
sim.piecewise.exp <- function(n,lambda,breaks){
  m <- length(breaks)
  haz <- lambda[-m]*diff(breaks)
  cumhaz <- cumsum(c(0,haz))
  e <- rexp(n,rate=1)
  I <- findInterval(e,cumhaz)
  times  <- breaks[I]+(e-cumhaz[I])/lambda[I]
  return(times)
}
