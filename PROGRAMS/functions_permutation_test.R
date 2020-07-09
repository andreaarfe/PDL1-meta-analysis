
# Functions to implement the design of Arfè et al.
# Manuscript: https://arxiv.org/abs/1902.00161
# Source: https://github.com/andreaarfe/Bayesian-optimal-tests-non-PH

# Prepares the data to compute the marginal likelihood
prep.data <- function(times,events,breaks){
  # event indicators
  ev <- rep(NA,times = length(times))
  I <- findInterval(times[events==1],
                    breaks,all.inside = FALSE,
                    left.open = TRUE, 
                    rightmost.closed = TRUE)
  ev[events==1] <- I
  ev <- factor(ev,levels = 1:length(breaks))
  
  # Person-time per interval
  delta <- c(diff(breaks),Inf) # length of intervals
  m <- t(sapply(times,function(t) pmin(delta,pmax(0,t-breaks))))
  m[is.na(m)] <- 0
  return(list(ev=ev,pt=m))
}

# Computes the log-marginal likelihood under the alternative hypothesis
# Assuming independent piecewise exponential distributions for the two arms
# and independent conjugate Gamma priors for the hazards
log_marg_lik <- function(ev,pt,arms,breaks,alpha0,beta0,alpha1,beta1){
  # Event per interval
  n0 <- table(ev[arms==0],useNA="no")
  n1 <- table(ev[arms==1],useNA="no")
  
  # Person-time per interval
  m0 <- rep(0,times=length(breaks))
  m1 <- rep(0,times=length(breaks))
  if(sum(arms==0)>0) m0 <- colSums(pt[arms==0,])
  if(sum(arms==1)>0) m1 <- colSums(pt[arms==1,])
  m0[is.na(m0)] <- 0
  m1[is.na(m1)] <- 0
  
  # Control arm
  LML.1.0 <- 0 
  if(sum(arms==0)>0){
    LML.1.0 <- sum(alpha0*log(beta0)-(alpha0+n0)*log(beta0+m0) 
                     + lgamma(alpha0+n0) - lgamma(alpha0))  
  } 
  # Experimental arm
  LML.1.1 <- 0 
  if(sum(arms==1)>0){
    LML.1.1 <- sum(alpha1*log(beta1)-(alpha1+n1)*log(beta1+m1)  
                   + lgamma(alpha1+n1) - lgamma(alpha1))       
  }
  # Log-marginal-likelihood under the alternative
  LML.1 <- LML.1.0 + LML.1.1
  #if(is.na(LML.1)) browser()
  return(LML.1)
}

# Implements the optimal permutation test 
# Outputs a list contatining the observed test statistics (tobs),
# the samples from the permutation distribution of the test statistics (tsim),
# the permutation p-value (pval), and the monte carlo standard error of the p-value (mcse)
perm_test <- function(times,events,arms,breaks,alpha0,beta0,alpha1,beta1,nsim=1e3){
  dt <- prep.data(times,events,breaks)
  tobs <- log_marg_lik(dt$ev,dt$pt,arms,breaks,alpha0,beta0,alpha1,beta1)
  n <- length(times)
  tsim <- vector(mode="numeric",length = nsim)
  for(i in 1:nsim){
    asim <- sample(arms,size=n,replace = FALSE)
    tsim[i] <- log_marg_lik(dt$ev,dt$pt,asim,breaks,alpha0,beta0,alpha1,beta1)
  }
  pval <- mean(tsim>=tobs)
  mcse <- sqrt(pval*(1-pval)/nsim)
  out <- list(tobs=tobs,pval=pval,mcse=mcse,tsim=tsim)
  return(out)
}

# Implements the optimal permutation test for the case of binary marker (Section 8 of the paper)
# Outputs a list contatining the observed test statistics (tobs),
# the samples from the permutation distribution of the test statistics (tsim),
# the permutation p-value (pval), and the monte carlo standard error of the p-value (mcse)
perm_test_marker <- function(times,events,arms,markers,breaks,
                             # parameters for z=0
                             alpha00,beta00,alpha10,beta10,
                             # parameters for z=1
                             alpha01,beta01,alpha11,beta11,
                             nsim=1e3){
  I1 <- which(markers == 1)
  I0 <- which(markers == 0)
  dt0 <- prep.data(times[I0],events[I0],breaks)
  dt1 <- prep.data(times[I1],events[I1],breaks)
  tobs <- log_marg_lik(dt0$ev,dt0$pt,arms[I0],breaks,alpha00,beta00,alpha10,beta10) +
    log_marg_lik(dt1$ev,dt1$pt,arms[I1],breaks,alpha01,beta01,alpha11,beta11)
  n0 <- length(I0)
  n1 <- length(I1)
  tsim <- vector(mode="numeric",length = nsim)
  for(i in 1:nsim){
    asim  <- arms
    asim[I0] <- sample(arms[I0],size=n0,replace = FALSE)
    asim[I1] <- sample(arms[I1],size=n1,replace = FALSE)
    tsim[i] <- log_marg_lik(dt0$ev,dt0$pt,asim[I0],breaks,alpha00,beta00,alpha10,beta10) +
      log_marg_lik(dt1$ev,dt1$pt,asim[I1],breaks,alpha01,beta01,alpha11,beta11)
  }
  pval <- mean(tsim>=tobs)
  mcse <- sqrt(pval*(1-pval)/nsim)
  out <- list(tobs=tobs,pval=pval,mcse=mcse,tsim=tsim)
  return(out)
}
