
library(isotone)
library(survival)
library(survRM2)
library(compiler)
enableJIT(3)
set.seed(34783)

# Loads the functions to implement the design of Arfè et al.
# Manuscript: https://arxiv.org/abs/1902.00161
# Source: https://github.com/andreaarfe/Bayesian-optimal-tests-non-PH
source("./PROGRAMS/functions_permutation_test.R")

# Loads data and MCMC output 
load("./DATA/analysis_dset.Rdata")
load("./DATA/fit_pem_type.Rdata")
NITER <- nrow(fit$pi)
NBURN <- NITER/2

#The RMST cut-point is chosen as the 80% quantile of the Kaplan-Meier curve in the chemotherapy arm obtained by aggregating all trials (no distiction between tumors).
tau <- quantile(
  survfit(Surv(time,event)~1, data=data[data$arm==0,]),
  probs = 0.80
)$quantile

# Prevalence of the pd-l1 classes
pi <- fit$pi
pi <- pi[(NBURN+1):NITER,] # discards burn-in
cl <- c("0%-1%", 
        "1%-5%", 
        "5%-10%", 
        "10%-50%", 
        "50%-80%", 
        "80%-100%")
colnames(pi) <- cl

# Delta-RMST (with monotonicity constraint) for NSCLC

# Piecewise exponential breakpoints.
breaks <- c(0,quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5))   

# Function to compute the RMST 
rmst <- function(tau,lambda,breaks){
  delta <- pmax(0,pmin(tau,breaks[2:length(breaks)]) - breaks[1:(length(breaks)-1)])
  delta <- c(delta, pmax(0, tau-breaks[length(breaks)]))
  H <- c(0,cumsum(lambda*delta))
  B <- (1-exp(-lambda*delta))/lambda
  RMST <- sum( exp(-H[1:(length(H)-1)])*B )
  return(RMST)
}

# Computes the delta-RMST
lambda0 <- fit$lambda[,1,,,2]
lambda1 <- fit$lambda[,2,,,2]
lambda0 <- lambda0[(NBURN+1):NITER,,]
lambda1 <- lambda1[(NBURN+1):NITER,,]
R <- dim(lambda0)[1]
C <- dim(lambda0)[2]
r0 <- matrix(0,nrow=R, ncol=C)
r1 <- matrix(0,nrow=R, ncol=C)
for(i in 1:R){
  for(j in 1:C){
    ll0 <- lambda0[i,j,]
    ll1 <- lambda1[i,j,]
    r0[i,j] <- rmst(tau,ll0,breaks)
    r1[i,j] <- rmst(tau,ll1,breaks)
  }
}
diff <- r1 - r0
V <- cov(diff)
Vinv <- solve(V)
ord <- cbind(1:5, 2:6)
monot2 <- t(sapply(1:nrow(diff), 
                   function(i) activeSet(ord, lfSolver, weights = Vinv, y = diff[i,])$x))
colnames(monot2) <- cl

### Sets the simulation parameters
rmst0 <- r0
diff <- monot2
rmst1 <- rmst0 + diff

# solve for the exponential parameters that will be used to generate data
pi.hat <- colMeans(pi)
rhat0 <- colMeans(rmst0)
rhat1 <- rhat0 + colMeans(diff) 
lambda0 <- sapply(1:length(pi.hat), 
                  function(c) uniroot(function(x) (1-exp(-x*tau))/x-rhat0[c], 
                                      interval=c(0.001,1))$root) 
lambda1 <- sapply(1:length(pi.hat), 
                  function(c) uniroot(function(x) (1-exp(-x*tau))/x-rhat1[c], 
                                      interval=c(0.001,1))$root) 

# Function to simulate a trial of given sample size
# The distribution of PD-L1 expression is assumed uniform within each class
sim.trial <- function(N=100, # sample size
                      pi=NULL,
                      lambda0=NULL,
                      lambda1=NULL,
                      randprob=0.5,
                      tcens=30,
                      cuts = c(0, 0.01, 0.05, 0.10, 0.50, 0.80, 1) # PD-L1 cut-points
){
  low  <- cuts[-length(cuts)]
  high <- cuts[-1]
  pdl1.class <- sample.int(length(pi), size=N, replace=TRUE, prob=pi)
  pdl1 <- (high[pdl1.class]-low[pdl1.class])*runif(N) + low[pdl1.class]
  arm  <- rbinom(N,1,randprob)
  I0 <- which(arm==0)
  I1 <- which(arm==1)
  times <- vector(mode="numeric", length=N)
  for(i in seq_along(arm)){
    lambda   <- ifelse(arm[i]==0, lambda0[pdl1.class[i]], lambda1[pdl1.class[i]])
    times[i] <- rexp(1, rate=lambda)
  }
  events <- as.numeric(times<=tcens)
  times <- pmin(times, tcens)
  return(data.frame(arm=arm,
                    time=times,
                    event=events,
                    pdl1=pdl1))
}


# meta-analytic estimates for the classes with PD-L1 <10% and >=10%
p1m <- pi[,1]+pi[,2]+pi[,3]
p1p <- 1-p1m
rmst01m <- rowSums(rmst0[,1:3]*pi[,1:3])/p1m
rmst11m <- rowSums(rmst1[,1:3]*pi[,1:3])/p1m
rmst01p <- (rowSums(rmst0 * pi) - p1m*rmst01m)/p1p
rmst11p <- (rowSums(rmst1 * pi) - p1m*rmst11m)/p1p

# solve for the exponential parameters to fix the prior distribution
lambda01m <- sapply(rmst01m,  
                    function(y) uniroot(function(x) ifelse(x>0, (1-exp(-x*tau))/x, tau)-y, 
                                        interval=c(0,1))$root)
lambda01p <- sapply(rmst01p,  
                    function(y) uniroot(function(x) ifelse(x>0, (1-exp(-x*tau))/x, tau)-y, 
                                        interval=c(0,1))$root)
lambda11m <- sapply(rmst11m,  
                    function(y) uniroot(function(x) ifelse(x>0, (1-exp(-x*tau))/x, tau)-y, 
                                        interval=c(0,1))$root)
lambda11p <- sapply(rmst11p,  
                    function(y) uniroot(function(x) ifelse(x>0, (1-exp(-x*tau))/x, tau)-y, 
                                        interval=c(0,1))$root)


# mean and variance for prior distributions
m01m <- mean(lambda01m)
m01p <- mean(lambda01p)
m11m <- mean(lambda11m)
m11p <- mean(lambda11p)
v01m <- var(lambda01m)
v01p <- var(lambda01p)
v11m <- var(lambda11m)
v11p <- var(lambda11p)

# Derives Gamma prior parameters (from mean/variance to shape/rate)
# mean = a/b, variance = a/b^2, a = shape, b = rate
# b = mean/variance, a = b * mean
b01m <- m01m/v01m
b01p <- m01p/v01p
b11m <- m11m/v11m
b11p <- m11p/v11p
a01m <- b01m * m01m
a01p <- b01p * m01p
a11m <- b11m * m11m
a11p <- b11p * m11p

# Checks the fit of the gamma approximation to the meta-analityic posteriors
# hist(lambda01m,xlim=c(0.05,0.1),freq=FALSE); curve(dgamma(x,shape=a01m, rate=b01m), xlim=c(0.05,0.1), add=TRUE)
# hist(lambda01p,xlim=c(0.05,0.1),freq=FALSE); curve(dgamma(x,shape=a01p, rate=b01p), xlim=c(0.05,0.1), add=TRUE)
# hist(lambda11m,xlim=c(0.05,0.1),freq=FALSE); curve(dgamma(x,shape=a11m, rate=b11m), xlim=c(0.05,0.1), add=TRUE)
# hist(lambda11p,xlim=c(0.04,0.06),freq=FALSE); curve(dgamma(x,shape=a11p, rate=b11p), xlim=c(0.01,0.06), add=TRUE)

# Number of iterations per simulation scenario
NSIM <- 10000

# design: all-comers design
pvals.1 <- matrix(0, nrow=NSIM, ncol=4)
colnames(pvals.1) <- c("naive_some",
                       "meta_some",
                       "naive_all",
                       "meta_all")

###############################################################
### Main analysis: scenario 1
###############################################################

for(iter in 1:NSIM){
  if(iter %% 100 == 0) print(iter)
  d <- sim.trial(N=500,pi.hat,lambda0,lambda1)
  dPOS <- d[d$pdl1>=0.01,]
  dNEG <- d[d$pdl1<0.01,]
  dPOS10 <- d[d$pdl1>=0.1,]
  dNEG10 <- d[d$pdl1<0.1,]
  p.meta.0 <- perm_test(times   = dNEG10$time,
                        events  = dNEG10$event,
                        arms    = dNEG10$arm,
                        breaks  = c(0,+Inf),
                        # parameters for PD-L1 negative (<10%)
                        alpha0 = a01m,
                        beta0  = b01m,
                        alpha1 = a11m,
                        beta1  = b11m)$pval
  p.meta.1 <- perm_test(times   = dPOS10$time,
                        events  = dPOS10$event,
                        arms    = dPOS10$arm,
                        breaks  = c(0,+Inf),
                        # Parameters for PD-L1 positive (>10%)
                        alpha0 = a01p,
                        beta0  = b01p,
                        alpha1 = a11p,
                        beta1  = b11p)$pval
  p1 <- with(dPOS,
             rmst2(time=time,
                   status=event,
                   arm=arm))$unadjusted.result[1,"p"]
  p0 <- with(dNEG,
             rmst2(time=time,
                   status=event,
                   arm=arm))$unadjusted.result[1,"p"]
  p.adj      <- p.adjust(c(p0,p1), method="holm")
  p.meta.adj <- p.adjust(c(p.meta.0,p.meta.1), method="holm")
  pvals.1[iter, "naive_some"] <- min(p.adj)
  #pvals.1[iter, "meta_some"]  <- p.meta_some
  pvals.1[iter, "meta_some"]  <- min(p.meta.adj)
  pvals.1[iter, "naive_all"]  <- max(p.adj)
  pvals.1[iter, "meta_all"]   <- max(p.meta.adj)
}

# Power estimates
sink(file="./RESULTS/meta_designs_power_2.txt")
print(colMeans(pvals.1<=0.05))
sink()

# saves the output
save(file="./DATA/meta_designs_sims_2.Rdata",
     list="pvals.1")

###############################################################
### Main analysis: scenario 2
###############################################################

# Simulations under the null hypothesis, i.e. assuming that
# ICIs have no effect regardless of PD-L1 levels

#for(iter in 1:NSIM){
#  if(iter %% 100 == 0) print(iter )
#  d <- sim.trial(N=500,pi.hat,lambda0,lambda0)
#  dPOS <- d[d$pdl1>=0.01,]
#  dNEG <- d[d$pdl1<0.01,]
#  dPOS10 <- d[d$pdl1>=0.1,]
#  dNEG10 <- d[d$pdl1<0.1,]
#  p.meta.0 <- perm_test(times   = dNEG10$time,
#                        events  = dNEG10$event,
#                        arms    = dNEG10$arm,
#                        breaks  = c(0,+Inf),
#                        # parameters for PD-L1 negative (<10%)
#                        alpha0 = a01m,
#                        beta0  = b01m,
#                        alpha1 = a11m,
#                        beta1  = b11m)$pval
#  p.meta.1 <- perm_test(times   = dPOS10$time,
#                        events  = dPOS10$event,
#                        arms    = dPOS10$arm,
#                        breaks  = c(0,+Inf),
#                        # Parameters for PD-L1 positive (>10%)
#                        alpha0 = a01p,
#                        beta0  = b01p,
#                        alpha1 = a11p,
#                        beta1  = b11p)$pval 
#  p1 <- with(dPOS,
#             rmst2(time=time,
#                   status=event,
#                   arm=arm))$unadjusted.result[1,"p"]
#  p0 <- with(dNEG,
#             rmst2(time=time,
#                   status=event,
#                   arm=arm))$unadjusted.result[1,"p"]
#  p.adj      <- p.adjust(c(p0,p1), method="holm")
#  p.meta.adj <- p.adjust(c(p.meta.0,p.meta.1), method="holm")
#  pvals.1[iter, "naive_some"] <- min(p.adj)
#  #pvals.1[iter, "meta_some"]  <- p.meta_some
#  pvals.1[iter, "meta_some"]  <- min(p.meta.adj)
#  pvals.1[iter, "naive_all"]  <- max(p.adj)
#  pvals.1[iter, "meta_all"]   <- max(p.meta.adj)
#}
#
## Power estimates
#sink(file="./RESULTS/meta_designs_power_2_null.txt")
#print(colMeans(pvals.1<=0.05))
#sink()
#
## saves the output
#save(file="./DATA/meta_designs_sims_2_null.Rdata",
#     list="pvals.1")

###############################################################
# Sensitivity analysis 1: different survival distribution
###############################################################

# Design tailored to meta-analytic results for NSCLC, but
# actual distribution is the one for other tumors (c.f. Fig. 1)

# Computes the parameters of the data-generating distribution 
# for other tumors

# lambda0 <- fit$lambda[,1,,,1]
# lambda1 <- fit$lambda[,2,,,1]
# lambda0 <- lambda0[(NBURN+1):NITER,,]
# lambda1 <- lambda1[(NBURN+1):NITER,,]
# R <- dim(lambda0)[1]
# C <- dim(lambda0)[2]
# r0 <- matrix(0,nrow=R, ncol=C)
# r1 <- matrix(0,nrow=R, ncol=C)
# for(i in 1:R){
#   for(j in 1:C){
#     ll0 <- lambda0[i,j,]
#     ll1 <- lambda1[i,j,]
#     r0[i,j] <- rmst(tau,ll0,breaks)
#     r1[i,j] <- rmst(tau,ll1,breaks)
#   }
# }
# diff <- r1 - r0
# V <- cov(diff)
# Vinv <- solve(V)
# ord <- cbind(1:5, 2:6)
# monot2 <- t(sapply(1:nrow(diff), 
#                    function(i) activeSet(ord, lfSolver, weights = Vinv, y = diff[i,])$x))
# colnames(monot2) <- cl
# 
# ### Sets the simulation parameters
# rmst0 <- r0
# diff <- monot2
# rmst1 <- rmst0 + diff
# 
# # solve for the exponential parameters that will be used to generate data
# pi.hat <- colMeans(pi)
# rhat0 <- colMeans(rmst0)
# rhat1 <- rhat0 + colMeans(diff) 
# lambda0 <- sapply(1:length(pi.hat), 
#                   function(c) uniroot(function(x) (1-exp(-x*tau))/x-rhat0[c], 
#                                       interval=c(0.001,1))$root) 
# lambda1 <- sapply(1:length(pi.hat), 
#                   function(c) uniroot(function(x) (1-exp(-x*tau))/x-rhat1[c], 
#                                       interval=c(0.001,1))$root) 
# 
# for(iter in 1:NSIM){
#   if(iter %% 100 == 0) print(iter)
#   d <- sim.trial(N=500,pi.hat,lambda0,lambda1)
#   dPOS <- d[d$pdl1>=0.01,]
#   dNEG <- d[d$pdl1<0.01,]
#   dPOS10 <- d[d$pdl1>=0.1,]
#   dNEG10 <- d[d$pdl1<0.1,]
#   p.meta.0 <- perm_test(times   = dNEG10$time,
#                         events  = dNEG10$event,
#                         arms    = dNEG10$arm,
#                         breaks  = c(0,+Inf),
#                         # parameters for PD-L1 negative (<10%)
#                         alpha0 = a01m,
#                         beta0  = b01m,
#                         alpha1 = a11m,
#                         beta1  = b11m)$pval
#   p.meta.1 <- perm_test(times   = dPOS10$time,
#                         events  = dPOS10$event,
#                         arms    = dPOS10$arm,
#                         breaks  = c(0,+Inf),
#                         # Parameters for PD-L1 positive (>10%)
#                         alpha0 = a01p,
#                         beta0  = b01p,
#                         alpha1 = a11p,
#                         beta1  = b11p)$pval
#   p1 <- with(dPOS,
#              rmst2(time=time,
#                    status=event,
#                    arm=arm))$unadjusted.result[1,"p"]
#   p0 <- with(dNEG,
#              rmst2(time=time,
#                    status=event,
#                    arm=arm))$unadjusted.result[1,"p"]
#   p.adj      <- p.adjust(c(p0,p1), method="holm")
#   p.meta.adj <- p.adjust(c(p.meta.0,p.meta.1), method="holm")
#   pvals.1[iter, "naive_some"] <- min(p.adj)
#   #pvals.1[iter, "meta_some"]  <- p.meta_some
#   pvals.1[iter, "meta_some"]  <- min(p.meta.adj)
#   pvals.1[iter, "naive_all"]  <- max(p.adj)
#   pvals.1[iter, "meta_all"]   <- max(p.meta.adj)
# }
# 
# # Power estimates
# sink(file="./RESULTS/meta_designs_power_2_sens1.txt")
# print(colMeans(pvals.1<=0.05))
# sink()

###############################################################
# Sensitivity analysis 1: different survival distribution
###############################################################

# Design tailored to meta-analytic that assume prevalence of
# PD-L1 classes does not vary with the type of assay, but 
# the true prevalences are different

# PD-L1 classes are generated assuming the same prevalences 
# estimated for assay 22C3 (Figure 2, panel b)

# # Load the posterior samples for 22C3
# load('./DATA/sens_ICH22C3.Rdata')
# 
# # Data-generating distribution
# pi <- fit.ICH22C3$pi
# pi <- pi[(NBURN+1):NITER,] # discards burn-in
# cl <- c("0%-1%", 
#         "1%-5%", 
#         "5%-10%", 
#         "10%-50%", 
#         "50%-80%", 
#         "80%-100%")
# colnames(pi) <- cl
# pi.hat <- colMeans(pi)
# 
# for(iter in 1:NSIM){
#   if(iter %% 100 == 0) print(iter)
#   d <- sim.trial(N=500,pi.hat,lambda0,lambda1)
#   dPOS <- d[d$pdl1>=0.01,]
#   dNEG <- d[d$pdl1<0.01,]
#   dPOS10 <- d[d$pdl1>=0.1,]
#   dNEG10 <- d[d$pdl1<0.1,]
#   p.meta.0 <- perm_test(times   = dNEG10$time,
#                         events  = dNEG10$event,
#                         arms    = dNEG10$arm,
#                         breaks  = c(0,+Inf),
#                         # parameters for PD-L1 negative (<10%)
#                         alpha0 = a01m,
#                         beta0  = b01m,
#                         alpha1 = a11m,
#                         beta1  = b11m)$pval
#   p.meta.1 <- perm_test(times   = dPOS10$time,
#                         events  = dPOS10$event,
#                         arms    = dPOS10$arm,
#                         breaks  = c(0,+Inf),
#                         # Parameters for PD-L1 positive (>10%)
#                         alpha0 = a01p,
#                         beta0  = b01p,
#                         alpha1 = a11p,
#                         beta1  = b11p)$pval
#   p1 <- with(dPOS,
#              rmst2(time=time,
#                    status=event,
#                    arm=arm))$unadjusted.result[1,"p"]
#   p0 <- with(dNEG,
#              rmst2(time=time,
#                    status=event,
#                    arm=arm))$unadjusted.result[1,"p"]
#   p.adj      <- p.adjust(c(p0,p1), method="holm")
#   p.meta.adj <- p.adjust(c(p.meta.0,p.meta.1), method="holm")
#   pvals.1[iter, "naive_some"] <- min(p.adj)
#   #pvals.1[iter, "meta_some"]  <- p.meta_some
#   pvals.1[iter, "meta_some"]  <- min(p.meta.adj)
#   pvals.1[iter, "naive_all"]  <- max(p.adj)
#   pvals.1[iter, "meta_all"]   <- max(p.meta.adj)
# }
# 
# # Power estimates
# sink(file="./RESULTS/meta_designs_power_2_sens2.txt")
# print(colMeans(pvals.1<=0.05))
# sink()

