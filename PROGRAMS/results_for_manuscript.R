
library(tidyr)
library(dplyr)
library(ggplot2)
library(survival)
library(gridExtra)
library(isotone)
source("./PROGRAMS/functions_piecewise_exponential.R")

# Output file with the results to be inserted as text in the manuscript 
outfile <- "./RESULTS/results_for_manuscript.txt"
cat("\n", file=outfile)

# Data
load("./DATA/analysis_dset.Rdata")
NSCLC <- as.numeric(data$trial %in% c("Checkmate 017",
                                      "Checkmate 057",
                                      #"Keynote 10",
                                      "POPLAR",
                                      "OAK",
                                      #"Checkmate 227", excluded because combination trial
                                      "Keynote 189",
                                      "JAVELIN Lung 200"))
data$tumor_type <- NSCLC

# Loads MCMC output and discards the burn-in iterations
load("./DATA/fit_pem_type.Rdata")
NITER <- nrow(fit$pi)
NBURN <- NITER/2
pi <- fit$pi
pi <- pi[(NBURN+1):NITER,] # discards burn-in
cl <- c("0%-1%", 
        "1%-5%", 
        "5%-10%", 
        "10%-50%", 
        "50%-80%", 
        "80%-100%")
colnames(pi) <- cl
pi.mat <- pi
pi <- as.data.frame(pi) %>% gather(PDL1,prob)
pi$PDL1 <- factor(pi$PDL1, levels=cl)

#The RMST cut-point is chosen as the 80% quantile of the Kaplan-Meier curve in the chemotherapy arm obtained by aggregating all trials (no distiction between tumors).
tau <- quantile(
  survfit(Surv(time,event)~1, data=data[data$arm==0,]),
  probs = 0.80
)$quantile

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

# Other tumors
lambda0 <- fit$lambda[,1,,,1]
lambda1 <- fit$lambda[,2,,,1]
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
monot <- t(sapply(1:nrow(diff), 
                  function(i) activeSet(ord, lfSolver, weights = Vinv, y = diff[i,])$x))
colnames(monot) <- cl
monot_summary <- as.data.frame(monot) %>% 
  gather(PDL1,delta) %>% 
  group_by(PDL1) %>%
  summarize(m=mean(delta),
            lcl=quantile(delta, probs = 0.025),
            ucl=quantile(delta, probs = 0.975)) 

# NSCLC
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
monot2_summary <- as.data.frame(monot2) %>% 
  gather(PDL1,delta) %>% 
  group_by(PDL1) %>%
  summarize(m=mean(delta),
            lcl=quantile(delta, probs = 0.025),
            ucl=quantile(delta, probs = 0.975)) 

#################################################################
# PD-L1 prevalences
#################################################################

# Prevalence of PD-L1 <=5%
p1 <- rowSums(pi.mat[,c("0%-1%","1%-5%")])
res1 <- c(mean(p1),quantile(p1, probs=c(0.025, 0.975)))
names <- c("Posterior mean", "95% Lower PI limit", "95% Upper PI limit")
cat("Prevalence of PD-L1 <=5%",
    paste(names, collapse="   "),
    paste(round(res1,2), collapse="   "),
    "",
    sep="\n",
    file=outfile,
    append=TRUE)

# Prevalence of 5%<PD-L1 <=50%
p2 <- rowSums(pi.mat[,c("5%-10%","10%-50%")])
res2 <- c(mean(p2),quantile(p2, probs=c(0.025, 0.975)))
names <- c("Posterior mean", "95% Lower PI limit", "95% Upper PI limit")
cat("Prevalence of 5%<PD-L1 <=50%",
    paste(names, collapse="   "),
    paste(round(res2,2), collapse="   "),
    "",
    sep="\n",
    append=TRUE,
    file=outfile)

# Prevalence of PD-L1 >50%
p3 <- 1-p1-p2
res3 <- c(mean(p3),quantile(p3, probs=c(0.025, 0.975)))
names <- c("Posterior mean", "95% Lower PI limit", "95% Upper PI limit")
cat("Prevalence of PD-L1 >50%",
    paste(names, collapse="   "),
    paste(round(res3,2), collapse="   "),
    "",
    sep="\n",
    append=TRUE,
    file=outfile)


#################################################################
# Summary statistics
#################################################################

# Total number of events and person-time at risk
res <- c(sum(data$event), sum(data$time)) 
names(res) <- c("Number of events", "Person-months of follow-up")
cat("Summary statistics",
    paste(names(res), collapse="   "),
    paste(round(res,2), collapse="   "),
    "",
    sep="\n",
    append=TRUE,
    file=outfile)

# Quintiles of event times (piecewise exponential cut-points)
# Piecewise exponential breakpoints.
breaks <- quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5)   
cat("Quintiles",
    paste(names(breaks), collapse="   "),
    paste(round(breaks,2), collapse="   "),
    "",
    sep="\n",
    append=TRUE,
    file=outfile)

#################################################################
# TREATMENT EFFECTS
#################################################################

# Summary of treatment effects

cat("Treatment effects, NSCLC",
    paste(c("Value",monot2_summary$PDL1), collapse="   "),
    paste(c("Estimate",round(monot2_summary$m,2)), collapse="   "),
    paste(c("LPL95",round(monot2_summary$lcl,2)), collapse="   "),
    paste(c("UPL95",round(monot2_summary$ucl,2)), collapse="   "),
    "",
    sep="\n",
    append=TRUE,
    file=outfile)

cat("Treatment effects, Other tumors",
    paste(c("Value",monot_summary$PDL1), collapse="   "),
    paste(c("Estimate",round(monot_summary$m,2)), collapse="   "),
    paste(c("LPL95",round(monot_summary$lcl,2)), collapse="   "),
    paste(c("UPL95",round(monot_summary$ucl,2)), collapse="   "),
    "",
    sep="\n",
    append=TRUE,
    file=outfile)






