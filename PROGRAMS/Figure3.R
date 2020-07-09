
library(tidyr)
library(dplyr)
library(ggplot2)
library(survival)
library(gridExtra)
library(isotone)
library(survRM2)
source("./PROGRAMS/functions_collapsed_gibbs_cpp_pem_types.R")
source("./PROGRAMS/functions_piecewise_exponential.R")
set.seed(45265)

# Plot of the distribution of PD-L1 and the class-specific treatment effects

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

# Prior specification
breaks <- c(0,quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5))
cuts   <- sort(unique(c(data$low,data$high)))
alpha  <- diff(cuts)/min(diff(cuts))
shape  <- matrix(10*log(2)/12,nrow=length(cuts)-1,ncol=length(breaks))
rate   <- matrix(10,nrow=length(cuts)-1,ncol=length(breaks))

# Censors at 0 all event times
data_cens <- data
data_cens$time <- 0
data_cens$event <- 0

# Fit the model
fit <- mix.exp.collapsed.cpp.pem.type(data    = data_cens, 
                                      cuts    = cuts,
                                      breaks  = breaks,
                                      alpha   = alpha,
                                      shape   = shape,
                                      rate    = rate,
                                      verbose = TRUE, 
                                      NITER   = 20000)

# Loads MCMC output and discards the burn-in iterations
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

### Estimation via weighted KM curves

#Computes the posterior estimate of the class membership probabilities for each subject.

# Function to estimate the RMST from a Kaplan-Meyer curve from a weighted 
# Kaplan-Meier curve.
# This is a modification of the "rmst1" function in the survRM2 package.
# survRM2::rmst1 is not visible, to find it run the command getAnywhere(rmst1).
# From the vignette of survRM2, the standard error of the RMST is computed
# using formulas provided in Miller, R. G. (1981). Survival Analysis. Wiley.
rmst.weights <- function(time, status, weights, tau) 
{
  ft = survfit(Surv(time, status) ~ 1, weights = weights)
  idx = ft$time <= tau
  wk.time = sort(c(ft$time[idx], tau))
  wk.surv = ft$surv[idx]
  wk.n.risk = ft$n.risk[idx]
  wk.n.event = ft$n.event[idx]
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  wk.var <- ifelse((wk.n.risk - wk.n.event) == 0, 0, wk.n.event/(wk.n.risk * 
                                                                   (wk.n.risk - wk.n.event)))
  wk.var = c(wk.var, 0)
  rmst.var = sum(cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se = sqrt(rmst.var)
  out <- list(est=rmst, se=rmst.se)
  return(out)
}

# Membership matrix
cuts <- sort(unique(c(data$low,data$high)))
X <- matrix(0,nrow=nrow(data),ncol=length(cuts)-1)
for(i in 1:nrow(X)){
  for(j in 1:ncol(X)){
    if(data$low[i]<=cuts[j] & cuts[j+1]<=data$high[i]) X[i,j] <- 1
  }
}

# Assigns subject-specific weights
P <- 0*X
for(i in 1:nrow(X)){
  temp <- t(X[i,]*t(pi.mat))
  temp <- temp/rowSums(temp)
  P[i,] <- colMeans(temp)
}

#The RMST cut-point is chosen as the 80% quantile of the Kaplan-Meier curve in the chemotherapy arm obtained by aggregating all trials (no distiction between tumors).
tau <- quantile(
  survfit(Surv(time,event)~1, data=data[data$arm==0,]),
  probs = 0.80
)$quantile

# Estimates the difference in RMST (with SE)
tr <- unique(data.frame(data$trial,data$tumor_type))
DiffRMST <- matrix(0,nrow=nrow(tr), ncol=ncol(X))
SE <- matrix(0,nrow=nrow(tr), ncol=ncol(X))
for(tt in seq_along(tr[,1])){
  for(c in 1:ncol(X)){
    I0 <- which(data$trial==tr[tt,1] & data$arm == 0)
    I1 <- which(data$trial==tr[tt,1] & data$arm == 1)
    rmst0 <- rmst.weights(data$time[I0],
                          data$event[I0],
                          P[I0,c],
                          tau)
    rmst1 <- rmst.weights(data$time[I1],
                          data$event[I1],
                          P[I1,c],
                          tau)
    DiffRMST[tt,c] <- rmst1$est - rmst0$est   # difference in rmst
    SE[tt,c] <- sqrt(rmst1$se^2 + rmst0$se^2) # standard error
  }
  #plot(survfit(Surv(time,event)~arm, data=data[I,], weights=P[I,6]),
  #     xlim=c(0,35)) 
}

# Plot of study specific estimates
colnames(DiffRMST) <- cl
d <- as.data.frame(DiffRMST)
d$trial <- tr$data.trial
d$type <- tr$data.tumor_type
colnames(SE) <- cl
e <- as.data.frame(SE)
e$trial <- tr$data.trial
dd <- d %>% gather(PDL1,est,-trial, -type)
ee <- e %>% gather(PDL1,se, -trial)
dd$se    <- ee$se
dd$lcl95 <- dd$est - 1.96 * dd$se
dd$ucl95 <- dd$est + 1.96 * dd$se
dd$type <- factor(dd$type)
levels(dd$type) <- c("Other tumors", "NSCLC")
dd$type <- relevel(dd$type, ref="NSCLC")
fig1 <-ggplot() +
  theme_minimal() +
  #theme(panel.grid.major = element_blank(), 
  #      panel.grid.minor = element_blank()) +
  scale_x_discrete(limits=cl) +
  geom_point(data=dd,aes(x=PDL1, y=est, col=type, size=1/se, shape=type),
             #shape=1,
             position=position_jitterdodge(dodge.width = 0.5)) +
  geom_smooth(data=dd, aes(x=PDL1, y=est, col=type, group=type),
              method="loess",lwd=1.2,
              position=position_dodge(width=0.5),
              se=FALSE,
              show.legend = FALSE) +
  ylab("Difference in RMST") +
  xlab("PD-L1 expression level") +
  labs(title = "Study-specific estimates") +
  scale_color_discrete(name="Tumor type") +
  scale_shape_manual(values=c(0,1),
                     name="Tumor type") +
  scale_size_continuous(name="1/Standard Error") + 
  geom_hline(yintercept=0, lty=2)
  
# Saves the plot
ggsave("./RESULTS/Figure_3.pdf",
       width = 6,
       height = 3.5,
       plot=fig1)


