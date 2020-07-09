
library(tidyr)
library(dplyr)
library(ggplot2)
library(survival)
library(meta)
library(survRM2)
library(gridExtra)

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

# study specific estimates
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
dd$est[dd$se==0] <- NA
dd$lcl95[dd$se==0] <- NA
dd$ucl95[dd$se==0] <- NA
dd$PDL1 <- factor(dd$PDL1, levels=cl)

# Forest plots

# NSCLC
out <- list()
for(i in seq_along(cl)) out[[i]] <- with(filter(dd,type=="NSCLC" & se!=0 & PDL1==cl[i]),metagen(est,se))
pval.Q <- sapply(seq_along(out), function(x) out[[x]]$pval.Q)
I2 <- sapply(seq_along(out), function(x) out[[x]]$I2)
I2.low <- sapply(seq_along(out), function(x) out[[x]]$lower.I2)
I2.high <- sapply(seq_along(out), function(x) out[[x]]$upper.I2)
het <- data.frame(PDL1=cl, 
                  label0 = paste0("PD-L1: ",cl),
                  label1=paste("Q-test p: ", 
                               format(round(pval.Q, 3), nsmall = 3)),
                  label2=paste("I-squared: ",
                               format(round(I2, 1), nsmall = 1),
                               sep=""))
het$label0 <- factor(het$label0, levels=paste0("PD-L1: ",cl))

p1 <- dd %>% 
  filter(type=="NSCLC") %>%
  left_join(het,by="PDL1") %>%
  ggplot(aes(x=est,y=trial)) +
  geom_point() +
  geom_errorbarh(aes(xmin=lcl95, xmax=ucl95),height=0.25) +
  theme_bw() +
  facet_grid(~label0 + label1 + label2) +
  labs(title="a) Study-specific estrimates by PD-L1 level - NSCLC") +
  xlab("Difference in RMST (months)") + 
  ylab("Trial") +
  theme(panel.spacing = unit(0.85, "lines"))
  

out.oth <- list()
for(i in seq_along(cl)) out.oth[[i]] <- with(filter(dd,type!="NSCLC" & se!=0 & PDL1==cl[i]),metagen(est,se))
pval.Q.oth <- sapply(seq_along(out), function(x) out.oth[[x]]$pval.Q)
I2.oth <- sapply(seq_along(out), function(x) out.oth[[x]]$I2)
I2.low.oth <- sapply(seq_along(out), function(x) out.oth[[x]]$lower.I2)
I2.high.oth <- sapply(seq_along(out), function(x) out.oth[[x]]$upper.I2)
het.oth <- data.frame(PDL1=cl,
                  label0 = paste0("PD-L1: ",cl),
                  label1=paste("Q-test p: ", 
                               format(round(pval.Q.oth, 3), nsmall = 3)),
                  label2=paste("I-squared: ",
                               format(round(I2.oth, 1), nsmall = 1),
                               sep=""))
het.oth$label0 <- factor(het.oth$label0, levels=paste0("PD-L1: ",cl))

levels(dd$trial)[6] <- "Checkmate 141"
p2 <- dd %>% 
  filter(type!="NSCLC") %>%
  left_join(het.oth,by="PDL1") %>%
  ggplot(aes(x=est,y=trial)) +
  geom_point() +
  geom_errorbarh(aes(xmin=lcl95, xmax=ucl95),height=0.25) +
  theme_bw() +
  facet_grid(~label0 + label1 + label2) +
  labs(title="b) Study-specific estrimates by PD-L1 level - Other tumors") +
  xlab("Difference in RMST (months)") + 
  ylab("Trial") +
  theme(panel.spacing = unit(0.85, "lines"))

fplot <- grid.arrange(p1,p2, nrow=2)

# Saves the plot
ggsave("./RESULTS/forestplot_heterogeneity.pdf",
       width = 10,
       height = 6,
       plot=fplot)
