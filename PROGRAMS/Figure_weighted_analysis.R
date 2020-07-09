
library(tidyr)
library(dplyr)
library(ggplot2)
library(survival)
library(gridExtra)
library(isotone)

# Plot of the distribution of PD-L1 and the class-specific treatment effects
source("./PROGRAMS/functions_piecewise_exponential.R")

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
load("./DATA/weighted_analysis.Rdata")
fit <- wfit
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

dist <- pi %>% group_by(PDL1) %>% summarise(height=mean(prob),
                                    lcl=quantile(prob, probs=0.025),
                                    ucl=quantile(prob, probs=0.975)) %>%
  ggplot(aes(x=PDL1,y=height)) + 
  geom_bar(stat="identity", #fill="#00BCD8"
           fill="cyan3") +
  geom_errorbar(aes(ymin=lcl,ymax=ucl), width=.5, lwd=0.5) +
  theme_minimal() + 
  ylab("Proportion of patients") +
  xlab("PD-L1 expression level") +
  labs(title="(a) Distribution of PD-L1 expression")

# Plot of the class-specific treatment effects, with Isotonic transofrmation

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
monot <- as.data.frame(monot) %>% 
  gather(PDL1,delta) %>% 
  group_by(PDL1) %>%
  summarize(m=mean(delta),
            lcl=quantile(delta, probs = 0.025),
            ucl=quantile(delta, probs = 0.975)) 
monot$type="Other tumors"

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
monot2 <- as.data.frame(monot2) %>% 
  gather(PDL1,delta) %>% 
  group_by(PDL1) %>%
  summarize(m=mean(delta),
            lcl=quantile(delta, probs = 0.025),
            ucl=quantile(delta, probs = 0.975)) 
monot2$type="NSCLC"

eff <- rbind(monot,monot2) %>% 
  ggplot(aes(x=PDL1, y=m, col=type, shape=type)) + 
  theme_minimal() + 
  scale_x_discrete(limits=cl) +
  geom_point(size=3, position=position_dodge(width=0.3)) + 
  geom_line(aes(x=PDL1, y=m, group=type),
            linetype="dotted",
            position=position_dodge(width = 0.3),
            lwd=1,
            show.legend = FALSE) +
  #ylim(0,5) + 
  geom_errorbar(aes(ymin=lcl,ymax=ucl, col=type), 
                width=0, lwd=1,
                position=position_dodge(width=0.3),
                show.legend = FALSE) + 
  xlab("PD-L1 expression level") + 
  ylab("Difference in RMST") +
  theme(legend.position = "bottom",
        legend.margin = margin(-5,-5,-5,-5)) +
  scale_color_discrete(name="Tumor type") +
  scale_shape_discrete(name="Tumor type") +
  labs(title="(b) Class-specific treatment effects") 

# Side by side plots
final <- grid.arrange(dist,eff,ncol=2)

# Saves the figure 
ggsave("./RESULTS/Figure_weighted_analysis.pdf",
       width = 9,
       height = 4.5,
       plot=final)
