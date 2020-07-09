
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
ll <- levels(data$trial)
ll[ll=="Checkmate 141 updated"] <- "Checkmate 141"
levels(data$trial) <- ll
data$tumor_type <- NSCLC
data$low[data$trial=="Keynote 189" & data$low==0.5] <- 0.01
data$high[data$trial=="Keynote 189" & data$high==0.5] <- 1

# Selects only trials that used 1% as a cut-off
d <- data %>% 
  group_by(trial,low,high) %>%
  filter(time>0) %>%
  filter( (low==0 & high==0.01) | (low==0.01 & high==1) )
d$trial <- as.factor(levels(d$trial)[as.numeric(d$trial)])

# Computes study- and PD-L1-specific (<1% or >1%) RMST estimates
studies <-  levels(d$trial)
res <- data.frame(trial="", 
                  strata="", 
                  est=NA, 
                  se=rep(NA,times=length(studies)*2),
                  stringsAsFactors=FALSE)
for(i in seq_along(studies)){
  s <- studies[i]
  less <- which(d$trial==s & d$low==0.00 & d$high==0.01)
  more <- which(d$trial==s & d$low==0.01 & d$high==1)
  rl_est <- NA
  rm_est <- NA
  rl_se <- NA
  rm_se <- NA
  if(length(less)>0){ 
    rl <- with(d[less,], rmst2(time,event,arm)) 
    rl_est <- rl$unadjusted.result[1,1]
    rl_se  <- (rl$unadjusted.result[1,3]-rl$unadjusted.result[1,2])/(2*1.96)
  }
  if(length(more)>0){ 
    rm <- with(d[more,], rmst2(time,event,arm)) 
    rm_est <- rm$unadjusted.result[1,1]
    rm_se  <- (rm$unadjusted.result[1,3]-rm$unadjusted.result[1,2])/(2*1.96)
  }
  for(j in -1:0){
    res[2*i+j,"trial"]  <- studies[i]
    res[2*i+j,"strata"] <- ifelse(j==-1,"PD-L1: 0%-1%","PD-L1: 1%-100%")
    res[2*i+j,"est"]    <- ifelse(j==-1, rl_est, rm_est)
    res[2*i+j,"se"]     <- ifelse(j==-1, rl_se, rm_se)
  }
}
res$trial  <- factor(res$trial)
res$strata <- factor(res$strata)
res$NSCLC  <- res$trial %in% c("Checkmate 017",
                               "Checkmate 057",
                               #"Keynote 10",
                               "POPLAR",
                               "OAK",
                               #"Checkmate 227", excluded because combination trial
                               "Keynote 189",
                               "JAVELIN Lung 200")
res$lcl95 <- res$est - 1.96 * res$se
res$ucl95 <- res$est + 1.96 * res$se

# Meta-analytic Q p-value for heterogeneity and I-squared statistics
res$Q  <- ""
res$I2 <- ""
ma1 <- with(res %>% filter(NSCLC==FALSE & strata=="PD-L1: 0%-1%"), metagen(TE=est, seTE=se))
ma2 <- with(res %>% filter(NSCLC==FALSE & strata=="PD-L1: 1%-100%"), metagen(TE=est, seTE=se))
ma3 <- with(res %>% filter(NSCLC==TRUE & strata=="PD-L1: 0%-1%"), metagen(TE=est, seTE=se))
ma4 <- with(res %>% filter(NSCLC==TRUE & strata=="PD-L1: 1%-100%"), metagen(TE=est, seTE=se))
res$Q[res$NSCLC==FALSE & res$strata=="PD-L1: 0%-1%"] <- paste0("Q-test p: ", round(ma1$pval.Q,3))
res$Q[res$NSCLC==FALSE & res$strata=="PD-L1: 1%-100%"] <- paste0("Q-test p: ", round(ma2$pval.Q,3))
res$Q[res$NSCLC==TRUE & res$strata=="PD-L1: 0%-1%"] <- paste0("Q-test p: ", round(ma3$pval.Q,3))
res$Q[res$NSCLC==TRUE & res$strata=="PD-L1: 1%-100%"] <- paste0("Q-test p: ", round(ma4$pval.Q,3))
res$I2[res$NSCLC==FALSE & res$strata=="PD-L1: 0%-1%"] <- paste0("I-squared: ", round(ma1$I2,1))
res$I2[res$NSCLC==FALSE & res$strata=="PD-L1: 1%-100%"] <- paste0("I-squared: ", round(ma2$I2,1))
res$I2[res$NSCLC==TRUE & res$strata=="PD-L1: 0%-1%"] <- paste0("I-squared: ", round(ma3$I2,1))
res$I2[res$NSCLC==TRUE & res$strata=="PD-L1: 1%-100%"] <- paste0("I-squared: ", round(ma4$I2,1))

# Plot
p1 <- res %>% 
  filter(NSCLC==TRUE) %>%
  ggplot(aes(x=est, y=trial)) +
  geom_point() + 
  facet_grid(~strata + Q + I2) + 
  geom_errorbarh(aes(xmin=lcl95, xmax=ucl95)) +
  theme_bw() +
  labs(title="a) Study-specific estrimates by PD-L1 level - NSCLC") +
  xlab("Difference in RMST (months)") + 
  ylab("Trial") +
  theme(panel.spacing = unit(0.85, "lines")) 

p2 <- res %>% 
  filter(NSCLC==FALSE) %>%
  ggplot(aes(x=est, y=trial)) +
  geom_point() + 
  facet_grid(~strata + Q + I2) + 
  geom_errorbarh(aes(xmin=lcl95, xmax=ucl95)) +
  labs(title="b) Study-specific estrimates by PD-L1 level - other tumors") +
  xlab("Difference in RMST (months)") + 
  ylab("Trial") +
  theme(panel.spacing = unit(0.85, "lines")) +
  theme_bw()

fullplot <- grid.arrange(p1, p2, nrow=2)

# Saves the output
ggsave("./RESULTS/forestplot_het_rawests.pdf",
       width = 6.5,
       height = 6,
       plot=fullplot)




















