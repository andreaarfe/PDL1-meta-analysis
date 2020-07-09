
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

### Stratified analysis: distribution of PD-L1 using TPS vs. CPS

load("./DATA/sens_IC.Rdata")
load("./DATA/sens_TC.Rdata")

NITER <- nrow(fit.IC$pi)
NBURN <- round(NITER/2)
pi.IC <- as.data.frame(fit.IC$pi[(1+NBURN):NITER,])
pi.TC <- as.data.frame(fit.TC$pi[(1+NBURN):NITER,])
pi.TC$assay <- "TPS"
pi.IC$assay <- "CPS"
pi <- rbind(pi.IC,pi.TC)
cl <- c("0%-1%", 
        "1%-5%", 
        "5%-10%", 
        "10%-50%", 
        "50%-80%", 
        "80%-100%")
names(pi) <- c(cl, "assay")
pi <- gather(pi, PDL1, prob, -assay)
pi$assay <- as.factor(pi$assay)
pi.ictc <- pi
rm(list=c("fit.IC","fit.TC","pi"))
gc()

p1 <- pi.ictc %>% 
  group_by(assay,PDL1) %>%
  summarise(m = mean(prob),
            lcl = quantile(prob, probs=0.025),
            ucl = quantile(prob, probs=0.975)) %>%
  ggplot(aes(x=PDL1, y=m, fill=assay)) +
  theme_minimal() + 
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.5, 
                position=position_dodge(0.85)) + 
  scale_x_discrete(limits=cl) +
  xlab("PD-L1 expression") +
  ylab("Proportion of patients") +
  labs(title="(a) PD-L1 distribution assessed by TPS or CPS assays") + 
  scale_fill_discrete(name="Assay type") + 
  theme(legend.position = "bottom",
        legend.margin = margin(-5,-5,-5,-5),
        plot.title = element_text(size=10))


### Stratified analysis: distribution of PD-L1 using NSCLC vs. others


load("./DATA/sens_NSCLC.Rdata")
load("./DATA/sens_notNSCLC.Rdata")

NITER <- nrow(fit.NSCLC$pi)
NBURN <- round(NITER/2)
pi.NSCLC <- as.data.frame(fit.NSCLC$pi[(1+NBURN):NITER,])
pi.notNSCLC <- as.data.frame(fit.notNSCLC$pi[(1+NBURN):NITER,])
pi.NSCLC$type <- "NSCLC"
pi.notNSCLC$type <- "Other Tumors"
pi <- rbind(pi.NSCLC,pi.notNSCLC)
cl <- c("0%-1%", 
        "1%-5%", 
        "5%-10%", 
        "10%-50%", 
        "50%-80%", 
        "80%-100%")
names(pi) <- c(cl, "type")
pi <- gather(pi, PDL1, prob, -type)
pi$type <- as.factor(pi$type)

p2 <- pi %>% 
  group_by(type,PDL1) %>%
  summarise(m = mean(prob),
            lcl = quantile(prob, probs=0.025),
            ucl = quantile(prob, probs=0.975)) %>%
  ggplot(aes(x=PDL1, y=m, fill=type)) +
  theme_minimal() + 
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.5, 
                position=position_dodge(0.85)) + 
  scale_x_discrete(limits=cl) +
  xlab("PD-L1 expression") +
  ylab("Proportion of patients") +
  labs(title="(c) PD-L1 distribution for NSCLC or other tumors") + 
  scale_fill_discrete(name="Tumor type") + 
  theme(legend.position = "bottom",
        legend.margin = margin(-5,-5,-5,-5),
        plot.title = element_text(size=10))


### Stratified analysis: distribution of PD-L1 by ICH type

load("./DATA/sens_ICH22C3.Rdata")
load("./DATA/sens_ICH288.Rdata")
load("./DATA/sens_ICH7310.Rdata")
load("./DATA/sens_ICHSP142.Rdata")

NITER <- nrow(fit.ICH22C3$pi)
NBURN <- round(NITER/2)
pi.ICH22C3 <- as.data.frame(fit.ICH22C3$pi[(1+NBURN):NITER,])
pi.ICH288 <- as.data.frame(fit.ICH288$pi[(1+NBURN):NITER,])
pi.ICH7310 <- as.data.frame(fit.ICH7310$pi[(1+NBURN):NITER,])
pi.ICHSP142 <- as.data.frame(fit.ICHSP142$pi[(1+NBURN):NITER,])
pi.ICH22C3$assay <- "22C3"
pi.ICH288$assay <- "28-8"
pi.ICH7310$assay <- "73-10"
pi.ICHSP142$assay <- "SP142"
pi <- rbind(pi.ICH22C3,
            pi.ICH288,
            pi.ICH7310,
            pi.ICHSP142
)
cl <- c("0%-1%", 
        "1%-5%", 
        "5%-10%", 
        "10%-50%", 
        "50%-80%", 
        "80%-100%")
names(pi) <- c(cl, "assay")
pi <- gather(pi, PDL1, prob, -assay)
pi$assay <- as.factor(pi$assay)
#rm(list=c("fit.IC","fit.TC","pi"))
#gc()

pich <- pi %>% 
  group_by(assay,PDL1) %>%
  summarise(m = mean(prob),
            lcl = quantile(prob, probs=0.025),
            ucl = quantile(prob, probs=0.975)) %>%
  ggplot(aes(x=PDL1, y=m, fill=assay)) +
  theme_minimal() + 
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.5, 
                position=position_dodge(0.85)) + 
  scale_x_discrete(limits=cl) +
  xlab("PD-L1 expression") +
  ylab("Proportion of patients") +
  labs(title="(b) PD-L1 distribution by IHC assay type") + 
  scale_fill_discrete(name="IHC assay type") + 
  theme(legend.position = "bottom",
        legend.margin = margin(-5,-5,-5,-5),
        plot.title = element_text(size=10))

# Saves the figure
final <- grid.arrange(p1,pich,p2,nrow=2)
ggsave("./RESULTS/Figure_2.pdf",
       plot = final,
       width = 9,
       height = 7)
