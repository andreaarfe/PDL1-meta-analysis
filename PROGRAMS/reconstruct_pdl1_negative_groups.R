
library(dplyr)

# reconstructs the PD-L1 negative data for those publications that did not
# report the corresponding Kaplan-Meier curve.
# The number of patients is obtained from the total sample size and the number
# of patients in the PD-L1 positive group. All survival times are censored at 0.

# List of PD-L1 groups to reconstruct
# IMvigor 211 <5%
# IMvigor 211 <1%
# JAVELIN LUNG 200 <50%
# JAVELIN LUNG 200 <1%
# JAVELIN LUNG 200 <80%
# KEYNOTE 10 <50%
# KEYNOTE 45 <10%
# KEYNOTE 61 <10%
# OAK <5%
# POPLAR <5%

# master dataset
load("./DATA/master_IPD.Rdata")

# Loads the sample sizes
d <- read.csv("./DATA/Publication_level_variables/pubs_variables.csv",
              colClasses = "character")
d$PMID <- as.factor(d$PMID)
d$SAMPLE.SIZE <- as.numeric(d$SAMPLE.SIZE)

# Extracts the sample sizes 
classes <- master %>% group_by(PMID,trial,low,high) %>% 
  count() %>%
  left_join(select(d,PMID,SAMPLE.SIZE),by="PMID")

# Recovers the groups and sample sizes to be recreated
I         <- which(classes$low>0 & classes$high==1.00)
neg       <- classes[I,]
neg$N_neg <- neg$SAMPLE.SIZE-neg$n
neg$high  <- neg$low
neg$low   <- 0
I1 <- with(neg,which(trial=="IMvigor 211" & high==0.05))
I1 <- c(I1, with(neg,which(trial=="IMvigor 211" & high==0.01)))
I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200" & high==0.50)))
# I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200 squamous" & high==0.01)))
# I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200 squamous" & high==0.80)))
# I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200 non-squamous" & high==0.01)))
# I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200 non-squamous" & high==0.80)))
I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200" & high==0.01)))
I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200" & high==0.80)))
I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200" & high==0.01)))
I1 <- c(I1, with(neg,which(trial=="JAVELIN Lung 200" & high==0.80)))
#I1 <- c(I1, with(neg,which(trial=="Keynote 10" & high==0.50)))
I1 <- c(I1, with(neg,which(trial=="Keynote 45" & high==0.10)))
I1 <- c(I1, with(neg,which(trial=="Keynote 61" & high==0.10)))
I1 <- c(I1, with(neg,which(trial=="OAK" & high==0.05)))
I1 <- c(I1, with(neg,which(trial=="POPLAR" & high==0.05)))
neg <- neg[I1,]

# Dataset with reconstructed observations
data <- master
for(i in 1:nrow(neg)){
  new_obs <- data.frame(PMID  = neg$PMID[i],
                        trial = neg$trial[i],
                        time  = rep(0, times=neg$N_neg[i]),
                        event = 0,
                        arm   = 0,
                        low   = neg$low[i],
                        high  = neg$high[i])
  data <- rbind(data, new_obs)
}

# Saves the dataset for the analyses
save(data,file="./DATA/analysis_dset.Rdata")
