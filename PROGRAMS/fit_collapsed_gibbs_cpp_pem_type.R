
load("./DATA/analysis_dset.Rdata")
source("./PROGRAMS/functions_collapsed_gibbs_cpp_pem_types.R")

# Identifies NSCLC trials
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

# Fit the model
fit <- mix.exp.collapsed.cpp.pem.type(data    = data, 
                                      cuts    = cuts,
                                      breaks  = breaks,
                                      alpha   = alpha,
                                      shape   = shape,
                                      rate    = rate,
                                      verbose = TRUE, 
                                      NITER   = 20000)

# Saves the output (no burn-in!)
save(fit, file="./DATA/fit_pem_type.Rdata")


