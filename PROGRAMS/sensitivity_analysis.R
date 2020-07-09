
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

# Identifies trials using assays based on infiltrating cells (IC)
IC <- as.numeric(data$trial %in% c("Keynote 61",
                                   "IMvigor 211",
                                   "IMvigor 211 high TMB",
                                   "Keynote 45",
                                   "OAK",
                                   "POPLAR"))

# Datasets for sensitivity analysis
data$tumor_type <- NSCLC
data.NSCLC <- data[NSCLC==1,]
data.NSCLC$arm <- as.factor(levels(data.NSCLC$arm)[data.NSCLC$arm])
data.notNSCLC <- data[NSCLC==0,]
data.notNSCLC$arm <- as.factor(levels(data.notNSCLC$arm)[data.notNSCLC$arm]) 
data.IC <- data[IC==1,]
data.IC$arm <- as.factor(levels(data.IC$arm)[data.IC$arm])
data.TC <- data[IC==0,]
data.TC$arm <- as.factor(levels(data.TC$arm)[data.TC$arm])

# Prior specification
breaks <- c(0,quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5))
cuts   <- sort(unique(c(data$low,data$high)))
alpha  <- diff(cuts)/min(diff(cuts))
shape  <- matrix(10*log(2)/12,nrow=length(cuts)-1,ncol=length(breaks))
rate   <- matrix(10,nrow=length(cuts)-1,ncol=length(breaks))

# Runs the Gibbs sampler
fit.NSCLC <- mix.exp.collapsed.cpp.pem.type(data    = data.NSCLC, 
                                            cuts    = cuts,
                                            breaks  = breaks,
                                            alpha   = alpha,
                                            shape   = shape,
                                            rate    = rate,
                                            verbose = FALSE, 
                                            NITER   = 20000)
save(fit.NSCLC, file="./DATA/sens_NSCLC.Rdata")
rm(list="fit.NSCLC")

fit.notNSCLC <- mix.exp.collapsed.cpp.pem.type(data    = data.notNSCLC, 
                                            cuts    = cuts,
                                            breaks  = breaks,
                                            alpha   = alpha,
                                            shape   = shape,
                                            rate    = rate,
                                            verbose = FALSE, 
                                            NITER   = 20000)
save(fit.notNSCLC, file="./DATA/sens_notNSCLC.Rdata")
rm(list="fit.notNSCLC")

fit.IC <- mix.exp.collapsed.cpp.pem.type(data    = data.IC, 
                                            cuts    = cuts,
                                            breaks  = breaks,
                                            alpha   = alpha,
                                            shape   = shape,
                                            rate    = rate,
                                            verbose = FALSE, 
                                            NITER   = 20000)
save(fit.IC, file="./DATA/sens_IC.Rdata")
rm(list="fit.IC")

fit.TC <- mix.exp.collapsed.cpp.pem.type(data    = data.TC, 
                                         cuts    = cuts,
                                         breaks  = breaks,
                                         alpha   = alpha,
                                         shape   = shape,
                                         rate    = rate,
                                         verbose = FALSE, 
                                         NITER   = 20000)
save(fit.TC, file="./DATA/sens_TC.Rdata")
rm(list="fit.TC")



