
load("./DATA/analysis_dset.Rdata")
source("./PROGRAMS/functions_collapsed_gibbs_cpp_pem_types.R")

# Identifies non-placebo-controlled trials
not_placebo <- as.numeric(!(data$trial %in% c("Keynote 189", "ATTRACTION-2")))

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

# Datasets for sensitivity analysis
data.nonplac <- data[not_placebo==1,]
data.nonplac$arm <- as.factor(levels(data.nonplac$arm)[data.nonplac$arm])

# Prior specification
breaks <- c(0,quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5))
cuts   <- sort(unique(c(data$low,data$high)))
alpha  <- diff(cuts)/min(diff(cuts))
shape  <- matrix(10*log(2)/12,nrow=length(cuts)-1,ncol=length(breaks))
rate   <- matrix(10,nrow=length(cuts)-1,ncol=length(breaks))

# Runs the Gibbs sampler
fit.nonplac <- mix.exp.collapsed.cpp.pem.type(data  = data.nonplac, 
                                            cuts    = cuts,
                                            breaks  = breaks,
                                            alpha   = alpha,
                                            shape   = shape,
                                            rate    = rate,
                                            verbose = TRUE, 
                                            NITER   = 20000)
save(fit.nonplac, file="./DATA/sens_not_placebo.Rdata")


