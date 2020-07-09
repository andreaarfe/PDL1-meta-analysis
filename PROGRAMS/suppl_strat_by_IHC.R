
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

# Identifies IHC 22C3 trials
ICH22C3 <- as.numeric(data$trial %in% c('Keynote 61',
																				'Keynote 189',
																				'Keynote 45'
																				))

# Identifies IHC 28-8  trials 
ICH288 <- as.numeric(data$trial %in% c('Checkmate 141 updated',
																			 'ATTRACTION-2'	,
																			 'Checkmate 057',
																			 'Checkmate 025',
																			 'Checkmate 017',
																			 'Checkmate 066'))

# Identifies IHC 73-10  trials 
ICH7310 <- as.numeric(data$trial %in% c('JAVELIN Lung 200',
																				'JAVELIN Gastric 300'
))

# Identifies IHC SP142  trials 
ICHSP142 <- as.numeric(data$trial %in% c('IMvigor 211',
																				 'OAK',
																				 'POPLAR'
))

# Datasets for sensitivity analysis
data$tumor_type <- NSCLC
data.ICH22C3  <- data[ICH22C3==1,]
data.ICH288   <- data[ICH288==1,]
data.ICH7310  <- data[ICH7310==1,]
data.ICHSP142 <- data[ICHSP142==1,]
data.ICH22C3$arm <- as.factor(levels(data.ICH22C3$arm)[data.ICH22C3$arm])
data.ICH288$arm <- as.factor(levels(data.ICH288$arm)[data.ICH288$arm])
data.ICH7310$arm <- as.factor(levels(data.ICH7310$arm)[data.ICH7310$arm])
data.ICHSP142$arm <- as.factor(levels(data.ICHSP142$arm)[data.ICHSP142$arm])

# Prior specification
breaks <- c(0,quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5))
cuts   <- sort(unique(c(data$low,data$high)))
alpha  <- diff(cuts)/min(diff(cuts))
shape  <- matrix(10*log(2)/12,nrow=length(cuts)-1,ncol=length(breaks))
rate   <- matrix(10,nrow=length(cuts)-1,ncol=length(breaks))

# Runs the Gibbs sampler
fit.ICH22C3  <- mix.exp.collapsed.cpp.pem.type(data    = data.ICH22C3 , 
                                            cuts    = cuts,
                                            breaks  = breaks,
                                            alpha   = alpha,
                                            shape   = shape,
                                            rate    = rate,
                                            verbose = TRUE, 
                                            NITER   = 20000)
save(fit.ICH22C3, file="./DATA/sens_ICH22C3.Rdata")
rm(list="fit.ICH22C3")

fit.ICH288  <- mix.exp.collapsed.cpp.pem.type(data    = data.ICH288 , 
																							 cuts    = cuts,
																							 breaks  = breaks,
																							 alpha   = alpha,
																							 shape   = shape,
																							 rate    = rate,
																							 verbose = TRUE, 
																							 NITER   = 20000)
save(fit.ICH288, file="./DATA/sens_ICH288.Rdata")
rm(list="fit.ICH288")

fit.ICH7310  <- mix.exp.collapsed.cpp.pem.type(data    = data.ICH7310 , 
																							 cuts    = cuts,
																							 breaks  = breaks,
																							 alpha   = alpha,
																							 shape   = shape,
																							 rate    = rate,
																							 verbose = TRUE, 
																							 NITER   = 20000)
save(fit.ICH7310, file="./DATA/sens_ICH7310.Rdata")
rm(list="fit.ICH7310")

fit.ICHSP142  <- mix.exp.collapsed.cpp.pem.type(data    = data.ICHSP142 , 
																							 cuts    = cuts,
																							 breaks  = breaks,
																							 alpha   = alpha,
																							 shape   = shape,
																							 rate    = rate,
																							 verbose = TRUE, 
																							 NITER   = 20000)
save(fit.ICHSP142, file="./DATA/sens_ICHSP142.Rdata")
rm(list="fit.ICHSP142")



