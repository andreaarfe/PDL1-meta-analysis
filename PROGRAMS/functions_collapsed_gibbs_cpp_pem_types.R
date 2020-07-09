
library(Rcpp)

sourceCpp("./PROGRAMS/cpp_code_pem_tumor_type.cpp")

# Function to estimate the exponential mixture model for treatment effects within
# each PD-L1 expression class
mix.exp.collapsed.cpp.pem.type <- function(data,
                              NITER = 1000,
                              cuts  = sort(unique(c(data$low,data$high))),
                              breaks = c(0,quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5)),     
                              alpha = rep(1,times=length(cuts)-1),
                              shape = matrix(1,nrow=length(cuts)-1,ncol=length(breaks)),
                              rate  = matrix(1,nrow=length(cuts)-1,ncol=length(breaks)),
                              weights = rep(1, times=nrow(data)),
                              verbose=FALSE){
  # Marix of PD-L1 class indicators
  X <- matrix(0,nrow=nrow(data),ncol=length(cuts)-1)
  for(i in 1:nrow(X)){
    for(j in 1:ncol(X)){
      if(data$low[i]<=cuts[j] & cuts[j+1]<=data$high[i]) X[i,j] <- 1
    }
  }
  # Events per interval
  EV <- matrix(0, nrow=nrow(data), ncol=length(breaks))
  for(i in 1:nrow(data)){
    if(data$event[i]==1){
      I <- findInterval(data$time[i],
                        breaks,
                        all.inside = FALSE,
                        left.open = TRUE, 
                        rightmost.closed = TRUE)
      EV[i,I] <- 1
    }
  }
  event <- EV
  # Person-time per interval
  PT <- matrix(0, nrow=nrow(data), ncol=length(breaks))
  delta <- c(diff(breaks),Inf)
  for(i in 1:nrow(data)){
    PT[i,] <- rowSums(sapply(data$time[i],function(t) pmin(delta,pmax(0,t-breaks))))
  }
  time <- PT
  # Applies individual-level weights
  time  <- time*weights
  event <- event*weights 
  # Initialization of class memberships
  k <- length(cuts)-1
  H <- length(breaks)
  class  <- sapply(1:nrow(X), function(i) sample.int(k,size=1,prob=X[i,]))
  Nclass <- as.vector(table(class))
  
  # Runs the collapsed data-augmentation Gibbs sampler
  trace <- cpp_gibbs_collapsed_pem_type(
                               ifelse(data$arm=="0",0,1),
                               time,
                               event,
                               data$tumor_type,
                               2,
                               k,
                               H,
                               class-1,
                               X,
                               Nclass,
                               alpha,
                               shape,
                               rate,
                               NITER,
                               verbose)
  return(list(cuts=cuts,
              breaks=breaks,
              pi=trace$pi,
              classes=trace$classes+1,
              lambda=array(trace$lambda, dim=c(NITER,2,k,H,2))))
}

# load("./DATA/analysis_dset.Rdata")
# 
# #Selects only NSCLC trials
# NSCLC <- as.numeric(data$trial %in% c("Checkmate 017",
#                                  "Checkmate 057",
#                                  "Keynote 10",
#                                  "POPLAR",
#                                  "OAK",
#                                  #"Checkmate 227", excluded because combination trial
#                                  "Keynote 189",
#                                  "JAVELIN Lung 200"))
# data$tumor_type <- NSCLC
# 
# 
# breaks <- c(0,quantile(data$time[data$arm==0 & data$event==1], probs=(1:4)/5))
# cuts   <- sort(unique(c(data$low,data$high)))
# alpha  <- diff(cuts)/min(diff(cuts))
# shape <- matrix(10*log(2)/12,nrow=length(cuts)-1,ncol=length(breaks))
# rate  <- matrix(10,nrow=length(cuts)-1,ncol=length(breaks))
# mix.exp.collapsed.cpp.pem.type(data, 
#                                alpha=alpha,
#                                rate=rate,
#                                shape=shape,
#                                verbose=TRUE, NITER=1000)
