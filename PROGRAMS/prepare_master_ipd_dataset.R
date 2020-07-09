
# imports the IPD data and prepares a master dataset
path <- "./DATA/IPD/"
infiles <- list.files(path=path, full.names = FALSE)
inpaths <- list.files(path=path, full.names = TRUE)
dsets <- list()
for(i in 1:length(inpaths)){
  dsets[[i]] <- read.csv(inpaths[i],stringsAsFactors = FALSE)
  dsets[[i]] <- dsets[[i]][,c("Trial","Time","Event","Arm","Expression","Units")]
  # Reads the PMID from the file name
  dsets[[i]]$PMID <- strsplit(infiles[i],"_")[[1]][2]
} 
master <- do.call("rbind",dsets) # binds all datasets into a single one

# Check that all times are in months
table(master$Units)
# Yes!

# Corrects some of the trial labels
master$Trial[master$Trial == "Checkmate 25"] <- "Checkmate 025"
master$Trial[master$Trial == "JAVELIN"]      <- "JAVELIN Gastric 300"
master$Trial[master$Trial == "JAVELIN Lung 200 squamous"] <- "JAVELIN Lung 200"
master$Trial[master$Trial == "JAVELIN Lung 200 non-squamous"] <- "JAVELIN Lung 200"


# recodes the cuts as high and low interval limits
labels <- unique(master$Expression)
lows   <- c(0.00, 0.05, 0.01, 0.00, 0.10, 0.00, 0.50, 0.80, 0.01)
highs  <- c(0.05, 1.00, 1.00, 0.01, 1.00, 0.10, 1.00, 1.00, 0.50)
#lows   <- c(0.01, 0.00, 0.05, 0.00, 0.10, 0.00, 0.50, 0.80, 0.01)
#highs  <- c(1.00, 0.01, 1.00, 0.05, 1.00, 0.10, 1.00, 1.00, 0.50)
names(lows)  <- labels
names(highs) <- labels
master$low  <- lows[master$Expression]
master$high <- highs[master$Expression]
names(master$low)   <- NULL
names(master$highs) <- NULL

# recodes string as factors
master$Trial      <- as.factor(master$Trial)
master$Expression <- as.factor(master$Expression)
#master$Outcome    <- as.factor(master$Outcome)
#master$Figure     <- as.factor(master$Figure)
master$Arm        <- as.factor(master$Arm)
#master$Disease    <- as.factor(master$Disease)
#master$Units      <- as.factor(toupper(master$Units))
master$PMID       <- as.factor(master$PMID) 

# prepares the data to fit the exponential misture model
data <- master
data$event <- master$Event
data$arm   <- master$Arm
data$time  <- master$Time
data$trial <- master$Trial
data <- data[,c("trial","PMID","time","event","arm","high","low")]

#  Codes the treatment and control arms in each trial
#  Experimental = any immunotherapy
#  Control = everything else
immunos <- c("atezolizumab",
             "avelumab",
             "nivolumab",
             #"nivolumab+ipilimumab",
             "pembrolizumab",
             "pembrolizumab-combo", # It's Keynote 189: combination of pembrolizumab and standard chemotherapy
             "pembrolizumab 10mg/kg",
             "pembrolizumab 2mg/kg")
data$arm <- as.factor(ifelse(data$arm %in% immunos,1,0))
master <- data

# saves some descriptives
# sink("./Results/descriptives_master.txt")
# summary(master)
# sink()

# saves the master dataset
save(list = "master", file="./DATA/master_IPD.Rdata")




