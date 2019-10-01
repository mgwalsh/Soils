# Stacked calibrations of soil compositional properties with Alpha ZnSe spectra
# M. Walsh, September 2019

# Required packages
# install.packages(c("devtools","caret","pls","glmnet","randomForest","gbm","Cubist","nnet","plyr","doParallel")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(pls)
  require(glmnet)
  require(randomForest)
  require(gbm)
  require(Cubist)
  require(nnet)
  require(plyr)
  require(doParallel)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Soils/blob/master/Alpha_recal_data.R
# ... or
# source_https <- function(url, ...) {
#  # load package
#  require(RCurl)
#  # parse and evaluate .R script
#  sapply(c(url, ...), function(u) {
#    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
#  })
# }
# source_https("https://github.com/mgwalsh/Soils/blob/master/Alpha_recal_data.R")
rm(list=setdiff(ls(), c("nbal"))) ## scrubs extraneous objects in memory)

# set randomization seed
seed <- 1385321
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(nbal$Fv, p = 9/10, list=F, times = 1)
cal <- nbal[ gsIndex,]
val <- nbal[-gsIndex,]

# GeoSurvey calibration labels
labs <- c("Fv") ## insert other labels (C,N,P ...) here!
lcal <- as.vector(t(cal[labs]))

# spectral calibration features
fcal <- cal[,15:1728]
fpca <- cal[,1729:1748] ## PCA variables

# PLS ---------------------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(seed)
tc <- trainControl(method="repeatedcv", number=10, repeats=3, allowParallel=T)
tg <- expand.grid(ncomp=seq(2,40, by=2)) ## model tuning steps

pls <- train(fcal, lcal,
             method = "pls",
             preProc = c("center", "scale"),
             tuneGrid = tg,
             trControl = tc)
print(pls)
stopCluster(mc)
fname <- paste("./Results/", labs, "_pls.rds", sep = "")
saveRDS(pls, fname)

# Elastic net <glmnet> ----------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(seed)
tc <- trainControl(method="cv", allowParallel=T)
# tg <- needs tuning

# model training
en <- train(fcal, lcal,
            method = "glmnet",
            family = "gaussian",
            tuneGrid = tg,
            trControl = tc)
print(en)
stopCluster(mc)
fname <- paste("./Results/", labs, "_en.rds", sep = "")
saveRDS(en, fname)

# Random forest <randomForest> --------------------------------------------
# Random forest with spectral PCA covariates
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(seed)
tc <- trainControl(method="cv", allowParallel=T)
tg <- expand.grid(mtry = seq(2,20, by=1)) ## model tuning

# model training
rf <- train(fpca, lcal,
            method = "rf",
            ntree = 501,
            tuneGrid = tg,
            trControl = tc)
print(rf)
stopCluster(mc)
fname <- paste("./Results/", labs, "_rf.rds", sep = "")
saveRDS(rf, fname)

# Generalized boosting <gbm> ----------------------------------------------
# Generalized boosting with spectral PCA variables
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(seed)
tc <- trainControl(method = "cv", allowParallel = T)
tg <- expand.grid(interaction.depth = seq(2,20, by=1), shrinkage = seq(0.02,0.1, by=0.02), n.trees = 501,
                  n.minobsinnode = 25) ## model tuning steps

gb <- train(fpca, lcal, 
            method = "gbm", 
            trControl = tc,
            tuneGrid = tg)
print(gb)
stopCluster(mc)
fname <- paste("./Results/", labs, "_gb.rds", sep = "")
saveRDS(gb, fname)

# Cubist <Cubist> ---------------------------------------------------------
# Cubist with spectral PCA variables
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(seed)
tc <- trainControl(method="repeatedcv", number=10, repeats=3, allowParallel = T)
# tg <- cubistControl() may need to be tuned

cu <- train(fpca, lcal, 
            method = "cubist", 
            trControl = tc)
print(cu)
stopCluster(mc)
fname <- paste("./Results/", labs, "_cu.rds", sep = "")
saveRDS(cu, fname)

# Neural net <nnet> -------------------------------------------------------
# nnet with spectral PCA variables
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(seed)
tc <- trainControl(method = "cv", allowParallel = T)
# tg <- needs tuning

nn <- train(fpca, lcal, 
            method = "nnet", 
            trControl = tc)
print(nn)
stopCluster(mc)
fname <- paste("./Results/", labs, "_nn.rds", sep = "")
saveRDS(nn, fname)

