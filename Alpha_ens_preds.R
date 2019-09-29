# Stacked calibrations of soil properties with Alpha ZnSe spectra
# M. Walsh, September 2019

# Required packages
# install.packages(c("devtools","caret","pls","randomForest","gbm","nnet","plyr","doParallel")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(pls)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(plyr)
  require(doParallel)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Cropland-Atlas/blob/master/ZM_GS19_data.R
rm(list=setdiff(ls(), c("nbal"))) ## scrubs extraneous objects in memory)
# nbal <- as.data.frame(nbal[complete.cases(nbal[ ,c(1:1727)]),]) ## removes incomplete cases

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(nbal$Fv, p = 9/10, list = F, times = 1)
cal <- nbal[ gsIndex,]
val <- nbal[-gsIndex,]

# GeoSurvey calibration labels
labs <- c("Fv") ## substitute other labels here
lcal <- as.vector(t(cal[labs]))

# raster calibration features
fcal <- cal[,14:1727]

# PLS ---------------------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
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

# Random forest <randomForest> --------------------------------------------
seed <- 12358
set.seed(seed)

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method="cv", allowParallel=T)
tg <- expand.grid(mtry = seq(10,100, by=10)) ## model tuning steps

# model training
rf1 <- train(fcal, lcal,
             preProc = c("center","scale"),
             method = "rf",
             ntree = 501,
             tuneGrid = tg,
             trControl = tc)
print(rf1)
stopCluster(mc)
fname <- paste("./Results/", labs, "_rf1.rds", sep = "")
saveRDS(rf1, fname)

# Random forest with spectral PCA covariates
seed <- 12358
set.seed(seed)

# PCA features
fpca <- cal[,1728:1747]

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method="cv", allowParallel=T)
tg <- expand.grid(mtry = seq(2,20, by=1)) ## model tuning steps

# model training
rf2 <- train(fpca, lcal,
             preProc = c("center","scale"),
             method = "rf",
             ntree = 501,
             tuneGrid = tg,
             trControl = tc)
print(rf2)
stopCluster(mc)
fname <- paste("./Results/", labs, "_rf2.rds", sep = "")
saveRDS(rf, fname)

# Generalized boosting <gbm> ----------------------------------------------
seed <- 12358
set.seed(seed)

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)
tg <- expand.grid(interaction.depth = seq(10,100, by=10), shrinkage = seq(0.02,0.1, by=0.02), n.trees = 501,
                  n.minobsinnode = 25) ## model tuning steps

gb <- train(fcal, lcal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg)
print(gb)
stopCluster(mc)
fname <- paste("./Results/", labs, "_gb.rds", sep = "")
saveRDS(gb, fname)
