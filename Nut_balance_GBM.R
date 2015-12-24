#' Boosted regression predictions of nutrient mass balance variables with HSTXT MIR data
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","doParallel","gbm")), dependencies=TRUE)
require(devtools)
require(caret)
require(doParallel)
require(gbm)
require(plyr)

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Soils/blob/master/Nut_balance_setup.R
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Soils/master/Nut_balance_setup.R"
# source_url(SourceURL)

# Mass balance variables
V1 <- nb60_cal$V1
V2 <- nb60_cal$V2
V3 <- nb60_cal$V3
V4 <- nb60_cal$V4
V5 <- nb60_cal$V5
V6 <- nb60_cal$V6
V7 <- nb60_cal$V7

# Spectral covariates
HSTXTc <- nb60_cal[c(8,24:3601)] ## Depth in profile plus HSTXT spectra
HSTXTv <- nb60_val[c(8,24:3601)] ## same for 12 randomly selected validation sites

# GBM models --------------------------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", number = 10, allowParallel = TRUE)

# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
V1.gbm <- train(HSTXTc, V1, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 5,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 10))
print(V1.gbm)
v1.imp <- varImp(V1.gbm)
plot(v1.imp, top=20)

# V2 = ilr [P,K,S,Ca,Mg | C,N]
V2.gbm <- train(HSTXTc, V2, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 5,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 10))
print(V2.gbm)
v2.imp <- varImp(V2.gbm)
plot(v2.imp, top=20)

# V3 = ilr [P,S | K,Ca,Mg]
V3.gbm <- train(HSTXTc, V3, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 5,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 10))
print(V3.gbm)
v3.imp <- varImp(V3.gbm)
plot(v3.imp, top=20)

# V4 = ilr [K | Ca,Mg]
V4.gbm <- train(HSTXTc, V4, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 5,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 10))
print(V4.gbm)
v4.imp <- varImp(V4.gbm)
plot(v4.imp, top=20)

# V5 = ilr [P | S]
V5.gbm <- train(HSTXTc, V5, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 5,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 10))
print(V5.gbm)
v5.imp <- varImp(V5.gbm)
plot(v5.imp, top=20)

# V6 = ilr [Ca | Mg]
V6.gbm <- train(HSTXTc, V6, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 5,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 10))
print(V6.gbm)
v6.imp <- varImp(V6.gbm)
plot(v6.imp, top=20)

# V7 = ilr [C | N]
V7.gbm <- train(HSTXTc, V7, 
                method = "gbm", 
                preProc = c("center", "scale"),
                trControl = tc,
                tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                       .interaction.depth = 5,
                                       .shrinkage = 0.1,
                                       .n.minobsinnode = 10))
print(V7.gbm)
v7.imp <- varImp(V7.gbm)
plot(v7.imp, top=20)

# Stop doParallel
stopCluster(mc)

# Test set predictions ----------------------------------------------------
V1_gbm <- predict(V1.gbm, HSTXTv)
V2_gbm <- predict(V2.gbm, HSTXTv)
V3_gbm <- predict(V3.gbm, HSTXTv)
V4_gbm <- predict(V4.gbm, HSTXTv)
V5_gbm <- predict(V5.gbm, HSTXTv)
V6_gbm <- predict(V6.gbm, HSTXTv)
V7_gbm <- predict(V7.gbm, HSTXTv)
pred <- cbind.data.frame(V1_gbm,V2_gbm,V3_gbm,V4_gbm,V5_gbm,V6_gbm,V7_gbm)
test <- nb60_val[c("SSN","V1","V2","V3","V4","V5","V6","V7")]
gbm_eval <- cbind(test, pred)

# Write data files --------------------------------------------------------
write.csv(gbm_eval, "GBM_pred.csv", row.names=F)
