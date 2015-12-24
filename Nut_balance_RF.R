#' Random Forest predictions of soil nutrient mass balance variables with HSTXT MIR data
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","doParallel,"randomForest")), dependencies=TRUE)
require(devtools)
require(caret)
require(doParallel)
require(randomForest)

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

# RFO models --------------------------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "oob", allowParallel = TRUE)
tg <- expand.grid(mtry=seq(20, 200, by=10))

# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
V1.rfo <- train(HSTXTc, V1,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(V1.rfo)
# v1.imp <- varImp(V1.rfo, useModel = FALSE) ## uncomment these if needed
# plot(v1.imp, top=20)

# V2 = ilr [P,K,S,Ca,Mg | C,N]
V2.rfo <- train(HSTXTc, V2,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(V2.rfo)
# v2.imp <- varImp(V2.rfo, useModel = FALSE)
# plot(v2.imp, top=20)

# V3 = ilr [P,S | K,Ca,Mg]
V3.rfo <- train(HSTXTc, V3,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(V3.rfo)
# v3.imp <- varImp(V3.rfo, useModel = FALSE)
# plot(v3.imp, top=20)

# V4 = ilr [K | Ca,Mg]
V4.rfo <- train(HSTXTc, V4,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(V4.rfo)
# v4.imp <- varImp(V4.rfo, useModel = FALSE)
# plot(v4.imp, top=20)

# V5 = ilr [P | S]
V5.rfo <- train(HSTXTc, V5,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(V5.rfo)
# v5.imp <- varImp(V5.rfo, useModel = FALSE)
# plot(v5.imp, top=20)

# V6 = ilr [Ca | Mg]
V6.rfo <- train(HSTXTc, V6,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(V6.rfo)
# v6.imp <- varImp(V6.rfo, useModel = FALSE)
# plot(v6.imp, top=20)

# V7 = ilr [C | N]
V7.rfo <- train(HSTXTc, V7,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(V7.rfo)
# v7.imp <- varImp(V7.rfo, useModel = FALSE)
# plot(v7.imp, top=20)

# Stop doParallel
stopCluster(mc)

# Test set predictions ----------------------------------------------------
V1_rfo <- predict(V1.rfo, HSTXTv)
V2_rfo <- predict(V2.rfo, HSTXTv)
V3_rfo <- predict(V3.rfo, HSTXTv)
V4_rfo <- predict(V4.rfo, HSTXTv)
V5_rfo <- predict(V5.rfo, HSTXTv)
V6_rfo <- predict(V6.rfo, HSTXTv)
V7_rfo <- predict(V7.rfo, HSTXTv)
pred <- cbind.data.frame(V1_rfo,V2_rfo,V3_rfo,V4_rfo,V5_rfo,V6_rfo,V7_rfo)
test <- nb60_val[c("SSN","V1","V2","V3","V4","V5","V6","V7")]
rfo_eval <- cbind(test, pred)

# Write data files --------------------------------------------------------
write.csv(rfo_eval, "RF_pred.csv", row.names=F)
