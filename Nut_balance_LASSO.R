#' LASSO regression predictions of nutrient mass balance variables with HSTXT MIR data
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# Required packages
# install.packages(c("caret","glmnet")), dependencies=TRUE)
require(devtools)
require(caret)
require(glmnet)

# Data setup --------------------------------------------------------------
SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Soils/master/Nut_balance_specdata.R"
source_url(SourceURL)

V1 <- nb60_cal$V1
V2 <- nb60_cal$V2
V3 <- nb60_cal$V3
V4 <- nb60_cal$V4
V5 <- nb60_cal$V5
V6 <- nb60_cal$V6
V7 <- nb60_cal$V7
HSTXTc <- nb60_cal[c(8,28:3605)] ## Depth in profile plus HSTXT spectra
HSTXTv <- nb60_val[c(8,28:3605)] ## same for 12 randomly selected validation sites

# LASSO models ------------------------------------------------------------
set.seed(1385321)

# Cross-validation setup
tc <- trainControl(method = "cv", number = 10)
  
# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
V1.las <- train(HSTXTc, V1,
                preProc = c("center", "scale"),
                method = "glmnet",
                family = "gaussian",
                tuneGrid = expand.grid(.alpha=1,.lambda=seq(0,0.1,by=0.01)),
                trControl = tc)
print(V1.las)
v1.imp <- varImp(V1.las)
plot(v1.imp, top=20)

# V2 = ilr [P,K,S,Ca,Mg | C,N]
V2.las <- train(HSTXTc, V2,
                preProc = c("center", "scale"),
                method = "glmnet",
                family = "gaussian",
                tuneGrid = expand.grid(.alpha=1,.lambda=seq(0,0.1,by=0.01)),
                trControl = tc)
print(V2.las)
v2.imp <- varImp(V2.las)
plot(v2.imp, top=20)

# V3 = ilr [P,S | K,Ca,Mg]
V3.las <- train(HSTXTc, V3,
                preProc = c("center", "scale"),
                method = "glmnet",
                family = "gaussian",
                tuneGrid = expand.grid(.alpha=1,.lambda=seq(0,0.1,by=0.01)),
                trControl = tc)
print(V3.las)
v3.imp <- varImp(V3.las)
plot(v3.imp, top=20)

# V4 = ilr [K | Ca,Mg]
V4.las <- train(HSTXTc, V4,
                preProc = c("center", "scale"),
                method = "glmnet",
                family = "gaussian",
                tuneGrid = expand.grid(.alpha=1,.lambda=seq(0,0.1,by=0.01)),
                trControl = tc)
print(V4.las)
v4.imp <- varImp(V4.las)
plot(v4.imp, top=20)

# V5 = ilr [P | S]
V5.las <- train(HSTXTc, V5,
                preProc = c("center", "scale"),
                method = "glmnet",
                family = "gaussian",
                tuneGrid = expand.grid(.alpha=1,.lambda=seq(0,0.1,by=0.01)),
                trControl = tc)
print(V5.las)
v5.imp <- varImp(V5.las)
plot(v5.imp, top=20)

# V6 = ilr [Ca | Mg]
V6.las <- train(HSTXTc, V6,
                preProc = c("center", "scale"),
                method = "glmnet",
                family = "gaussian",
                tuneGrid = expand.grid(.alpha=1,.lambda=seq(0,0.1,by=0.01)),
                trControl = tc)
print(V6.las)
v6.imp <- varImp(V6.las)
plot(v6.imp, top=20)

# V7 = ilr [C | N]
V7.las <- train(HSTXTc, V7,
                preProc = c("center", "scale"),
                method = "glmnet",
                family = "gaussian",
                tuneGrid = expand.grid(.alpha=1,.lambda=seq(0,0.1,by=0.01)),
                trControl = tc)
print(V7.las)
v7.imp <- varImp(V7.las)
plot(v7.imp, top=20)

# Test set predictions ----------------------------------------------------
V1p <- predict(V1.las, HSTXTv)
V2p <- predict(V2.las, HSTXTv)
V3p <- predict(V3.las, HSTXTv)
V4p <- predict(V4.las, HSTXTv)
V5p <- predict(V5.las, HSTXTv)
V6p <- predict(V6.las, HSTXTv)
V7p <- predict(V7.las, HSTXTv)
rreg <- cbind.data.frame(V1p,V2p,V3p,V4p,V5p,V6p,V7p)
test <- nb60_val[c("SSN","V1","V2","V3","V4","V5","V6","V7")]
pred <- cbind(test, rreg)

# Write data files --------------------------------------------------------
write.csv(pred, "LASSO_pred.csv", row.names=F)




