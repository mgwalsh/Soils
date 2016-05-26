#' PLS regression predictions of nutrient mass balance variables with HSTXT MIR,
#' depth in profile and GeoSurvey data
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# Required packages
# install.packages(c("devtools","caret","doParallel","pls")), dependencies=TRUE)
suppressPackageStartupMessages({
require(devtools)
require(compositions)
require(caret)
require(doParallel)
require(pls)
})

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

# Covariates
HSTXTc <- nb60_cal[c(8,25:26,27:3604)] ## Depth, CP, WP and HSTXT spectra
HSTXTv <- nb60_val[c(8,25:26,27:3604)] ## same for 12 randomly selected validation sites

# PLS models --------------------------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", number = 10, repeats = 3, allowParallel = TRUE)

# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
V1.pls <- train(HSTXTc, V1,
                preProc = c("center", "scale"),
                method = "pls",
                tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                trControl = tc)
print(V1.pls)
v1.imp <- varImp(V1.pls)
plot(v1.imp, top=20)

# V2 = ilr [P,K,S,Ca,Mg | C,N]
V2.pls <- train(HSTXTc, V2,
                preProc = c("center", "scale"),
                method = "pls",
                tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                trControl = tc)
print(V2.pls)
v2.imp <- varImp(V2.pls)
plot(v2.imp, top=20)

# V3 = ilr [P,S | K,Ca,Mg]
V3.pls <- train(HSTXTc, V3,
                preProc = c("center", "scale"),
                method = "pls",
                tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                trControl = tc)
print(V3.pls)
v3.imp <- varImp(V3.pls)
plot(v3.imp, top=20)

# V4 = ilr [K | Ca,Mg]
V4.pls <- train(HSTXTc, V4,
                preProc = c("center", "scale"),
                method = "pls",
                tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                trControl = tc)
print(V4.pls)
v4.imp <- varImp(V4.pls)
plot(v4.imp, top=20)

# V5 = ilr [P | S]
V5.pls <- train(HSTXTc, V5,
                preProc = c("center", "scale"),
                method = "pls",
                tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                trControl = tc)
print(V5.pls)
v5.imp <- varImp(V5.pls)
plot(v5.imp, top=20)

# V6 = ilr [Ca | Mg]
V6.pls <- train(HSTXTc, V6,
                preProc = c("center", "scale"),
                method = "pls",
                tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                trControl = tc)
print(V6.pls)
v6.imp <- varImp(V6.pls)
plot(v6.imp, top=20)

# V7 = ilr [C | N]
V7.pls <- train(HSTXTc, V7,
                preProc = c("center", "scale"),
                method = "pls",
                tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                trControl = tc)
print(V7.pls)
v7.imp <- varImp(V7.pls)
plot(v7.imp, top=20)

# Stop doParallel
stopCluster(mc)

# Test set predictions ----------------------------------------------------
V1_pls <- predict(V1.pls, HSTXTv)
V2_pls <- predict(V2.pls, HSTXTv)
V3_pls <- predict(V3.pls, HSTXTv)
V4_pls <- predict(V4.pls, HSTXTv)
V5_pls <- predict(V5.pls, HSTXTv)
V6_pls <- predict(V6.pls, HSTXTv)
V7_pls <- predict(V7.pls, HSTXTv)
pred <- cbind.data.frame(V1_pls,V2_pls,V3_pls,V4_pls,V5_pls,V6_pls,V7_pls)
test <- nb60_val[c("SSN","V1","V2","V3","V4","V5","V6","V7")]
pls_eval <- cbind(test, pred)

# Validation plots --------------------------------------------------------
# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
plot(V1~V1_pls, pls_eval, xlim=c(-19, -11), ylim=c(-19, -11), xlab="V1 predicted", ylab="V1 measured")
abline(c(0,1), col="red", lwd=2)

# V2 = ilr [P,K,S,Ca,Mg | C,N]
plot(V2~V2_pls, pls_eval, xlim=c(-20, -6), ylim=c(-20, -6), xlab="V2 predicted", ylab="V2 measured")
abline(c(0,1), col="red", lwd=2)

# Write data files --------------------------------------------------------
write.csv(pls_eval, "PLS_pred.csv", row.names=F)
