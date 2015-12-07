#' PLS regression predictions of nutrient mass balance variables with HSTXT MIR data
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# Required packages
# install.packages(c("caret","glmnet")), dependencies=TRUE)
require(devtools)
require(caret)
require(pls)

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


# PLS models --------------------------------------------------------------
set.seed(1385321)

# Cross-validation setup
tc <- trainControl(method = "cv", number = 10)

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

# Test set predictions ----------------------------------------------------
V1p <- predict(V1.pls, HSTXTv)
V2p <- predict(V2.pls, HSTXTv)
V3p <- predict(V3.pls, HSTXTv)
V4p <- predict(V4.pls, HSTXTv)
V5p <- predict(V5.pls, HSTXTv)
V6p <- predict(V6.pls, HSTXTv)
V7p <- predict(V7.pls, HSTXTv)
pred <- cbind.data.frame(V1p,V2p,V3p,V4p,V5p,V6p,V7p)
test <- nb60_val[c("SSN","V1","V2","V3","V4","V5","V6","V7")]
eval <- cbind(test, pred)

# Write data files --------------------------------------------------------
write.csv(eval, "PLS_pred.csv", row.names=F)
