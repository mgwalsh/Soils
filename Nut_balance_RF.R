#' Random Forest predictions of nutrient mass balance variables with HSTXT MIR data
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# Required packages
# install.packages(c("caret","randomForest")), dependencies=TRUE)
require(devtools)
require(caret)
require(randomForest)

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

# Random Forest models ----------------------------------------------------
set.seed(1385321)

# Cross-validation setup
tc <- trainControl(method = "oob")

# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
V1.rfo <- train(HSTXTc, V1,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = expand.grid(mtry=seq(70,110,by=5)),
                trControl = tc)
print(V1.rfo)

# V2 = ilr [P,K,S,Ca,Mg | C,N]
V2.rfo <- train(HSTXTc, V2,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = expand.grid(mtry=seq(70,110,by=5)),
                trControl = tc)
print(V2.rfo)

# V3 = ilr [P,S | K,Ca,Mg]
V3.rfo <- train(HSTXTc, V3,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = expand.grid(mtry=seq(70,110,by=5)),
                trControl = tc)
print(V3.rfo)

# V4 = ilr [K | Ca,Mg]
V4.rfo <- train(HSTXTc, V4,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = expand.grid(mtry=seq(70,110,by=5)),
                trControl = tc)
print(V4.rfo)

# V5 = ilr [P | S]
V5.rfo <- train(HSTXTc, V5,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = expand.grid(mtry=seq(70,110,by=5)),
                trControl = tc)
print(V5.rfo)

# V6 = ilr [Ca | Mg]
V6.rfo <- train(HSTXTc, V6,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = expand.grid(mtry=seq(70,110,by=5)),
                trControl = tc)
print(V6.rfo)

# V7 = ilr [C | N]
V7.rfo <- train(HSTXTc, V7,
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = expand.grid(mtry=seq(70,110,by=5)),
                trControl = tc)
print(V7.rfo)

# Test set predictions ----------------------------------------------------
V1p <- predict(V1.rfo, HSTXTv)
V2p <- predict(V2.rfo, HSTXTv)
V3p <- predict(V3.rfo, HSTXTv)
V4p <- predict(V4.rfo, HSTXTv)
V5p <- predict(V5.rfo, HSTXTv)
V6p <- predict(V6.rfo, HSTXTv)
V7p <- predict(V7.rfo, HSTXTv)
pred <- cbind.data.frame(V1p,V2p,V3p,V4p,V5p,V6p,V7p)
test <- nb60_val[c("SSN","V1","V2","V3","V4","V5","V6","V7")]
eval <- cbind(test, pred)

# Write data files --------------------------------------------------------
write.csv(eval, "RF_pred.csv", row.names=F)
