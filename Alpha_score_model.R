# Alpha-MIR scoring rules for soil remediation decisions
# M. Walsh, October 2019

# install.packages(c("downloader","arm","dismo"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(arm)
  require(dismo)
})
rm(list = ls())

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("alpha_data", showWarnings=F)
setwd("./alpha_data")

# Download
download("https://osf.io/wyrup?raw=1", "alpha_long_preds.csv.zip", mode="wb")
unzip("alpha_long_preds.csv.zip", overwrite=T)
rals <- read.table("alpha_long_preds.csv", header=T, sep=",") ## RAL and MIR predictions

# Scoring rules -----------------------------------------------------------
y50 <- glmer(c50~mir+(1+mir|cpart), family=binomial(link="logit"), data=rals)
summary(y50)
rals$rscore <- fitted(y50)
write.csv(rals, "./Results/rals_2019.csv", row.names=F)

# Receiver-operator characteristics ---------------------------------------
p <- rals[ which(rals$cpart=="C" & rals$c50==1), ] ## substitute other properties here (N,P,K ...)
p <- p[,6]
a <- rals[ which(rals$cpart=="C" & rals$c50==0), ] ## substitute other properties here
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on validation set
plot(e, 'ROC') ## plot ROC curve

