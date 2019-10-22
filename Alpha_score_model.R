# Alpha-MIR scoring rules for soil remediation decisions
# M. Walsh, October 2019

# install.packages(c("downloader","arm"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(arm)
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
c50 <- glmer(c50~mir+(1+mir|cpart), family=binomial(link="logit"), data=rals)
summary(c50)
rals$rscore <- fitted(c50)
write.csv(rals, "./Results/rals_2019.csv", row.names=F)

