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
# Fv
p <- rals[ which(rals$cpart=="Fv" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Fv" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# C
p <- rals[ which(rals$cpart=="C" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="C" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# Ca
p <- rals[ which(rals$cpart=="Ca" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Ca" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# N
p <- rals[ which(rals$cpart=="N" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="N" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# Mg
p <- rals[ which(rals$cpart=="Mg" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Mg" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# K
p <- rals[ which(rals$cpart=="K" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="K" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# Fe
p <- rals[ which(rals$cpart=="Fe" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Fe" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# Mn
p <- rals[ which(rals$cpart=="Mn" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Mn" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# Na
p <- rals[ which(rals$cpart=="Na" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Na" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# S
p <- rals[ which(rals$cpart=="S" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="S" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# P
p <- rals[ which(rals$cpart=="P" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="P" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# Cu
p <- rals[ which(rals$cpart=="Cu" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Cu" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

# Zn
p <- rals[ which(rals$cpart=="Zn" & rals$c50==1), ]
p <- p[,6]
a <- rals[ which(rals$cpart=="Zn" & rals$c50==0), ]
a <- a[,6]
e <- evaluate(p=p, a=a) ## calculate ROC's on test set
plot(e, 'ROC') ## plot ROC curve

