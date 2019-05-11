# CEC imbalance data setup
# CEC & LDPSA data from 60 sentinel sites
# M. Walsh, May 2019

# install.packages(c("downloader","caret"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(caret)
})

# Data setup --------------------------------------------------------------
# Create a data folder in  your current working directory
dir.create("CEC_data", showWarnings=F)
setwd("./CEC_data")

# Download
download("https://www.dropbox.com/s/p9ql8dqxjg6x4pe/CEC_imbalances.zip?raw=1", "CEC_imbalances.zip", mode="wb")
unzip("CEC_imbalances.zip", overwrite=T)
pro <- read.table("profiles.csv", header=T, sep=",") ## profile locations and site names
cec <- read.table("cec.csv", header=T, sep=",") ## CEC data
soc <- read.table("soc.csv", header=T, sep=",") ## SOC data
psa <- read.table("psa.csv", header=T, sep=",") ## LDPSA data
samp <- merge(pro, cec, by="sid")
samp <- merge(samp, soc, by="sid")
samp <- merge(samp, psa, by="sid")
samp$w1k <- 7.594*(0.0034+0.0387*exp(-0.5*(log(samp$w1dg)+1.533/0.7671)^2))
samp$c4k <- 7.594*(0.0034+0.0387*exp(-0.5*(log(samp$c4dg)+1.533/0.7671)^2))
