# Alpha ZnSe spectrometer recalibration data setup
# pH, EC, Hp, C, N and M3 data from 60 AfSIS1 sentinel sites and TZ & GH AfSIS2 sites
# M. Walsh, September 2019

# install.packages(c("downloader","compositions","archetypes","RColorBrewer"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(compositions)
  require(archetypes)
})
rm(list = ls())

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("alpha_data", showWarnings=F)
setwd("./alpha_data")

# Download
download("https://osf.io/ut7ya?raw=1", "alpha_ref_data.zip", mode="wb")
unzip("alpha_ref_data.zip", overwrite=T)
wet <- read.table("wet.csv", header=T, sep=",") ## pH, EC, Hp, C, N & M3 data
vars <- c("SSID","C","N","P","K","S","Ca","Mg","Na","Fe","Mn","Cu","Zn")
wet <- na.omit(wet[vars])
alpha <- read.table("alpha.csv", header=T, sep=",") ## alpha ZnSe spectral data

# Compositional data analysis setup ---------------------------------------
vars <- c("C","N","P","K","S","Ca","Mg","Na","Fe","Mn","Cu","Zn")
wet$Fv <- 1000000-rowSums(wet[vars]) ## calculates "fill value" (Fv), in mg/kg soil

# Centered log ratio (clr) transform
vars <- c("C","N","P","K","S","Ca","Mg","Na","Fe","Mn","Cu","Zn","Fv")
nbal <- wet[vars]
nbal <- as.data.frame(clr(nbal)) ## centered log ratio (clr) transform
nbal <- cbind(wet$SSID, nbal)
colnames(nbal)[colnames(nbal)=="wet$SSID"] <- "SSID"

# merge with spectral data
nbal <- merge(nbal, alpha, by="SSID")
write.csv(nbal, "nbal_2019.csv", row.names=F)

