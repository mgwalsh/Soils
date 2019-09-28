# Alpha ZnSe spectrometer recalibration data setup
# pH, EC, Hp, C, N and M3 data from 60 AfSIS1 sentinel sites and TZ & GH AfSIS2 sites
# M. Walsh, September 2019

# install.packages(c("downloader","compositions","archetypes","RColorBrewer"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(compositions)
  require(archetypes)
  require(RColorBrewer)
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
alpha <- read.table("alpha.csv", header=T, sep=",") ## alpha ZnSe spectral data

# Compositional data analysis setup ---------------------------------------
vars <- c("C","N","P","S","K","Ca","Mg","Na","Fe","Mn","Zn","Cu")
nbal <- na.omit(wet[vars])
nbal$Fv <- 1000000-rowSums(nbal[vars]) ## calculates "fill value" (Fv), in mg/kg soil

