# Alpha ZnSe spectrometer recalibration data setup
# pH, EC, Hp, C, N and M3 data from 60 AfSIS1 sentinel sites and TZ & GH AfSIS2 sites
# M. Walsh, September 2019

# install.packages(c("downloader","compositions"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(compositions)
})
rm(list = ls())

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("alpha_data", showWarnings=F)
setwd("./alpha_data")
dir.create("Results", showWarnings = F)

# Download
download("https://osf.io/ut7ya?raw=1", "alpha_ref_data.zip", mode="wb")
unzip("alpha_ref_data.zip", overwrite=T)
wet <- read.table("wet.csv", header=T, sep=",") ## pH, EC, Hp, C, N & M3 data
vars <- c("SSID","pH","Hp","Ca","Mg")
lreq <- na.omit(wet[vars])
lreq$Hpa <- ifelse(lreq$pH >= 7.0, 0.0, lreq$Hp*10) ## adjusts Hp to meq/100 gm for lime requirement calcs
lreq$camg <- lreq$Ca/lreq$Mg ## calculates Ca:Mg ratio based on Mehlich-3 measurements
vars <- c("SSID","C","N","P","K","S","Ca","Mg","Na","Fe","Mn","Cu","Zn")
wet <- na.omit(wet[vars])
alpha <- read.table("alpha.csv", header=T, sep=",") ## Alpha ZnSe spectral data

# Hp vs pH plot
par(pty="s")
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(I(Hp*10)~pH, xlab="pH (in water)", ylab="Hp (meq/kg)", cex.lab=1.3, 
     xlim=c(4,10), ylim=c(0,50), lreq)
abline(h=2.0, col="red", lwd=2)
abline(v=6.5, col="blue", lwd=2)

# Compositional data analysis setup ---------------------------------------
vars <- c("C","N","P","K","S","Ca","Mg","Na","Fe","Mn","Cu","Zn")
wet$Fv <- 1000000-rowSums(wet[vars]) ## calculates "fill value" (Fv), in mg/kg soil

# Centered log ratio (clr) transform
vars <- c("C","N","P","K","S","Ca","Mg","Na","Fe","Mn","Cu","Zn","Fv")
nbal <- wet[vars]
nbal <- as.data.frame(clr(nbal)) ## centered log ratio (clr) transform
nbal <- cbind(wet$SSID, nbal)
colnames(nbal)[colnames(nbal)=="wet$SSID"] <- "SSID"

# Alpha principal components ----------------------------------------------
alpha.pca <- prcomp(alpha[,2:1715], center=T, scale=T)
screeplot(alpha.pca)
pcas <- predict(alpha.pca, alpha)
pcas <- pcas[,1:20]

# Merge & write files -----------------------------------------------------
alpha <- cbind(alpha, pcas)
lreq <- merge(lreq, alpha, by="SSID")
nbal <- merge(nbal, alpha, by="SSID")
write.csv(lreq, "./Results/lreq_2019.csv", row.names=F)
write.csv(nbal, "./Results/nbal_2019.csv", row.names=F)

