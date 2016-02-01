#' Geochemical enrichment/depletion factor models with AfSIS-1 XRF data:
#' XRF data from 60 sentinel sites
#' M. Walsh, January 2016

# install.packages(c("downloader","quantreg","MASS","colorRamps","RColorBrewer"), dependencies=T)
require(downloader)
require(quantreg)
require(MASS)
require(colorRamps)
require(RColorBrewer)

# Data setup --------------------------------------------------------------
# Create a data folder in  your current working directory
dir.create("XRF_data", showWarnings=F)
setwd("./XRF_data")

# Download
download("https://www.dropbox.com/s/6mr5ubca0jrhk25/XRF60.zip?dl=0", "XRF60.zip", mode="wb")
unzip("XRF60.zip", overwrite=T)
prof <- read.table("Profiles.csv", header=T, sep=",") ## profile locations and site names
samp <- read.table("Samples.csv", header=T, sep=",") ## sample ID's and depths
xrfd <- read.table("XRF.csv", header=T, sep=",") ## XRF data
samp <- merge(prof, samp, by="PID")
xrfd <- merge(samp, xrfd, by="SSN")

# Setup topsoil / subsoil contrast
top <- subset(xrfd, Depth==10, select=c(PID,Site,Lat,Lon,Al,P,K,S)) 
sub <- subset(xrfd, Depth==35, select=c(PID,Al,P,K,S))
del <- merge(top, sub, by="PID") ## merge by profile ID
colnames(del) <- c("PID","Site","Lat","Lon","Alt","Pt","Kt","St","Als","Ps","Ks","Ss")
del$PALt <- del$Pt/del$Alt
del$PALs <- del$Ps/del$Als
del$KALt <- del$Kt/del$Alt
del$KALs <- del$Ks/del$Als
del$SALt <- del$St/del$Alt
del$SALs <- del$Ss/del$Als

# Parallel coordinates plot of the raw data
xdata <- xrfd[,8:48]
xdata$Depth <- as.factor(xdata$Depth)
k <- adjustcolor(brewer.pal(3, "Set1")[xdata$Depth], alpha=0.8)
parcoord(xdata[,2:21], col = k)
parcoord(xdata[,22:41], col = k)

# Quantile regression benchmarks ------------------------------------------
# Predict P-quantile
tau <- 0.5 ## set quantile level
PAL.rq <- rq(log(PALt)~log(PALs), tau = tau, data=del)
del$PALp <- exp(predict(PAL.rq, del))
del$dPAL <- (del$PALt/del$PALp-1)*100

# Predict K-quantile
tau <- 0.5 ## set quantile level
KAL.rq <- rq(log(KALt)~log(KALs), tau = tau, data=del)
del$KALp <- exp(predict(KAL.rq, del))
del$dKAL <- (del$KALt/del$KALp-1)*100

# Predict S-quantile
tau <- 0.5 ## set quantile level
SAL.rq <- rq(log(SALt)~log(SALs), tau = tau, data=del)
del$SALp <- exp(predict(SAL.rq, del))
del$dSAL <- (del$SALt/del$SALp-1)*100
