#' Geochemical weathering index model setup with AfSIS-1 XRF data
#' XRF data from 60 sentinel sites
#' M. Walsh, January 2016

# install.packages(c("downloader","MASS","colorRamps","RColorBrewer","compositions"), dependencies=T)
require(downloader)
require(MASS)
require(colorRamps)
require(RColorBrewer)
require(compositions)

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

# Parallel coordinates plot of all of the XRF data
xdata <- xrfd[,8:48]
xdata$Depth <- as.factor(xdata$Depth)
k <- adjustcolor(brewer.pal(3, "Set1")[xdata$Depth], alpha=0.8)
parcoord(xdata[,2:21], col = k)
parcoord(xdata[,22:41], col = k)

# Calculate equivalent oxide wt% units and weathering indices
attach(xrfd)
xrfd$Al2O3 <- (Al*1.8895)/10000
xrfd$CaO   <- (Ca*1.2046)/10000
xrfd$Na2O  <- (Na*1.3992)/10000
xrfd$K2O   <- (K*1.3480)/10000
xrfd$CIA <- (Al*1.8895)/(Al*1.8895+K*1.2046+(Ca*1.3992-(P*2.2916*10/3))+Na*1.3480)*100
xrfd$CIW <- (Al*1.8895)/(Al*1.8895+(Ca*1.3992-(P*2.2916*10/3))+Na*1.3480)*100
detach(xrfd)
plot(CIW ~ CIA, data=xrfd, xlim=c(0,100), ylim=c(0,100))
abline(c(0,1), lwd=2, col="red")

# Compositional analysis setup
vars <- c("SSN","Site","Lat","Lon","Depth","Al2O3","CaO","Na2O","K2O","CIA","CIW")
wi60 <- na.omit(xrfd[vars])
cpart <- c("Al2O3","CaO","Na2O","K2O")

# Sequential binary partion & isometric log ratio (ilr) transform
cdata <- acomp(wi60[cpart])
bpart <- t(matrix(c(1,-1,-1,-1,
                    0,-1,-1, 1,
                    0, 1,-1, 0), ncol=4, nrow=3, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart, type="lines") ## mass balance mobile graph				
idata <- as.data.frame(ilr(cdata, V=bpart))
wi60 <- cbind(wi60, idata)
