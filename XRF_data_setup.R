#' Geochemical weathering index model setup with AfSIS-1 XRF data
#' XRF data from 60 sentinel sites
#' M. Walsh & J. Chen, January 2016

# install.packages(c("downloader","MASS","colorRamps","RColorBrewer","compositions","archetypes"), dependencies=T)
require(downloader)
require(MASS)
require(colorRamps)
require(RColorBrewer)
require(compositions)
require(archetypes)

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

# Calculate equivalent oxide wt% units and weathering index
attach(xrfd)
xrfd$wA  <- (Al*1.8895)/10000
xrfd$wCN <- (Ca*1.2046)/10000 + (Na*1.3992)/10000
xrfd$wK  <- (K*1.3480)/10000
detach(xrfd)

# Compositional analysis setup
vars <- c("SSN","Site","Lat","Lon","Depth","wA","wCN","wK")
wi60 <- na.omit(xrfd[vars])
cpart <- c("wCN","wK","wA")

# Sequential binary partion & isometric log ratio (ilr) transform
cdata <- acomp(wi60[cpart])
colnames(cdata) <- c("CN","K","A")

# Binary partition
bpart <- t(matrix(c(-1,-1, 1,
                    -1, 1, 0), ncol=3, nrow=2, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart, type="lines") ## compositional balance mobile graph				
idata <- as.data.frame(ilr(cdata, V=bpart))
wi60 <- cbind(wi60, idata)

# Identification of archetypes / endmembers -------------------------------
set.seed(85321)
wi60.arc <- stepArchetypes(data=wi60[,9:10], k=1:4, nrep=5, verbose=T)
screeplot(wi60.arc)

# Select no. of archetypes, screeplot suggests 3 endmembers
wi60.arc3 <- bestModel(wi60.arc[[3]])
idata <- as.data.frame(parameters(wi60.arc3))

# Inverse ilr transform 
ilrInv_nonorthonormal <- function(idata, V){
  n <- dim(idata)[2]
  LT <- (cbind(diag(n), rep(-1, n)))%*%V
  cdata_tmp <- as.matrix(idata)%*%solve(LT)
  cdata_tmp <- cbind(cdata_tmp, -rowSums(cdata_tmp))
  orig_data <- clrInv(cdata_tmp)
  return(orig_data)
}
arch <- ilrInv_nonorthonormal(idata, bpart)
arch <- as.data.frame(arch)
colnames(arch) <- c("CN","K","A")
carc <- acomp(arch)

# A-CN-K ternary plot
plot(cdata, cex=0.5, col="grey", center=T)
crust <- c(15.4,3.59+3.57,2.8) ## average composition of upper crust reference point
crust <- acomp(crust)
plot(crust, cex=1.3, col="red", add=T)
plot(carc, cex=1.3, col="blue", add=T)

# Plot of endmember compositions
mypal <- brewer.pal(3,"Blues")
barplot(carc, xlim=c(0,5), ylab=list("Proportion", cex=1), col=mypal, legend.text=T)

