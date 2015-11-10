#' Soil nutrient mass balances with AfSIS-1 data:
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, Oct. 2015

# install.packages(c("downloader","compositions","arm","rgdal"), dependencies=T)
require(downloader)
require(compositions)
require(arm)
require(rgdal)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("NB60_data", showWarnings=F)
setwd("./NB60_data")

# Download
download("https://www.dropbox.com/s/k9pti8a4fvaxjlm/Nutbal60.zip?dl=0", "Nutbal60.zip", mode="wb")
unzip("Nutbal60.zip", overwrite=T)
prof <- read.table("Profiles.csv", header=T, sep=",") ## profile locations and site names
samp <- read.table("Samples.csv", header=T, sep=",") ## sample data
dat <- merge(prof, samp, by="PID")

# Compositional analysis setup
vars <- c("SSN","Site","Lat","Lon","Depth","C","N","P","K","S","Ca","Mg")
nb60 <- na.omit(dat[vars])
fpart <- c("C","N","P","K","S","Ca","Mg") ## all values in mg/kg
nb60$Fv <- 1000000-rowSums(nb60[fpart]) ## calculates "fill value" (Fv), in mg/kg soil
cpart <- c("C","N","P","K","S","Ca","Mg","Fv")

# Sequential binary partion & integrated log ratio (ilr) transform
cdata <- acomp(nb60[cpart])
bpart <- t(matrix(c( 1, 1, 1, 1, 1, 1, 1,-1,
                    -1,-1, 1, 1, 1, 1, 1, 0,
                     0, 0, 1,-1, 1,-1,-1, 0,
                     0, 0, 0, 1, 0,-1,-1, 0,
                     0, 0, 1, 0,-1, 0, 0, 0,
                     0, 0, 0, 0, 0, 1,-1, 0,
                     1,-1, 0, 0, 0, 0, 0, 0), ncol=8, nrow=7, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart) ## mass balance mobile graph				
idata <- as.data.frame(ilr(cdata, V=bpart))
nb60 <- cbind(nb60, idata)

# Topsoil / subsoil contrast ----------------------------------------------
tsc.glmer <- glmer(factor(Depth)~V1+V2+V4+V5+(1|Site), family="binomial"(link=logit), data=nb60)
summary(tsc.glmer)
tsc.ranef <- ranef(tsc.glmer)
tsc.se <- se.coef(tsc.glmer)
coefplot(tsc.ranef$Site[,1], tsc.se$Site[,1], varnames=rownames(tsc.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="SFI")
fix <- fixef(tsc.glmer) ## extract mean effects

# Proposed soil nutrient balance / fertility index (SFI), based on topsoil/subsoil contrast 
attach(nb60)
nb60$SFI <- (V1*fix[2]+V2*fix[3]+V4*fix[4]+V5*fix[5]+fix[1])*-1
detach(nb60)

# Topsoil / subsoil (SFI) contrast ecdf plot
top <- subset(nb60, Depth==10, select=c(V1,V2,V3,V4,V5,V6,V7,SFI))
quantile(top$SFI) ## value above the 50% topsoil SFI quantile ~ high fertility soils
sub <- subset(nb60, Depth==35, select=c(SFI))
plot(ecdf(top$SFI), main="", xlab="SFI", ylab="Cum. proportion of observations", xlim=c(-4,4), verticals=T, lty=1, lwd=2, col="red", do.points=F)
abline(0.5,0, lty=2, col="grey")
plot(ecdf(sub$SFI), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)

# Add grid coordinates ----------------------------------------------------
# Project to Africa LAEA from LonLat
nb60.laea <- as.data.frame(project(cbind(nb60$Lon, nb60$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(nb60.laea) <- c("x","y")
nb60 <- cbind(nb60.laea, nb60)

# Generate AfSIS 100m resolution grid cell ID's (GID)
res.pixel <- 100 ## set pixel resolution in m
xgid <- ceiling(abs(nb60$x)/res.pixel)
ygid <- ceiling(abs(nb60$y)/res.pixel)
gidx <- ifelse(nb60$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(nb60$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
nb60.gid <- cbind(GID, nb60)

# Write data files --------------------------------------------------------
write.csv(nb60.gid, "nb60.csv", row.names=F)

# nb60 locations in LonLat
nb60_loc <- subset(nb60.gid, Depth==10, select=c(Lat,Lon))
write.csv(nb60_loc, "nb60_loc.csv", row.names=F)
