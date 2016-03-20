#' Soil nutrient mass balance set-up with AfSIS1 data:
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# install.packages(c("downloader","compositions","colorRamps","RColorBrewer","MASS","rgdal"), dependencies=T)
suppressPackageStartupMessages({
require(downloader)
require(compositions)
require(MASS)
require(colorRamps)
require(RColorBrewer)
require(rgdal)
})

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("NB60_data", showWarnings=F)
setwd("./NB60_data")

# Download
download("https://www.dropbox.com/s/k9pti8a4fvaxjlm/Nutbal60.zip?dl=0", "Nutbal60.zip", mode="wb")
unzip("Nutbal60.zip", overwrite=T)
prof <- read.table("Profiles.csv", header=T, sep=",") ## profile locations and site names
samp <- read.table("Samples.csv", header=T, sep=",") ## sample data
geos <- read.table("nb60_GS.csv", header=T, sep=",") ## GeoSurvey data
dat <- merge(prof, samp, by="PID")

# load Mehlich-3 Al,B,Cu,Fe,Mn & Zn data
download("https://www.dropbox.com/s/bivkvxrjno8fo67/nb60_micro.csv?dl=0", "nb60_micro.csv", mode="wb")
mic <- read.table("nb60_micro.csv", header=T, sep=",")
dat <- merge(dat, mic, by="SSN")

# Parallel coordinates plot of the raw data
xvars <- c("B","C","N","Mg","P","S","K","Ca","Mn","Fe","Cu","Zn","Depth")
xdata <- dat[xvars]
xdata$Depth <- as.factor(xdata$Depth)
k <- adjustcolor(brewer.pal(3, "Set1")[xdata$Depth], alpha=.4)
parcoord(xdata[,1:12], col = k)

# Compositional data analysis setup
vars <- c("SSN","Site","Lat","Lon","Depth","B","C","N","Mg","P","S","K","Ca","Mn","Fe","Cu","Zn")
qdat <- na.omit(dat[vars])

# Fill value / subcomposition calculations
cvar <- c("B","C","N","Mg","P","S","K","Ca","Mn","Fe","Cu","Zn") ## all values in mg/kg
ivar <- c("C","N","P","K","S","Ca","Mg") ## all values in mg/kg
qdat$Fc <- 1000000-rowSums(qdat[cvar]) ## calculates "fill value" (Fc), in mg/kg soil
qdat$Fi <- 1000000-rowSums(qdat[ivar]) ## calculates "fill value" (Fi), in mg/kg soil

# Log ratio transforms and sequential binary partion
# clr transform
cvar <- c("B","C","N","Mg","P","S","K","Ca","Mn","Fe","Cu","Zn","Fc") ## all values in mg/kg
cdat <- acomp(qdat[cvar])
clrt <- as.data.frame(clr(cdat)) ## centered log ratio (clr) transform

# ilr transform
ivar <- c("C","N","P","K","S","Ca","Mg","Fi") ## all values in mg/kg 
idat <- acomp(qdat[ivar])
bpart <- t(matrix(c( 1, 1, 1, 1, 1, 1, 1,-1,
                    -1,-1, 1, 1, 1, 1, 1, 0,
                     0, 0, 1,-1, 1,-1,-1, 0,
                     0, 0, 0, 1, 0,-1,-1, 0,
                     0, 0, 1, 0,-1, 0, 0, 0,
                     0, 0, 0, 0, 0, 1,-1, 0,
                     1,-1, 0, 0, 0, 0, 0, 0), ncol=8, nrow=7, byrow=T))
CoDaDendrogram(X=acomp(idat), signary=bpart, type="lines") ## mass balance mobile graph				
ilrt <- as.data.frame(ilr(idat, V=bpart)) ## isometric log ratio (ilr) transform

# Assemble nurtrient composition dataframe
cdat <- cbind(clrt, ilrt)
vars <- c("SSN","Site","Lat","Lon","Depth")
nb60 <- cbind(qdat[vars], cdat)

# Project to Africa LAEA from LonLat
nb60.laea <- as.data.frame(project(cbind(nb60$Lon, nb60$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(nb60.laea) <- c("x","y")
nb60 <- cbind(nb60.laea, nb60)

# Add AfSIS grid coordinates ----------------------------------------------
# Generate AfSIS 100m resolution grid cell ID's (GID)
res.pixel <- 100 ## set pixel resolution in m
xgid <- ceiling(abs(nb60$x)/res.pixel)
ygid <- ceiling(abs(nb60$y)/res.pixel)
gidx <- ifelse(nb60$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(nb60$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
nb60 <- cbind(GID, nb60)

# Project GeoSurvey observations to Africa LAEA from LonLat
geos.laea <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.laea) <- c("x","y")
geos <- cbind(geos.laea, geos)

# Generate GeoSurvey grid cell ID's (GID) & add in GeoSurvey results
xgid <- ceiling(abs(geos$x)/res.pixel)
ygid <- ceiling(abs(geos$y)/res.pixel)
gidx <- ifelse(geos$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(geos$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
geos <- cbind(GID, geos)
geos <- geos[c(1,6:8)]

# Merge nb60 w Geosurvey observations by GID
nb60 <- merge(nb60, geos, by="GID")

# Train/Test set partition ------------------------------------------------
sites <- names(table(nb60$Site))
set.seed(1385321)
train <- sample(sites, 0.8*length(sites)) ## sample 80% of sites for calibration
nb60_cal <- nb60[ nb60$Site%in%train, ] ## calibration data
nb60_val <- nb60[!nb60$Site%in%train, ] ## validation data

# Load and merge HSTXT MIR spectra ----------------------------------------
download("https://www.dropbox.com/s/6fvipxqlmg704g3/hstxt_MIR.csv.zip?dl=0", "hstxt_MIR.csv.zip", mode="wb")
unzip("hstxt_MIR.csv.zip", overwrite=T)
mir <- read.table("hstxt_MIR.csv", header=F, sep=",", stringsAsFactors=F)
mir <- as.data.frame(mir)
names(mir) <- c("SSN", paste("m", signif(seq(7497.964, 599.76, length.out=3578), 6), sep=""))
nb60_cal <- merge(nb60_cal, mir, by="SSN")
nb60_cal <- na.omit(nb60_cal)
nb60_val <- merge(nb60_val, mir, by="SSN")
nb60_val <- na.omit(nb60_val)

# Write data files --------------------------------------------------------
write.csv(nb60, "nb60_dat.csv", row.names=F)
write.csv(nb60_cal, "nb60_cal.csv", row.names=F)
write.csv(nb60_val, "nb60_val.csv", row.names=F)
