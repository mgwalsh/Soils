#' Soil nutrient mass balances with AfSIS-1 data:
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, December 2015

# install.packages(c("downloader","compositions","rgdal"), dependencies=T)
require(downloader)
require(compositions)
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
geos <- read.table("nb60_GS.csv", header=T, sep=",") ## GeoSurvey data
dat <- merge(prof, samp, by="PID")

# Compositional analysis setup
vars <- c("SSN","Site","Lat","Lon","Depth","C","N","P","K","S","Ca","Mg")
nb60 <- na.omit(dat[vars])
fpart <- c("C","N","P","K","S","Ca","Mg") ## all values in mg/kg
nb60$Fv <- 1000000-rowSums(nb60[fpart]) ## calculates "fill value" (Fv), in mg/kg soil
cpart <- c("C","N","P","K","S","Ca","Mg","Fv")

# Sequential binary partion & isometric log ratio (ilr) transform
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
nb60 <- cbind(GID, nb60)

# Train/Test set partition ------------------------------------------------
sites <- names(table(nb60$Site))
set.seed(1385321)
train <- sample(sites, 0.8*length(sites))
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
write.csv(nb60_cal, "nb60_cal.csv", row.names=F)
write.csv(nb60_val, "nb60_val.csv", row.names=F)
