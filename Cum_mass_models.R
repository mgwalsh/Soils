#' Example cumulative soil and soil organic carbon mass analyses with multilevel models
#' Data extracted from http://www.isric.org/data/africa-soil-profiles-database-version-01-1
#' M. Walsh, June 2015

# Load packages ------------------------------------------------------------
# install.packages(c("downloader","proj4","arm")), dependencies=TRUE)
require(downloader)
require(proj4)
require(arm)

#+ Data download -----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("AFSOC_data", showWarnings=F)
dat_dir <- "./AFSOC_data"

# AFSOC data
download("https://www.dropbox.com/s/pa0v8xg6vyuqzn4/AFSOC.zip?dl=0", "./AFSOC_data/AFSOC.zip", mode="wb")
unzip("./AFSOC_data/AFSOC.zip", exdir="./AFSOC_data", overwrite=T)
profile <- read.table(paste(dat_dir, "/Profiles.csv", sep=""), header=T, sep=",")
profile <- na.omit(profile)
samples <- read.table(paste(dat_dir, "/Horizons.csv", sep=""), header=T, sep=",")
samples <- na.omit(samples)

#+ Generate coordinate reference and GID's ---------------------------------
# Project profile coords to Africa LAEA from LonLat
profile.laea <- as.data.frame(project(cbind(profile$Lon, profile$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(profile.laea) <- c("x","y")
profile <- cbind(profile, profile.laea)

# Generate AfSIS grid cell ID's (GID)
# these are included to decluster profiles measured within a distance of 1 km of each other
res.pixel <- 1000
xgid <- ceiling(abs(profile$x)/res.pixel)
ygid <- ceiling(abs(profile$y)/res.pixel)
gidx <- ifelse(profile$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(profile$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
profile.gid <- cbind(profile, GID)
cmdat <- merge(profile.gid, samples, by="PID")

#+ Multilevel models ------------------------------------------------------
# Random intercept model for cumulative soil mass (CSM) with depth (Bot) in profile
# note Profile ID's (PID) variance components are nested within Grid Cell ID's (GID)
csmd.lmer <- lmer(log(CSM)~log(Bot)+(1|GID)+(1|PID:GID), data=cmdat)
display(csmd.lmer)
plot(log(CSM)~fitted(csmd.lmer), cmdat)

# Random intercept model for cumulative soil organic carbon mass (CCM) with depth (Bot) in profile
ccmd.lmer <- lmer(log(CCM)~log(Bot)+(1|GID)+(1|PID:GID), data=cmdat)
display(ccmd.lmer)
plot(log(CCM)~fitted(ccmd.lmer), cmdat)

# Random intercept model for cumulative soil organic carbon mass (CCM) 
# with cumulative soil mass (CSM) in profile
ccmm.lmer <- lmer(log(CCM)~log(CSM)+(1|GID)+(1|PID:GID), data=cmdat)
display(ccmm.lmer)
plot(log(CCM)~fitted(ccmm.lmer), cmdat)
