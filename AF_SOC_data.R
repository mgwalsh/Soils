# Africa-wide SOC data setup 
# M. Walsh, September 2017

# Required packages
# install.packages(c("downloader","rgdal","raster","caret")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(caret)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("AF_SOC", showWarnings=F)
setwd("./AF_SOC")

# download SOC data
download("https://www.dropbox.com/s/ubc00ukm5k433cb/SOC.csv.zip?raw=1", "SOC.csv.zip", mode="wb")
unzip("SOC.csv.zip", overwrite=T)
soc <- read.table("SOC.csv", header=T, sep=",")
soc <- soc[!(soc$lon==0 & soc$lat==0),] ## delete non-georeferenced profiles

# download Africa Gtifs and stack in raster (note this is a big 1Gb+ download)
download("https://www.dropbox.com/s/8kw4jitwp1n1bmc/AF_test_grids.zip?raw=1", "AF_test_grids.zip", mode="wb")
unzip("AF_test_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Data setup ---------------------------------------------------------------
# project SOC coords to grid CRS
soc.proj <- as.data.frame(project(cbind(soc$lon, soc$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(soc.proj) <- c("x","y")
soc <- cbind(soc, soc.proj)
coordinates(soc) <- ~x+y
projection(soc) <- projection(grids)

# extract gridded variables at profile locations
socgrid <- extract(grids, soc)
soc <- as.data.frame(cbind(soc, socgrid))
write.csv(soc, "socdat.csv", row.names = FALSE)
