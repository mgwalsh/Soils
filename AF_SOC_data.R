# Africa-wide SOC data setup and national estimates
# M. Walsh, September 2017

# Required packages
# install.packages(c("downloader","rgdal","raster","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(leaflet)
  require(htmlwidgets)
  require(arm)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("AF_SOC", showWarnings=F)
setwd("./AF_SOC")

# download SOC profile data (courtesy ISRIC)
download("https://osf.io/u8wng?raw=1", "socdat.zip", mode="wb")
unzip("socdat.zip", overwrite=T)
prof <- read.table("Profiles.csv", header=T, sep=",") ## profile locations & collection year
prof <- prof[complete.cases(prof[ ,3:4]),] ## delete non-georeferenced profiles
samp <- read.table("Samples.csv", header=T, sep=",") ## sample data

# download Africa Gtifs and stack in raster (note this is a big 450Mb+ download)
download("https://www.dropbox.com/s/8kw4jitwp1n1bmc/AF_test_grids.zip?raw=1", "AF_test_grids.zip", mode="wb")
unzip("AF_test_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)
plot(grids, axes = FALSE)

# Data setup ---------------------------------------------------------------
# project SOC coords to grid CRS
prof.proj <- as.data.frame(project(cbind(prof$lon, prof$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(prof.proj) <- c("x","y")
prof <- cbind(prof, prof.proj)
coordinates(prof) <- ~x+y
projection(prof) <- projection(grids)

# extract gridded variables at profile locations
socgrid <- extract(grids, prof)
prof <- as.data.frame(cbind(prof, socgrid))
soc <- merge(prof, samp, by="pid") ## merge samples by profile ID (pid)
# socc <- soc[complete.cases(soc[ ,c(17,19:20)]),] ## delete cases with incomplete Sand, pH or SOC measurements
write.csv(soc, "socdat.csv", row.names = F)

# Profile location map widget ---------------------------------------------
# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(prof$lon, prof$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'SOC_locs.html', selfcontained = T)



