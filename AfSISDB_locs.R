# AfSISDB publicly available soil sample locations 
# M. Walsh, February 2019

# Required packages
# install.packages(c("downloader","rgdal","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(leaflet)
  require(htmlwidgets)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("AFDB", showWarnings=F)
setwd("./AFDB")

# download profile location data 
download("https://www.dropbox.com/s/mavd0mts8ycpta5/afsisdb_locs.csv?raw=1", "afsisdb_locs.csv", mode="wb")
samp <- read.table("afsisdb_locs.csv", header=T, sep=",") ## sample locations

# Profile location map widget ---------------------------------------------
# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(samp$lon, samp$lat, clusterOptions = markerClusterOptions())
w ## plot widget 

# save widget
saveWidget(w, 'AFDB_soil_samples.html', selfcontained = T)
