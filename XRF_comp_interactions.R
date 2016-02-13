#' Compositional affinities between essential plant nutrients in soils
#' XRF data from 60 sentinel sites
#' M. Walsh, Februaray 2016

# install.packages(c("downloader","colorRamps","RColorBrewer","compositions","circlize"), dependencies=T)
require(downloader)
require(colorRamps)
require(RColorBrewer)
require(compositions)
require(circlize)

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

# Extract essential nutrients
vars <- c("Na","Mg","P","S","Cl","K","Ca","Mn","Fe","Co","Ni","Cu","Zn","Mo")
enut <- xrfd[vars]

# Calculate compositional correlations
edat <- as.data.frame(clr(acomp(enut)))
ecor <- (cor(edat))^2 ## affinities = sum of compositional coefficients of determination (R^2)
ecor[lower.tri(ecor, diag=F)] <- 0
ecor <- ifelse(ecor < 0.1 | ecor == 1.0, 0, ecor)

# Plot chord diagram of affinities ----------------------------------------
set.seed(12358)
circos.par(gap.degree = 2)
chordDiagram(ecor, directional = F, annotationTrack = "grid",
             preAllocateTracks = list(list(track.height = 0.05),
                                      list(track.height = 0.05)))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), mean(ylim), sector.index, facing = "clockwise", niceFacing = T)
}, bg.border = NA)
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  circos.axis("bottom", major.tick.percentage = 0.2, labels.cex = 0.6)
}, bg.border = NA)
circos.clear()
