#' Multivariate soil nutrient balance benchmarks & archetype classification
#' with AfSIS-1 data: C,N and Mehlich-3 P,K,S,Ca & Mg from 60 sentinel sites
#' M. Walsh, Oct. 2015

# install.packages(c("rgdal","compositions","archetypes","arm"), dependencies=T)
require(downloader)
require(rgdal)
require(compositions)
require(archetypes)
require(arm)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("NB60_data", showWarnings=F)
setwd("./NB60_data")

# Data download
download("https://www.dropbox.com/s/wxmd9cx5m9h5b4r/Nutbal_60.csv.zip?dl=0", "Nutbal_60.csv.zip", mode="wb")
unzip("Nutbal_60.csv.zip", overwrite=T)
dat <- read.table("Nutbal_60.csv", header=T, sep=",")

# Compositional data analysis setup
vars <- c("Site","Lat","Lon","Depth","C","N","P","K","S","Ca","Mg")
nb60 <- na.omit(dat[vars])
nuts <- c("C","N","P","K","S","Ca","Mg")
nb60$Fv <- 1000000-rowSums(nb60[nuts]) ## calculates "fill value" (Fv), in mg/kg soil

# Sequential binary partion & integrated log ratio (ilr) transform
fpart <- c("C","N","P","K","S","Ca","Mg","Fv")
cdata <- acomp(nb60[fpart])
bpart <- t(matrix(c( 1, 1, 1, 1, 1, 1, 1,-1,
                    -1,-1, 1, 1, 1, 1, 1, 0,
                     0, 0, 1,-1, 1,-1,-1, 0,
                     0, 0, 0, 1, 0,-1,-1, 0,
                     0, 0, 1, 0,-1, 0, 0, 0,
                     0, 0, 0, 0, 0, 1,-1, 0,
                     1,-1, 0, 0, 0, 0, 0, 0), ncol=8, nrow=7, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart)  				
idata <- as.data.frame(ilr(cdata, V=bpart))
nb60.comp <- cbind(nb60, idata)
write.csv(nb60.comp, "nb60_comp.csv", row.name = F)







