#' Multivariate soil nutrient benchmarks & archetype classification,
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

# Download
download("https://www.dropbox.com/s/wxmd9cx5m9h5b4r/Nutbal_60.csv.zip?dl=0", "Nutbal_60.csv.zip", mode="wb")
unzip("Nutbal_60.csv.zip", overwrite=T)
dat <- read.table("Nutbal_60.csv", header=T, sep=",")

# Compositional analysis setup
vars <- c("Site","Lat","Lon","Depth","C","N","P","K","S","Ca","Mg")
nb60 <- na.omit(dat[vars])
fpart <- c("C","N","P","K","S","Ca","Mg")
nb60$Fv <- 1000000-rowSums(nb60[nuts]) ## calculates "fill value" (Fv), in mg/kg soil

# Sequential binary partion & integrated log ratio (ilr) transform
cdata <- acomp(nb60[fpart])
bpart <- t(matrix(c(-1,-1, 1, 1, 1, 1, 1,
                     0, 0, 1,-1, 1,-1,-1,
                     0, 0, 0, 1, 0,-1,-1,
                     0, 0, 1, 0,-1, 0, 0,
                     0, 0, 0, 0, 0, 1,-1,
                     1,-1, 0, 0, 0, 0, 0), ncol=7, nrow=6, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart)  				
idata <- as.data.frame(ilr(cdata, V=bpart))

# Identify compositional archetypes ---------------------------------------
set.seed(1235813)
idata.arc <- stepArchetypes(data=idata, k=1:10, nrep=10, verbose=T)
screeplot(idata.arc)

# Select no. of archetypes (screeplot suggests k=7 archetypes)
idata.arc7 <- bestModel(idata.arc[[7]])
invarc <- ilrInv(parameters(idata.arc7), V=bpart)

# Classify ilr's by "dominant archetypes" (DA's)
arc <- predict(idata.arc7, idata)
DA <- apply(arc, 1, which.max)
nb60 <- cbind(nb60, idata, DA)
write.csv(nb60, "nb60_comp.csv", row.names=F)








