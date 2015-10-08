#' Soil nutrient mass balances with AfSIS-1 data:
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, Oct. 2015

# install.packages(c("downloader","compositions","arm"), dependencies=T)
require(downloader)
require(compositions)
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
CoDaDendrogram(X=acomp(cdata), signary=bpart)  				
idata <- as.data.frame(ilr(cdata, V=bpart))
nb60 <- cbind(nb60, idata)
write.csv(nb60, "nb60_comp.csv", row.names=F)

# Site level summaries / expected values ----------------------------------
# V1 = ilr[C,N,P,K,S,Ca,Mg|Fv]
V1.lmer <- lmer(V1~I(Depth/100)+(1|Site), data=nb60)
summary(V1.lmer)

# Plot of site-level random effects and standard errors
V1.ranef <- ranef(V1.lmer)
V1.se <- se.coef(V1.lmer)
coefplot(V1.ranef$Site[,1], V1.se$Site[,1], varnames=rownames(V1.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [C,N,P,K,S,Ca,Mg | Fv]")

# V2 = ilr[P,K,S,Ca,Mg|C,N]
V2.lmer <- lmer(V2~I(Depth/100)+(1|Site), data=nb60)
summary(V2.lmer)
V2.ranef <- ranef(V2.lmer)
V2.se <- se.coef(V2.lmer)
coefplot(V2.ranef$Site[,1], V2.se$Site[,1], varnames=rownames(V2.ranef$Site), xlim=c(-7,7), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [P,K,S,Ca,Mg | C,N]")

# V3 = ilr[P,S|K,Ca,Mg]
V3.lmer <- lmer(V3~I(Depth/100)+(1|Site), data=nb60)
summary(V3.lmer)
V3.ranef <- ranef(V3.lmer)
V3.se <- se.coef(V3.lmer)
coefplot(V3.ranef$Site[,1], V3.se$Site[,1], varnames=rownames(V3.ranef$Site), xlim=c(-5,5), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [P,S | K,Ca,Mg]")

# V4 = ilr[K|Ca,Mg]
V4.lmer <- lmer(V4~I(Depth/100)+(1|Site), data=nb60)
summary(V4.lmer)
V4.ranef <- ranef(V4.lmer)
V4.se <- se.coef(V4.lmer)
coefplot(V4.ranef$Site[,1], V4.se$Site[,1], varnames=rownames(V4.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [K | Ca,Mg]")

# V5 = ilr[P|S]
V5.lmer <- lmer(V5~I(Depth/100)+(1|Site), data=nb60)
summary(V5.lmer)
V5.ranef <- ranef(V5.lmer)
V5.se <- se.coef(V5.lmer)
coefplot(V5.ranef$Site[,1], V5.se$Site[,1], varnames=rownames(V5.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [P | S]")

# V6 = ilr[Ca|Mg]
V6.lmer <- lmer(V6~I(Depth/100)+(1|Site), data=nb60)
summary(V6.lmer)
V6.ranef <- ranef(V6.lmer)
V6.se <- se.coef(V6.lmer)
coefplot(V6.ranef$Site[,1], V6.se$Site[,1], varnames=rownames(V6.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [Ca | Mg]")

# V7 = ilr[C|N]
V7.lmer <- lmer(V7~I(Depth/100)+(1|Site), data=nb60)
summary(V7.lmer)
V7.ranef <- ranef(V7.lmer)
V7.se <- se.coef(V7.lmer)
coefplot(V7.ranef$Site[,1], V7.se$Site[,1], varnames=rownames(V7.ranef$Site), xlim=c(-0.6,0.6), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [C | N]")
