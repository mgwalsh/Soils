#' Exploratory soil nutrient mass balance benchmarks with AfSIS-1 data:
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh & J. Chen Nov. 2015

# install.packages(c("devtools"), dependencies=T)
require(devtools)

# Data setup --------------------------------------------------------------
SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Soils/master/Nut_balance_setup.R"
source_url(SourceURL)

# Project GeoSurvey data to Africa LAEA from LonLat
geos.laea <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.laea) <- c("x","y")
geos <- cbind(geos.laea, geos)

# Generate AfSIS 100m resolution grid cell ID's (GID)
res.pixel <- 100 ## set pixel resolution in m
xgid <- ceiling(abs(geos$x)/res.pixel)
ygid <- ceiling(abs(geos$y)/res.pixel)
gidx <- ifelse(geos$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(geos$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
geos.gid <- cbind(GID, geos)

# merge GeoSurvey observations with nutrient mass balance samples
vars <- c("GID","BP","CP","WP")
geos.gid <- geos.gid[vars]
nb60 <- merge(nb60.gid, geos.gid, by="GID")

# Site level summaries ----------------------------------------------------
# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
v1.lmer <- lmer(V1~I(Depth/100)+CP+WP+I(Depth/100)*CP+I(Depth/100)*WP+(1|Site), data=nb60)
summary(v1.lmer)
v1.ranef <- ranef(v1.lmer)
v1.se <- se.coef(v1.lmer)
coefplot(v1.ranef$Site[,1], v1.se$Site[,1], varnames=rownames(v1.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[C,N,P,K,S,Ca,Mg | Fv]")

# V2 = ilr [P,K,S,Ca,Mg | C,N]
v2.lmer <- lmer(V2~I(Depth/100)+CP+WP+I(Depth/100)*CP+I(Depth/100)*WP+(1|Site), data=nb60)
summary(v2.lmer)
v2.ranef <- ranef(v2.lmer)
v2.se <- se.coef(v2.lmer)
coefplot(v2.ranef$Site[,1], v2.se$Site[,1], varnames=rownames(v2.ranef$Site), xlim=c(-10,10), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[P,K,S,Ca,Mg | C,N]")

# V3 = ilr [P,S | K,Ca,Mg]
v3.lmer <- lmer(V3~I(Depth/100)+CP+WP+I(Depth/100)*CP+I(Depth/100)*WP+(1|Site), data=nb60)
summary(v3.lmer)
v3.ranef <- ranef(v3.lmer)
v3.se <- se.coef(v3.lmer)
coefplot(v3.ranef$Site[,1], v3.se$Site[,1], varnames=rownames(v3.ranef$Site), xlim=c(-6,6), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[P,S | K,Ca,Mg]")

# V4 = ilr [K | Ca,Mg]
v4.lmer <- lmer(V4~I(Depth/100)+CP+WP+I(Depth/100)*CP+I(Depth/100)*WP+(1|Site), data=nb60)
summary(v4.lmer)
v4.ranef <- ranef(v4.lmer)
v4.se <- se.coef(v4.lmer)
coefplot(v4.ranef$Site[,1], v4.se$Site[,1], varnames=rownames(v4.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[K | Ca,Mg]")

# V5 = ilr [P | S]
v5.lmer <- lmer(V5~I(Depth/100)+CP+WP+I(Depth/100)*CP+I(Depth/100)*WP+(1|Site), data=nb60)
summary(v5.lmer)
v5.ranef <- ranef(v5.lmer)
v5.se <- se.coef(v5.lmer)
coefplot(v5.ranef$Site[,1], v5.se$Site[,1], varnames=rownames(v5.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[P | S]")

# V6 = ilr [Ca | Mg]
v6.lmer <- lmer(V6~I(Depth/100)+CP+WP+I(Depth/100)*CP+I(Depth/100)*WP+(1|Site), data=nb60)
summary(v6.lmer)
v6.ranef <- ranef(v6.lmer)
v6.se <- se.coef(v6.lmer)
coefplot(v6.ranef$Site[,1], v6.se$Site[,1], varnames=rownames(v6.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[Ca | Mg]")

# V7 = ilr [C | N]
v7.lmer <- lmer(V7~I(Depth/100)+CP+WP+I(Depth/100)*CP+I(Depth/100)*WP+(1|Site), data=nb60)
summary(v7.lmer)
v7.ranef <- ranef(v7.lmer)
v7.se <- se.coef(v7.lmer)
coefplot(v7.ranef$Site[,1], v7.se$Site[,1], varnames=rownames(v7.ranef$Site), xlim=c(-1,1), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[C | N]")

# Topsoil / subsoil contrast ----------------------------------------------
tsc.glmer <- glmer(factor(Depth)~V1+V2+V4+V5+(1|Site), family="binomial"(link=logit), data=nb60)
summary(tsc.glmer)
tsc.ranef <- ranef(tsc.glmer)
tsc.se <- se.coef(tsc.glmer)
coefplot(tsc.ranef$Site[,1], tsc.se$Site[,1], varnames=rownames(tsc.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="SFI")
fix <- fixef(tsc.glmer) ## extract mean effects

# Proposed soil nutrient balance / fertility index (SFI), based on topsoil/subsoil contrast 
attach(nb60)
nb60$SFI <- (V1*fix[2]+V2*fix[3]+V4*fix[4]+V5*fix[5]+fix[1])*-1
detach(nb60)

# Topsoil / subsoil (SFI) contrast ecdf plot
top <- subset(nb60, Depth==10, select=c(SFI))
quantile(top$SFI) ## value above the 50% topsoil SFI quantile ~ high fertility soils
sub <- subset(nb60, Depth==35, select=c(SFI))
plot(ecdf(top$SFI), main="", xlab="SFI", ylab="Cum. proportion of observations", xlim=c(-4,4), verticals=T, lty=1, lwd=2, col="red", do.points=F)
abline(0.5,0, lty=2, col="grey")
plot(ecdf(sub$SFI), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)

# Train/Test set partition ------------------------------------------------
sites <- rownames(tsc.ranef$Site)
set.seed(5321)
train <- sample(sites, 0.8*length(sites))
nb60_cal <- nb60[ nb60$Site%in%train, ] ## calibration data
nb60_val <- nb60[!nb60$Site%in%train, ] ## validation data