#' Exploratory soil nutrient mass balance benchmarks with AfSIS-1 data:
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, Nov. 2015

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
top <- subset(nb60, Depth==10, select=c(V1,V2,V3,V4,V5,V6,V7,SFI))
quantile(top$SFI) ## value above the 50% topsoil SFI quantile ~ high fertility soils
sub <- subset(nb60, Depth==35, select=c(SFI))
plot(ecdf(top$SFI), main="", xlab="SFI", ylab="Cum. proportion of observations", xlim=c(-4,4), verticals=T, lty=1, lwd=2, col="red", do.points=F)
abline(0.5,0, lty=2, col="grey")
plot(ecdf(sub$SFI), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)


