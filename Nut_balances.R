#' Extractable soil nutrient mass balance benchmarks with AfSIS-1 data:
#' C,N, Mehlich-3 extractable elements, XRF A-CN-K data and Laser Diffraction data from 60 sentinel sites
#' M. Walsh, January 2016

# install.packages(c("devtools","downloader","arm","quantreg","circlize","RColorBrewer"), dependencies=T)
suppressPackageStartupMessages({
require(devtools)
require(downloader)
require(arm)
require(quantreg)
require(circlize)
require(RColorBrewer)
})

# Data setup --------------------------------------------------------------
# Run this first
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Soils/master/Nut_balance_setup.R"
# source_url(SourceURL)

# load XRF A-CN-K reference data (see: https://github.com/mgwalsh/Soils/blob/master/XRF_data_setup.R)
download("https://www.dropbox.com/s/i73ce1l3k8lsufl/wi60.csv?dl=0", "wi60.csv", mode="wb")
wi60 <- read.table("wi60.csv", header=T, sep=",")
nb60 <- merge(nb60, wi60[c(1,6:14)], by="SSN")

# load Laser Diffraction data (see: https://github.com/mgwalsh/LDPSA/blob/master/LDPSA_starter.R)
download("https://www.dropbox.com/s/kwwyk9wgl0i7yog/LDPSA_comp.csv.zip?dl=0", "LDPSA_comp.csv.zip", mode="wb")
unzip("LDPSA_comp.csv.zip", overwrite=T)
ldsp <- read.table("LDPSA_comp.csv", header=T, sep=",")
ldsp <- subset(ldsp, TRT=="c4", select=c(SSN, Sand, Silt, Clay, V1, V2))
colnames(ldsp) <- c("SSN","Sand","Silt","Clay","tV1","tV2")
nb60 <- merge(nb60, ldsp, by="SSN")

# Plot chord diagram of affinities ------------------------------------------
# Extract variables
vars <- c("Fc","B","C","N","Mg","P","S","K","Ca","Mn","Fe","Cu","Zn")
mnus <- nb60[vars]

# Calculate compositional affinities
mcors <- (cor(mnus))^2 ## affinities = sum of compositional coefficients of determination
mcors[lower.tri(mcors, diag=F)] <- 0
mcors <- ifelse(mcors < 0.05 | mcors == 1.0, 0, mcors)

# Chord diagram
set.seed(11)
circos.par(gap.degree = 3)
chordDiagram(mcors, directional = F, annotationTrack = "grid",
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

# Quantile enrichment/ depletion factor (EDF) models ----------------------
# Data setup
vars <- c("SSN","Site","Depth","CP","wV1","tV1","N","P","K","V1","V2","V3") 
edfs <- nb60[vars]
tauQL <- 0.5 ## set dependent variable quantile reference level

# clr [N] | Depth, GeoSurvey, A-CN-K, LDPSA
N.rq <- rq(N~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edfs) 
summary(N.rq, se="boot", bsmethod="xy")

# N-quantile regression coefficients plot
est <- as.vector(coef(summary(N.rq, se="boot", bsmethod="xy"))[2:6,1])
ses <- as.vector(coef(summary(N.rq, se="boot", bsmethod="xy"))[2:6,2])
names <- c("log(Depth)","Cropland","CWI","PWI","CWI*PWI")
coefplot(est, ses, varnames=names, cex.var=1, cex.pts=1.5, col.pts="blue", CI=2, xlim=c(-0.6,0.6), main="", mar=c(0.2,4,3,1))

# clr [P] | Depth, GeoSurvey, A-CN-K, LDPSA
P.rq <- rq(P~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edfs) 
summary(P.rq, se="boot", bsmethod="xy")

# P-quantile regression coefficients plot
est <- as.vector(coef(summary(P.rq, se="boot", bsmethod="xy"))[2:6,1])
ses <- as.vector(coef(summary(P.rq, se="boot", bsmethod="xy"))[2:6,2])
coefplot(est, ses, cex.pts=1.5, offset=0.1, CI=2, main="", add=T)

# clr [K] | Depth, GeoSurvey, A-CN-K, LDPSA
K.rq <- rq(K~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edfs) 
summary(K.rq, se="boot", bsmethod="xy")

# K-quantile regression coefficients plot
est <- as.vector(coef(summary(K.rq, se="boot", bsmethod="xy"))[2:6,1])
ses <- as.vector(coef(summary(K.rq, se="boot", bsmethod="xy"))[2:6,2])
coefplot(est, ses, cex.pts=1.5, col.pts="red", offset=0.2, CI=2, main="", add=T)

# V1 = ilr [C,N,P,K,Ca,Mg,S | Fv] | Depth, GeoSurvey, A-CN-K, LDPSA
V1.rq <- rq(V1~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edfs) 
summary(V1.rq, se="boot", bsmethod="xy")

# V1-quantile regression coefficient plot
est <- as.vector(coef(summary(V1.rq, se="boot", bsmethod="xy"))[2:6,1])
ses <- as.vector(coef(summary(V1.rq, se="boot", bsmethod="xy"))[2:6,2])
names <- c("log(Depth)","Cropland","CWI","PWI","CWI*PWI")
coefplot(est, ses, varnames=names, cex.var=1, cex.pts=1.5, col.pts="blue", CI=2, xlim=c(-1.5,1.5), main="", mar=c(0.2,4,3,1))

# V2 = ilr [P,K,Ca,Mg,S | C,N] | Depth, GeoSurvey, A-CN-K, LDPSA
V2.rq <- rq(V2~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edfs) 
summary(V2.rq, se="boot", bsmethod="xy")

# V2-quantile regression coefficients plot
est <- as.vector(coef(summary(V2.rq, se="boot", bsmethod="xy"))[2:6,1])
ses <- as.vector(coef(summary(V2.rq, se="boot", bsmethod="xy"))[2:6,2])
coefplot(est, ses, cex.pts=1.5, offset=0.1, CI=2, main="", add=T)

# V3 = ilr [P,S | K,Ca,Mg] | Depth, GeoSurvey, A-CN-K, LDPSA
V3.rq <- rq(V3~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edfs) 
summary(V3.rq, se="boot", bsmethod="xy")

# V3-quantile regression coefficients plot
est <- as.vector(coef(summary(V3.rq))[2:6,1])
ses <- as.vector(coef(summary(V3.rq))[2:6,2])
coefplot(est, ses, cex.pts=1.5, col.pts="red", offset=0.2, CI=2, main="", add=T)

# Site-level enrichment-depletion factors (EDFs) --------------------------
# N EDF REML estimates
edfs$NEDF <- predict(N.rq, edfs)
N.lmer <- lmer(I(N-NEDF)~1+(1|Site), edfs)
summary(N.lmer)
ran <- ranef(N.lmer)
ses <- se.coef(N.lmer)
nam <- rownames(ran$Site) 
coefplot(ran$Site[31:60,1], ses$Site[31:60,1], varnames=nam[31:60], cex.var=0.8, col.pts="blue", xlim=c(-1.5,1.5), CI=2, main="")

# P EDF REML estimates
edfs$PEDF <- predict(P.rq, edfs)
P.lmer <- lmer(I(P-PEDF)~1+(1|Site), edfs)
summary(P.lmer)
ran <- ranef(P.lmer)
ses <- se.coef(P.lmer)
nam <- rownames(ran$Site) 
coefplot(ran$Site[31:60,1], ses$Site[31:60,1], offset=0.2, CI=2, main="", add=T)

# K EDF REML estimates
edfs$KEDF <- predict(K.rq, edfs)
K.lmer <- lmer(I(K-KEDF)~1+(1|Site), edfs)
summary(K.lmer)
ran <- ranef(K.lmer)
ses <- se.coef(K.lmer)
nam <- rownames(ran$Site) 
coefplot(ran$Site[31:60,1], ses$Site[31:60,1], col.pts="red", offset=0.4, CI=2, main="", add=T)
