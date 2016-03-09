#' Soil nutrient mass balance benchmarks with AfSIS-1 data:
#' C,N, Mehlich-3 extractable elements, XRF A-CN-K data and Laser Diffraction data from 60 sentinel sites
#' M. Walsh, January 2016

# install.packages(c("devtools","downloader","compositions","arm","quantreg","circlize","RColorBrewer"), dependencies=T)
suppressPackageStartupMessages({
require(devtools)
require(downloader)
require(compositions)
require(arm)
require(quantreg)
require(circlize)
require(RColorBrewer)
})

# Data setup --------------------------------------------------------------
# Run this first
SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Soils/master/Nut_balance_setup.R"
source_url(SourceURL)

# load Mehlich-3 Al,B,Cu,Fe,Mn & Zn data
download("https://www.dropbox.com/s/bivkvxrjno8fo67/nb60_micro.csv?dl=0", "nb60_micro.csv", mode="wb")
mic <- read.table("nb60_micro.csv", header=T, sep=",")
nb60 <- merge(nb60, mic, by="SSN")

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
vars <- c("Fv","C","N","P","K","Ca","Mg","S")
varl <- c("C","N","P","K","Ca","Mg","S","Mn","B","Cu","Zn","Fe")
mnus <- nb60[vars]
mnul <- nb60[varl]

# Calculate compositional affinities
mdats <- as.data.frame(clr(acomp(mnus)))
mcors <- (cor(mnus))^2 ## affinities = sum of compositional coefficients of determination (R^2)
mcors[lower.tri(mcors, diag=F)] <- 0
mcors <- ifelse(mcors < 0.05 | mcors == 1.0, 0, mcors)

mdatl <- as.data.frame(clr(acomp(mnul)))
mcorl <- (cor(mdatl))^2 ## affinities = sum of compositional coefficients of determination (R^2)
mcorl[lower.tri(mcorl, diag=F)] <- 0
mcorl <- ifelse(mcorl < 0.1 | mcorl == 1.0, 0, mcorl)

# Chord diagrams
# Macro nutrient affinities
set.seed(1235813)
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

# Macro+micro nutrient affinities
set.seed(123)
circos.par(gap.degree = 3)
chordDiagram(mcorl, directional = F, annotationTrack = "grid",
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
edf <- cbind(nb60[c(1,5,8,17:25,36:37,41,45:46)], mdats)
tauQL <- 0.5 ## set dependent variable quantile reference level

# clr [N] | Depth, GeoSurvey, A-CN-K, LDPSA
N.rq <- rq(N~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edf) 
summary(N.rq)

# N-quantile regression coefficients plot
est <- as.vector(coef(summary(N.rq))[2:6,1])
ses <- as.vector(coef(summary(N.rq))[2:6,2])
names <- c("log(Depth)","Cropland","CWI","PWI","CWI*PWI")
coefplot(est, ses, varnames=names, cex.var=1, cex.pts=1.5, col.pts="blue", CI=2, xlim=c(-0.6,0.6), main="")

# clr [P] | Depth, GeoSurvey, A-CN-K, LDPSA
P.rq <- rq(P~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edf) 
summary(P.rq)

# P-quantile regression coefficients plot
est <- as.vector(coef(summary(P.rq))[2:6,1])
ses <- as.vector(coef(summary(P.rq))[2:6,2])
coefplot(est, ses, cex.pts=1.5, offset=0.1, CI=2, main="", add=T)

# clr [K] | Depth, GeoSurvey, A-CN-K, LDPSA
K.rq <- rq(K~log(Depth)+CP+wV1*tV1, tau = tauQL, data=edf) 
summary(K.rq)

# K-quantile regression coefficients plot
est <- as.vector(coef(summary(K.rq))[2:6,1])
ses <- as.vector(coef(summary(K.rq))[2:6,2])
coefplot(est, ses, cex.pts=1.5, col.pts="red", offset=0.2, CI=2, main="", add=T)

# P site-level summaries
edf$PQL <- predict(P.rq, edf)
PQL.lmer <- lmer(I(P-PQL)~1+(1|Site), edf)
summary(PQL.lmer)
PQL.coef <- coef(PQL.lmer)
PQL.se <- se.coef(PQL.lmer)
coefplot(PQL.coef$Site[,1], PQL.se$Site[,1], varnames=rownames(PQL.coef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# K site-level summaries
edf$KQL <- predict(K.rq, edf)
KQL.lmer <- lmer(I(K-KQL)~1+(1|Site), edf)
summary(KQL.lmer)
KQL.coef <- coef(KQL.lmer)
KQL.se <- se.coef(KQL.lmer)
coefplot(KQL.coef$Site[,1], KQL.se$Site[,1], varnames=rownames(KQL.coef$Site), xlim=c(-1.5,1.5), CI=2, cex.var=0.6, cex.pts=0.9, main="")


# ilr [C,N,P,K,Ca,Mg,S | Fv] | Depth, GeoSurvey, A-CN-K, LDPSA
tauQL <- 0.5 ## set dependent variable quantile reference level
V1.rq <- rq(V1~I(Depth/10)+wV1*tV1, tau = tauQL, data=edf) 
summary(V1.rq)

# V1 quantile regression coefficient plot
est <- as.vector(coef(summary(V1.rq))[2:5,1])
ses <- as.vector(coef(summary(V1.rq))[2:5,2])
names <- c("Depth","CWI","PWI","CWI*PWI")
coefplot(est, ses, varnames=names, CI=2, main="")

# V1 site-level summaries
edf$V1QL <- predict(V1.rq, edf)
V1QL.lmer <- lmer(I(V1-V1QL)~1+(1|Site), edf)
summary(V1QL.lmer)
V1QL.coef <- coef(V1QL.lmer)
V1QL.se <- se.coef(V1QL.lmer)
coefplot(V1QL.coef$Site[,1], V1QL.se$Site[,1], varnames=rownames(V1QL.coef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# ilr [P,K,Ca,Mg,S | C,N] | Depth, GeoSurvey, A-CN-K, LDPSA
tauQL <- 0.5 ## set dependent variable quantile reference level
V2.rq <- rq(V2~I(Depth/100)+wV1*tV1, tau = tauQL, data=edf) 
summary(V2.rq)
edf$V2QL <- predict(V2.rq, edf)

# V2 site-level summaries
V2QL.lmer <- lmer(I(V2-V2QL)~1+(1|Site), edf)
summary(V2QL.lmer)
V2QL.coef <- coef(V2QL.lmer)
V2QL.se <- se.coef(V2QL.lmer)
coefplot(V2QL.coef$Site[,1], V2QL.se$Site[,1], varnames=rownames(V2QL.coef$Site), xlim=c(-6,6), CI=2, cex.var=0.6, cex.pts=0.9, main="")
