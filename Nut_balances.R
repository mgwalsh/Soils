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

# load Laser Diffraction data (see: https://github.com/mgwalsh/LDPSA/blob/master/XRF_data_setup.R)
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
edf <- cbind(nb60[c(1,5,8,17:25,36:37,41)], mdats)

# clr [P] | Depth, GeoSurvey, A-CN-K composition
tauQL <- 0.5 ## set dependent variable quantile reference level
P.rq <- rq(P~I(Depth/100)+BP*CP+wV1*wV2, tau = tauQL, data=edf) 
summary(P.rq)
edf$PQL <- predict(P.rq, edf)

# PQL site-level summaries
PQL.lmer <- lmer(I(P-PQL)~1+(1|Site), edf)
summary(PQL.lmer)
PQL.coef <- coef(PQL.lmer)
PQL.se <- se.coef(PQL.lmer)
coefplot(PQL.coef$Site[,1], PQL.se$Site[,1], varnames=rownames(PQL.coef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# clr [K] | Depth, GeoSurvey, A-CN-K composition
tauQL <- 0.5 ## set dependent variable quantile reference level
K.rq <- rq(K~I(Depth/100)+BP*CP+wV1*wV2, tau = tauQL, data=edf) 
summary(K.rq)
edf$KQL <- predict(K.rq, edf)

# KQL site-level summaries
KQL.lmer <- lmer(I(K-KQL)~1+(1|Site), edf)
summary(KQL.lmer)
KQL.coef <- coef(KQL.lmer)
KQL.se <- se.coef(KQL.lmer)
coefplot(KQL.coef$Site[,1], KQL.se$Site[,1], varnames=rownames(KQL.coef$Site), xlim=c(-1.5,1.5), CI=2, cex.var=0.6, cex.pts=0.9, main="")

# ilr [C,N,P,K,Ca,Mg,S | Fv] | Depth, GeoSurvey, A-CN-K composition
tauQL <- 0.5 ## set dependent variable quantile reference level
V1.rq <- rq(V1~I(Depth/100)+BP*CP+wV1*wV2, tau = tauQL, data=edf) 
summary(V1.rq)

# ilr [P,K,Ca,Mg,S | C,N] | Depth, GeoSurvey, A-CN-K composition
tauQL <- 0.5 ## set dependent variable quantile reference level
V2.rq <- rq(V2~I(Depth/100)+BP*CP+wV1*wV2, tau = tauQL, data=edf) 
summary(V2.rq)
