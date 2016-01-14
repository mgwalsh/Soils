#' Soil nutrient mass balance benchmarks with AfSIS-1 data:
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites
#' M. Walsh, January 2016

# install.packages(c("devtools","arm","quantreg"), dependencies=T)
require(devtools)
require(arm)
require(quantreg)

# Data setup --------------------------------------------------------------
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Soils/master/Nut_balance_setup.R"
# source_url(SourceURL)

# Enrichment factor models ------------------------------------------------
# load XRF Aluminum reference data
download("https://www.dropbox.com/s/0j6hl56pwy2yiwt/Al_XRF.csv?dl=0", "Al_XRF.csv", mode="wb")
xrf <- read.table("Al_XRF.csv", header=T, sep=",")
nb60 <- merge(nb60, xrf, by="SSN")

# P | Al
nb60$PAL <- nb60$P/nb60$Al
PAL.rq <- rq(PAL~I(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(PAL.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

tau <- 0.25 ## set quantile reference level
PAL.rq <- rq(PAL~I(Depth/100)+CP+WP, tau = tau, data=nb60) 
nb60$PALp <- predict(PAL.rq, nb60)
nb60$PALd <- log(nb60$PAL/nb60$PALp)
nb60$PALq <- ifelse(nb60$PALp > nb60$PAL, 1, 0) ## predict above/below quantile value
prop.table(table(nb60$PALq))
hist(nb60$PALd)

# K | Al
nb60$KAL <- nb60$K/nb60$Al
KAL.rq <- rq(KAL~I(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(KAL.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

tau <- 0.25 ## set quantile reference level
KAL.rq <- rq(KAL~I(Depth/100)+CP+WP, tau = tau, data=nb60) 
nb60$KALp <- predict(KAL.rq, nb60)
nb60$KALd <- log(nb60$KAL/nb60$KALp)
nb60$KALq <- ifelse(nb60$KALp > nb60$KAL, 1, 0) ## predict above/below quantile value
prop.table(table(nb60$KALq))
hist(nb60$KALd)

# S | Al
nb60$SAL <- nb60$S/nb60$Al
SAL.rq <- rq(SAL~I(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(SAL.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

tau <- 0.25 ## set quantile reference level
SAL.rq <- rq(SAL~I(Depth/100)+CP+WP, tau = tau, data=nb60) 
nb60$SALp <- predict(SAL.rq, nb60)
nb60$SALd <- log(nb60$SAL/nb60$SALp)
nb60$SALq <- ifelse(nb60$SALp > nb60$SAL, 1, 0) ## predict above/below quantile value
prop.table(table(nb60$SALq))
hist(nb60$SALd)

# Site-level EDF summaries ------------------------------------------------
# PAL = P | Al
PAL.lmer <- lmer(PALd~1+(1|Site), data=nb60)
summary(PAL.lmer)
PAL.coef <- coef(PAL.lmer)
PAL.se <- se.coef(PAL.lmer)
coefplot(PAL.coef$Site[,1], PAL.se$Site[,1], varnames=rownames(PAL.coef$Site), xlim=c(-2,4), CI=2, cex.var=0.6, cex.pts=1.0, main="")

# KAL = K | Al
KAL.lmer <- lmer(KALd~1+(1|Site), data=nb60)
summary(KAL.lmer)
KAL.coef <- coef(KAL.lmer)
KAL.se <- se.coef(KAL.lmer)
coefplot(KAL.coef$Site[,1], KAL.se$Site[,1], varnames=rownames(KAL.coef$Site), xlim=c(-2,4), CI=2, cex.var=0.6, cex.pts=1.0, main="")

#  SAL = S | Al
SAL.lmer <- lmer(SALd~1+(1|Site), data=nb60)
summary(SAL.lmer)
SAL.coef <- coef(SAL.lmer)
SAL.se <- se.coef(SAL.lmer)
coefplot(SAL.coef$Site[,1], SAL.se$Site[,1], varnames=rownames(SAL.coef$Site), xlim=c(-2,3), CI=2, cex.var=0.6, cex.pts=1.0, main="")

# Quantile regression -----------------------------------------------------
# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
V1.rq <- rq(V1~log(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(V1.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

# V2 = ilr [P,K,S,Ca,Mg | C,N]
V2.rq <- rq(V2~log(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(V2.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

# V3 = ilr [P,S | K,Ca,Mg]
V3.rq <- rq(V3~log(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(V3.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

# V4 = ilr [K | Ca,Mg]
V4.rq <- rq(V4~log(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(V4.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

# V5 = ilr [P | S]
V5.rq <- rq(V5~log(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(V5.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

# V6 = ilr [Ca | Mg]
V6.rq <- rq(V6~log(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(V6.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots

# V7 = ilr [C | N]
V7.rq <- rq(V7~log(Depth/100)+CP+WP, tau = seq(0.05, 0.95, by = 0.05), data=nb60)
plot(summary(V7.rq), main = c("Intercept","Depth","Cropland","Woody cover")) ## Coefficient plots
