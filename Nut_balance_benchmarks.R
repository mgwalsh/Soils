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

# Site-level main effects summaries ---------------------------------------
# V1 = ilr [C,N,P,K,S,Ca,Mg | Fv]
v1.lmer <- lmer(V1~log(Depth/100)+CP+WP+(1|Site), data=nb60)
summary(v1.lmer)
v1.ranef <- ranef(v1.lmer)
v1.se <- se.coef(v1.lmer)
coefplot(v1.ranef$Site[,1], v1.se$Site[,1], varnames=rownames(v1.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr [C,N,P,K,S,Ca,Mg | Fv]")

# V2 = ilr [P,K,S,Ca,Mg | C,N]
v2.lmer <- lmer(V2~log(Depth/100)+CP+WP+(1|Site), data=nb60)
summary(v2.lmer)
v2.ranef <- ranef(v2.lmer)
v2.se <- se.coef(v2.lmer)
coefplot(v2.ranef$Site[,1], v2.se$Site[,1], varnames=rownames(v2.ranef$Site), xlim=c(-10,10), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[P,K,S,Ca,Mg | C,N]")

# V3 = ilr [P,S | K,Ca,Mg]
v3.lmer <- lmer(V3~log(Depth/100)+CP+WP+(1|Site), data=nb60)
summary(v3.lmer)
v3.ranef <- ranef(v3.lmer)
v3.se <- se.coef(v3.lmer)
coefplot(v3.ranef$Site[,1], v3.se$Site[,1], varnames=rownames(v3.ranef$Site), xlim=c(-6,6), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[P,S | K,Ca,Mg]")

# V4 = ilr [K | Ca,Mg]
v4.lmer <- lmer(V4~I(Depth/100)+CP+WP+(1|Site), data=nb60)
summary(v4.lmer)
v4.ranef <- ranef(v4.lmer)
v4.se <- se.coef(v4.lmer)
coefplot(v4.ranef$Site[,1], v4.se$Site[,1], varnames=rownames(v4.ranef$Site), xlim=c(-3,3), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[K | Ca,Mg]")

# V5 = ilr [P | S]
v5.lmer <- lmer(V5~I(Depth/100)+CP+WP+(1|Site), data=nb60)
summary(v5.lmer)
v5.ranef <- ranef(v5.lmer)
v5.se <- se.coef(v5.lmer)
coefplot(v5.ranef$Site[,1], v5.se$Site[,1], varnames=rownames(v5.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[P | S]")

# V6 = ilr [Ca | Mg]
v6.lmer <- lmer(V6~I(Depth/100)+CP+WP+(1|Site), data=nb60)
summary(v6.lmer)
v6.ranef <- ranef(v6.lmer)
v6.se <- se.coef(v6.lmer)
coefplot(v6.ranef$Site[,1], v6.se$Site[,1], varnames=rownames(v6.ranef$Site), xlim=c(-2,2), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[Ca | Mg]")

# V7 = ilr [C | N]
v7.lmer <- lmer(V7~I(Depth/100)+CP+WP+(1|Site), data=nb60)
summary(v7.lmer)
v7.ranef <- ranef(v7.lmer)
v7.se <- se.coef(v7.lmer)
coefplot(v7.ranef$Site[,1], v7.se$Site[,1], varnames=rownames(v7.ranef$Site), xlim=c(-1,1), CI=2, cex.var=0.6, cex.pts=1.0, main="ilr[C | N]")

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
