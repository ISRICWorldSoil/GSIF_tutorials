# title         : TEX_AfSIS.R
# purpose       : Mapping soil textures for Africa (20 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl, G.B.M. Heuvelink and B. Kempen
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles "afsp.geo.rda" and WorldGrids maps
# outputs       : Predictions of silt sand and clay at 6 depths;
# remarks 1     : ;


#install.packages("GSIF", repos="http://R-Forge.R-project.org")
library(GSIF)
library(maptools)
library(plotKML)
library(gstat)
data(R_pal)
#load(file="TEX_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## Load input data
load("../samples/afsp.geo.rda") # Points
load("afgrid20.rda") # Covariates
load("../soilmasks/admin.af.rda")
proj4string(admin.af) <- get("ref_CRS", envir = plotKML.opts)
admin.af <- spTransform(admin.af, CRS(af.csy))
admin.af <- as(admin.af, "SpatialLines")

## transformation function:
logits = function(x){log((x/100)/(1-x/100))}
## fix values outside the psychical range:
sel = afsp.geo@data$methodid=="SLTPPT"|afsp.geo@data$methodid=="SNDPPT"|afsp.geo@data$methodid=="CLYPPT"
afsp.geo@data[sel,"observedValue"] <- ifelse(afsp.geo@data[sel,"observedValue"]==0, 1, ifelse(afsp.geo@data[sel,"observedValue"]==100, 99, afsp.geo@data[sel,"observedValue"]))

## define regression models:
formulaString.SLT <- as.formula(paste('logits(SLTPPT) ~ DEMSRE0a + SLPSRT0a + TDMMOD0a + TDSMOD0a + TNMMOD0a + TNSMOD0a + EVMMOD0a + EVSMOD0a + PREGSM0a +', paste("G0", c(1:7,9), "ESA0a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA0a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG0x + TWISRE0x + L3POBI0a'))
formulaString.SND <- as.formula(paste('logits(SNDPPT) ~ DEMSRE0a + SLPSRT0a + TDMMOD0a + TDSMOD0a + TNMMOD0a + TNSMOD0a + EVMMOD0a + EVSMOD0a + PREGSM0a +', paste("G0", c(1:7,9), "ESA0a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA0a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG0x + TWISRE0x + L3POBI0a'))
formulaString.CLY <- as.formula(paste('logits(CLYPPT) ~ DEMSRE0a + SLPSRT0a + TDMMOD0a + TDSMOD0a + TNMMOD0a + TNSMOD0a + EVMMOD0a + EVSMOD0a + PREGSM0a +', paste("G0", c(1:7,9), "ESA0a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA0a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG0x + TWISRE0x + L3POBI0a'))

## fit the model (takes ca 3 x 2 mins):
SLT.m <- fit.gstatModel(observations=afsp.geo, formulaString.SLT, covariates=afgrid20, dimensions="3D")
SND.m <- fit.gstatModel(observations=afsp.geo, formulaString.SND, covariates=afgrid20, dimensions="3D")
CLY.m <- fit.gstatModel(observations=afsp.geo, formulaString.CLY, covariates=afgrid20, dimensions="3D")

## plot variogram:
pnts1 <- SpatialPointsDataFrame(SLT.m@sp, SLT.m@regModel$model)
pnts1$residual <- residuals(SLT.m@regModel)
v1 = SLT.m@vgmModel
class(v1) = c("variogramModel", "data.frame")
plot(variogram(residual~1, pnts1[runif(length(pnts1$residual))<.2,]), v1)

pnts2 <- SpatialPointsDataFrame(SND.m@sp, SND.m@regModel$model)
pnts2$residual <- residuals(SND.m@regModel)
## Manually fix vgm parameters?
SND.m@vgmModel$psill = c(1.2, 0.4)
v2 = SND.m@vgmModel
class(v2) = c("variogramModel", "data.frame")
plot(variogram(residual~1, pnts2[runif(length(pnts2$residual))<.2,]), v2)

pnts3 <- SpatialPointsDataFrame(CLY.m@sp, CLY.m@regModel$model)
pnts3$residual <- residuals(CLY.m@regModel)
## Manually fix vgm parameters?
CLY.m@vgmModel$psill = c(1.05, 0.25)
v3 = CLY.m@vgmModel
class(v3) = c("variogramModel", "data.frame")
plot(variogram(residual~1, pnts3[runif(length(pnts3$residual))<.2,]), v3)

## Remove classes not represented otherwise the predict.glm reports an error of type "factor 'SGEUSG0x' has new levels"!
fix.c = levels(afgrid20$SGEUSG0x)[!(levels(afgrid20$SGEUSG0x) %in% levels(SLT.m@regModel$model$SGEUSG0x))]
for(j in fix.c){
  afgrid20$SGEUSG0x[afgrid20$SGEUSG0x == j] <- "pCm"
}
summary(afgrid20$SGEUSG0x)
## Prepare new locations and make predictions:
new3D <- sp3D(afgrid20)

## Predict at all 6 depths!
SLT.l <- lapply(new3D, FUN=function(x){predict(SLT.m, x, nmin=60, nmax=90, nfold=0)})
SND.l <- lapply(new3D, FUN=function(x){predict(SND.m, x, nmin=60, nmax=90, nfold=0)})
CLY.l <- lapply(new3D, FUN=function(x){predict(CLY.m, x, nmin=60, nmax=90, nfold=0)})
## TH: takes ca 3 x 3 minutes!

## back-transformation function:
invlogit = function(x){exp(x)/(1+exp(x))*100}
## for the back-transformation for the mean value see Diggle and Ribeiro, 2007, p. 148:
invlogit.m = function(x, v){((1+exp(-x))^(-1)-.5*v*exp(-x)*(1-exp(-x))*(1+exp(-x))^(-3) )*100}

for(i in 1:length(SLT.l)){
   SLT.l[[i]]@predicted$SLTPPT.t <- round(invlogit.m(SLT.l[[i]]@predicted$SLTPPT, SLT.l[[i]]@predicted$var1.var))
   SLT.l[[i]]@predicted$SLTPPT_L <- round(invlogit(SLT.l[[i]]@predicted$SLTPPT - 1.645*sqrt(SLT.l[[i]]@predicted$var1.var)))
   SLT.l[[i]]@predicted$SLTPPT_U <- round(invlogit(SLT.l[[i]]@predicted$SLTPPT + 1.645*sqrt(SLT.l[[i]]@predicted$var1.var)))
   SND.l[[i]]@predicted$SNDPPT.t <- round(invlogit.m(SND.l[[i]]@predicted$SNDPPT, SND.l[[i]]@predicted$var1.var))
   SND.l[[i]]@predicted$SNDPPT_L <- round(invlogit(SND.l[[i]]@predicted$SNDPPT - 1.645*sqrt(SND.l[[i]]@predicted$var1.var)))
   SND.l[[i]]@predicted$SNDPPT_U <- round(invlogit(SND.l[[i]]@predicted$SNDPPT + 1.645*sqrt(SND.l[[i]]@predicted$var1.var)))
   CLY.l[[i]]@predicted$CLYPPT.t <- round(invlogit.m(CLY.l[[i]]@predicted$CLYPPT, CLY.l[[i]]@predicted$var1.var))
   CLY.l[[i]]@predicted$CLYPPT_L <- round(invlogit(CLY.l[[i]]@predicted$CLYPPT - 1.645*sqrt(CLY.l[[i]]@predicted$var1.var)))
   CLY.l[[i]]@predicted$CLYPPT_U <- round(invlogit(CLY.l[[i]]@predicted$CLYPPT + 1.645*sqrt(CLY.l[[i]]@predicted$var1.var)))
   ## standardize values so they sum up to 100%:
   sums <- rowSums(cbind(v1=SLT.l[[i]]@predicted$SLTPPT.t, v2=SND.l[[i]]@predicted$SNDPPT.t, v3=CLY.l[[i]]@predicted$CLYPPT.t))
   SLT.l[[i]]@predicted$SLTPPT.t <- round(SLT.l[[i]]@predicted$SLTPPT.t / sums * 100, 0)
   SND.l[[i]]@predicted$SNDPPT.t <- round(SND.l[[i]]@predicted$SNDPPT.t / sums * 100, 0)
   CLY.l[[i]]@predicted$CLYPPT.t <- round(CLY.l[[i]]@predicted$CLYPPT.t / sums * 100, 0)
}

## write to a GIS format:
for(i in 1:length(SLT.l)){
  writeGDAL(SLT.l[[i]]@predicted["SLTPPT.t"], paste("SLTPPT", "_sd", i, "_M.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(SLT.l[[i]]@predicted["SLTPPT_U"], paste("SLTPPT", "_sd", i, "_U.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(SLT.l[[i]]@predicted["SLTPPT_L"], paste("SLTPPT", "_sd", i, "_L.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(SND.l[[i]]@predicted["SNDPPT.t"], paste("SNDPPT", "_sd", i, "_M.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(SND.l[[i]]@predicted["SNDPPT_U"], paste("SNDPPT", "_sd", i, "_U.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(SND.l[[i]]@predicted["SNDPPT_L"], paste("SNDPPT", "_sd", i, "_L.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(CLY.l[[i]]@predicted["CLYPPT.t"], paste("CLYPPT", "_sd", i, "_M.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(CLY.l[[i]]@predicted["CLYPPT_U"], paste("CLYPPT", "_sd", i, "_U.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(CLY.l[[i]]@predicted["CLYPPT_L"], paste("CLYPPT", "_sd", i, "_L.tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
}

## compress all produced maps:
system("7za a SLTPPT_glmrk.tif.7z SLTPPT_sd*.tif")
system("7za a SNDPPT_glmrk.tif.7z SNDPPT_sd*.tif")
system("7za a CLYPPT_glmrk.tif.7z CLYPPT_sd*.tif")

## Plot in R:
r <- c(0, 80)
rx <- rev(as.character(round(c(round(r[1], 1), NA, round(mean(r), 1), NA, round(r[2], 1)), 2)))
png(file = "Fig_predicted_Africa_CLYPPT_6_depths.png", res = 150, width = 1300, height = 1300)
windows(width = 12, height = 12)
dev.off() 
par(mfrow=c(2,3), mar=c(.0,.0,2.0,.0), oma=c(0,0,0,0))
image(raster(CLY.l[[1]]@predicted["CLYPPT.t"]), col=R_pal[[8]], main='0-5 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CLY.l[[2]]@predicted["CLYPPT.t"]), col=R_pal[[8]], main='5-15 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CLY.l[[3]]@predicted["CLYPPT.t"]), col=R_pal[[8]], main='15-30 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CLY.l[[4]]@predicted["CLYPPT.t"]), col=R_pal[[8]], main='30-60 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CLY.l[[5]]@predicted["CLYPPT.t"]), col=R_pal[[8]], main='60-100 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CLY.l[[6]]@predicted["CLYPPT.t"]), col=R_pal[[8]], main='100-200 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
dev.off()

r <- c(15, 100)
rx <- rev(as.character(round(c(round(r[1], 1), NA, round(mean(r), 1), NA, round(r[2], 1)), 2)))
png(file = "Fig_predicted_Africa_SNDPPT_6_depths.png", res = 150, width = 1300, height = 1300)
windows(width = 12, height = 12)
dev.off() 
par(mfrow=c(2,3), mar=c(.0,.0,2.0,.0), oma=c(0,0,0,0))
image(raster(SND.l[[1]]@predicted["SNDPPT.t"]), col=R_pal[[8]], main='0-5 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(SND.l[[2]]@predicted["SNDPPT.t"]), col=R_pal[[8]], main='5-15 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(SND.l[[3]]@predicted["SNDPPT.t"]), col=R_pal[[8]], main='15-30 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(SND.l[[4]]@predicted["SNDPPT.t"]), col=R_pal[[8]], main='30-60 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(SND.l[[5]]@predicted["SNDPPT.t"]), col=R_pal[[8]], main='60-100 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(SND.l[[6]]@predicted["SNDPPT.t"]), col=R_pal[[8]], main='100-200 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
dev.off()

save.image(file="TEX_AfSIS.RData", compress="xz")

# end of script;