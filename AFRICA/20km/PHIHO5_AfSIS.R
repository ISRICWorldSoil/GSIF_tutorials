# title         : PHIHO5_AfSIS.R
# purpose       : Mapping soil pH for Africa (20 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl, G.B.M. Heuvelink and B. Kempen
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles "afsp.geo.rda" and WorldGrids maps
# outputs       : Predictions of soil pH at 6 depths;
# remarks 1     : ;


#install.packages("GSIF", repos="http://R-Forge.R-project.org")
library(GSIF) 
library(maptools)
library(plotKML)
library(gstat)
data(R_pal)
#load(file="PHIHO5_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## Load input data
load("../samples/afsp.geo.rda") # Points
load("afgrid20.rda") # Covariates
load("../soilmasks/admin.af.rda")
proj4string(admin.af) <- get("ref_CRS", envir = plotKML.opts)
admin.af <- spTransform(admin.af, CRS(af.csy))
admin.af <- as(admin.af, "SpatialLines")

## Normal gaussian model:
formulaString.PHI <- as.formula(paste('PHIHO5 ~ DEMSRE0a + SLPSRT0a + TDMMOD0a + TDSMOD0a + TNMMOD0a + TNSMOD0a + EVMMOD0a + EVSMOD0a + PREGSM0a +', paste("G0", c(1:7,9), "ESA0a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA0a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG0x + TWISRE0x + L3POBI0a'))
PHI.m <- fit.gstatModel(observations=afsp.geo, formulaString.PHI, covariates=afgrid20, dimensions = "3D")
summary(PHI.m@regModel)
## plot variogram:
pnts <- SpatialPointsDataFrame(PHI.m@sp, PHI.m@regModel$model)
pnts$residual <- residuals(PHI.m@regModel)
## Manually fix vgm parameters?
PHI.m@vgmModel$psill = c(0.65, 0.25)
v = PHI.m@vgmModel
class(v) = c("variogramModel", "data.frame")
## model plot:
png(file = "Fig_GLMK_model_Africa_PHIHO5.png", res = 150, width = 1300, height = 650)
windows(width = 12, height = 6)
dev.off() 
par(mfrow=c(1,2), mar=c(4.5,4.5,.8,.8))
fv = fitted.values(PHI.m@regModel)
Ov = PHI.m@regModel$model[,1]
plot(y=fv, x=Ov, pch=21, asp=1, col="red", xlab='observed', col.main = rgb(0.99,0.99,0.99), ylab='predicted', xlim=range(c(fv, Ov), na.rm=TRUE), ylim=range(c(fv, Ov), na.rm=TRUE))
pnts.s = pnts[runif(length(pnts@data[,1]))<.2, c("PHIHO5", "residual")]
vv = variogram(PHIHO5~1, pnts.s)
vvres = variogram(residual~1, pnts.s)
plot(x=vvres$dist, y=vvres$gamma, pch="+", col="red", xlab='distance', ylab='gamma', ylim = c(0, max(vv$gamma)))
vline <- variogramLine(v, maxdist=max(vvres$dist), n=length(vvres$dist))
lines(x=vline$dist, y=vline$gamma, lwd=1)
lines(x=c(0, vv$dist), y=rep(var(PHI.m@regModel$model[,1], na.rm=TRUE), length(vv$dist)+1), lty=2)
dev.off()
## close to pure nugget effect

## Remove classes not represented otherwise the predict.glm reports an error of type "factor 'SGEUSG0x' has new levels"!
fix.c = levels(afgrid20$SGEUSG0x)[!(levels(afgrid20$SGEUSG0x) %in% levels(PHI.m@regModel$model$SGEUSG0x))]
for(j in fix.c){
  afgrid20$SGEUSG0x[afgrid20$SGEUSG0x == j] <- "pCm"
}
summary(afgrid20$SGEUSG0x)
## Prepare new locations and make predictions:
new3D <- sp3D(afgrid20)

## Predict at all 6 depths!
PHI.l <- lapply(new3D, FUN=function(x){predict(PHI.m, x, nmin=60, nmax=90, nfold=0)})
## TH: takes ca 3 minutes!

## write to a GIS format:
for(i in 1:length(PHI.l)){
  PHI.l[[i]]@predicted$PHIHO5 <- round(PHI.l[[i]]@predicted$PHIHO5, 1)
  PHI.l[[i]]@predicted$PHIHO5_L <- round(PHI.l[[i]]@predicted$PHIHO5 - 1.645*sqrt(PHI.l[[i]]@predicted$var1.var), 1)
  PHI.l[[i]]@predicted$PHIHO5_U <- round(PHI.l[[i]]@predicted$PHIHO5 + 1.645*sqrt(PHI.l[[i]]@predicted$var1.var), 1)
  writeGDAL(PHI.l[[i]]@predicted["PHIHO5"], paste("PHIHO5", "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-99999)
  writeGDAL(PHI.l[[i]]@predicted["PHIHO5_L"], paste("PHIHO5", "_sd", i, "_L.tif", sep=""), "GTiFF", mvFlag=-99999)
  writeGDAL(PHI.l[[i]]@predicted["PHIHO5_U"], paste("PHIHO5", "_sd", i, "_U.tif", sep=""), "GTiFF", mvFlag=-99999)
}

## compress all produced maps:
system("7za a PHIHO5_glmrk.tif.7z PHIHO5_sd*.tif")

## Plot in R:
r <- c(4.2, 7.8)
rx <- rev(as.character(round(c(round(r[1], 1), NA, round(mean(r), 1), NA, round(r[2], 1)), 2)))
png(file = "Fig_predicted_Africa_PHIHO5_6_depths.png", res = 150, width = 1300, height = 1300)
windows(width = 12, height = 12)
dev.off() 
par(mfrow=c(2,3), mar=c(.0,.0,2.0,.0), oma=c(0,0,0,0))
image(raster(PHI.l[[1]]@predicted["PHIHO5"]), col=rev(R_pal[[1]]), main='0-5 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=R_pal[[1]][c(1,5,10,15,20)], horiz=FALSE, bty="n")
image(raster(PHI.l[[2]]@predicted["PHIHO5"]), col=rev(R_pal[[1]]), main='5-15 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=R_pal[[1]][c(1,5,10,15,20)], horiz=FALSE, bty="n")
image(raster(PHI.l[[3]]@predicted["PHIHO5"]), col=rev(R_pal[[1]]), main='15-30 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=R_pal[[1]][c(1,5,10,15,20)], horiz=FALSE, bty="n")
image(raster(PHI.l[[4]]@predicted["PHIHO5"]), col=rev(R_pal[[1]]), main='30-60 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=R_pal[[1]][c(1,5,10,15,20)], horiz=FALSE, bty="n")
image(raster(PHI.l[[5]]@predicted["PHIHO5"]), col=rev(R_pal[[1]]), main='60-100 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=R_pal[[1]][c(1,5,10,15,20)], horiz=FALSE, bty="n")
image(raster(PHI.l[[6]]@predicted["PHIHO5"]), col=rev(R_pal[[1]]), main='100-200 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=R_pal[[1]][c(1,5,10,15,20)], horiz=FALSE, bty="n")
dev.off()

## Export to Google Earth:
plotKML(PHI.l[[1]]@observed["PHIHO5"], file="PHI_observed_sd1.kml", points_names=paste(signif(PHI.l[[1]]@observed$PHIHO5, 3)), altitude=PHI.l[[1]]@observed@coords[,3]*1000+350, kmz=TRUE, colour_scale=rev(R_pal[[1]]), z.lim=r)
pl = grid2poly(PHI.l[[1]]@predicted["PHIHO5"])  
## Operation not recommended on a standard PC!!
kml(pl, file="PHIHO5_sd1.kml", colour=PHIHO5, colour_scale=rev(R_pal[[1]]), z.lim=r, kmz=TRUE)
unlink("PHIHO5_sd1.kml")
unlink("PHI_observed_sd1.kml")

save.image(file="PHIHO5_AfSIS.RData")

# end of script;