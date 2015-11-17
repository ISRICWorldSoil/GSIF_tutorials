# title         : ORCDRC_AfSIS.R
# purpose       : Mapping soil organic carbon for Africa (20 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl, G.B.M. Heuvelink and B. Kempen
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles "afsp.rda" and WorldGrids maps
# outputs       : Predictions of pH, clay content and organic carbon at 6 depths;
# remarks 1     : This script takes ca 15 mins;


#install.packages("GSIF", repos="http://R-Forge.R-project.org")
library(GSIF)
library(maptools)
library(plotKML)
data(SAGA_pal)
library(gstat)
library(aqp)
#load(file="ORCDRC_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## Load input data
## Points:
load("../samples/afsp.geo.rda")
## Covariates:
load("afgrid20.rda")
names(afgrid20)
## country borders:
load("../soilmasks/admin.af.rda")
proj4string(admin.af) <- get("ref_CRS", envir = plotKML.opts)
admin.af <- spTransform(admin.af, CRS(af.csy))
admin.af <- as(admin.af, "SpatialLines")

## Log-model:
formulaString.ORC <- as.formula(paste('log1p(ORCDRC) ~ DEMSRE0a + SLPSRT0a + TDMMOD0a + TDSMOD0a + TNMMOD0a + TNSMOD0a + EVMMOD0a + EVSMOD0a + PREGSM0a +', paste("G0", c(1:7,9), "ESA0a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA0a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG0x + TWISRE0x + L3POBI0a'))
## fit model:
ORC.m <- fit.gstatModel(observations=afsp.geo, formulaString.ORC, covariates=afgrid20, dimensions = "3D")
summary(ORC.m@regModel)

## variogram:
pnts <- SpatialPointsDataFrame(ORC.m@sp, ORC.m@regModel$model)
pnts$residual <- residuals(ORC.m@regModel)
## Manually fix vgm parameters?
ORC.m@vgmModel$psill = c(0.32, 0.1)
v = ORC.m@vgmModel
class(v) = c("variogramModel", "data.frame")
## model plot:
png(file = "Fig_GLMK_model_Africa_ORCDRC.png", res = 150, width = 1300, height = 650)
windows(width = 12, height = 6)
dev.off() 
par(mfrow=c(1,2), mar=c(4.5,4.5,.8,.8))
fv = fitted.values(ORC.m@regModel)
Ov = ORC.m@regModel$model[,1]
plot(y=fv, x=Ov, pch=21, asp=1, col="red", xlab='observed', col.main = rgb(0.99,0.99,0.99), ylab='predicted', xlim=range(c(fv, Ov), na.rm=TRUE), ylim=range(c(fv, Ov), na.rm=TRUE))
pnts.s = pnts[runif(length(pnts@data[,1]))<.2, c("log1p(ORCDRC)", "residual")]
names(pnts.s)[1] = "ORCDRC"
vv = variogram(ORCDRC~1, pnts.s)
vvres = variogram(residual~1, pnts.s)
plot(x=vvres$dist, y=vvres$gamma, pch="+", col="red", xlab='distance', ylab='gamma', ylim = c(0, max(vv$gamma)))
vline <- variogramLine(v, maxdist=max(vvres$dist), n=length(vvres$dist))
lines(x=vline$dist, y=vline$gamma, lwd=1)
lines(x=c(0, vv$dist), y=rep(var(ORC.m@regModel$model[,1], na.rm=TRUE), length(vv$dist)+1), lty=2)
dev.off()

## Remove classes not represented otherwise the predict.glm reports an error of type "factor 'SGEUSG0x' has new levels"!
fix.c = levels(afgrid20$SGEUSG0x)[!(levels(afgrid20$SGEUSG0x) %in% levels(ORC.m@regModel$model$SGEUSG0x))]
for(j in fix.c){
  afgrid20$SGEUSG0x[afgrid20$SGEUSG0x == j] <- "pCm"
}
## Prepare new locations and make predictions:
new3D <- sp3D(afgrid20)
## Predict at all 6 depths!
ORC.l <- lapply(new3D, FUN=function(x){predict(ORC.m, x, nmin=60, nmax=90, nfold=0)})
## TH: takes ca 5 minutes!

## back-transform and write to a GIS format (mean value and 90% probability range):
for(i in 1:length(ORC.l)){
  ORC.l[[i]]@predicted$ORCDRC.t <- expm1(ORC.l[[i]]@predicted$var1.pred + ORC.l[[i]]@predicted$var1.var/2)
  ORC.l[[i]]@predicted$ORCDRC.t_L <- expm1(ORC.l[[i]]@predicted$var1.pred - 1.645*sqrt(ORC.l[[i]]@predicted$var1.var))  
  ORC.l[[i]]@predicted$ORCDRC.t_U <- expm1(ORC.l[[i]]@predicted$var1.pred + 1.645*sqrt(ORC.l[[i]]@predicted$var1.var))
  writeGDAL(ORC.l[[i]]@predicted["ORCDRC.t"], paste("ORCDRC", "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
  writeGDAL(ORC.l[[i]]@predicted["ORCDRC.t_L"], paste("ORCDRC", "_sd", i, "_L.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
  writeGDAL(ORC.l[[i]]@predicted["ORCDRC.t_U"], paste("ORCDRC", "_sd", i, "_U.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
  #writeGDAL(ORC.l[[i]]@predicted["var1.var"], paste("ORCDRC", "_sd", i, "_se.tif", sep=""), "GTiFF", mvFlag=-99999)
}

## compress all produced maps:
system("7za a ORCDRC_glmrk.tif.7z ORCDRC_sd*.tif")

## Plot in R:
data(SAGA_pal)
r <- c(0, 60)
rx <- rev(as.character(round(c(round(r[1], 0), NA, round(expm1(mean(log1p(r))), 0), NA, round(r[2], 0)), 2)))
png(file = "Fig_predicted_Africa_ORCDRC_6_depths.png", res = 150, width = 1300, height = 1300)
windows(width = 12, height = 6)
dev.off() 
par(mfrow=c(2,3), mar=c(.0,.0,2.0,.0), oma=c(0,0,0,0))
image(log1p(raster(ORC.l[[1]]@predicted["ORCDRC.t"])), col=SAGA_pal[[1]], main='0-5 cm', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(log1p(raster(ORC.l[[2]]@predicted["ORCDRC.t"])), col=SAGA_pal[[1]], main='5-15 cm', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(log1p(raster(ORC.l[[3]]@predicted["ORCDRC.t"])), col=SAGA_pal[[1]], main='15-30 cm', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(log1p(raster(ORC.l[[4]]@predicted["ORCDRC.t"])), col=SAGA_pal[[1]], main='30-60 cm', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(log1p(raster(ORC.l[[5]]@predicted["ORCDRC.t"])), col=SAGA_pal[[1]], main='60-100 cm', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(log1p(raster(ORC.l[[6]]@predicted["ORCDRC.t"])), col=SAGA_pal[[1]], main='100-200 cm', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
dev.off()

## plot the model errors:
ORC.l[[1]]@predicted$glm.var <- stats::predict.glm(ORC.m@regModel, newdata=new3D[[1]], type="response", se.fit = TRUE, na.action = na.pass)$se.fit
rv = c(0.011, 0.1)
rl = rev(as.character(c(signif(rv[1], 3), NA, signif(mean(rv), 3), NA, signif(rv[2], 3))))
png(file = "Fig_mapping_error_Africa_ORCDRC_1sd.png", res = 150, width = 1200, height = 600)
windows(width = 12, height = 6)
dev.off() 
par(mfrow=c(1,2), mar=c(.0,.0,2.0,.0), oma=c(0,0,0,0))
image(log1p(raster(ORC.l[[1]]@predicted["ORCDRC.t"])), col=SAGA_pal[[1]], cex.main=.8, main='predictions (0-5 cm)', axes = FALSE, xlab="", ylab="", zlim=log1p(r), asp=1)
lines(admin.af, col="black")
points(pnts, pch="+", cex=.3)
legend("bottomleft", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n", cex=.8)
image(raster(ORC.l[[1]]@predicted["glm.var"]), col=SAGA_pal[[10]], cex.main=.8, main='prediction variance (0-5 cm)', axes = FALSE, xlab="", ylab="", zlim=rv, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rl, fill=rev(SAGA_pal[[10]][c(1,5,10,15,20)]), horiz=FALSE, bty="n", cex=.8)
dev.off()

## Export to Google Earth:
ORC.l[[1]]@observed$ORCDRC.t <- expm1(ORC.l[[1]]@observed$ORCDRC)
plotKML(ORC.l[[1]]@observed["ORCDRC"], folder.name="ORCDRC", file="ORC_observed_sd1.kml", points_names=paste(signif(ORC.l[[1]]@observed$ORCDRC.t, 3)), altitude=ORC.l[[1]]@observed@coords[,3]*1000+350, kmz=TRUE, colour_scale=SAGA_pal[[1]], z.lim=log1p(c(0,60)))
pl = grid2poly(ORC.l[[1]]@predicted["ORCDRC.t"])  
## Operation not recommended on a standard PC!!
kml(pl, folder.name="predicted", file="ORCDRC_sd1.kml", colour=log1p(ORCDRC.t), colour_scale = SAGA_pal[[1]], z.lim=log1p(c(0,60)), kmz=TRUE)
unlink("ORCDRC_sd1.kml")
unlink("ORC_observed_sd1.kml")

save.image(file="ORCDRC_AfSIS.RData", compress="xz")

# end of script; 