# title         : BLD_AfSIS.R
# purpose       : Mapping bulk density in kg per cubic meter (20 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl, G.B.M. Heuvelink and B. Kempen
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles "afsp.geo.rda" and WorldGrids maps
# outputs       : Predictions of bulk density at 6 depths;
# remarks 1     : Accuracy of predictions is probably low;


#install.packages("GSIF", repos="http://R-Forge.R-project.org")
library(GSIF) 
library(maptools)
library(plotKML)
library(gstat)
data(R_pal)
#load(file="BLD_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## Load input data
load("../samples/afsp.geo.rda") # Points
load("afgrid20.rda") # Covariates
load("../soilmasks/admin.af.rda")
proj4string(admin.af) <- get("ref_CRS", envir = plotKML.opts)
admin.af <- spTransform(admin.af, CRS(af.csy))
admin.af <- as(admin.af, "SpatialLines")
    
## Gaussian model:
formulaString.BLD <- as.formula(paste('BLD ~ DEMSRE0a + SLPSRT0a + TDMMOD0a + TDSMOD0a + TNMMOD0a + TNSMOD0a + EVMMOD0a + EVSMOD0a + PREGSM0a +', paste("G0", c(1:7,9), "ESA0a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA0a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG0x + TWISRE0x + L3POBI0a'))
BLD.m <- fit.gstatModel(observations=afsp.geo, formulaString.BLD, covariates=afgrid20, dimensions = "3D")
summary(BLD.m@regModel)
## plot variogram:
pnts <- SpatialPointsDataFrame(BLD.m@sp, BLD.m@regModel$model)
pnts$residual <- residuals(BLD.m@regModel)
## Manually fix vgm parameters?
BLD.m@vgmModel$psill = c(.029, .021)
v = BLD.m@vgmModel
class(v) = c("variogramModel", "data.frame")
plot(variogram(residual~1, pnts[runif(length(pnts$residual))<.2,]), v)
## close to pure nugget effect!!

## Remove classes not represented otherwise the predict.glm reports an error of type "factor 'SGEUSG0x' has new levels"!
fix.c = levels(afgrid20$SGEUSG0x)[!(levels(afgrid20$SGEUSG0x) %in% levels(BLD.m@regModel$model$SGEUSG0x))]
for(j in fix.c){
  afgrid20$SGEUSG0x[afgrid20$SGEUSG0x == j] <- "pCm"
}
summary(afgrid20$SGEUSG0x)
## Prepare new locations and make predictions:
new3D <- sp3D(afgrid20)

## Predict at all 6 depths!
BLD.l <- lapply(new3D, FUN=function(x){predict(BLD.m, x, nmin=60, nmax=90, nfold=0)})
## TH: takes ca 3 minutes!

## write to a GIS format:
for(i in 1:length(BLD.l)){
  BLD.l[[i]]@predicted$BLD <- round(BLD.l[[i]]@predicted$BLD, 1)
  BLD.l[[i]]@predicted$BLD_L <- round(BLD.l[[i]]@predicted$BLD - 1.645*sqrt(BLD.l[[i]]@predicted$var1.var), 1)
  BLD.l[[i]]@predicted$BLD_U <- round(BLD.l[[i]]@predicted$BLD + 1.645*sqrt(BLD.l[[i]]@predicted$var1.var), 1) 
  writeGDAL(BLD.l[[i]]@predicted["BLD"], paste("BLD", "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-99999)
  writeGDAL(BLD.l[[i]]@predicted["BLD_L"], paste("BLD", "_sd", i, "_L.tif", sep=""), "GTiFF", mvFlag=-99999)
  writeGDAL(BLD.l[[i]]@predicted["BLD_U"], paste("BLD", "_sd", i, "_U.tif", sep=""), "GTiFF", mvFlag=-99999)  
}

## compress all produced maps:
system("7za a BLD_glmrk.tif.7z BLD_sd*.tif")

## Plot in R:
r <- c(1, 1.8)
rx <- rev(as.character(round(c(round(r[1], 1), NA, round(mean(r), 1), NA, round(r[2], 1)), 2)))
png(file = "Fig_predicted_Africa_BLD_6_depths.png", res = 150, width = 1300, height = 1300)
windows(width = 12, height = 12)
dev.off() 
par(mfrow=c(2,3), mar=c(.0,.0,2.0,.0), oma=c(0,0,0,0))
image(raster(BLD.l[[1]]@predicted["BLD"]), col=R_pal[[1]], main='0-5 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(BLD.l[[2]]@predicted["BLD"]), col=R_pal[[1]], main='5-15 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(BLD.l[[3]]@predicted["BLD"]), col=R_pal[[1]], main='15-30 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(BLD.l[[4]]@predicted["BLD"]), col=R_pal[[1]], main='30-60 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(BLD.l[[5]]@predicted["BLD"]), col=R_pal[[1]], main='60-100 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(BLD.l[[6]]@predicted["BLD"]), col=R_pal[[1]], main='100-200 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
dev.off()

save.image(file="BLD_AfSIS.RData", compress="xz")

# end of script;