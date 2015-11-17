# title         : CRFVOL_AfSIS.R
# purpose       : Mapping volume percentage of coarse fragments (> 2 mm) (20 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl, G.B.M. Heuvelink and B. Kempen
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles "afsp.geo.rda" and WorldGrids maps
# outputs       : Predictions of volume percentage of coarse fragments at 6 depths;
# remarks 1     : The maps produced here are of critically low quality;


#install.packages("GSIF", repos="http://R-Forge.R-project.org")
library(GSIF)
library(maptools)
library(plotKML)
data(R_pal)
library(gstat)
library(aqp)
#load(file="CRFVOL_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"

## Load input data
load("../samples/afsp.geo.rda") # Points
load("afgrid20.rda") # Covariates
load("../soilmasks/admin.af.rda")
proj4string(admin.af) <- get("ref_CRS", envir = plotKML.opts)
admin.af <- spTransform(admin.af, CRS(af.csy))
admin.af <- as(admin.af, "SpatialLines")

## TH: this is highly skewed variable hence it probably requires a zero-inflated model [http://www.ats.ucla.edu/stat/r/dae/zipoisson.htm]!  
## Poisson regression:
formulaString.CRF <- as.formula(paste('CRFVOL ~ DEMSRE0a + SLPSRT0a + TDMMOD0a + TDSMOD0a + TNMMOD0a + TNSMOD0a + EVMMOD0a + EVSMOD0a + PREGSM0a +', paste("G0", c(1:7,9), "ESA0a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA0a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG0x + TWISRE0x + L3POBI0a'))
CRF.m <- fit.gstatModel(observations=afsp.geo, formulaString.CRF, covariates=afgrid20, family=poisson, dimensions = "3D")
## ... model fitting is difficult
summary(CRF.m@regModel)
## plot variogram:
pnts <- SpatialPointsDataFrame(CRF.m@sp, CRF.m@regModel$model)
pnts$residual <- residuals(CRF.m@regModel)
## Manually fix vgm parameters?
CRF.m@vgmModel$psill = c(16, 4)
v = CRF.m@vgmModel
class(v) = c("variogramModel", "data.frame")
plot(variogram(residual~1, pnts[runif(length(pnts$residual))<.2,]), v)
## close to pure nugget effect!!

## Remove classes not represented otherwise the predict.glm reports an error of type "factor 'SGEUSG0x' has new levels"!
fix.c = levels(afgrid20$SGEUSG0x)[!(levels(afgrid20$SGEUSG0x) %in% levels(CRF.m@regModel$model$SGEUSG0x))]
for(j in fix.c){
  afgrid20$SGEUSG0x[afgrid20$SGEUSG0x == j] <- "pCm"
}
## Prepare new locations and make predictions:
new3D <- sp3D(afgrid20)

## Predict at all 6 depths!
CRF.l <- lapply(new3D, FUN=function(x){predict(CRF.m, x, nmin=60, nmax=90, nfold=0)})
## TH: takes ca 3 minutes!

## write to a GIS format:
for(i in 1:length(CRF.l)){
  ## TH: how to derive upper and lower limits for a Poisson model?
  writeGDAL(CRF.l[[i]]@predicted["CRFVOL"], paste("CRFVOL", "_sd", i, ".tif", sep=""), "GTiFF", type = "Byte", mvFlag=255)
  writeGDAL(CRF.l[[i]]@predicted["var1.var"], paste("CRFVOL", "_sd", i, "_var.tif", sep=""), "GTiFF", mvFlag=-99999)
}

## compress all produced maps:
system("7za a CRFVOL_glmrk.tif.7z CRFVOL_sd*.tif")

## Plot in R:
r <- c(0, 25)
rx <- rev(as.character(round(c(round(r[1], 1), NA, round(mean(r), 1), NA, round(r[2], 1)), 2)))
png(file = "Fig_predicted_Africa_CRFVOL_6_depths.png", res = 150, width = 1300, height = 1300)
windows(width = 12, height = 12)
dev.off() 
par(mfrow=c(2,3), mar=c(.0,.0,2.0,.0), oma=c(0,0,0,0))
image(raster(CRF.l[[1]]@predicted["CRFVOL"]), col=R_pal[[8]], main='0-5 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CRF.l[[2]]@predicted["CRFVOL"]), col=R_pal[[8]], main='5-15 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CRF.l[[3]]@predicted["CRFVOL"]), col=R_pal[[8]], main='15-30 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CRF.l[[4]]@predicted["CRFVOL"]), col=R_pal[[8]], main='30-60 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CRF.l[[5]]@predicted["CRFVOL"]), col=R_pal[[8]], main='60-100 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
image(raster(CRF.l[[6]]@predicted["CRFVOL"]), col=R_pal[[8]], main='100-200 cm', axes = FALSE, xlab="", ylab="", zlim=r, asp=1)
lines(admin.af, col="black")
legend("bottomleft", rx, fill=rev(R_pal[[8]][c(1,5,10,15,20)]), horiz=FALSE, bty="n")
dev.off()

save.image(file="CRFVOL_AfSIS.RData", compress="xz")

# end of script;