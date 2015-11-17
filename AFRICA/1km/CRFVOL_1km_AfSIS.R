# title         : CRFVOL_1km_AfSIS.R
# purpose       : Mapping volume percentage of coarse fragments (> 2 mm) (1 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl, G.B.M. Heuvelink and B. Kempen
# last update   : In Wageningen, NL, March 2013.
# inputs        : Africa soil profiles "afsp.rda" and WorldGrids maps
# outputs       : Predictions of pH at 6 depths;
# remarks 1     : This is LARGE data! Not recommended for a PC without at least 16GB of RAM and multicore processor. The script is extra long because at many places we were forced to run things in loops i.e. tile by tile;

setwd("G:\\AFSIS\\1km\\")
#system("R CMD INSTALL GSIF_0.3-0.tar.gz")
#detach(package:GSIF)
library(GSIF)
library(maptools)
library(plotKML)
data(SAGA_pal)
library(gstat)
library(aqp)
#load(file="CRFVOL_1km_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
fw.path <- utils::readRegistry("SOFTWARE\\WOW6432NODE\\FWTools")$Install_Dir
gdalwarp <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdalbuildvrt <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))

## Load input data
## Points:
load("../samples/afsp.geo.rda")
## Covariates:
load("afgrid1km.rda")
## 5 km mask:
load("../soilmasks/afgrid5km.rda")

## TH: this is highly skewed variable hence it probably requires a zero-inflated model [http://www.ats.ucla.edu/stat/r/dae/zipoisson.htm]!  

## Fit model:
ov.afsp <- overlay(x=afgrid1km, y=afsp.geo, method="CRFVOL", var.type = "numeric")
names(ov.afsp)[which(names(ov.afsp)=="observedValue")] = "CRFVOL"
## Log-model:
formulaString.CRF <- as.formula(paste('CRFVOL ~ DEMSRE3a + SLPSRT3a + TDMMOD3a + TDSMOD3a + TNMMOD3a + TNSMOD3a + EVMMOD3a + EVSMOD3a + PREGSM1a +', paste("G0", c(1:7,9), "ESA3a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA3a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG3x + TWISRE3a + L3POBI3a'))
CRF.m <- fit.regModel(formulaString.CRF, rmatrix=ov.afsp, afgrid1km, method="GLM", family=poisson, stepwise=TRUE, dimensions = "3D")
## takes >3 minutes;
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

## Remove classes not represented otherwise the predict.glm reports an error of type "factor 'SGEUSG0x' has new levels"!
fix.c = levels(afgrid1km$SGEUSG3x)[!(levels(afgrid1km$SGEUSG3x) %in% levels(CRF.m@regModel$model$SGEUSG3x))]
for(j in fix.c){
  afgrid1km$SGEUSG3x[afgrid1km$SGEUSG3x == j] <- "pCm"
}
## takes ca 3 mins!
gc()

## 3D points:
observed <- SpatialPointsDataFrame(CRF.m@sp, data=data.frame(observedValue=CRF.m@regModel$model[,1]))
observed@data[,"residuals"] <- residuals(CRF.m@regModel)
#observed[1:3,]
## remove duplicates:
if(length(zerodist(observed))>0){
  observed <- remove.duplicates(observed)
}

## output grid:
rkgrid1km <- afgrid1km["SMKMOD3x"]
# number of tiles:
L = 10
sel <- c(seq(1, nrow(afgrid1km), by=round(nrow(afgrid1km)/L)), nrow(afgrid1km))
gc()

## predict values at each depth:
system.time(for(i in 1:6){
  if(is.na(file.info(paste("CRFVOL", "_sd", i, "_M.tif", sep=""))$size)){
  afgrid1km@data[,"altitude"] <- rep(get("stdepths", envir = GSIF.opts)[i], nrow(afgrid1km))
  gc()
  ## regression part (in tiles):
  for(j in 1:L){  ## takes about 10 mins and > 8GB!
    rp <- stats::predict.glm(CRF.m@regModel, newdata=afgrid1km@data[sel[j]:sel[j+1],], type="response", se.fit = TRUE, na.action = na.pass)
    gc()
    ## copy values:
    rkgrid1km@data[sel[j]:sel[j+1], paste("CRFVOL", "_sd", i, "_glm", sep="")] <- rp$fit
    rkgrid1km@data[sel[j]:sel[j+1], paste("CRFVOL", "_sd", i, "_var", sep="")] <- rp$se.fit^2
    gc()
  }
  
  ## prediction locations for OK:  
  new3D5km <- sp3D(afgrid5km, stdepths = get("stdepths", envir = GSIF.opts)[i], stsize = get("stsize", envir = GSIF.opts)[i])[[1]]
  suppressWarnings(proj4string(new3D5km) <- observed@proj4string)
  ## subset points to speed up kriging:
  subset.observed = observed[observed@coords[,3] > new3D5km@bbox[3,1]-.15 & observed@coords[,3] < new3D5km@bbox[3,2]+.15,]
  rkgrid5km <- gstat::krige(residuals~1, locations=subset.observed, newdata=new3D5km, model = v, nmin = 80, nmax = 100, debug.level = -1, block = 1000)
  ## down-scale to 1 km:
  rk <- gdalwarp(rkgrid5km, GridTopology=afgrid1km@grid, pixsize=afgrid1km@grid@cellsize[1], resampling_method="cubicspline", tmp.file=TRUE)
  ## convert to SpatialPixels (necessary!):
  rk <- SpatialPixelsDataFrame(afgrid1km@coords[,1:2], data=rk@data[afgrid1km@grid.index,], proj4string=afgrid1km@proj4string, grid=afgrid1km@grid)
  
  ## sum the trend and residual part:
  #gc()
  #rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_var", sep="")] <- rkgrid1km@data[,paste("CRFVOL", "_sd", i, "_var", sep="")]+ rk$var1.var
  gc()
  rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_M", sep="")] <- round(rkgrid1km@data[,paste("CRFVOL", "_sd", i, "_glm", sep="")] + rk$var1.pred)
  gc()
  rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_L", sep="")] <- round(rkgrid1km@data[,paste("CRFVOL", "_sd", i, "_glm", sep="")] + rk$var1.pred - 1.645*sqrt(rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_var", sep="")]))
  rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_L", sep="")] <- ifelse(rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_L", sep="")] < 0, 0, rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_L", sep="")])
  gc()
  rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_U", sep="")] <- round(rkgrid1km@data[,paste("CRFVOL", "_sd", i, "_glm", sep="")] + rk$var1.pred + 1.645*sqrt(rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_var", sep="")]))
  rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_U", sep="")] <- ifelse(rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_L", sep="")] >99, 100, rkgrid1km@data[, paste("CRFVOL", "_sd", i, "_U", sep="")])
  
  ## write to GeoTiff:
  writeGDAL(rkgrid1km[paste("CRFVOL", "_sd", i, "_M", sep="")], paste("CRFVOL", "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-99999)
  gc()
  writeGDAL(rkgrid1km[paste("CRFVOL", "_sd", i, "_L", sep="")], paste("CRFVOL", "_sd", i, "_L.tif", sep=""), "GTiFF", mvFlag=-99999)
  gc()
  writeGDAL(rkgrid1km[paste("CRFVOL", "_sd", i, "_U", sep="")], paste("CRFVOL", "_sd", i, "_U.tif", sep=""), "GTiFF", mvFlag=-99999)
  gc()
  #writeGDAL(rkgrid1km[paste("CRFVOL", "_sd", i, "_var", sep="")], paste("CRFVOL", "_sd", i, "_var.tif", sep=""), "GTiFF", mvFlag=-99999)
  #gc()
  }
})
## time elapsed: 147 minutes!!

## compress all produced maps:
system("7za a CRFVOL_1km_glmrk.tif.7z CRFVOL_sd*.tif")
## add a readme file:
#system("7za a CRFVOL_1km_glmrk.tif.7z README_CRFVOL.txt")

save.image(file="CRFVOL_AfSIS.RData")

# end of script; 