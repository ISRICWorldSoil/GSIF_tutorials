# title         : TEX_1km_AfSIS.R
# purpose       : Mapping soil textures for Africa (1 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.isric.org/]
# producer      : Prepared by T. Hengl, G.B.M. Heuvelink and B. Kempen
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles "afsp.geo.rda" and WorldGrids maps
# outputs       : Predictions of silt sand and clay at 6 depths;
# remarks 1     : This is LARGE data! Not recommended for a PC without at least 16GB of RAM and multicore processor. The script is extra long because at many places we were forced to run things in loops i.e. tile by tile;


#setwd("E:\\gsif\\AFRICA\\1km")
#install.packages("GSIF", repos=c("http://R-Forge.R-project.org"))
library(GSIF)
library(maptools)
library(gstat)
library(aqp)
library(splines)

#load(file="TEX_1km_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
fw.path <- utils::readRegistry("SOFTWARE\\WOW6432NODE\\FWTools")$Install_Dir
gdalwarp <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdalbuildvrt <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))

## Load input data
## Points:
if(is.na(file.info("../samples/afsp.geo.rda")$size)){
  download.file("http://gsif.isric.org/lib/exe/fetch.php?media=afsp.geo.rda", "../samples/afsp.geo.rda")
}
load("../samples/afsp.geo.rda")
## Covariates:
if(is.na(file.info("afgrid1km.rda")$size)){
  download.file("http://worldgrids.org/rda/afgrid1km.rda", "afgrid1km.rda")
}
load("afgrid1km.rda")
## 5 km mask:
if(is.na(file.info("../soilmasks/afgrid5km.rda")$size)){
  download.file("http://gsif.isric.org/lib/exe/fetch.php?media=afgrid5km.rda", "../soilmasks/afgrid5km.rda")
}
load("../soilmasks/afgrid5km.rda")

logits = function(x){log((x/100)/(1-x/100))}
## fix values outside the psychical range:
sel = afsp.geo@data$methodid=="SLTPPT"|afsp.geo@data$methodid=="SNDPPT"|afsp.geo@data$methodid=="CLYPPT"
afsp.geo@data[sel,"observedValue"] <- ifelse(afsp.geo@data[sel,"observedValue"]==0, 1, ifelse(afsp.geo@data[sel,"observedValue"]==100, 99, afsp.geo@data[sel,"observedValue"]))

## define regression models:
formulaString.SLT <- as.formula(paste('logits(SLTPPT) ~ DEMSRE3a + SLPSRT3a + TDMMOD3a + TDSMOD3a + TNMMOD3a + TNSMOD3a + EVMMOD3a + EVSMOD3a + PREGSM1a +', paste("G0", c(1:7,9), "ESA3a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA3a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG3x + TWISRE3a + L3POBI3a'))
formulaString.SND <- as.formula(paste('logits(SNDPPT) ~ DEMSRE3a + SLPSRT3a + TDMMOD3a + TDSMOD3a + TNMMOD3a + TNSMOD3a + EVMMOD3a + EVSMOD3a + PREGSM1a +', paste("G0", c(1:7,9), "ESA3a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA3a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG3x + TWISRE3a + L3POBI3a'))
formulaString.CLY <- as.formula(paste('logits(CLYPPT) ~ DEMSRE3a + SLPSRT3a + TDMMOD3a + TDSMOD3a + TNMMOD3a + TNSMOD3a + EVMMOD3a + EVSMOD3a + PREGSM1a +', paste("G0", c(1:7,9), "ESA3a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA3a", sep="", collapse="+"), '+ ns(altitude, df=2) + SGEUSG3x + TWISRE3a + L3POBI3a'))

## fit the models (takes ca 3 x 2 mins):
SLT.m <- fit.gstatModel(observations=afsp.geo, formulaString.SLT, covariates=afgrid1km, dimensions="3D")
gc()
SND.m <- fit.gstatModel(observations=afsp.geo, formulaString.SND, covariates=afgrid1km, dimensions="3D")
gc()
CLY.m <- fit.gstatModel(observations=afsp.geo, formulaString.CLY, covariates=afgrid1km, dimensions="3D")
gc()

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
CLY.m@vgmModel$psill = c(1.05, 0.3)
v3 = CLY.m@vgmModel
class(v3) = c("variogramModel", "data.frame")
plot(variogram(residual~1, pnts3[runif(length(pnts3$residual))<.2,]), v3)

## Remove classes not represented otherwise the predict.glm reports an error of type "factor 'SGEUSG0x' has new levels"!
fix.c = levels(afgrid1km$SGEUSG3x)[!(levels(afgrid1km$SGEUSG3x) %in% levels(SLT.m@regModel$model$SGEUSG3x))]
for(j in fix.c){
  afgrid1km$SGEUSG3x[afgrid1km$SGEUSG3x == j] <- "pCm"
}
## takes ca 3 mins!
gc()

m.lst <- list(SLT.m, SND.m, CLY.m)
n.lst <- c("SLTPPT", "SNDPPT", "CLYPPT")
v.lst <- list(v1, v2, v3)

## back-transformation function:
invlogit = function(x){exp(x)/(1+exp(x))*100}
## for the back-transformation for the mean value see Diggle and Ribeiro, 2007, p. 148:
invlogit.m = function(x, v){((1+exp(-x))^(-1)-.5*v*exp(-x)*(1-exp(-x))*(1+exp(-x))^(-3) )*100}

for(k in 1:length(m.lst)){
  ## 3D points:
  observed <- SpatialPointsDataFrame(m.lst[[k]]@sp, data=data.frame(observedValue=m.lst[[k]]@regModel$model[,1]))
  observed@data[,"residuals"] <- residuals(m.lst[[k]]@regModel)
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
  for(i in 1:6){
    if(is.na(file.info(paste(n.lst[k], "_sd", i, "_M.tif", sep=""))$size)){
    afgrid1km@data[,"altitude"] <- rep(get("stdepths", envir = GSIF.opts)[i], nrow(afgrid1km))
    gc()
    ## regression part (in tiles):
    for(j in 1:L){  ## takes about 10 mins and > 8GB!
      rp <- stats::predict.glm(m.lst[[k]]@regModel, newdata=afgrid1km@data[sel[j]:sel[j+1],], type="response", se.fit = TRUE, na.action = na.pass)
      gc()
      ## copy values:
      rkgrid1km@data[sel[j]:sel[j+1], paste(n.lst[k], "_sd", i, "_glm", sep="")] <- rp$fit
      rkgrid1km@data[sel[j]:sel[j+1], paste(n.lst[k], "_sd", i, "_var", sep="")] <- rp$se.fit^2
      gc()
    }
  
    ## prediction locations for OK:  
    new3D5km <- sp3D(afgrid5km, stdepths = get("stdepths", envir = GSIF.opts)[i], stsize = get("stsize", envir = GSIF.opts)[i])[[1]]
    suppressWarnings(proj4string(new3D5km) <- observed@proj4string)    
    ## subset points to speed up kriging:
    subset.observed = observed[observed@coords[,3] > new3D5km@bbox[3,1]-.15 & observed@coords[,3] < new3D5km@bbox[3,2]+.15,]
    rkgrid5km <- gstat::krige(residuals~1, locations=subset.observed, newdata=new3D5km, model = v.lst[[k]], nmin = 80, nmax = 100, debug.level = -1)
    ## down-scale to 1 km:
    rk <- gdalwarp(rkgrid5km, GridTopology=afgrid1km@grid, pixsize=afgrid1km@grid@cellsize[1], resampling_method="cubicspline", tmp.file=TRUE)
    ## convert to SpatialPixels (necessary!):
    rk <- SpatialPixelsDataFrame(afgrid1km@coords[,1:2], data=rk@data[afgrid1km@grid.index,], proj4string=afgrid1km@proj4string, grid=afgrid1km@grid)
  
    ## sum the trend and residual part:
    gc()
    rkgrid1km@data[, paste(n.lst[k], "_sd", i, "_var", sep="")] <- rkgrid1km@data[,paste(n.lst[k], "_sd", i, "_var", sep="")]+ rk$var1.var
    gc()
    rkgrid1km@data[, paste(n.lst[k], "_sd", i, "_M", sep="")] <- round(invlogit.m(rkgrid1km@data[,paste(n.lst[k], "_sd", i, "_glm", sep="")]+ rk$var1.pred, rkgrid1km@data[, paste(n.lst[k], "_sd", i, "_var", sep="")]))
    gc()
    rkgrid1km@data[, paste(n.lst[k], "_sd", i, "_L", sep="")] <- round(invlogit(rkgrid1km@data[,paste(n.lst[k], "_sd", i, "_glm", sep="")]+ rk$var1.pred - 1.645*sqrt(rkgrid1km@data[, paste(n.lst[k], "_sd", i, "_var", sep="")])))
    gc()
    rkgrid1km@data[, paste(n.lst[k], "_sd", i, "_U", sep="")] <- round(invlogit(rkgrid1km@data[,paste(n.lst[k], "_sd", i, "_glm", sep="")]+ rk$var1.pred + 1.645*sqrt(rkgrid1km@data[, paste(n.lst[k], "_sd", i, "_var", sep="")])))
  
    ## write to GeoTiff:
    writeGDAL(rkgrid1km[paste(n.lst[k], "_sd", i, "_M", sep="")], paste(n.lst[k], "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
    gc()
    writeGDAL(rkgrid1km[paste(n.lst[k], "_sd", i, "_L", sep="")], paste(n.lst[k], "_sd", i, "_L.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
    gc()
    writeGDAL(rkgrid1km[paste(n.lst[k], "_sd", i, "_U", sep="")], paste(n.lst[k], "_sd", i, "_U.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
    gc()
  }
  }
}
## time elapsed: 3 x 147 minutes!!

## ensure that the values for M sum to 100%:
for(i in 1:6){
  gc()
  x = readGDAL(paste(n.lst[1], "_sd", i, "_M.tif", sep=""))
  x@data[,2] <- readGDAL(paste(n.lst[2], "_sd", i, "_M.tif", sep=""))$band1
  x@data[,3] <- readGDAL(paste(n.lst[3], "_sd", i, "_M.tif", sep=""))$band1
  names(x) <- n.lst
  gc()
  sums <- rowSums(x@data)
  xr = range(sums, na.rm=TRUE)
  gc()
  if(xr[1]<100|xr[2]>100){
    x$SLTPPT <- round(x$SLTPPT / sums * 100, 0)
    x$SNDPPT <- round(x$SNDPPT / sums * 100, 0)
    x$CLYPPT <- round(x$CLYPPT / sums * 100, 0)    
    writeGDAL(x[n.lst[1]], paste(n.lst[1], "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
    writeGDAL(x[n.lst[2]], paste(n.lst[2], "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
    writeGDAL(x[n.lst[3]], paste(n.lst[3], "_sd", i, "_M.tif", sep=""), "GTiFF", mvFlag=-9999, type="Int16")
    gc()
  }
}

## compress all produced maps:
system("7za a SNDPPT_1km_glmrk.tif.7z SNDPPT_sd*.tif")
system("7za a SLTPPT_1km_glmrk.tif.7z SLTPPT_sd*.tif")
system("7za a CLYPPT_1km_glmrk.tif.7z CLYPPT_sd*.tif")
## add a readme file:
#system("7za a SNDPPT_1km_glmrk.tif.7z README_SNDPPT.txt")
## save the model:
save(SLT.m, file="SLT.m.rda", compress="xz")
save(SND.m, file="SND.m.rda", compress="xz")
save(CLY.m, file="CLY.m.rda", compress="xz")

#save.image(file="TEX_1km_AfSIS.RData")

# end of script;