#########################################
#                                       #
# 3D+T demo of the cookfarm data        #
#                                       #
# The demo script has been prepared by: #
# Tomislav Hengl                        #
# ISRIC --- World Soil Information      #
# tom.hengl@wur.nl                      #
#                                       #
#########################################

library(rgdal)
library(sp)
library(spacetime)
library(aqp)
library(splines)
library(randomForest)
library(plyr)
library(plotKML)
## load the data set:
data(cookfarm)

## gridded data:
grid10m <- cookfarm$grids
gridded(grid10m) <- ~x+y
proj4string(grid10m) <- CRS(cookfarm$proj4string)
spplot(grid10m["DEM"], col.regions=SAGA_pal[[1]])

## soil profiles:
profs <- cookfarm$profiles
levels(cookfarm$profiles$HZDUSD)
## Bt horizon:
sel.Bt <- grep("Bt", profs$HZDUSD, ignore.case=FALSE, fixed=FALSE)
profs$Bt <- 0
profs$Bt[sel.Bt] <- 1
depths(profs) <- SOURCEID ~ UHDICM + LHDICM
site(profs) <- ~ TAXSUSDA + Easting + Northing
coordinates(profs) <- ~Easting + Northing
proj4string(profs) <- CRS(cookfarm$proj4string)
profs.geo <- GSIF::as.geosamples(profs)

## fit models for Bt horizon and soil pH:
m.Bt <- GSIF::fit.gstatModel(profs.geo, Bt~DEM+TWI+MUSYM+Cook_fall_ECa
   +Cook_spr_ECa+ns(altitude, df = 4), grid10m, fit.family = binomial(logit))
plot(m.Bt)
## fix vgm manually:
m.Bt@vgmModel$psill <- c(0.37,0.09)
m.Bt@vgmModel$range <- c(0,60)
plot(m.Bt)
m.PHI <- GSIF::fit.gstatModel(profs.geo, PHIHOX~DEM+TWI+MUSYM+Cook_fall_ECa
    +Cook_spr_ECa+ns(altitude, df = 4), grid10m)
plot(m.PHI)

## prepare 3D locations:
s <- c(-.3,-.6,-.9,-1.2,-1.5)
new3D <- GSIF::sp3D(grid10m, stdepths=s)
## predict using 3D regression-kriging:
rk.Bt <- lapply(new3D, function(x){
  predict(m.Bt, predictionLocations=x, nfold=0, nmax=50)
  })
rk.PHI <- lapply(new3D, function(x){
  predict(m.PHI, predictionLocations=x, nfold=0, nmax=50)
  })
## visualize 3D pH maps in Google Earth:
z0 = mean(grid10m$DEM, na.rm=TRUE)
## export grids:
for(j in 1:length(rk.PHI)){
  kml( slot(rk.PHI[[j]], "predicted")["PHIHOX"],
     folder.name=paste("cookfarm_PHIHOX_", j, sep=""),
     file=paste("cookfarm_PHIHOX_", j, ".kml", sep=""),
     colour=PHIHOX, colour_scale=R_pal[["pH_pal"]],
     z.lim=c(5,7.2),
     raster_name=paste("cookfarm_PHIHOX_", j, ".png", sep=""),
     altitude=z0+400+(s[j]*200) )
}
## Copy predictions:
for(j in 1:length(rk.Bt)){
  grid10m@data[,paste0("Bt_Port", j)] <- slot(rk.Bt[[j]],
  "predicted")@data[,"Bt"]
}
for(j in 1:length(rk.PHI)){
  grid10m@data[,paste0("PHI_Port", j)] <- slot(rk.PHI[[j]],
    "predicted")@data[,"PHIHOX"]
}
## predicted Bt horizon and soil pH:
spplot(grid10m[rev(paste0("Bt_Port", 1:5))],
   col.regions=SAGA_pal[["SG_COLORS_YELLOW_RED"]])
spplot(grid10m[rev(paste0("PHI_Port", 1:5))],
   col.regions=R_pal[["pH_pal"]])

## prepare spacetime regression-matrix:
d.xy <- data.frame(profs@sp)
d.xy$SOURCEID <- profs@site$SOURCEID
## Cumulative rainfall:
cookfarm$weather$Precip_cum <- ave(cookfarm$weather$Precip_wrcc,
   rev(cumsum(rev(cookfarm$weather$Precip_wrcc)==0)), FUN=cumsum)
## focus on volumetric water content:
t <- cookfarm$readings[,c("SOURCEID","Date",paste0("Port",1:5,"VW"))]
VW.df <- plyr::join(data.frame(
   SOURCEID=rep(t$SOURCEID, 5), Date=rep(t$Date, 5),
   VW=signif(as.vector(unlist(t[,paste0("Port",1:5,"VW")])), 3),
   altitude=as.vector(unlist(sapply(s, function(x){
       rep(x, length(t$SOURCEID))
       })))), d.xy, by="SOURCEID", type="left")
## Create a space time object:
stVW <- STIDF(SpatialPoints(VW.df[,c("Easting","Northing","altitude")],
   proj4string=grid10m@proj4string),
   time=VW.df$Date-.5, data=VW.df[,c("SOURCEID","VW")],
   endTime=as.POSIXct(VW.df$Date+.5))

## prepare spacetime covariates
## 3D static covariates:
grid10m.3D <- data.frame(DEM=rep(grid10m$DEM, 5),
   TWI=rep(grid10m$TWI, 5), NDRE.M=rep(grid10m$NDRE.M, 5),
   NDRE.Sd=rep(grid10m$NDRE.Sd, 5),
   Bt=as.vector(unlist(grid10m@data[,paste0("Bt_Port",1:5)])),
   PHI=as.vector(unlist(grid10m@data[,paste0("PHI_Port",1:5)])),
   Easting=rep(grid10m@coords[,1], 5), Northing=rep(grid10m@coords[,2], 5),
   altitude=as.vector(unlist(sapply(s, function(x){
      rep(x, length(grid10m$DEM))
      })))
)
coordinates(grid10m.3D) <- ~ Easting + Northing + altitude
gridded(grid10m.3D) = TRUE
proj4string(grid10m.3D) <- grid10m@proj4string
sp.grid <- as(grid10m.3D, "SpatialPixels")
begin <- min(c(cookfarm$readings$Date, cookfarm$readings$Date), na.rm=TRUE)
endTime <- as.POSIXct(max(c(cookfarm$readings$Date,
    cookfarm$readings$Date), na.rm=TRUE))
grid10m.st0 <- STFDF(sp.grid, time=begin, data=grid10m.3D@data,
   endTime=endTime)
## Overlay in space and time (thanx to spacetime package):
ov.s0 <- over(stVW, grid10m.st0)
ov.d <- plyr::join(data.frame(Date=time(stVW)), cookfarm$weather, type="left")
## Prepare the st regression matrix:
rm.lst <- do.call(cbind, list(stVW@data, as.data.frame(stVW@sp@coords),
   ov.s0, ov.d))
## Estimate cumulative days:
rm.lst$cday <- floor(unclass(stVW@endTime)/86400-.5)
rm.lst <- rm.lst[complete.cases(rm.lst),]
rm.lst$cdayt <- cos((rm.lst$cday-min(rm.lst$cday))*pi/180)

## Fit space-time regression models:
fm <- VW ~ DEM+TWI+NDRE.M+NDRE.Sd+Bt+PHI+Precip_cum+MaxT_wrcc+MinT_wrcc+cdayt
fit.LM <- lm(fm, rm.lst)
summary(fit.LM)
## This is a LARGE data set and fitting random forest could take hours...
subs <- runif(nrow(rm.lst))<.02
fit.RF <- randomForest(fm, rm.lst[subs,], importance=TRUE)
varImpPlot(fit.RF, main="Water content (volumetric)", bg="blue", pch=22)

## predict values using spacetime RF model:
dates.lst <- seq(as.Date("2012-07-01"),as.Date("2012-07-30"), by=1)
xy <- attr(grid10m.3D@coords, "dimnames")[[2]]
for(j in 1:length(dates.lst)){
  date <- dates.lst[j]
  newD <- cbind(data.frame(grid10m.3D),
     data.frame(Date=rep(as.Date(date), length(grid10m.3D))))
  newD <- plyr::join(newD, cookfarm$weather, type="left")
  newD$cday <- floor(unclass(newD$Date)/86400-.5)
  newD$cdayt <- cos((newD$cday-min(rm.lst$cday))*pi/180)
  outn <- paste("VW", gsub("-", "_", date), sep="_")
  newD[,outn] <- predict(fit.RF, newD)
  for(k in 1:length(s)){
    tif.out <- paste0(outn, "_Port_", k, ".tif")
    outP <- newD[newD$altitude==s[k], c(xy,outn)]
    gridded(outP) <- as.formula(paste("~ ", xy[1], "+", xy[2]))
    proj4string(outP) <- grid10m@proj4string
    writeGDAL(outP[outn], tif.out, "GTiff", mvFlag=-99999)
  }
}