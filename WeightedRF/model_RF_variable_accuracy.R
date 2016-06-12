## Spatial prediction of clay content using measurements of variable accuracy (i.e. variable measurement error)
## Dataset from the NRCSS (http://ncsslabdatamart.sc.egov.usda.gov/) / NASIS soil profile DB.
## Tom.Hengl@isric.org

library(plyr)
library(rgdal)
library(quantregForest)
library(GSIF)
library(raster)
library(plotKML)
plotKML.env(convert="convert", show.env=FALSE)
options(bitmapType='cairo')

## load point data:
carson <- read.csv(file="carson_CLYPPT.csv")
## load covariates:
covs1km <- readRDS("carson_covs1km.rds")
coordinates(carson) <- ~X+Y
proj4string(carson) = covs1km@proj4string
writeOGR(carson, "carson_CLYPPT.shp", "carson_CLYPPT", "ESRI Shapefile")
## subset to speed up model fitting:
carson.sp <- sample.grid(carson, cell.size=c(5000,5000), n=2)
plot(raster(covs1km["DEMMRG5_1km"]))
points(carson.sp[[1]], pch="+")

## focus on mapping clay content:
summary(carson$CLYPPT)
## clay content has been determined with variable accuracy:
hist(carson$CLYPPT.sd)
## based on the CLYPPT.sd we can determine the weights as:
carson.sp[[1]]$w <- round(max(carson.sp[[1]]$CLYPPT.sd, na.rm=TRUE)^2/carson.sp[[1]]$CLYPPT.sd^2)
carson.sp[[1]]@data[1:5,]

## regression model:
fm <- as.formula(paste("CLYPPT ~", paste(names(covs1km), collapse="+")))
ov2 <- over(carson.sp[[1]], covs1km)
ov2 <- cbind(data.frame(carson.sp[[1]][c("CLYPPT","w")]), ov2)
ov2 <- plyr::rename(ov2, replace=c("X"="s1","Y"="s2"))
## multiply observations based on the measurement accuracy:
ov1 <- ov2[rep(row.names(ov2), ov2$w),]
## much larger data set now

## fit weighted RF:
mw1 <- fit.regModel(fm, ov1, covs1km, method="quantregForest")
pr.mw1 <- predict(mw1, covs1km, zmin=0, zmax=100, vgmmodel=NULL)
#plot(pr.mw1)
plot(stack(pr.mw1@predicted), col=SAGA_pal[[1]])
writeGDAL(pr.mw1@predicted["CLYPPT"], "mw1_pred.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")

## compare with the model without weights:
mw2 <- fit.regModel(fm, ov2, covs1km, method="quantregForest")
pr.mw2 <- predict(mw2, covs1km, zmin=0, zmax=100, vgmmodel=NULL)
#plot(pr.mw2)
plot(stack(pr.mw2@predicted), col=SAGA_pal[[1]])
writeGDAL(pr.mw2@predicted["CLYPPT"], "mw2_pred.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")

## plot predictions and errors next to each other:
par(mfrow=c(2,2))
par(mai=c(0.2,0,0.5,0), oma=c(0,0,0,0), xaxs='i', yaxs='i')
image(raster(pr.mw1@predicted[2]), col=SAGA_pal[[1]], main="Clay % (weighted)", asp=1, zlim=c(0,65), axes=FALSE, xlab="", ylab="")
image(raster(pr.mw2@predicted[2]), col=SAGA_pal[[1]], main="Clay %", asp=1, zlim=c(0,65), axes=FALSE, xlab="", ylab="")
image(raster(pr.mw1@predicted[1]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], main="Error (weighted)", asp=1, zlim=c(0,830), axes=FALSE, xlab="", ylab="")
points(carson.sp[[1]], pch="+")
image(raster(pr.mw2@predicted[1]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], main="Error", asp=1, zlim=c(0,830), axes=FALSE, xlab="", ylab="")
points(carson.sp[[1]], pch="+")

mean(pr.mw1@predicted@data[,1], na.rm=TRUE)
mean(pr.mw2@predicted@data[,1], na.rm=TRUE)
kml(pr.mw1@predicted[2], colour="CLYPPT", colour_scale=SAGA_pal[[1]])

## Conclusions:
## 1. not a large difference between the two (weighted / no weights) - slightly different covariates become more important
## 2. the errors for weighted predictions seem to be slightly smaller on the end (20%)
