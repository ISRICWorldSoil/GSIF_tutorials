## Afsoilmapper demo
## By: Tom.Hengl@isric.org and milan.kili11@gmail.com

library(spatstat)
library(rgdal)
library(raster)
library(sp)
library(gstat)
library(quantregForest)
library(ranger)
library(plotKML)
#install.packages("GSIF", repos="http://R-Forge.R-project.org")
library(GSIF) ## v0.5-2
library(RSAGA)
library(leaflet)

gdal.dir <- shortPathName("C:/Program files/GDAL")
gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe")

## Input parameters:
target <- readOGR("samples.gml", "samples")
## The GML driver does not support coordinate reference systems
proj4st <- "+proj=merc +lon_0=0 +lat_ts=6.64456744726493 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4string(target) <- proj4st
tvar <- "N"
target <- target[tvar]
pixel.size <- 100
t.dens <- 0.05
sigma <- 0.3*sqrt(diff(target@bbox[1,])*diff(target@bbox[2,]))
target.min <- min(target@data[,1])
target.max <- max(target@data[,1])
#method <- "quantregForest"
method <- "ranger"
grid2poly <- FALSE
dens.maps <- TRUE
vgmFun <- "Nug"
vgmmodel <- NULL
pal <- SAGA_pal[[1]]

## Standard covariates (could be geotifs, virtual mosaicks)
y <- c("H:\\AFSIS\\100m\\SRTMGL3_AF_100m.tif", "H:\\AFSIS\\100m\\GLC2010_100m.tif", "H:\\AFSIS\\100m\\PHIHOX_M_sl2_250m.tif")
factors = c(FALSE, TRUE, FALSE)
## Prepare covariates:
covariates <- makePixels(as(target, "SpatialPoints"), y, factors, pixel.size = pixel.size, sigma=sigma, t.dens=t.dens)
if(dens.maps==TRUE){
  grid.dist <- buffer.dist(target, covariates[1], cut(target@data[,1], breaks=seq(target.min, target.max, length=6)))
  covariates@data <- cbind(covariates@data, grid.dist@data)
}
## Fit model:
fm <- as.formula(paste(names(target)[1], "~", paste(names(covariates), collapse="+")))
m <- fit.gstatModel(target, fm, covariates, method=method, vgmFun=vgmFun)
p <- predict.gstatModel(m, covariates, vgmmodel=vgmmodel)

## Outputs:
## plot in GE
plotKML(p, file.name=paste0(names(target),".kml"), colour_scale=pal, grid2poly=grid2poly, obj.summary=FALSE)
## download GeoTiff

## Plot output in a browser:
r = raster(p@predicted[names(target)])
#plot(r)
pal <- colorNumeric(R_pal[["pH_pal"]], values(r), na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(r, colors=pal, opacity=0.6) %>%
  addLegend(pal=pal, values=values(r), title=names(target))
