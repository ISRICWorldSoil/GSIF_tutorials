# title         : spatial_predictors_meuse.R
# purpose       : Comparison of spatial predictors implemented in gstat (includes idea for the Best Combined Spatial Predictor)
# reference     : Hengl, T. 2009. A Practical Guide to Geostatistical Mapping, 2nd Edt. University of Amsterdam, www.lulu.com, 291 p. ISBN: ISBN 978-90-9024981-0
# producer      : T. Hengl
# address       : In Wageningen
# inputs        : meuse dataset (comes with gstat);
# outputs       : predictions of zinc, assessment of the mapping accuracy;
# remarks 1     : to get the maps right in GE, you will need to fix the epsg:28992 manually in "library/rgdal/proj/epsg"
# remarks 2     : live article at [http://spatial-analyst.net/wiki/index.php?title=Best_Combined_Spatial_Predictor]

library(gstat)
library(maptools)
library(rgdal)
library(lattice)
library(GSIF)
library(plotKML)
data(SAGA_pal)
#trellis.par.set(sp.theme())
library(RCurl)
urlExists = url.exists("http://spatialreference.org/ref/sr-org/6781/proj4/")
if(urlExists){ 
  nl.rd <- CRS(getURL("http://spatialreference.org/ref/sr-org/6781/proj4/"))
} else {
  nl.rd <- CRS("+init=epsg:28992")
}

## load data:
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- nl.rd
## load grids:
data(meuse.grid)
meuse.grid$X <- meuse.grid$x
meuse.grid$Y <- meuse.grid$y
meuse.grid$soil <- as.factor(meuse.grid$soil)
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE
fullgrid(meuse.grid) <- TRUE
proj4string(meuse.grid) <- nl.rd
# get values from grids to points:
meuse.ov <- over(meuse, meuse.grid)
meuse.ov <- cbind(as.data.frame(meuse), meuse.ov)
names(meuse.ov)

##------------------------------------------------------------
## Generate predictions
##------------------------------------------------------------

## inverse distances:
zinc.id <- krige(zinc~1, meuse, meuse.grid)
## regression on coordinates (trend):
zinc.tr <- krige(zinc~1, meuse, meuse.grid, degree=2, nmax=20)
zinc.tr$zinc.pred <- ifelse(zinc.tr$var1.pred<0, 0, zinc.tr$var1.pred)
## linear regression:
zinc.lm <- krige(log1p(zinc)~sqrt(dist), meuse, meuse.grid)
zinc.lm$zinc.pred <- expm1(zinc.lm$var1.pred) 
## ordinary kriging: 
vt.fit <- fit.variogram(variogram(log1p(zinc)~1, meuse), vgm(1, "Exp", 300, 1))
plot(variogram(log1p(zinc)~1, meuse), vt.fit)
zinc.ok <- krige(log1p(zinc)~1, meuse, meuse.grid, vt.fit)
zinc.ok$zinc.pred <- expm1(zinc.ok$var1.pred) 
## normalized error:
zinc.ok$svar <- zinc.ok$var1.var/var(log1p(meuse$zinc))
## GLM prediction (all possible predictors):
m.glm <- glm(zinc~sqrt(dist)+soil+ffreq, meuse.ov, family=gaussian(link=log))
summary(m.glm)
p.glm <- predict(m.glm, newdata=meuse.grid, type="response", se.fit=TRUE)
zinc.glm <- as(meuse.grid["soil"], "SpatialPointsDataFrame")
zinc.glm$var1.pred <- p.glm$fit
zinc.glm$var1.var <- p.glm$se.fit^2 + p.glm[["residual.scale"]]^2
zinc.glm$svar <- zinc.glm$var1.var/(m.glm$null.deviance/m.glm$df.null)
summary(zinc.glm$svar)
gridded(zinc.glm) <- TRUE
fullgrid(zinc.glm) <- TRUE
## regression-kriging:
vr.fit <- fit.variogram(variogram(log1p(zinc)~sqrt(dist)+soil+ffreq, meuse), vgm(1, "Exp", 300, 1))
plot(variogram(log1p(zinc)~sqrt(dist)+soil+ffreq, meuse), vr.fit)
zinc.rk <- krige(log1p(zinc)~sqrt(dist)+soil+ffreq, meuse, meuse.grid, vr.fit)
zinc.rk$zinc.pred <- expm1(zinc.rk$var1.pred) 

## plot all predictions next to each other:
meuse.grid$zinc.id <- zinc.id$var1.pred
meuse.grid$zinc.tr <- zinc.tr$zinc.pred
meuse.grid$zinc.lm <- zinc.lm$zinc.pred
meuse.grid$zinc.ok <- zinc.ok$zinc.pred
meuse.grid$zinc.glm <- zinc.glm$var1.pred
meuse.grid$zinc.rk <- zinc.rk$zinc.pred

spplot(meuse.grid[c("zinc.ok", "zinc.glm", "zinc.rk", "zinc.id", "zinc.tr", "zinc.lm")], col.regions=SAGA_pal[[1]])
## grey(rev(seq(0,1,0.025)))

## In the GSIF package we run everything at once:
mod = log1p(zinc)~sqrt(dist)+soil+ffreq
zinc.m <- fit.gstatModel(formulaString=mod, meuse, as(meuse.grid, "SpatialPixelsDataFrame"))
zinc.rk2 <- predict(zinc.m, as(meuse.grid, "SpatialPixelsDataFrame"))
zinc.rk2
zinc.rk2@observed@data$zinc <- signif(zinc.rk2@observed@data$zinc, 3)
plotKML(zinc.rk2, png.width = gridparameters(meuse.grid)[1,"cells.dim"]*5, png.height = gridparameters(meuse.grid)[2,"cells.dim"]*5)

##------------------------------------------------------------
## Spatial interpolation comparison
##------------------------------------------------------------

## Comparison of the techniques (PCA):

pc.preds <- prcomp(~zinc.id+zinc.tr+zinc.lm+zinc.ok+zinc.glm+zinc.rk, scale=F, meuse.grid)
biplot(pc.preds, cex=0.9, xlabs=rep(".", length(pc.preds$x[,1])), xlim=c(-0.06,0.025), ylim=c(-0.04,0.03))
## this shows that ID and OK, and LM and GLM are highly correlated (as expected); RK lies in-between (as expected); all predictors are positively correlated (of course);

## Mapping accuracy (cross-validation):
zinc.ok.cv <- krige.cv(log1p(zinc)~1, meuse, vt.fit, nfold=5)
zinc.rk.cv <- krige.cv(log1p(zinc)~sqrt(dist)+soil+ffreq, meuse, vr.fit, nfold=5)
## portion of variance exaplained by OK:
(1-var(zinc.ok.cv$residual, na.rm=T)/var(log1p(meuse$zinc), na.rm=T))*100
## portion of variance exaplained by RK:
(1-var(zinc.rk.cv$residual, na.rm=T)/var(log1p(meuse$zinc), na.rm=T))*100
## statistical comparison:
t.test(zinc.ok.cv$residual, zinc.rk.cv$residual) 
## the methods do not really differ significantly ?

## computing time:
system.time(krige(log1p(zinc)~1, meuse, meuse.grid, vt.fit))[3]
system.time(krige(log1p(zinc)~sqrt(dist)+soil+ffreq, meuse, meuse.grid, vr.fit))[3]
## RK takes a bit more time than OK;

##------------------------------------------------------------
## Best Combined Spatial Prediction (OK and GLM):
##------------------------------------------------------------

meuse.grid$zinc.BCSP <- (zinc.glm$var1.pred/zinc.glm$svar+zinc.ok$zinc.pred/zinc.ok$svar)/(1/zinc.glm$svar+1/zinc.ok$svar)
## plots next to each other:
spplot(meuse.grid[c("zinc.ok", "zinc.BCSP", "zinc.glm")], col.regions=SAGA_pal[[1]], sp.layout=list("sp.points", meuse, col="black", cex=.5))
BSCP.pol <- grid2poly(meuse.grid["zinc.BCSP"])
kml(BSCP.pol, colour=zinc.BCSP, colour_scale=SAGA_pal[[1]], kmz=TRUE)
kml(meuse, colour=zinc, labels=zinc, balloon=TRUE, kmz=TRUE)

## Comparison BCSP and RK
## derive residuals in response scale:
pnts <- m.glm$data
coordinates(pnts) <- ~ x+y
proj4string(pnts) <- nl.rd
pnts$glm.res <-  m.glm$data$zinc - fitted.values(m.glm)
vrglm.fit <- fit.variogram(variogram(glm.res~1, pnts), vgm(1, "Exp", 300, 1))
plot(variogram(glm.res~1, pnts), vrglm.fit)
zinc.okglm <- krige(glm.res~1, pnts, meuse.grid, vrglm.fit)
## sum the two GLM predictions and residuals:
meuse.grid$zinc.rkglm <- zinc.glm$var1.pred + zinc.okglm$var1.pred 
meuse.grid$zinc.rkglm <- ifelse(meuse.grid$zinc.rkglm<=0, 0, meuse.grid$zinc.rkglm)

## compare the two approaches:
scatter.smooth(log1p(meuse.grid$zinc.rkglm), log1p(meuse.grid$zinc.BCSP), col="grey")
spplot(meuse.grid[c("zinc.rkglm", "zinc.BCSP")], col.regions=SAGA_pal[[1]], sp.layout=list("sp.points", meuse, col="black", cex=.5))
## The BSCP gives more emphasis on the spatial auto-correlation of zinc (i.e. on OK predictions). This is possibly because the auto-correlation in zinc reflects partially auto-correlation in predictors. 

## export maps to ASCII format:
meuse.grid$mask <- ifelse(is.na(meuse.grid$soil), NA, 1)
meuse.grid$sqrtdist <- sqrt(meuse.grid$dist)
write.asciigrid(meuse.grid["mask"], "mask_map.asc", na.value=-1)
write.asciigrid(meuse.grid["sqrtdist"], "sqrtdist.asc", na.value=-1)
writeOGR(meuse, "meuse.shp", "meuse", "ESRI Shapefile")

##------------------------------------------------------------
## Two-stage sampling:
##------------------------------------------------------------

## stage 1 - transects:
t1 <- Line(matrix(c(178504,179517,330218,329776), ncol=2))
t2 <- Line(matrix(c(179129,180455,331115,330402), ncol=2))
t3 <- Line(matrix(c(179585,180666,331945,331591), ncol=2))
t4 <- Line(matrix(c(180570,181114,332692,332420), ncol=2))
t5 <- Line(matrix(c(181094,181447,333630,333121), ncol=2))
transects <- SpatialLines(list(Lines(list(t1, t2, t3, t4, t5), ID=c("t1", "t2", "t3", "t4", "t5"))), CRS("+init=epsg:28992"))
pts <- spsample(transects, n=50, type="regular")
pts.ov <- overlay(meuse.grid, pts)
transects.sp <- list("sp.points", pts.ov, col="black")
spplot(meuse.grid["dist"], sp.layout=transects.sp)

## study area of interest:
write.asciigrid(meuse.grid["mask"], "meuse_mask.asc")
rsaga.esri.to.sgrd(in.grids="meuse_mask.asc", out.sgrd="meuse_mask.sgrd", in.path=getwd())
## raster to polygon conversion;
rsaga.geoprocessor(lib="shapes_grid", module=6, param=list(GRID="meuse_mask.sgrd", SHAPES="meuse_mask.shp", CLASS_ALL=1))
mask <- readShapePoly("meuse_mask.shp", proj4string=CRS("+init=epsg:28992"), force_ring=T)


## stage 2 - optimized sample:
library(intamapInteractive)
## (spatial simulated annealing) and criterion "mukv"
observations <- pts.ov[c("zinc.rk", "dist")]
names(observations) <- c("zinc", "dist")
attr(observations@coords, "dimnames")[[2]] <- c("x", "y")
predGrid <- data.frame(meuse.grid["dist"])
predGrid$dist <- NULL
coordinates(predGrid) <- ~x+y

## OK model:
observations <- pts.ov["zinc.rk"]
names(observations) <- "zinc"
proj4string(mask) <- CRS(as.character(NA)); proj4string(observations) <- CRS(as.character(NA)); proj4string(predGrid) <- CRS(as.character(NA))
## compute the Mukv of the current network
calculateMukv(observations, predGrid, model=vt.fit, nmax=60)  # formulaString=log1p(zinc)~1
windows()
## add 100 more points (SSA method):
optim.ok <- optimizeNetwork(observations, predGrid, candidates=mask, method="ssa", action="add", nDiff=100, model=vt.fit, criterion="MUKV", plot=FALSE, nr_iterations=50, nmax=40)
## compute the Mukv of the extended network:
calculateMukv(optim.ok, predGrid, model=vt.fit, nmax=60)

plt1 <- spplot(meuse.grid["mask"], col.regions=grey(rev(seq(0.3,1,0.025))), sp.layout=list(transects.sp, list("sp.points", optim.ok[-c(1:49),], pch=21, col="black", cex=1)), main="Optimized network (SSA)")

## RK model:
## in intamapInteractive still difficult to optimize samples using covariates
library(spatstat)
proj4string(observations) <- meuse.grid@proj4string
zinc.rk2 <- krige(log1p(zinc)~sqrt(dist), observations, meuse.grid, vr.fit)
zinc.rk2$weight <- zinc.rk2$var1.var^2
dens.weight <- as.im(as.image.SpatialGridDataFrame(zinc.rk2["weight"]))
new.rk <- rpoint(100, f=dens.weight)
new.rk <- as.SpatialPoints.ppp(new.rk)
plt2 <- spplot(meuse.grid["dist"], col.regions=grey(rev(seq(0.3,1,0.025))), sp.layout=list(transects.sp, list("sp.points", new.rk, pch=21, col="black", cex=1)), main="New sampling locations (RK)")

print(plt1, split=c(1,1,2,1), more=T)
print(plt2, split=c(2,1,2,1), more=F)

## to be continued!

## end of script