#########################################
#                                       #
# 3D+T demo of the cookfarm data        #
#                                       #
# The demo script has been prepared by: #
# Benedikt Graeler                      #
# Institute for Geoinformatics          #
# ben.graeler@uni-muenster.de           #
# ifgi.de/graeler                       #
#                                       #
#########################################

data(cookfarm)

# time <- proc.time() # the full script takes 10-15 minutes
## krigeST
library(spacetime)
library(gstat)
library(xts)
library(sp)

# add coordinates to the readings data.frame
matchIDs <- match(cookfarm$readings$SOURCEID, cookfarm$profiles$SOURCEID)
measC <- cookfarm$readings[,c("Date","Port1C","Port2C","Port3C","Port4C","Port5C")]
measC$Easting <-  cookfarm$profiles$Easting[matchIDs]
measC$Northing <-  cookfarm$profiles$Northing[matchIDs]
measC$DoY <- as.numeric(format(as.POSIXct(measC$Date), format="%j"))

# deseasonalize per depth
deseason <- function(data, var, par) {
  data[,var] - (par[1]+par[2]*sin(((par[3]+data[,"DoY"])/365)*2*pi))
}

varVec <- c("Port1C","Port2C","Port3C","Port4C","Port5C")
opt <- matrix(NA,length(varVec), 3)
for (dVar in 1:length(varVec)) {
  opt[dVar,] <- optim(c(0,1,0), 
                      function(par) {
                        sqrt(mean(deseason(measC, varVec[dVar], par)^2, na.rm = T))
                      })$par
  measC[,paste("resid",varVec[dVar],sep="")] <- deseason(measC, varVec[dVar], opt[dVar,])  
}

# re-order data.frame in a long format
measC <- data.frame(Easting=rep(measC$Easting, 5),
                    Northing=rep(measC$Northing, 5),
                    altitude=rep(1:5*-0.3, each=nrow(measC)),
                    Date=rep(measC$Date,5),
                    DoY=rep(measC$DoY, 5),
                    C=c(measC$Port1C, 
                        measC$Port2C,
                        measC$Port3C,
                        measC$Port4C,
                        measC$Port5C),
                    resid=c(measC$residPort1C, 
                            measC$residPort2C,
                            measC$residPort3C,
                            measC$residPort4C,
                            measC$residPort5C))

# abuse time as depth axis and encode depth levels as days
measC$dd <- -measC$altitude/3*10*24*3600 # depth level -> days in secs

# move the 3D data set by 2 km every day avoiding any overlap between days
measC$EastingShift <- measC$Easting - (as.numeric(measC$Date)-14976)*2000 

# create a spacetime data structure
measCSt <- STIDF(SpatialPoints(measC[,c("EastingShift","Northing")]),
                 as.POSIXct(measC$dd, tz = "GMT", origin = "2005-01-01"),
                 measC[,c("resid","C")])
measCSt <- as(measCSt,"STFDF")

# calculate the horizointal-vertical sample variogram
svgmCDepthAsTime <- variogramST(resid~1, measCSt, 
                                tlags=0:4, assumeRegular = T, 
                                width=75, cutoff=375, na.omit=F)
svgmCDepthAsTime$timelag <- svgmCDepthAsTime$timelag*0.3
plot(svgmCDepthAsTime, wireframe=T)

# fit a variogram surface, assuming a metric model
fvgmCDepthAsTime <- fit.StVariogram(svgmCDepthAsTime, 
                                    vgmST("metric",
                                          joint = vgm(0.02, "Sph", 150, 0.02),
                                          stAni=100))
attr(fvgmCDepthAsTime, "optim.output")$value

plot(svgmCDepthAsTime, fvgmCDepthAsTime, wireframe=T, all=T,
     ylab="depth diff [m]", scales=list(arrows=F))


fvgmCDepthAsTime$stAni # 1 m in depth corresponds to ~486 m in the plain

# build a re-scaled spatio-temporal data set with 3D spatial component
measC$scaleAlt <- measC$altitude * fvgmCDepthAsTime$stAni

measCSt <- STIDF(SpatialPoints(measC[,c("Easting","Northing","scaleAlt")]),
                 as.POSIXct(measC$Date, tz = "GMT", origin = "2005-01-01"),
                 measC[,c("resid", "C")])
measCSt <- as(measCSt, "STFDF")

# 3D space-time sample variogram
svgmC3DT <- variogramST(resid~1, measCSt, 
                        tlags=0:9, assumeRegular=T, 
                        width=75, cutoff=375)

plot(svgmC3DT, wireframe=T, all=T, zlim=c(0,1.5),
     xlab = "3D distance [m]", scales=list(arrows=F),
     col.regions=bpy.colors())

# fit a variogram
fvgmC3DT <- fit.StVariogram(svgmC3DT, 
                            vgmST("sumMetric",
                                  space = vgm(1, "Exp", 50, 0),
                                  time = vgm(1, "Exp", 9, 0),
                                  joint = vgm(1, "Sph", 20, 0),
                                  stAni=0.3),
                            lower=c(0,0,0,
                                    0,0,0,
                                    0,0,0,
                                    0),
                            method="L-BFGS-B")
attr(fvgmC3DT, "optim.output")$value # 0.031

plot(svgmC3DT, fvgmC3DT, wireframe=T, all=T, 
     scales=list(arrows=F), xlab = "3D distance [m]")

# predict depth level wise for July:
dLevel <- 1 # 1..5

## gridded data:
grid10m <- cookfarm$grids
gridded(grid10m) <- ~x+y
proj4string(grid10m) <- CRS(cookfarm$proj4string)

predSTF <- as.data.frame(cbind(coordinates(grid10m),
                               c(1:5*(-0.3)*fvgmCDepthAsTime$stAni)[dLevel]))
colnames(predSTF) <- c("Easting", "Northing", "altitude")
coordinates(predSTF) <- ~Easting+Northing+altitude
predSTF <- STF(predSTF, measCSt@time[547:577])

# drop NAs
measCSt <- as(measCSt, "STSDF")
noNA <- !is.na(measCSt@data$resid)
measCSt@index <- measCSt@index[noNA,]
measCSt@data <- measCSt@data[noNA,]

# do the prediction for the first 7 days of July 2012
preds <- NULL
pars <- opt[dLevel,]
for(day in 1:7) { # day <- 1
  cat(paste("day:", day), "\n")
  pred <- krigeST(resid~1, measCSt[,(541+day):(551+day)], predSTF[,day,drop=F],
                  fvgmC3DT)$var1.pred
  doy <- as.numeric(format(index(predSTF[,day,drop=F]@time), format="%j"))
  pred <- data.frame(residPred=pred, 
                     pred=pred + (pars[1]+pars[2]*sin(((pars[3]+doy)/365)*2*pi)))
  preds <- rbind(preds, pred)
}

stCpred <- STFDF(SpatialPixels(SpatialPoints(coordinates(predSTF@sp)[,1:2])),
                 predSTF@time[1:7], preds)

stplot(stCpred[,,"pred"])
# proc.time() - time