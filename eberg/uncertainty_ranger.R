## Some examples with deriving uncertainty (for spatial prediction) using the ranger package
## tom.hengl@isric.org

list.of.packages = c("GSIF", "plotKML", "entropy", "plyr", "parallel", "ranger", "doParallel", "doMC", "rgdal", "caret")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(plotKML)
library(GSIF)
library(entropy)
library(ranger)
library(rgdal)
library(caret)
library(parallel)
library(doParallel)
library(doMC)
library(plyr)
cpus = 4
registerDoMC(cpus)

## Uncertainty factor type variable:
data(eberg)
data(eberg_grid)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
eberg_spc <- spc(eberg_grid, ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6)
eberg_grid@data <- cbind(eberg_grid@data, eberg_spc@predicted@data)
## overlay and create a regression-matrix:
ov <- over(eberg, eberg_grid)
m <- cbind(ov, eberg@data)
## clean-up target variable:
xg = summary(m$TAXGRSC, maxsum=length(levels(m$TAXGRSC)))
str(xg)
selg.levs = attr(xg, "names")[xg > 5]
m$soiltype <- m$TAXGRSC
m$soiltype[which(!m$TAXGRSC %in% selg.levs)] <- NA
m$soiltype <- droplevels(m$soiltype)
m$soiltype <- as.factor(m$soiltype)
## fit model:
fm1 = as.formula(paste("soiltype ~", paste0("PC", 1:10, collapse = "+")))
TAXGRSC.rf <- ranger(fm1, data = m[complete.cases(m[,all.vars(fm1)]),], write.forest = TRUE, probability = TRUE)
## predict values:
eberg_soiltype = eberg_grid[1]
eberg_soiltype@data <- data.frame(predict(TAXGRSC.rf, eberg_grid@data, na.action = na.pass)$predictions)
str(eberg_soiltype@data)
## uncertainty can be represented using the Shannon Entropy:
entropy.empirical(c(0,0,0,0,1))
entropy.empirical(c(1/5,1/5,1/5,1/5,1/5))
entropy.empirical(unlist(eberg_soiltype@data[1,]))
## apply function in parallel:
eberg_soiltype$SE = unlist(alply(eberg_soiltype@data, 1, .fun=function(x){entropy.empirical(unlist(x))}, .parallel = TRUE))
hist(eberg_soiltype$SE)
## divide by maximum possible uncertainty:
eberg_soiltype$SE = round(eberg_soiltype$SE/entropy.empirical(rep(1/length(levels(m$soiltype)),length(levels(m$soiltype))))*100)
## plot:
plot(raster(eberg_soiltype["SE"]), col=SAGA_pal[[1]])
points(eberg[!is.na(eberg$TAXGRSC),], pch="+")

## Uncertainty numeric variable:
fm2 = as.formula(paste("SNDMHT_D ~", paste0("PC", 1:10, collapse = "+")))
mS = m[complete.cases(m[,all.vars(fm2)]),all.vars(fm2)]
SNDMHT.rf <- ranger(fm2, data = mS, write.forest = TRUE)
SNDMHT.rf
## Sample CV errors based on the caret approach:
fitControl <- trainControl(method="repeatedcv", number=5, repeats=1)
mFit <- train(fm2, data=mS, method="ranger", trControl=fitControl)
mS$resid = abs(resid(mFit))
mS$predicted = mFit$finalModel$predictions
xyplot(predicted~SNDMHT_D, mS, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), xlab="measured", ylab="predicted (machine learning)")
## scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE))
## model absolute residuals as function of covariates:
var.fm2 = as.formula(paste("resid ~", paste0("PC", 1:10, collapse = "+")))
var.SNDMHT.rf <- ranger(var.fm2, data = mS, write.forest = TRUE, mtry = mFit$finalModel$mtry)
var.SNDMHT.rf
## predict uncertainty:
eberg_SNDMHT = eberg_grid[1]
eberg_SNDMHT$SNDMHT_D <- predict(mFit, eberg_grid@data, na.action = na.pass)
plot(raster(eberg_SNDMHT["SNDMHT_D"]), col=SAGA_pal[[1]])
eberg_SNDMHT$var.SNDMHT_D <- predict(var.SNDMHT.rf, eberg_grid@data, na.action = na.pass)$predictions
summary(eberg_SNDMHT$var.SNDMHT_D)
eberg_SNDMHT$var.SNDMHT_Df = ifelse(eberg_SNDMHT$var.SNDMHT_D>12, 12, eberg_SNDMHT$var.SNDMHT_D)
## uncertainty is function of location / covariates, but also of values (higher values - higher uncertainty usually):
plot(raster(eberg_SNDMHT["var.SNDMHT_Df"]), col=SAGA_pal[[1]])
points(eberg[!is.na(eberg$SNDMHT_D),], pch="+")
