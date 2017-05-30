## Some examples with deriving uncertainty (for spatial prediction) using the ranger package
## tom.hengl@isric.org

list.of.packages = c("GSIF", "plotKML", "entropy", "plyr", "parallel", "ranger", "raster", "doParallel", "doMC", "rgdal", "caret", "dismo", "snowfall")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(plotKML)
library(GSIF)
library(entropy)
library(ranger)
library(rgdal)
library(dismo)
library(raster)
library(caret)
library(parallel)
library(doParallel)
library(doMC)
library(plyr)
library(snowfall)
cpus = 8
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
mD = m[complete.cases(m[,all.vars(fm1)]),all.vars(fm1)]
TAXGRSC.rf <- ranger(fm1, data = mD, write.forest = TRUE, probability = TRUE)
## predict values:
eberg_soiltype = eberg_grid[1]
eberg_soiltype@data <- data.frame(predict(TAXGRSC.rf, eberg_grid@data, na.action = na.pass)$predictions)
str(eberg_soiltype@data)
## uncertainty can be represented using the Shannon Entropy (entropy package); some examples:
entropy.empirical(c(0,0,0,0,1)) ## smallest entropy
entropy.empirical(c(1/5,1/5,1/5,1/5,1/5)) ## highest entropy
## scaled entropy using log-base as number of classes:
p=c(0,0.8,0.2,0,0)
sum(-p*log(p,base=length(p)), na.rm=TRUE)
entropy.empirical(p)/entropy.empirical(rep(1/length(p),length(p)))
## same result!
entropy.empirical(unlist(eberg_soiltype@data[1,])) ## example

## apply entropy function in parallel:
eberg_soiltype$SE = unlist(alply(eberg_soiltype@data, 1, .fun=function(x){entropy.empirical(unlist(x))}, .parallel = TRUE))
hist(eberg_soiltype$SE, col="grey")
## divide by maximum possible uncertainty:
eberg_soiltype$SE = round(eberg_soiltype$SE/entropy.empirical(rep(1/length(levels(m$soiltype)),length(levels(m$soiltype))))*100)
## plot:
plot(raster(eberg_soiltype["SE"]), col=SAGA_pal[[1]])
points(eberg[!is.na(eberg$TAXGRSC),], pch="+")
## high uncertainty correlates with extrapolation areas

## Uncertainty for numeric variables:
fm2 = as.formula(paste("SNDMHT_D ~", paste0("PC", 1:10, collapse = "+")))
mS = m[complete.cases(m[,all.vars(fm2)]),all.vars(fm2)]
SNDMHT.rf <- ranger(fm2, data = mS, write.forest = TRUE)
SNDMHT.rf
sqrt(SNDMHT.rf$prediction.error) ## mean error = +/-15

## Estimate CV residuals with re-fitting ("repeated CV"):
predict_cv_resid = function(formulaString, data, nfold){
  sel <- dismo::kfold(data, k=nfold)
  out = list(NULL)
  for(j in 1:nfold){
    s.train <- data[!sel==j,]
    s.test <- data[sel==j,]
    m <- ranger(formulaString, data=s.train, write.forest=TRUE)
    pred <- predict(m, s.test, na.action = na.pass)$predictions
    obs.pred <- as.data.frame(list(s.test[,all.vars(formulaString)[1]], pred))
    names(obs.pred) = c("Observed", "Predicted")
    obs.pred[,"ID"] <- row.names(s.test)
    obs.pred$fold = j
    out[[j]] = obs.pred
  }
  out <- plyr::rbind.fill(out)
  out <- out[order(as.numeric(out$ID)),]
  return(out)
}
resid.SNDMHT = predict_cv_resid(fm2, data=mS, nfold=5)
xyplot(Predicted~Observed, resid.SNDMHT, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), xlab="measured", ylab="predicted (ranger)")
mS$resid = abs(resid.SNDMHT$Observed - resid.SNDMHT$Predicted)

## model absolute residuals as function of covariates:
var.fm2 = as.formula(paste("resid ~", paste0("PC", 1:10, collapse = "+")))
var.SNDMHT.rf <- ranger(var.fm2, data = mS, write.forest = TRUE)
var.SNDMHT.rf
## predict uncertainty:
eberg_SNDMHT = eberg_grid[1]
eberg_SNDMHT$SNDMHT_D <- predict(SNDMHT.rf, eberg_grid@data, na.action = na.pass)
plot(raster(eberg_SNDMHT["SNDMHT_D"]), col=SAGA_pal[[1]])
## prediction error:
eberg_SNDMHT$var.SNDMHT_D <- predict(var.SNDMHT.rf, eberg_grid@data, na.action = na.pass)$predictions
summary(eberg_SNDMHT$var.SNDMHT_D)
## mean/median is a bit smaller than what we get with ranger OOB - probably because in most of the study area SND content is low
eberg_SNDMHT$var.SNDMHT_Df = ifelse(eberg_SNDMHT$var.SNDMHT_D>20, 20, eberg_SNDMHT$var.SNDMHT_D)
plot(raster(eberg_SNDMHT["var.SNDMHT_Df"]), col=SAGA_pal[[1]])
points(eberg[!is.na(eberg$SNDMHT_D),], pch="+")
## uncertainty is function of location / covariates, but also of values (higher values - higher uncertainty usually):

## Since ranger v0.7.2 errors can be extracted directly from the package based on the method of Wager et al. (2014) 
## https://github.com/imbs-hl/ranger/issues/136
SNDMHT.rfu <- ranger(fm2, data = mS, num.trees = 85, write.forest = TRUE, keep.inbag = TRUE)
eberg_SNDMHT = eberg_grid[1]
x = predict(SNDMHT.rfu, eberg_grid@data, type = "se") ## this is somewhat computationally intensive if the number of trees is high!
eberg_SNDMHT$SNDMHT_D <- x$predictions
eberg_SNDMHT$sd.SNDMHT_D <- x$se ## standard errors i.e. absolute errors +/-
plot(raster(eberg_SNDMHT["SNDMHT_D"]), col=SAGA_pal[[1]])
plot(raster(eberg_SNDMHT["sd.SNDMHT_D"]), col=SAGA_pal[[1]]) ## in some areas errors are quite high!
points(eberg[!is.na(eberg$SNDMHT_D),], pch="+")
## compare numbers:
summary(x$se)
SNDMHT.rfu
## median se about 15 hence good match!