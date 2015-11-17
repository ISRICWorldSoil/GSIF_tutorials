## Testing machine learning methods for prediction of soil classes and soil properties
## tom.hengl@isric.org

#library(caretEnsemble)
library(plotKML)
library(sp)
library(randomForest)
library(nnet)
library(e1071)
library(GSIF)

set.seed(42)
data(eberg)
data(eberg_grid)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
eberg_spc <- spc(eberg_grid, ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6)
eberg_grid@data <- cbind(eberg_grid@data, eberg_spc@predicted@data)
## overlay:
ov <- over(eberg, eberg_grid)
m <- cbind(ov, eberg@data)

# ------------------------------------------------------------
# SOIL TYPES / FACTORS
# ------------------------------------------------------------

## clean-up target variable:
summary(m$TAXGRSC)
m$soiltype <- NA
for(i in levels(m$TAXGRSC)){
  sel <- grep(pattern=i, m$TAXGRSC)
  if(length(sel)>5){
    m$soiltype[sel] <- i
  }
}
m <- m[complete.cases(m[,1:(ncol(eberg_grid)+2)]),]
m$soiltype <- as.factor(m$soiltype)

## subsample
s <- sample.int(nrow(m), 500)
## CV error for RF
TAXGRSC.rf <- randomForest(x=m[-s,1:4], y=m$soiltype[-s], xtest=m[s,1:4], ytest=m$soiltype[s])
TAXGRSC.rf$test$confusion[,"class.error"]

## Fit 3 independent machine learning models:
TAXGRSC.rf <- randomForest(x=m[,1:4], y=m$soiltype)
## prediction success per class:
TAXGRSC.mn <- multinom(soiltype~PRMGEO6+DEMSRT6+TWISRT6+TIRAST6, m)
TAXGRSC.svm <- svm(soiltype~PRMGEO6+DEMSRT6+TWISRT6+TIRAST6, m, probability = TRUE, cross=5)
TAXGRSC.svm$tot.accuracy
## Make ensemble predictions:
probs1 <- predict(TAXGRSC.mn, eberg_grid@data, type="probs", na.action = na.pass) 
probs2 <- predict(TAXGRSC.rf, eberg_grid@data, type="prob", na.action = na.pass)
probs3 <- predict(TAXGRSC.svm, eberg_grid@data, probability=TRUE, na.action = na.pass)
#probs3 <- predict(ens, eberg_grid@data, probability=TRUE, na.action = na.pass)
lt <- list(probs1, probs2, attr(probs3, "probabilities"))
## Simple average:
probs <- Reduce("+", lt) / length(lt)
eberg_soiltype = eberg_grid
eberg_soiltype@data <- data.frame(probs)
## check that probs sum up to 1:
ch <- rowSums(eberg_soiltype@data)
summary(ch)
plot(raster::stack(eberg_soiltype), col=SAGA_pal[[1]], zlim=c(0,1))
#plotKML(eberg_soiltype["Rendzina"], colour_scale=SAGA_pal[[1]], z.lim=c(0,1))
#plotKML(eberg_soiltype["Braunerde"], colour_scale=SAGA_pal[[1]], z.lim=c(0,1))

# ------------------------------------------------------------
# SOIL PROPERTIES / NUMERIC
# ------------------------------------------------------------

names(eberg)
## randomForest:
library(h2o)
localH2O = h2o.init()
eberg.hex = as.h2o(m, conn = h2o.getConnection(), destination_frame = "eberg.hex")
eberg.grid = as.h2o(eberg_grid@data, conn = h2o.getConnection(), destination_frame = "eberg.grid")
View(eberg.hex@mutable$col_names)
names(m)[22]
RF.m = h2o.randomForest(y = 22, x = 6:16, training_frame = eberg.hex, ntree = 50, depth = 100, importance=T)
eberg_grid$RFx <- as.data.frame(h2o.predict(RF.m, eberg.grid, na.action=na.pass))$predict
rf.VI = RF.m@model$variable_importances
print(rf.VI)
plot(raster(eberg_grid["RFx"]))
#plotKML(eberg_grid["RFx"], colour_scale=SAGA_pal[[10]], z.lim=c(0,100))
shape = "http://maps.google.com/mapfiles/kml/paddle/wht-blank.png"
kml(eberg, colour=SNDMHT_A, shape=shape, labels=SNDMHT_A, colour_scale=SAGA_pal[[10]])
kml_View("eberg.kml")
## remove objects:
#h2o.rm(eberg_grid, conn = h2o.getConnection())

