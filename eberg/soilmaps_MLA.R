## Testing machine learning methods for prediction of soil classes and soil properties
## tom.hengl@isric.org

list.of.packages <- c("nnet", "plyr", "ROCR", "randomForest", "plyr", "parallel", "psych", "mda", "Cubist", "h2o", "dismo", "grDevices", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "mxnet", "doParallel", "caret", "bartMachine")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#library(devtools)
#devtools::install_github('dmlc/mxnet/R-package')
#library(caretEnsemble)
library(plotKML)
library(sp)
library(randomForest)
library(nnet)
library(e1071)
library(GSIF)
library(plyr)
library(raster)
library(caret)
library(Cubist)
library(GSIF)
library(xgboost)
library(snowfall)
library(mxnet)
#options(java.parameters = "-Xmx5g")
library(bartMachine)
set_bart_machine_num_cores(4)

set.seed(42)
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
dim(m)

# ------------------------------------------------------------
# SOIL TYPES / FACTORS
# ------------------------------------------------------------

## clean-up target variable:
xg = summary(m$TAXGRSC, maxsum=length(levels(m$TAXGRSC)))
str(xg)
selg.levs = attr(xg, "names")[xg > 5]
m$soiltype <- m$TAXGRSC
m$soiltype[which(!m$TAXGRSC %in% selg.levs)] <- NA
m$soiltype <- droplevels(m$soiltype)
str(summary(m$soiltype, maxsum=length(levels(m$soiltype))))

## use only complete cases:
m <- m[complete.cases(m[,1:(ncol(eberg_grid)+2)]),]
m$soiltype <- as.factor(m$soiltype)

## subsample
s <- sample.int(nrow(m), 500)
## CV error for RF
TAXGRSC.rf <- randomForest(x=m[-s,paste0("PC",1:10)], y=m$soiltype[-s], xtest=m[s,paste0("PC",1:10)], ytest=m$soiltype[s])
TAXGRSC.rf$test$confusion[,"class.error"]

## Fit 4 independent machine learning models
## Random Forest
TAXGRSC.rf <- randomForest(x=m[,paste0("PC",1:10)], y=m$soiltype, keep.forest = TRUE)
## prediction success per class:
TAXGRSC.rf
## Multinomial Logistic Regression
fm = as.formula(paste("soiltype~", paste(paste0("PC",1:10), collapse="+")))
TAXGRSC.mn <- multinom(fm, m)
## Support Vector Machines
TAXGRSC.svm <- svm(fm, m, probability=TRUE, cross=5) ## Takes few secs
TAXGRSC.svm$tot.accuracy
## Deep Neural Net
train.x = data.matrix(m[,all.vars(fm)[-1]])
## mxnet requires that the classes are coded as 0, 1, 2 ... p 
train.yf = as.factor(as.vector(m[,all.vars(fm)[1]]))
train.y = as.integer(train.yf)-1
TAXGRSC.nn <- mx.mlp(train.x, train.y, hidden_node=ncol(train.x), out_node=length(levels(m$soiltype)), out_activation="softmax", num.round=10, array.batch.size=15, learning.rate=0.1, momentum=0.9, eval.metric=mx.metric.accuracy, array.layout="rowmajor")

## Make ensemble predictions:
probs1 <- predict(TAXGRSC.mn, eberg_grid@data, type="probs", na.action = na.pass)
probs2 <- predict(TAXGRSC.rf, eberg_grid@data, type="prob", na.action = na.pass)
probs3 <- attr(predict(TAXGRSC.svm, eberg_grid@data, probability=TRUE, na.action = na.pass), "probabilities")
#probs3 <- predict(ens, eberg_grid@data, probability=TRUE, na.action = na.pass)
probs4 <- t(predict(TAXGRSC.nn, data.matrix(eberg_grid@data[,all.vars(fm)[-1]]), array.layout="rowmajor"))
attr(probs4, "dimnames")[[2]] = levels(train.yf)

leg <- levels(m$soiltype)
lt <- list(probs1[,leg], probs2[,leg], probs3[,leg], probs4[,leg])
## Simple average (an alternative would be to use the "caretEnsemble" package):
probs <- Reduce("+", lt) / length(lt)
eberg_soiltype = eberg_grid
#eberg_soiltype@data <- data.frame(probs4)
eberg_soiltype@data <- data.frame(probs)
## check that probs sum up to 1:
ch <- rowSums(eberg_soiltype@data)
summary(ch)
plot(raster::stack(eberg_soiltype), col=SAGA_pal[[1]], zlim=c(0,1))
#plotKML(eberg["TAXGRSC"])
setwd("./img")
kml(eberg_soiltype["Rendzina"], colour=Rendzina, colour_scale=SAGA_pal[[1]], z.lim=c(0,1))
kml(eberg_soiltype["Braunerde"], colour=Braunerde, colour_scale=SAGA_pal[[1]], z.lim=c(0,1))
#plotKML(eberg_soiltype["Regosol"], colour_scale=SAGA_pal[[1]], z.lim=c(0,1))

# ------------------------------------------------------------
# SOIL PROPERTIES / NUMERIC (h2o)
# ------------------------------------------------------------

names(eberg)
## randomForest for predicting Sand content (top soil):
library(h2o)
localH2O = h2o.init(nthreads = -1)
eberg.hex = as.h2o(m, destination_frame = "eberg.hex")
eberg.grid = as.h2o(eberg_grid@data, destination_frame = "eberg.grid")
names(m)[22] 
RF.m = h2o.randomForest(y = which(names(m)=="SNDMHT_A"), x = which(names(m) %in% paste0("PC",1:10)), training_frame = eberg.hex, ntree = 50)
RF.m
plot(m$SNDMHT_A, as.data.frame(h2o.predict(RF.m, eberg.hex, na.action=na.pass))$predict, asp=1, xlab="measured", ylab="predicted")
rf.VI = RF.m@model$variable_importances
print(rf.VI)
eberg_grid$RFx <- as.data.frame(h2o.predict(RF.m, eberg.grid, na.action=na.pass))$predict
plot(raster(eberg_grid["RFx"]), col=SAGA_pal[[1]], zlim=c(10,90))
points(eberg, pch="+")
rf.R2 = RF.m@model$training_metrics@metrics$r2
rf.R2

## Deep learning
DL.m = h2o.deeplearning(y = which(names(m)=="SNDMHT_A"), x = which(names(m) %in% paste0("PC",1:10)), training_frame = eberg.hex)
DL.m
eberg_grid$DLx <- as.data.frame(h2o.predict(DL.m, eberg.grid, na.action=na.pass))$predict
plot(raster(eberg_grid["DLx"]), col=SAGA_pal[[1]], zlim=c(10,90))
points(eberg, pch="+")
dl.R2 = DL.m@model$training_metrics@metrics$r2

## merge (weighted average):
eberg_grid$SNDMHT_A <- rowSums(cbind(eberg_grid$RFx*rf.R2, eberg_grid$DLx*dl.R2), na.rm=TRUE)/(rf.R2+dl.R2)
plot(raster(eberg_grid["SNDMHT_A"]), col=SAGA_pal[[1]], zlim=c(10,90))
points(eberg, pch="+")
plotKML(eberg_grid["SNDMHT_A"], colour_scale=SAGA_pal[[1]], z.lim=c(0,100))
shape = "http://maps.google.com/mapfiles/kml/paddle/wht-blank.png"
kml(eberg, colour=SNDMHT_A, shape=shape, labels=SNDMHT_A, colour_scale=SAGA_pal[[10]])
kml_View("eberg.kml")
## remove objects:
#h2o.rm(eberg_grid)
h2o.shutdown()

# ------------------------------------------------------------
# SOIL PROPERTIES / NUMERIC (3D)
# ------------------------------------------------------------

## Soil organic carbon in 3D
## Model fitting using CARET package [http://topepo.github.io/caret/]
## under Unix we would also run:
#library(doMC)
#registerDoMC(cores = 8)

data(edgeroi)
edgeroi.spc <- join(edgeroi$sites, edgeroi$horizons, type='inner')
coordinates(edgeroi.spc) <- ~ LONGDA94 + LATGDA94
proj4string(edgeroi.spc) <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
edgeroi.spc <- spTransform(edgeroi.spc, CRS("+init=epsg:28355"))
## load the 250 m grids:
con <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids.rda")
load(con)
gridded(edgeroi.grids) <- ~x+y
proj4string(edgeroi.grids) <- CRS("+init=epsg:28355")
## overlay points and grids:
ov2 <- over(edgeroi.spc, edgeroi.grids)
m2 <- cbind(ov2, edgeroi.spc@data)
m2$DEPTH <- m2$UHDICM + (m2$LHDICM - m2$UHDICM)/2
summary(m2$DEPTH)
formulaStringP2 = ORCDRC ~ DEMSRT5+TWISRT5+PMTGEO5+EV1MOD5+EV2MOD5+EV3MOD5+DEPTH
mP2 <- m2[complete.cases(m2[,all.vars(formulaStringP2)]),]

## process parallelized by default:
ctrl <- trainControl(method="repeatedcv", number=5, repeats=1)
## First, fine-tune parameters using smaller number of samples 
sel <- sample.int(nrow(mP2), 500)
tr.ORCDRC.rf <- train(formulaStringP2, data=mP2[sel,], method = "rf", trControl = ctrl, tuneLength = 3)
tr.ORCDRC.rf
## fit the final model using fine-tuned parameters
## http://stats.stackexchange.com/questions/23763/is-there-a-way-to-disable-the-parameter-tuning-grid-feature-in-caret
ORCDRC.rf <- train(formulaStringP2, data=mP2, method = "rf", tuneGrid=data.frame(mtry=7), trControl=trainControl(method="none"))
w1 = 100*max(tr.ORCDRC.rf$results$Rsquared)
varImpPlot(ORCDRC.rf$finalModel)
## the same using the "Cubist" package:
tr.ORCDRC.cb <- train(formulaStringP2, data=mP2[sel,], method = "cubist", trControl = ctrl, tuneLength = 3)
tr.ORCDRC.cb
ORCDRC.cb <- train(formulaStringP2, data=mP2, method = "cubist", tuneGrid=data.frame(committees = 1, neighbors = 0), trControl=trainControl(method="none"))
w2 = 100*max(tr.ORCDRC.cb$results$Rsquared)
## "XGBoost" package:
ORCDRC.gb <- train(formulaStringP2, data=mP2, method = "xgbTree", trControl=ctrl)
w3 = 100*max(ORCDRC.gb$results$Rsquared)
## "bartMachine" package (Bayesian Additive Regression Trees):
tuneGrid.bM <- expand.grid(num_trees=c(20,80), k=2, alpha=0.95, beta=2, nu=3)
tr.ORCDRC.bM <- train(formulaStringP2, data=mP2[sel,], method = "bartMachine", trControl=ctrl, tuneGrid=tuneGrid.bM)
w4 = 100*max(tr.ORCDRC.bM$results$Rsquared)
ORCDRC.bM <- train(formulaStringP2, data=mP2, method = "bartMachine", tuneGrid=tr.ORCDRC.bM$finalModel$tuneValue, trControl=trainControl(method="none"))

## Ensemble prediction:
edgeroi.grids$DEPTH = 5
edgeroi.grids$Random_forest <- predict(ORCDRC.rf, edgeroi.grids@data, na.action = na.pass) 
edgeroi.grids$Cubist <- predict(ORCDRC.cb, edgeroi.grids@data, na.action = na.pass)
edgeroi.grids$XGBoost <- predict(ORCDRC.gb, edgeroi.grids@data, na.action = na.pass)
## Takes few minutes...
## See computing times in: https://arxiv.org/abs/1312.2171
edgeroi.grids$bartMachine <- predict(ORCDRC.bM, edgeroi.grids@data, na.action = na.pass)
## final prediction as weighted average:
edgeroi.grids$ORCDRC_5cm <- (edgeroi.grids$Random_forest*w1+edgeroi.grids$Cubist*w2+edgeroi.grids$XGBoost*w3+edgeroi.grids$bartMachine*w4)/(w1+w2+w3+w4)
plot(log1p(stack(edgeroi.grids[c("Random_forest","Cubist","XGBoost","bartMachine","ORCDRC_5cm")])), col=SAGA_pal[[1]], zlim=log1p(c(5,65)))
## Visual comparison:
plot(log1p(raster(edgeroi.grids["bartMachine"])), col=SAGA_pal[[1]], zlim=log1p(c(5,65)))
points(edgeroi.spc, pch="+", col="black")
plot(log1p(raster(edgeroi.grids["Random_forest"])), col=SAGA_pal[[1]], zlim=log1p(c(5,65)))
points(edgeroi.spc, pch="+", col="black")
#plotKML(edgeroi.grids["ORCDRC_5cm"], colour_scale=SAGA_pal[[1]], z.lim=c(5,45))

## aggregate observed values for 0-5 cm depths:
m2$cl <- cut(m2$DEPTH, breaks=c(0,5,15,30,60,100,200))
ORCDRC.mv <- ddply(m2, .(cl,SOURCEID), summarize, aggregated = mean(ORCDRC))
ORCDRC.mv <- join(ORCDRC.mv, as.data.frame(edgeroi.spc)[,c("SOURCEID","LONGDA94","LATGDA94")], type="left", match = "first")
ORCDRC.mv.s <- ORCDRC.mv[ORCDRC.mv$cl=="(0,5]",]
ORCDRC.mv.s <- ORCDRC.mv.s[!is.na(ORCDRC.mv.s$LONGDA94),]
coordinates(ORCDRC.mv.s) <- ~ LONGDA94 + LATGDA94
proj4string(ORCDRC.mv.s) <- CRS("+init=epsg:28355") 
#plotKML(ORCDRC.mv.s["aggregated"])

## Test accuracy of mapping ORCDRC:
source_https <- function(url, ...) {
  require(RCurl)
  cat(getURL(url, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), file = basename(url))
  source(basename(url))
}
source_https("https://raw.githubusercontent.com/ISRICWorldSoil/SoilGrids250m/master/grids/cv/cv_functions.R")

test.ORC <- cv_numeric(formulaStringP2, rmatrix=mP2, nfold=5, idcol="SOURCEID", Log=TRUE)
str(test.ORC)
## Plot CV results (use log-scale):
plt0 <- xyplot(test.ORC[[1]]$Observed~test.ORC[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="measured", xlab="predicted (machine learning)")
plt0
