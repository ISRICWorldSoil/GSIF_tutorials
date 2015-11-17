# title         : RK_vs_RandomForestK.R
# purpose       : Comparison of regression-kriging and RandomForest-Kriging using Meuse data set
# reference     : Hengl, T. 2009. A Practical Guide to Geostatistical Mapping, 2nd Edt. University of Amsterdam, www.lulu.com, 291 p. ISBN: ISBN 978-90-9024981-0
# producer      : T. Hengl
# address       : In Wageningen
# inputs        : meuse dataset (comes with gstat);
# outputs       : predictions of organic matter %, assessment of the mapping accuracy;
# remarks 1     : random forest package [http://cran.r-project.org/web/packages/randomForest/index.html]

library(sp)
library(GSIF)
library(randomForest)
library(plotKML)
library(gstat)
library(raster)
data(SAGA_pal)

## load the Meuse data set:
demo(meuse, echo=FALSE)
## regression-kriging:
omm <- fit.gstatModel(meuse, log1p(om)~dist+ffreq, meuse.grid)
om.rk <- predict(omm, meuse.grid)
om.rk
v = omm@vgmModel
class(v) = c("variogramModel", "data.frame")
## compare with universal kriging:
om.uk <- krige(log1p(om)~dist+ffreq, meuse[-omm@regModel$na.action,], meuse.grid, model=v)
str(om.uk)
## TH: conclusion -> there is a systematic difference (RK simple is more optimistic than UK var), but the RK var simple vs UK var are on a straight line;

## fit randomForest model:
omm2 <- fit.gstatModel(meuse, log1p(om)~dist+ffreq, meuse.grid, method="randomForest")
om.rfk <- predict(omm2, meuse.grid)

## compare predictions:
meuse.grid$om.rk <- om.rk@predicted$om
meuse.grid$om.rfk <- om.rfk@predicted$om
spplot(meuse.grid[c("om.rk", "om.rfk")], col.regions=SAGA_pal[[1]], names.attr=c("regression-kriging", "randomForest-kriging"), sp.layout=list("sp.points", meuse, pch="+", cex=1.5, col="black"))
                            
## overlay points and fit a GLM:
meuse.ov <- over(meuse, meuse.grid)
meuse.ov <- cbind(meuse.ov, meuse@data)
## compare with GLM only:
omm0 <- glm(log1p(om)~dist+ffreq, meuse.ov, family=gaussian())
om.glm <- predict.glm(omm0, meuse.grid, se.fit=TRUE)
str(om.glm)
meuse.grid$om.glmvar <- om.glm[["se.fit"]]^2+om.glm[["residual.scale"]]^2
summary(meuse.grid$om.glmvar)
summary(om.uk$var1.var)
## TH: GLM var is now much higher becase the residual scale is high:
om.glm[["residual.scale"]]^2
par(mfrow=c(1,2))
plot(om.uk$var1.var, om.glm$se.fit^2+om.glm[["residual.scale"]]^2, xlab="Universal kriging variance", ylab="GLM variance", asp=1, xlim=c(0,.15), ylim=c(0,.15))
abline(a=0, b=1)
## the difference:
plot(om.uk$var1.var, om.rk@predicted$var1.var, xlab="Universal kriging variance", ylab="Regression-kriging (simple) variance", asp=1, xlim=c(0,.15), ylim=c(0,.15))
abline(a=0, b=1)

## fit a random forest model (step-by-step)
meuse.ov <- meuse.ov[-omm0$na.action,]
omm3 <- randomForest(log1p(om)~dist+ffreq, meuse.ov)
om.rf <- predict(omm3, meuse.grid)
str(om.rf)
## TH: This does not say anything about the prediction uncertainty

## map the error for random forest suggestion #1 by Forrest Stevens forrest@ufl.edu:
om.rf.all <- predict(omm3, meuse.grid, predict.all=TRUE)
meuse.grid$om.rfvar <- apply(om.rf.all$individual, MARGIN=1, var)
## plot the output:
spplot(meuse.grid[c("om.glmvar", "om.rfvar")], col.regions=SAGA_pal[[12]], names.attr=c("regression variance", "randomForest variance"), sp.layout=list("sp.points", meuse, pch="+", cex=1.5, col="black"))
## TH: I think this only gives us idea of the diversity of the trees available.

## map the error for random forest suggestion #2 by Forrest Stevens forrest@ufl.edu:
library(quantregForest)
ommX <- quantregForest(y=log1p(meuse.ov$om), x=meuse.ov[,c("dist", "ffreq")])
## let us assume normal distribution (-1 / 1 s.d.), which is probably incorrect:
om.rf.quant <- predict(ommX, meuse.grid@data[,c("dist", "ffreq")], quantiles=c((1-.682)/2, 1-(1-.682)/2))
## TH: this takes few seconds...
str(om.rf.quant)
meuse.grid$om.rfvar2 <- ((om.rf.quant[,1] - om.rf.quant[,2])/2)^2
spplot(meuse.grid[c("om.glmvar", "om.rfvar2")], col.regions=SAGA_pal[[12]], names.attr=c("regression variance", "randomForest variance"), sp.layout=list("sp.points", meuse, pch="+", cex=1.5, col="black"))
## TH: this looks conviencing - GLM makes standard error while for random forest the error is more diverse

## randomForest-kriging is also directly available via the GSIF package:
omm4 <- fit.gstatModel(meuse, log1p(om)~dist+ffreq, meuse.grid, method="quantregForest")
om.rfk1 <- predict(omm4, meuse.grid)
omm5 <- fit.gstatModel(meuse, log1p(om)~dist+ffreq, meuse.grid, method="quantregForest")
om.rfk2 <- predict(omm5, meuse.grid)
## TH: this takes again few seconds...
meuse.grid$om.rfkvar <- om.rfk1@predicted$var1.var 
meuse.grid$om.ukvar <- om.uk$var1.var
spplot(meuse.grid[c("om.ukvar", "om.rfkvar")], col.regions=SAGA_pal[[12]], names.attr=c("universal kriging variance", "randomForest-kriging variance"), sp.layout=list("sp.points", meuse, pch="+", cex=1.5, col="black"))

## plot together everything:
pre.r <- range(c(om.uk$var1.pred, om.rfk1@predicted$var1.pred, om.rfk1@predicted$var1.pred), na.rm=TRUE)
var.r <- range(c(om.uk$var1.var, om.rfk@predicted$var1.var, om.rfk@predicted$var1.var), na.rm=TRUE)
tvar.rk <- 1-var(om.rk@validation$residual, na.rm=T)/var(om.rk@validation$observed, na.rm=T)
tvar.rfk1 <- 1-var(om.rfk1@validation$residual, na.rm=T)/var(om.rfk1@validation$observed, na.rm=T)
tvar.rfk2 <- 1-var(om.rfk2@validation$residual, na.rm=T)/var(om.rfk2@validation$observed, na.rm=T)
rxV <- rev(as.character(round(seq(pre.r[1], pre.r[2], length.out=5), 2)))
rvV <- rev(as.character(round(seq(var.r[1], var.r[2], length.out=5), 4)))

par(mfrow=c(2,3), mar=c(0,0,2.5,0), oma=c(0,0,0,0))
image(raster(om.uk["var1.pred"]), main="regression-kriging", asp=1, axes=FALSE, xlab="", ylab="", col=SAGA_pal[[1]], zlim=pre.r)
points(meuse, pch="+", cex=.8)
legend("topleft", rxV, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n", y.intersp=1.5)
image(raster(om.rfk1@predicted["var1.pred"]), main="randomForest-kriging #1", asp=1, axes=FALSE, xlab="", ylab="", col=SAGA_pal[[1]], zlim=pre.r)
points(meuse, pch="+", cex=.8)
image(raster(om.rfk2@predicted["var1.pred"]), main="randomForest-kriging #2", asp=1, axes=FALSE, xlab="", ylab="", col=SAGA_pal[[1]], zlim=pre.r)
points(meuse, pch="+", cex=.8)
image(raster(om.uk["var1.var"]), main=paste("variance (explained: ", signif(tvar.rk*100, 2),"%)", sep=""), asp=1, axes=FALSE, xlab="", ylab="", col=SAGA_pal[[12]], zlim=var.r)
points(meuse, pch="+", cex=.8)
legend("topleft", rvV, fill=rev(SAGA_pal[[12]][c(1,5,10,15,20)]), horiz=FALSE, bty="n", y.intersp=1.5)
image(raster(om.rfk1@predicted["var1.var"]), main=paste("variance (explained: ", signif(tvar.rfk1*100, 2),"%)", sep=""), asp=1, axes=FALSE, xlab="", ylab="", col=SAGA_pal[[12]], zlim=var.r)
points(meuse, pch="+", cex=.8)
image(raster(om.rfk2@predicted["var1.var"]), main=paste("variance (explained: ", signif(tvar.rfk1*100, 2),"%)", sep=""), asp=1, axes=FALSE, xlab="", ylab="", col=SAGA_pal[[12]], zlim=var.r)
points(meuse, pch="+", cex=.8)

## end of script;