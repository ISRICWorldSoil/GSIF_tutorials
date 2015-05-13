## Confidence limits of estimating OM (spatially)
library(gstat)
library(sp)
library(GSIF)
library(plotKML)
library(raster)

## load the Meuse data set:
demo(meuse, echo=FALSE)
sd(meuse$om, na.rm=TRUE)
## regression-kriging:
omm <- fit.gstatModel(meuse, log1p(om)~dist+soil, meuse.grid)
om.rk <- predict(omm, meuse.grid)
## confidence limits:
pt1 <- data.frame(x=179390, y=330820)
coordinates(pt1) <- ~x+y
proj4string(pt1) = proj4string(meuse.grid)
pt1.om <- over(pt1, om.rk@predicted["om"])
pt1.om.sd <- over(pt1, om.rk@predicted["var1.var"])
expm1(pt1.om-1.645*sqrt(pt1.om.sd)); expm1(pt1.om+1.645*sqrt(pt1.om.sd))
(expm1(pt1.om+sqrt(pt1.om.sd)) - expm1(pt1.om-sqrt(pt1.om.sd)))/2

om.rk
## amount of variation explained -> ca 50%
meuse.grid$om.rk <- expm1(om.rk@predicted$om)
## Generate simulations:
om.rksim <- predict(omm, meuse.grid, nsim=5)
ov <- as(om.rksim@realizations, "SpatialGridDataFrame")
meuse.grid$om.sim1 <- expm1(ov@data[,1][meuse.grid@grid.index])

## confidence limits (whole plots):
par(mfrow=c(1,2))
boxplot(om~ffreq, omm@regModel$data, col="grey", xlab="Flooding frequency classes", ylab="Organic matter in %", main="Sampled (N = 153)", ylim=c(0,20))
boxplot(om.sim1~ffreq, meuse.grid, col="grey", xlab="Flooding frequency classes", ylab="Organic matter in %", main="Predicted (spatial simulations)", ylim=c(0,20))
## Fig_confidence_limits.png --> confidence limits are about the same width (grey boxes) which should be the case
## The difference is that the plot on the right shows estimated (based on the gstatModel) confidence limits for the whole area / population.

## standard deviation (per ffreq):
lapply(levels(meuse.grid$ffreq), function(x){sapply(subset(meuse.grid@data, ffreq==x, select=om.sim1), sd, na.rm=TRUE)})
## 3.6% for class "1", 3.0 for class "2" and 1.9% for class "3"

## confidence limits (mean per ffreq):
## error of the mean (http://www.cyclismo.org/tutorial/R/confidence.html)
sd.om <- qt(0.975, df=length(meuse$om)-1)*sd(meuse$om, na.rm=TRUE)/sqrt(length(meuse$om))
sd.om  ## 0.54%
## Estimating the confidence limits has been solved in most of regression models
## So the easiest thing to do is to fit a lm or glm and then evaluate the results:
omm0 <- lm(om~ffreq-1, omm@regModel$data)
om.r <- predict(omm0, meuse.grid, se.fit=TRUE)
meuse.grid$se.fit <- om.r$se.fit
signif(mean(meuse.grid$se.fit, na.rm=TRUE), 3)
## error of the mean = 0.48% --> almost the same number as "omm0" (the difference is that regression estimates the mean error for the whole population!);
aggregate(meuse.grid$se.fit, by=list(meuse.grid$ffreq), mean, na.rm=TRUE)
## while the actual total error is:
aggregate(sqrt(meuse.grid$se.fit^2+om.r$residual.scale^2), by=list(meuse.grid$ffreq), mean, na.rm=TRUE)
## which shows that the variance within units is about 3.3 and it is the same all over.

## For the gstat model, the estimated means are:
aggregate(meuse.grid$om.rk, by=list(meuse.grid$ffreq), mean, na.rm=TRUE)
## the error of estimating individual values is:
aggregate(meuse.grid$om.sim1, by=list(meuse.grid$ffreq), sd, na.rm=TRUE)
## again, the standard errors a bit smaller in 'meuse.grid$om.sim1' than in 'lm(om~ffreq)' because part of variation has been explained by the gstatModel! This is especially clear for ffreq class "3".

## end of script;