## PROOF1: Every universal kriging in gstat is basically OLS-based

library(gstat)
demo(meuse, echo=FALSE)
rgm <- lm(log1p(om)~dist, meuse)
meuse <- meuse[-rgm$na.action,]
meuse$res <- resid(rgm)
m <- vgm(0, "Exp", 166.0288, 0.1329353)
uk <- krige(log1p(om)~dist, meuse, m, newdata=meuse.grid[2000,])$var1.pred
rk <- predict(rgm, meuse.grid[2000,]) + krige(res~1, meuse, m, newdata=meuse.grid[2000,])$var1.pred
uk
rk
## exactly the same
## The literature suggests that the predictions should be based on Generalized Least Squares estimation!

## PROOF2: GLM-kriging is often a better choice than TransGaussian kriging
## Comparison GLM vs log-regression for predicting ORC

library(GSIF)
library(sp)
demo(meuse, echo=FALSE)
## fit a regression-tree:
omm1 <- fit.gstatModel(meuse, log1p(om)~dist+ffreq, meuse.grid)
omm2 <- fit.gstatModel(meuse, om~dist+ffreq, meuse.grid, fit.family = gaussian(log))
## compare the two models:
cv.1 <- validate(omm1)
sqrt(mean((expm1(cv.1$validation$var1.pred+cv.1$validation$var1.var/2)-expm1(cv.1$validation$observed))^2))
cv.2 <- validate(omm2)
sqrt(mean((cv.2$validation$var1.pred-cv.2$validation$observed)^2))

## plot regression lines next to each other:
par(mfrow=c(1,2))
plot(y=expm1(cv.1$validation$var1.pred+cv.1$validation$var1.var/2), x=expm1(cv.1$validation$observed), asp=1, main="trans-Gaussian (log-transformed)", ylab="Predicted", xlab="Observed", xlim=c(0,18), ylim=c(0,18))
abline(a=0, b=1, lw=2, col="black")
plot(y=cv.2$validation$var1.pred, x=cv.2$validation$observed, asp=1, main="GLM-kriging (log-link)", ylab="Predicted", xlab="Observed", xlim=c(0,18), ylim=c(0,18))
abline(a=0, b=1, lw=2, col="black")
## Summary: GLM-kriging is slighly better in predicting higher values