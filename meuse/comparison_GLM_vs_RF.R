## Comparison GLM vs RandomForest
library(gstat)
library(randomForest)
library(quantregForest)
library(scales)
library(psych)
library(intamap)

demo(meuse, echo=FALSE)
meuse$value = log(meuse$zinc)
output = interpolate(meuse, meuse.grid)
plot(output)

## compare GLM vs RF
m <- glm(zinc~log1p(dist)+ffreq, meuse, family=gaussian(link=log))
plot(m$fitted.values~m$y, asp=1)
abline(0,1)
rf <- quantregForest(x=meuse@data[,c("dist","ffreq")], y=meuse$zinc)
plot(rf$predicted~rf$y, asp=1)
abline(0,1)

meuse.grid$glm.zinc <- predict(m, meuse.grid@data, type="response")
meuse.grid$rf.zinc <- predict(rf, meuse.grid@data)[,2]
## plot next to each other:
spplot(meuse.grid[c("glm.zinc","rf.zinc")], col.regions=SAGA_pal[[1]])
with(meuse.grid@data, psych::scatter.hist(glm.zinc, rf.zinc, xlab="GLM", ylab="randomForest", pch=19, title="", col=scales::alpha("lightblue", 0.4), cex=1.5))
