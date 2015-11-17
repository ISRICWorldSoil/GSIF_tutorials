
library(spatstat)
library(maptools)
library(GSIF)
data(USDA.TT)

load("D:/SPDB/sprofs.rda")
tdf <- sprofs$horizons[,c("CLYPPT", "SLTPPT", "SNDPPT")]
## remove missing values:
tdf <- tdf[!is.na(tdf$SNDPPT)&!is.na(tdf$SLTPPT)&!is.na(tdf$CLYPPT),] 
## subset to 15%:
#tdf <- tdf[runif(nrow(tdf))<.15,]
## the fractions have to sum to 100%:
tdf$Sum = rowSums(tdf)
for(i in c("CLYPPT", "SLTPPT", "SNDPPT")) { tdf[,i] <- tdf[,i]/tdf$Sum * 100 }
tdf <- tdf[!(tdf$Sum>105|tdf$Sum<95),]
names(tdf)[1:3] <- c("CLAY", "SILT", "SAND")
geo <- TT.geo.get(class.sys = "USDA.TT")
## 2D probabilty density on an x-y grid using global data:
tex.img <- TT.kde2d(geo, tri.data = tdf, n=100)
tex.img$x <- tex.img$x/100; tex.img$y <- tex.img$y/100
tex.im <- as.im(tex.img)
image(tex.im, main="2D probabilty density image for texture fractions (USDA)", col=SAGA_pal[[3]][2:20])
lines(as(USDA.TT, "SpatialLines"))
tex.grid <- as.SpatialGridDataFrame.im(tex.im)
tex.grid$TEXMHT <- over(tex.grid, USDA.TT)[,"ID"]
## rename:
USDA.TT.im <- as.data.frame(tex.grid)

USDA.TT.im$v <- signif(USDA.TT.im$v, 3)
USDA.TT.im <- USDA.TT.im[!is.na(USDA.TT.im$TEXMHT)&!is.na(USDA.TT.im$v),]
save(USDA.TT.im, file="USDA.TT.im.rda", compress="xz")
