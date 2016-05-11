## Some examples of importing, converting and analyzing soil data in R;
## On-line version: [http://gsif.isric.org/doku.php?id=wiki:soil_data]
## by Tom.Hengl@isric.org

# load packages:
library(maptools)
library(rgdal)
library(GSIF)
library(soiltexture)
library(aqp)
library(plyr)
library(plotKML)

#--------------------------------------------------
# Example: converting texture-by-hand classes to fractions
#--------------------------------------------------

## plot the triangle:
TT.plot(class.sys = "USDA.TT")
TT.classes.tbl(class.sys="USDA.TT", collapse=", ")

## convert texture-by-feel classes to fractions:
vert <- TT.vertices.tbl(class.sys = "USDA.TT")
## convert texture fractions to coordinates:
vert$x <- 1-vert$SAND+(vert$SAND-(1-vert$SILT))*0.5
vert$y <- vert$CLAY*sin(pi/3)
## plot(vert$x, vert$y)
#vert.xy <- vert
#coordinates(vert.xy) <- ~x+y
#writeOGR(vert.xy, "vert.shp", "vert", "ESRI Shapefile")

USDA.TT <- data.frame(TT.classes.tbl(class.sys = "USDA.TT", collapse = ", "))
## get the polygon breaks for a class
TT.pnt <- as.list(rep(NA, length(USDA.TT$name)))
poly.lst <- as.list(rep(NA, length(USDA.TT$name)))
for(i in 1:length(USDA.TT$name)){ 
  ## strip the vertices numbers:
  TT.pnt[[i]] <- as.integer(strsplit(unclass(paste(USDA.TT[i, "points"])), ", ")[[1]])
  ## create a list of polygons:
  poly.lst[[i]] <- vert[TT.pnt[[i]],c("x","y")]
  ## add extra point:
  poly.lst[[i]] <- Polygons(list(Polygon(rbind(poly.lst[[i]], poly.lst[[i]][1,]))), ID=i)
}
## convert texture triangle to a spatial object:
poly.sp <- SpatialPolygons(poly.lst, proj4string=CRS(as.character(NA)))
poly.USDA.TT <- SpatialPolygonsDataFrame(poly.sp, data.frame(ID=USDA.TT$name), match.ID=FALSE)
spplot(poly.USDA.TT)
#save(list("get.TF.from.XY", "poly.USDA.TT"), file="poly_USDA_TT.RData")
#writeOGR(poly.USDA.TT, "poly_USDA_TT.shp", "poly_USDA_TT", "ESRI Shapefile")

## function to get TF from XY coordinates in the texture triangle:
get.TF.from.XY <- function(df, xcoord, ycoord) {
    df$CLAY <- df[,ycoord]/sin(pi/3)
    df$SAND <- (2 - df$CLAY - 2 * df[,xcoord]) * 0.5
    df$SILT <- 1 - (df$SAND + df$CLAY)
    return(df)
}


## estimate fractions for some class:
USDA.TT.cnt <- data.frame(t(sapply(slot(poly.USDA.TT, "polygons"), slot, "labpt")))
USDA.TT.cnt$name <- poly.USDA.TT$ID
USDA.TT.cnt <- get.TF.from.XY(USDA.TT.cnt, "X1", "X2")
USDA.TT.cnt[,c("SAND","SILT","CLAY")] <- signif(USDA.TT.cnt[,c("SAND","SILT","CLAY")], 2)
USDA.TT.cnt[,c("name","SAND","SILT","CLAY")]

## estimate uncertainty for some class using simulations:
sim.Cl <- data.frame(spsample(poly.USDA.TT[poly.USDA.TT$ID=="clay",], type="random", n=100))
sim.Cl <- get.TF.from.XY(sim.Cl, "x", "y")
## average error:
sd(sim.Cl$SAND); sd(sim.Cl$SILT); sd(sim.Cl$CLAY)

## plot Africa Soil Profile Data in texture triangle:
require(GSIF)
data(afsp)
tdf <- afsp$horizons[,c("CLYPPT", "SLTPPT", "SNDPPT")]
## remove missing values:
tdf <- tdf[!is.na(tdf$SNDPPT)&!is.na(tdf$SLTPPT)&!is.na(tdf$CLYPPT),] 
## subset to 15%:
tdf <- tdf[runif(nrow(tdf))<.15,]
## the fractions have to sum to 100%:
tdf$Sum = rowSums(tdf)
for(i in c("CLYPPT", "SLTPPT", "SNDPPT")) { tdf[,i] <- tdf[,i]/tdf$Sum * 100 }
names(tdf)[1:3] <- c("CLAY", "SILT", "SAND")
TT.plot(class.sys = "USDA.TT", tri.data = tdf, grid.show = FALSE, pch="+", cex=.4, col="red")


#--------------------------------------------------
# Example with Munsell color codes
#--------------------------------------------------

## load color table for all Munsell values:
load(file("http://gsif.isric.org/lib/exe/fetch.php?media=munsell_rgb.rdata"))
str(munsell.rgb)
library(colorspace)
## pick a random colour:
munsell.rgb[round(runif(1)*2350, 0),]
## convert to HSI values:
as(RGB(R=munsell.rgb[1007,"R"]/255, G=munsell.rgb[1007,"G"]/255, B=munsell.rgb[1007,"B"]/255), "HSV")
aqp::munsell2rgb(the_hue = "10B", the_value = 2, the_chroma = 12)
col2kml("#003A7CFF")
col2rgb("#003A7CFF")

## make a plot of RGB values of first soil horizon:
data(afsp)
paste(afsp$horizons[1:10,"MCOMNS"])
mcol <- join(afsp$horizons[,c("SOURCEID","MCOMNS","UHDICM","LHDICM")], afsp$sites[,c("SOURCEID","LONWGS84","LATWGS84")], type='inner') 
mcol <- mcol[!is.na(mcol$MCOMNS),]
str(mcol)
## reformat values:
mcol$Munsell <- sub(" ", "", sub("/", "_", mcol$MCOMNS))
hue.lst <- expand.grid(c("2.5", "5", "7.5", "10"), c("YR","GY","BG","YE","YN","YY","R","Y","B","G"))
hue.lst$mhue <- paste(hue.lst$Var1, hue.lst$Var2, sep="") 
for(j in hue.lst$mhue[1:28]){ 
  mcol$Munsell <- sub(j, paste(j, "_", sep=""), mcol$Munsell, fixed=TRUE) 
}
mcol$depth <- mcol$UHDICM + (mcol$LHDICM-mcol$UHDICM)/2 
## merge munsell colors and RGB values:
mcol.RGB <- merge(mcol, munsell.rgb, by="Munsell")
str(mcol.RGB)
mcol.RGB <- mcol.RGB[!is.na(mcol.RGB$R),]
mcol.RGB$Rc <- round(mcol.RGB$R/255, 3)
mcol.RGB$Gc <- round(mcol.RGB$G/255, 3)
mcol.RGB$Bc <- round(mcol.RGB$B/255, 3)
mcol.RGB$col <- rgb(mcol.RGB$Rc, mcol.RGB$Gc, mcol.RGB$Bc)
mcol.RGB <- mcol.RGB[mcol.RGB$depth>0 & mcol.RGB$depth<30 & !is.na(mcol.RGB$col),]
coordinates(mcol.RGB) <- ~ LONWGS84+LATWGS84
## plot whole of Africa:
load(file("http://gsif.isric.org/lib/exe/fetch.php?media=admin.af.rda"))
#load("admin.af.rda")
proj4string(admin.af) <- "+proj=longlat +datum=WGS84"
country <- as(admin.af, "SpatialLines")
par(mar=c(.0,.0,.0,.0))
plot(country, col="darkgrey", asp=1)
points(mcol.RGB, pch=21, bg=mcol.RGB$col, col="black")
dev.off()

## plot a single profile in Google Earth
## sample profile from Nigeria:
library(fossil)
lon = 3.90; lat = 7.50; id = "ISRIC:NG0017"; FAO1988 = "LXp" 
top = c(0, 18, 36, 65, 87, 127) 
bottom = c(18, 36, 65, 87, 127, 181)
ORCDRC = c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
hue = c("7.5YR", "7.5YR", "2.5YR", "5YR", "5YR", "10YR")
value = c(3, 4, 5, 5, 5, 7); chroma = c(2, 4, 6, 8, 4, 3)
## prepare a SoilProfileCollection:
prof1 <- join(data.frame(id, top, bottom, ORCDRC, hue, value, chroma), 
   data.frame(id, lon, lat, FAO1988), type='inner')
prof1$soil_color <- with(prof1, munsell2rgb(hue, value, chroma))
depths(prof1) <- id ~ top + bottom
site(prof1) <- ~ lon + lat + FAO1988 
coordinates(prof1) <- ~ lon + lat
proj4string(prof1) <- CRS("+proj=longlat +datum=WGS84")
prof1
## Not run: 
plotKML(prof1, var.name="ORCDRC", color.name="soil_color")

#--------------------------------------------------
# Using MLA's to fit PTFs
#--------------------------------------------------

library(randomForestSRC)
library(ggRandomForests)
library(ggplot2)
library(scales)

## ISRIC WISE data set:
load(file("http://gsif.isric.org/lib/exe/fetch.php?media=sprops.wise.rda"))
str(SPROPS.WISE)
plot(SPROPS.WISE$LONWGS84, SPROPS.WISE$LATWGS84, pch="+")

## Pedo-transfer function for bulk density:
#s.n <- sample.int(nrow(SPROPS.WISE), 10000)
#rfsrc_BD <- rfsrc(BLD ~ ORCDRC + PHIHOX + SNDPPT + CLYPPT + CRFVOL + DEPTH, data=SPROPS.WISE[s.n,])
rfsrc_BD <- rfsrc(BLD ~ ORCDRC + PHIHOX + SNDPPT + CLYPPT + CRFVOL + DEPTH, data=SPROPS.WISE)
rfsrc_BD
plot.variable(rfsrc_BD, xvar.names=c("ORCDRC","PHIHOX","SNDPPT","CLYPPT","CRFVOL","DEPTH"))
#plot.variable(rfsrc_BD, xvar.names=c("ORCDRC","PHIHOX","SNDPPT","CLYPPT","CRFVOL","DEPTH"), partial=TRUE, smooth.lines=TRUE, sorted=TRUE, npts=50)
## Test prediction:
predict(rfsrc_BD, data.frame(ORCDRC=1.2, PHIHOX=7.6, SNDPPT=45, CLYPPT=12, CRFVOL=0, DEPTH=20))$predicted
predict(rfsrc_BD, data.frame(ORCDRC=150, PHIHOX=4.6, SNDPPT=25, CLYPPT=35, CRFVOL=0, DEPTH=20))$predicted


## Pedo-transfer function for soil classification:
load(file("http://gsif.isric.org/lib/exe/fetch.php?media=wise_tax.rda"))
str(WISE_tax)
## Legend:
leg <- read.csv(file("http://gsif.isric.org/lib/exe/fetch.php?media=taxousda_greatgroups.csv"))
str(leg)

## add few numeric soil properties:
x.PHIHOX <- aggregate(SPROPS.WISE$PHIHOX, by=list(SPROPS.WISE$SOURCEID), FUN=mean, na.rm=TRUE); names(x.PHIHOX)[1] = "SOURCEID"
x.CLYPPT <- aggregate(SPROPS.WISE$CLYPPT, by=list(SPROPS.WISE$SOURCEID), FUN=mean, na.rm=TRUE); names(x.CLYPPT)[1] = "SOURCEID"
WISE_tax$PHIHOX <- join(WISE_tax, x.PHIHOX, type="left")$x
WISE_tax$CLYPPT <- join(WISE_tax, x.CLYPPT, type="left")$x

## Model to translate soil classes from one system to the other
WISE_tax.sites <- WISE_tax[complete.cases(WISE_tax[,c("TAXNWRB","PHIHOX","CLYPPT","TAXOUSDA")]),]
## clean up names:
WISE_tax.sites$TAXOUSDA.f <- NA
for(j in leg$Suborder){
  sel <- grep(j, WISE_tax.sites$TAXOUSDA, ignore.case=TRUE)
  WISE_tax.sites$TAXOUSDA.f[sel] = j
}
WISE_tax.sites$TAXOUSDA.f <- as.factor(WISE_tax.sites$TAXOUSDA.f)
WISE_tax.sites$TAXNWRB <- as.factor(paste(WISE_tax.sites$TAXNWRB))
TAXNUSDA.rf <- rfsrc(TAXOUSDA.f ~ TAXNWRB + PHIHOX + CLYPPT, data=WISE_tax.sites)
TAXNUSDA.rf
## test it:
newdata = data.frame(TAXNWRB=factor("Calcaric Cambisol", levels=levels(WISE_tax.sites$TAXNWRB)), PHIHOX=7.8, CLYPPT=12)
x <- data.frame(predict(TAXNUSDA.rf, newdata, type="prob")$predicted)
x[,order(1/x)[1:2]]
