## INTRO TO R!
## http://www.geostat-course.org/Baby_steps_R
## maintained by: T. Hengl (tom.hengl@wur.nl)

# See also: 
## 1. Quick-R [http://www.statmethods.net];
## 2. R by example [http://www.mayin.org/ajayshah/KB/R/index.html];
## 3. Statistics with R on youtube [http://www.youtube.com/watch?v=n2RZp8WSIgI]

## R: session management
## Your R objects are stored in a workspace
## To list the objects in your workspace:  
ls()

## To remove ALL objects in your workspace:
rm(list=ls()) # or use Remove all objects in the Misc menu

## To save your workspace to a file, you may type  
save.image() # or use Save Workspace. in the File menu

## Assignment "<-" used to indicate assignment
x <- c(1,2,3,4,5,6,7)
x <- c(1:7)
x <- 1:7
## note: as of version 1.4 "=" is also a valid assignment operator

## R as a calculator:
5 + (6 + 7) * pi^2
log(exp(1))
log(1000, 10)
sin(pi/3)^2 + cos(pi/3)^2
Sin(pi/3)^2 + cos(pi/3)^2
## R is case sensitive! 
## Error: couldn't find function "Sin"

log2(32)
sqrt(2)
seq(0, 5, length=6)
# plots in R:
plot(sin(seq(0, 2*pi, length=100)))

## Basic (atomic) data types
## Logical
x -> y
x <- T
x <- T; y <- F
x; y
## Numerical
a <- 5; b <- sqrt(2)
a; b
## Character
a <- "1"; b <- 1
a; b
a <- "character"
b <- "a"; c <- a
a; b; c

## Vectors, Matrices, Arrays
## Vector: Ordered collection of data of the same data type
x <- c(5.2, 1.7, 6.3) 
log(x) 
y <- 1:8
z <- seq(from=1, to=1.4, by = 0.1)
y + z
length(y)
mean(y + z)
## Matrix: Rectangular table of data of the same type
m <- matrix(1:12, 4, byrow = T); m 
y <- -1:2
m.new <- m + y
t(m.new)
dim(m)
dim(t(m.new))

## R is designed to handle statistical data and therefore predestined to deal with missing values / Numbers that are "not available"
x <- c(1, 2, 3, NA)
x + 3
## "Not a number"
log(c(0, 1, 2))
0/0

## Subsetting
## It is often necessary to extract a subset of a vector or matrix
## R offers a couple of neat ways to do that
x <- c("a", "b", "c", "d", "e", "f", "g", "h")
x <- letters[1:8]
x[1]
x[3:5]
x[-(3:5)]
x[c(T, F, T, F, T, F, T, F)]
x[x <= "d"]
m[,2]
m[3,]

## Lists are very flexible
## Lists: Session
Empl <- list(employee="Anna", spouse="Fred", children=3, child.ages=c(4,7,9))
str(Empl)
Empl[[4]]
Empl$child.a
Empl[4]   # a sublist consisting of the 4th component of Empl
names(Empl) <- letters[1:4]
Empl <- c(Empl, service=8)
unlist(Empl)  # converts it to a vector. Mixed types will be converted to character, giving a character vector.

## with matrices
x.mat <- matrix(c(3, -1, 2, 2, 0, 3), ncol = 2)
x.mat	   
dimnames(x.mat) <- list(c("L1","L2","L3"), c("R1","R2"))
x.mat

## subsetting:
cw <- chickwts
## question: where did this object "chickwts" come from?

## Let's look at it:
str(cw)
cw[3,2]
View(cw)
## Labels in data frames
labels(cw)
## you can edit the data using fix:
fix(cw)
## but this is not a good idea! better edit manually:
x <- which(cw$weight==136)
cw[x,"weight"] <- 140

## Control structures and functions
## Loops in R
for(i in 1:length(x)) {
    x[i] <- rnorm(1)
}
x

## Functions:
add <- function(a,b) 
{ result = a+b
  return(result) }
add(2,5)

## General Form of Functions
# function(arguments) {
#       expression
#     }

## GETTING HELP:
## If you are in doubt...
help(predict)
## 'predict' is a generic function for predictions from the results
##  of various model fitting functions. 
help(predict.lm)
##  'predict.lm' produces predicted values, obtained by evaluating the regression function in the frame 'newdata' 
data(cars)
plot(cars)
str(cars)
cars$fit <- signif(predict(lm(dist~speed, cars), cars), 3)

## see also:
system("open http://cran.r-project.org/doc/contrib/Short-refcard.pdf")

## Calling Conventions for Functions
## Arguments may be specified in the same order in which they occur in function definition, in which case the values are supplied in order.
## Arguments may be specified as name=value, when the order in which the arguments appear is irrelevant. 
## Above two rules can be mixed.

x1 <- rnorm(100)
y1 <- rnorm(100)
t.test(x1, y1, var.equal=F, conf.level=.99)
t.test(var.equal=F, conf.level=.99, x1, y1)

## Plots
x <- 1:10
y <- 2*x + rnorm(10,0,1)
plot(x,y,type="p") # Try l,b,s,o,h,n

## Interactive Graphics Functions
locator(x, type="p") 
## Waits for the user to select locations on the current plot using the left mouse button. This continues until n (default 500) points have been selected.

identify(x, y, labels = seq_along(x)) 
## Allow the user to highlight any of the points defined by x and y. 

text(x,y,"Hey")
## Write text at coordinate x,y.

## Plots for Multivariate Data
pairs(stack.x)
x <- 1:20/20
y <- 1:20/20
z <- outer(x,y,function(a,b){cos(10*a*b)/(1+a*b^2)})
contour(x,y,z)
persp(x,y,z)
image(x,y,z)

## Histogram:
par(mfrow=c(1,2))   # set up multiple plots
simdata <- rchisq(100,8)
hist(simdata)  # default number of bins
hist(simdata, breaks=2)  # etc,4,20

## Working with spatial data
## install and load additional R packages:
require(sp); require(rgdal); require(RSAGA)

## install SAGA GIS (http://sourceforge.net/projects/saga-gis/) and GDAL (http://gdal.org)

## from table data to SpatialPointsDataFrame:
data(meuse)
str(meuse)
coordinates(meuse) <- ~x+y
str(meuse, max.level=2)
## plot this data:
spplot(meuse[1], axes = TRUE)

## We can create spatial objects from scratch! For example a DEM:
dem <- expand.grid(x=seq(100,600,100), y=seq(100,600,100), KEEP.OUT.ATTRS=FALSE)
dem$Z <- as.vector(c(23, 24, 34, 38, 45, 51, 24, 20, 20, 28, 18, 49, 22, 20, 19, 14, 38, 45, 19, 15 ,13, 21, 23, 25, 14, 11, 18, 11, 18, 19, 10, 16, 23, 16, 9, 6))
str(dem)
## 3D Real-Time Visualization Device System for R:
library(rgl)
rgl.open()
terrain3d(x=seq(100,600,100), y=seq(100,600,100), z=dem$Z*5, back="lines")
## convert to a sp object:
gridded(dem) <- ~x+y
dem <- as(dem, "SpatialGridDataFrame")
str(dem)
spplot(dem[1],  main="DEM", scales = list(draw = FALSE), col.regions=topo.colors(25))
## export to GIS format:
writeGDAL(dem, "dem6.sdat", "SAGA", mvFlag = -99999)

## Controlling SAGA from R:
rsaga.env()
## SAGA Module Library Documentation (http://www.saga-gis.org/saga_module_doc/)
rsaga.geoprocessor(lib="ta_channels", module=5, param=list(DEM="dem6.sgrd", DIRECTION="channels.sgrd", CONNECTION="route.sgrd", SEGMENTS="channels.shp", BASINS="tmp.sgrd", THRESHOLD=1), check.module.exists = FALSE, warn=FALSE) 
## import the derived streams:
dem$route <- readGDAL("route.sdat")$band1
channels <- readOGR("channels.shp", "channels")
spplot(dem[2],col.regions=rev(gray(0:20/20)), main="Flow connectivity", sp.layout=list("sp.lines",channels,col="red"))

## Controlling GDAL from R
gdal.path <- "C:\\Program Files\\GDAL"
## the path depends on your Operating system!
library(gdalUtils)
gdal_setInstallation(gdal.path)
## The version number:
getOption("gdalUtils_gdalPath")[[1]]$version

## downloading data from web (MODIS Leaf Area Index):
library(RCurl)
MOD15A2 <- "ftp://ladsftp.nascom.nasa.gov/allData/5/MOD15A2/2015/121/"
download.file(paste0(MOD15A2, "MOD15A2.A2015121.h18v03.005.2015133080613.hdf"), "MOD15A2.A2015121.h18v03.005.2015133080613.hdf") 
GDALinfo("MOD15A2.A2015121.h18v03.005.2015133080613.hdf")
## unfortunatelly, rgdal does not recognize this format
gdalinfo("MOD15A2.A2015121.h18v03.005.2015133080613.hdf")

## define the local coordinate system:
NL.prj <- "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs +towgs84=565.237,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812"
## specify the working dir using the Unix convention:
workd <- shortPathName(normalizePath(getwd()))
outn <- "LAI_2015_121_NL.tif"
unlink(outn)
gdalwarp(paste0('HDF4_EOS:EOS_GRID:\"', workd, '\\MOD15A2.A2015121.h18v03.005.2015133080613.hdf\":MOD_Grid_MOD15A2:Lai_1km'), outn, t_srs=NL.prj, r="near", te=c(0, 300000, 280000, 625000),  tr=c(500, 500))
LAI_2012_04_22_NL <- readGDAL("LAI_2015_121_NL.tif")
library(raster)
plot(raster(LAI_2012_04_22_NL))
## description of the dataset is at:
## [https://lpdaac.usgs.gov/products/modis_products_table/]

# end of script;