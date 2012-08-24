# Purpose        : Initial settings;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : for more info see [http://cran.r-project.org/doc/manuals/R-exts.html];


################## NEW GSIF CLASSES ##############

## A new class for models fitted in gstat:
setClass("gstatModel", representation(regModel = "glm", vgmModel = "data.frame", sp = "SpatialPoints"), validity = function(object) {
    cn = c("model", "psill", "range", "kappa", "ang1", "ang2", "ang3", "anis1", "anis2")
    if(any(!(names(object@vgmModel) %in% cn)))
      return(paste("Expecting only column names:", paste(cn, collapse=", ")))
    if(!all(cn %in% names(object@vgmModel))){
      x <- cn[!(cn %in% names(object@vgmModel))]
      return(paste("Missing column names:", paste(x, collapse=", "))) 
      }
})

### GSIF soil property maps class:
setClass("GlobalSoilMap", representation (varname = 'character', TimeSpan.begin = 'POSIXct', TimeSpan.end = 'POSIXct', sd1 = 'SpatialPixelsDataFrame', sd2 = 'SpatialPixelsDataFrame', sd3 = 'SpatialPixelsDataFrame', sd4 = 'SpatialPixelsDataFrame', sd5 = 'SpatialPixelsDataFrame', sd6 = 'SpatialPixelsDataFrame'), validity = function(object){
   soilvars = read.csv(system.file("soilvars.csv", package="GSIF"))
   if(object@varname %in% soilvars$varname)
      return(paste("'property'", object@property, "not specified in the Soil Reference Library.", "See", system.file("soilvars.csv", package="GSIF"), "for more details."))
   if(object@TimeSpan.begin > object@TimeSpan.end)
      return("'TimeSpan.begin' must indicate time before or equal to 'TimeSpan.end'") 
   if(ncol(object@sd1)<2|ncol(object@sd2)<2|ncol(object@sd3)<2|ncol(object@sd4)<2|ncol(object@sd5)<2|ncol(object@sd6)<2)
      return("Object in slot 'sd' with at least two realizations (or predictions and variances) required")
   # check the projection system:
   require(plotKML)
   if(!all(check_projection(object@sd1)|check_projection(object@sd2)|check_projection(object@sd3)|check_projection(object@sd4)|check_projection(object@sd5)|check_projection(object@sd6))){
      ref_CRS = get("ref_CRS", envir = GSIF.opts)
      return(paste("The GlobalSoilMap object requires grids to be projected in the", ref_CRS, "projection"))
   }
   # check the target resolution:
   grd.lst <- get("cellsize", envir = GSIF.opts)
   if(!any(object@sd1@grid@cellsize %in% grd.lst)|!any(object@sd2@grid@cellsize %in% grd.lst)|!any(object@sd3@grid@cellsize %in% grd.lst)|!any(object@sd4@grid@cellsize %in% grd.lst)|!any(object@sd5@grid@cellsize %in% grd.lst)|!any(object@sd6@grid@cellsize %in% grd.lst))
      return(paste("Recommended grid cell size does not correspond to one of the following:", paste(signif(grd.lst, 4), collapse=", "))) 
   # check the bounding boxes:
   if(!(any(object@sd1@bbox %in% as.list(object@sd2@bbox, object@sd3@bbox, object@sd4@bbox, object@sd5@bbox, object@sd6@bbox))))
      return("The bounding box of all 'sd' slots is not standard") 
})


## georecord class:
setClass("geosamples", representation (registry = 'character', methods = 'data.frame', data = 'data.frame'), validity = function(object) {
   cnames <- c("observationid", "sampleid", "longitude", "latitude", "locationError", "TimeSpan.begin", "TimeSpan.end", "altitude", "altitudeMode", "sampleArea", "sampleThickness", "observedValue", "methodid", "measurementError")
   if(any(!(names(object@data) %in% cnames)))
      return(paste("Expecting only column names:", paste(cnames, collapse=", ")))
   mnames <- c("methodid", "description", "units", "detectionLimit")
   if(any(!(names(object@methods) %in% mnames)))
      return(paste("Expecting only column names:", paste(mnames, collapse=", ")))
   if(any(!(levels(as.factor(paste(object@methods$methodid))) %in% levels(as.factor(paste(object@data$methodid))))))
      return("'methodid' levels in the methods table and data table do not match")
   if(!any(class(object@data$TimeSpan.begin) %in% "POSIXct") | !any(class(object@data$TimeSpan.end) %in% "POSIXct")) {
      return("'TimeSpan.begin' and 'TimeSpan.end' of class 'POSIXct' required")
      } 
      else {
      sel <- !is.na(object@data$TimeSpan.begin)&!is.na(object@data$TimeSpan.end)  
      if(any(object@data$TimeSpan.begin[sel] > object@data$TimeSpan.end[sel]))
        return("'TimeSpan.begin' must indicate time before or equal to 'TimeSpan.end'")      
      }
   if(any(object@data$measurementError[!is.na(object@data$measurementError)] < 0))
       return("'measurementError' must be positive numbers")
   if(any(object@data$sampleArea[!is.na(object@data$sampleArea)] < 0))
       return("'sampleArea' must be positive numbers")
   if(any(object@data$sampleThickness[!is.na(object@data$sampleThickness)] < 0))
       return("'sampleThickness' must be positive numbers")
   # test if it is a longlat object:
   if(any(object@data$longitude>180|object@data$longitude< -180|object@data$latitude< -90|object@data$latitude> 90))
      return("longitude and latitude values in the range -180 to 180 and -90 to 90 required") 
})

## WPS class
setClass("WPS", representation (server = 'list', inRastername = 'character'), validity = function(object) {
   cnames <- c("URI", "service.name", "version", "request", "identifier")
   if(any(!(names(object@server) %in% cnames)))
      return(paste("Expecting only column names:", paste(cnames, collapse=", ")))
   # check if URI exists:
   uri = paste(paste(object@server$URI, "?", sep=""), object@server$version, object@server$service, "request=GetCapabilities", sep="&") 
   require(RCurl)
   try(z <- getURI(uri, .opts=curlOptions(header=TRUE, nobody=TRUE, transfertext=TRUE, failonerror=FALSE)))
   if(!length(x <- grep(z, pattern="404 Not Found"))==0)
      return("Server error: 404 Not Found")
})

## SpatialComponents class
setClass("SpatialComponents", representation (predicted = "SpatialPixelsDataFrame", pca = "list"), validity = function(object) {
   cnames <- attr(object@pca$rotation, "dimnames")[[1]]
   pnames <- attr(object@pca$rotation, "dimnames")[[2]]
   if(!length(object@pca$sdev)==length(cnames)|!length(object@pca$sdev)==length(pnames))
      return("Number of components of the 'sdev' and 'rotation' objects do not match")
   # check if column names match:
   if(!all(pnames %in% names(object@predicted)))
      return("Column names in the 'predicted' slot and 'pca' slots do not match")
})

## SpatialMemberships class
setClass("SpatialMemberships", representation (predicted = "SpatialPixelsDataFrame", model = "list", mu = "SpatialPixelsDataFrame", class.c = "matrix", class.sd = "matrix", confusion = "matrix"), validity = function(object) {
   # check if column names match:
   if(!all(names(object@mu) %in% levels(object@predicted@data[,1])))
      return("Class names in the 'predicted' and 'mu' slots do not match")
   # check if the row names in the class.sd, class.c match:
   if(!all(row.names(object@class.c) %in% levels(object@predicted@data[,1])))
      return("Row names in the 'class.c' slot and 'predicted' slots do not match")
   if(!all(row.names(object@class.sd) %in% levels(object@predicted@data[,1])))
      return("Row names in the 'class.sd' slot and 'predicted' slots do not match")
   if(ncol(object@mu@data)<2)
      return("A minimum of two membership maps required")   
   # check if all mu's sum to 1 (plus minus 1%):
   if(!all(rowSums(object@mu@data)>.99&rowSums(object@mu@data)<1.01))
      return("Some rows in the 'mu' slot do not sum up to 1")
   # check if the confusion matrix has kappa > 0
   if(length(object@confusion)==0|attr(object@confusion, "error")==0)
      return("Not possible to derive confusion table or no significant match detected")
})



################## generic functions ##############

if(!isGeneric("getID")){
  setGeneric("getID", function(obj, ...){standardGeneric("getID")})
}

if(!isGeneric("as.data.frame")){
  setGeneric("as.data.frame", function(x, row.names = NULL, optional = FALSE, ...){standardGeneric("as.data.frame")})
}

if(!isGeneric("predict")){
  setGeneric("predict", function(object, ...){standardGeneric("predict")})
}

if(!isGeneric("mpspline")){
  setGeneric("mpspline", function(obj, ...){standardGeneric("mpspline")})
}

if(!isGeneric("as.geosamples")){
  setGeneric("as.geosamples", function(obj, ...){standardGeneric("as.geosamples")})
}

if(!isGeneric("getProcess")){
  setGeneric("getProcess", function(x, ...){standardGeneric("getProcess")})
}

if(!isGeneric("describe")){
  setGeneric("describe", function(x, ...){standardGeneric("describe")})
}

if(!isGeneric("merge")){
  setGeneric("merge", function(x, y, ...){standardGeneric("merge")})
}

if(!isGeneric("subset")){
  setGeneric("subset", function(x, ...){standardGeneric("subset")})
}

if(!isGeneric("spc")){
  setGeneric("spc", function(obj, formulaString, ...){standardGeneric("spc")})
}

if(!isGeneric("make.3Dgrid")){
  setGeneric("make.3Dgrid", function(obj, ...){standardGeneric("make.3Dgrid")})
}

if (!isGeneric("fit.gstatModel")){
  setGeneric("fit.gstatModel", function(observations, formulaString, covariates, ...){standardGeneric("fit.gstatModel")})
}

if (!isGeneric("spmultinom")){
  setGeneric("spmultinom", function(formulaString, observations, covariates, ...){standardGeneric("spmultinom")})
}

if (!isGeneric("validate")){
  setGeneric("validate", function(obj, ...){standardGeneric("validate")})
}

if (!isGeneric("spfkm")){
  setGeneric("spfkm", function(formulaString, observations, covariates, ...){standardGeneric("spfkm")})
}

if (!isGeneric("sp3D")){
  setGeneric("sp3D", function(obj, ...){standardGeneric("sp3D")})
}

if (!isGeneric("write.data")){
  setGeneric("write.data", function(obj, ...){standardGeneric("write.data")})
}

if (!isGeneric("gdalwarp")){
  setGeneric("gdalwarp", function(obj, ...){standardGeneric("gdalwarp")})
}

if (!isGeneric("MaxEnt")){
  setGeneric("MaxEnt", function(occurrences, covariates, ...){standardGeneric("MaxEnt")})
}

################## STANDARD ENVIRONMENTS ##############

## setup the plotKML environment:
GSIF.opts <- new.env(hash=TRUE)

## Standard settings:
GSIF.env <- function(
    wps.server = "http://wps.worldgrids.org",
    ref_CRS = "+proj=longlat +datum=WGS84",
    NAflag = -99999,
    license_url = "http://creativecommons.org/licenses/by/3.0/",
    project_url = "http://gsif.r-forge.r-project.org/",
    stdepths = c(-2.5, -10, -22.5, -45, -80, -150)/100,
    stsize = c(5, 10, 15, 30, 40, 100)/100,
    cellsize = rev(c(6/120, 3/120, 1/120, 1/240, 1/600, 1/1200, 1/3600)),
    show.env = TRUE
    ){
    
    assign("wps.server", wps.server, envir=GSIF.opts)
    assign("ref_CRS", ref_CRS, envir=GSIF.opts)
    assign("NAflag", NAflag, envir=GSIF.opts)
    assign("license_url", license_url, envir=GSIF.opts)
    assign("project_url", project_url, envir=GSIF.opts)
    assign("stdepths", stdepths, envir=GSIF.opts)
    assign("stsize", stsize, envir=GSIF.opts)
    assign("cellsize", cellsize, envir=GSIF.opts)
    
    GSIF.opts <- list(wps.server, ref_CRS, NAflag, license_url, project_url, stdepths, stsize, cellsize)
    names(GSIF.opts) <- c("location of the WPS", "referent CRS", "NA flag value", "lisence URL", "project home", "standard depths", "standard thicknesses", "grid cell size")
    
    if(show.env){  return(GSIF.opts)  }
 
}

# load GSIF.opts with some basic information
GSIF.env(show.env = FALSE)

# end of script;
