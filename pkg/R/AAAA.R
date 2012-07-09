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
      return(paste("Expecting only column names:", cn))
    if(!all(cn %in% names(object@vgmModel))){
      x <- cn[!(cn %in% names(object@vgmModel))]
      return(paste("Missing column names:", x)) 
      }
})

### GSIF soil property maps class:
setClass("GlobalSoilMap", representation (varname = 'character', sd1 = 'SpatialPixelsDataFrame', sd2 = 'SpatialPixelsDataFrame', sd3 = 'SpatialPixelsDataFrame', sd4 = 'SpatialPixelsDataFrame', sd5 = 'SpatialPixelsDataFrame', sd6 = 'SpatialPixelsDataFrame', model = 'list', validation = 'SpatialPointsDataFrame'), validity <- function(obj) {
   # Check column names:
   soilvars = read.csv(system.file("soilvars.csv", package="GSIF"))
   if(obj@varname %in% soilvars$varname)
      warning(paste("'property'", obj@property, "not specified in the Soil Reference Library.", "See", system.file("soilvars.csv", package="GSIF"), "for more details."))
   if(ncol(obj@sd1)<2)
      return("Object in slot 'sd' with at least two realizations required")
   if(ncol(obj@sd1)<5)
      warning("Using <5 realizations can result in artifacts")
   # check the projection system:
   require(plotKML)
   if(check_projection(obj@sd1)){
      ref_CRS = get("ref_CRS", envir = plotKML.opts)
      return(paste("The GlobalSoilMap object requires grids to be projected in the", ref_CRS, "projection"))
   }
   # check the target resolution:
   grd.lst <- c(6/120,3/120,1/120,1/240,1/600,3/3600)
   if(!(any(eberg_grid@grid@cellsize) %in% grd.lst))
      warning(paste("Recommended grid cell size does not correspond to one of the following:", signif(grd.lst, 4))) 
   # check if validation slot is complete:
   ov <- extract(brick(obj@sd1), obj@validation)
   if(length(ov)==0)
      return("'sd1' and 'validation' do not overlap spatially")
})

## georecord class:
setClass("geosamples", representation (registry = 'character', methods = 'data.frame', data = 'data.frame'), validity <- function(obj) {
   cnames <- c("observationid", "sampleid", "longitude", "latitude", "locationError", "TimeSpan.begin", "TimeSpan.end", "altitude", "altitudeMode", "volume", "observedValue", "methodid", "measurementError")
   if(any(!(names(obj@data) %in% cnames)))
      return(paste("Expecting only column names:", cnames))
   mnames <- c("methodid", "description", "units", "detectionLimit")
   if(any(!(names(obj@methods) %in% mnames)))
      return(paste("Expecting only column names:", mnames))
   if(!length(levels(as.factor(paste(obj@methods$methodid))))==length(as.factor(paste(obj@data$methodid))))
      return("'methodid' levels in the methods table and data table do not match")
   if(!is.na(obj@data$TimeSpan.begin)){
      if(!any(class(obj@data$TimeSpan.begin) %in% "POSIXct") | !any(class(obj@data$TimeSpan.end) %in% "POSIXct")){
        return("'TimeSpan.begin' and 'TimeSpan.end' of class 'POSIXct' required")
      } 
      else {  
        if((obj@data$TimeSpan.end - obj@data$TimeSpan.begin)<0)
        return("'TimeSpan.end' must have positive time difference from 'TimeSpan.begin'")
      }
   }
   if(!is.na(obj@data$measurementError)){ 
      if(obj@data$measurementError < 0)
       return("'measurementError' must be positive numbers")
   }
   if(!is.na(obj@data$volume)){ 
      if(obj@data$volume < 0)
       return("'volume' must be positive numbers")
   }
   if(!(obj@data$dimension == 0|obj@data$dimension == 1|obj@data$dimension == 2|obj@data$dimension == 3))
      return("'dimension' must be an integer number: 0, 1, 2, or 3")
   # test if it is a longlat object:
   if(any(obj@data$longitude>180|obj@data$longitude< -180|obj@data$latitude< -90|obj@data$latitude> 90))
      return("longitude and latitude values in the range -180 to 180 and -90 to 90 required") 
})

## WPS class
setClass("WPS", representation (server = 'list', inRastername = 'character'), validity <- function(obj) {
   cnames <- c("URI", "service.name", "version", "request", "identifier")
   if(any(!(names(obj@server) %in% cnames)))
      return(paste("Expecting only column names:", cnames))
   # check if URI exists:
   uri = paste(paste(obj@server$URI, "?", sep=""), obj@server$version, obj@server$service, "request=GetCapabilities", sep="&") 
   require(RCurl)
   try(z <- getURI(uri, .opts=curlOptions(header=TRUE, nobody=TRUE, transfertext=TRUE, failonerror=FALSE)))
   if(!length(x <- grep(z, pattern="404 Not Found"))==0)
      return("Server error: 404 Not Found")
})

## SpatialComponents class
setClass("SpatialComponents", representation (predicted = "SpatialPixelsDataFrame", pca = "list"), validity <- function(obj) {
   cnames <- attr(obj@pca$rotation, "dimnames")[[1]]
   pnames <- attr(obj@pca$rotation, "dimnames")[[2]]
   if(!length(obj@pca$sdev)==length(cnames)|!length(obj@pca$sdev)==length(pnames))
      return("Number of components of the 'sdev' and 'rotation' objects do not match")
   # check if column names match:
   if(!all(pnames %in% names(obj@predicted)))
      return("Column names in the 'predicted' slot and 'pca' slots do not match")
})

## SpatialMemberships class
setClass("SpatialMemberships", representation (predicted = "SpatialPixelsDataFrame", model = "list", mu = "SpatialPixelsDataFrame", class.c = "matrix", class.sd = "matrix", confusion = "matrix"), validity <- function(obj) {
   # check if column names match:
   if(!all(levels(obj@predicted@data[,1]) %in% names(obj@mu)))
      return("Class names in the 'predicted' and 'mu' slots do not match")
   # check if the row names in the class.sd, class.c match:
   if(!all(levels(obj@predicted@data[,1]) %in% row.names(obj@class.c)))
      return("Row names in the 'class.c' slot and 'predicted' slots do not match")
   if(!all(levels(obj@predicted@data[,1]) %in% row.names(obj@class.sd)))
      return("Row names in the 'class.sd' slot and 'predicted' slots do not match")
   if(ncol(obj@mu)<2)
      return("A minimum of two membership maps required")   
   # check if all mu's sum to 1:
   if(!all(rowSums(obj@mu)==1))
      return("Some rows in the 'mu' slot do not sum up to 1")
   # check if the confusion matrix has kappa > 0
   if(length(obj@confusion)==0|attr(obj@confusion, "error")==0)
      return("Not possible to derive confusion table or no significant match detected")
})



################## generic functions ##############

if(!isGeneric("getID")){
  setGeneric("getID", function(obj, ...){standardGeneric("getID")})
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

if(!isGeneric("as.data.frame.default")) {
	setGeneric("as.data.frame.default", function(x, ...){standardGeneric("as.data.frame.default")})
}	

if(!isGeneric("describe")){
  setGeneric("describe", function(x, ...){standardGeneric("describe")})
}

if(!isGeneric("spc")){
  setGeneric("spc", function(obj, formulaString, ...){standardGeneric("spc")})
}

if (!isGeneric("fit.gstatModel")){
  setGeneric("fit.gstatModel", function(observations, formulaString, covariates, ...){standardGeneric("fit.gstatModel")})
}

if (!isGeneric("spmultinom")){
  setGeneric("spmultinom", function(formulaString, rmatrix, newdata, ...){standardGeneric("spmultinom")})
}

if (!isGeneric("spfkm")){
  setGeneric("spfkm", function(formulaString, observations, covariates, ...){standardGeneric("spfkm")})
}


################## STANDARD ENVIRONMENTS ##############

## setup the plotKML environment:
GSIF.opts <- new.env(hash=TRUE)

## Standard settings:
GSIF.env <- function(
    wps.server = '',
    ref_CRS,
    NAflag,
    license_url,
    project_url,
    show.env = TRUE
    ){
    
    # require(plotKML)
    if(missing(wps.server)) { wps.server <- "http://wps.worldgrids.org" }
    if(missing(ref_CRS)) { ref_CRS <- "+proj=longlat +datum=WGS84" }
    if(missing(NAflag)) { NAflag <- -99999 }
    if(missing(license_url)) { license_url <- "http://creativecommons.org/licenses/by/3.0/" }
    if(missing(project_url)) { project_url <- "http://gsif.r-forge.r-project.org/" }
 
    assign("wps.server", wps.server, envir=GSIF.opts)
    assign("ref_CRS", ref_CRS, envir=GSIF.opts)
    assign("NAflag", NAflag, envir=GSIF.opts)
    assign("license_url", license_url, envir=GSIF.opts)
    assign("project_url", project_url, envir=GSIF.opts)
    
    GSIF.opts <- list(wps.server, ref_CRS, NAflag, license_url, project_url)
    names(GSIF.opts) <- c("location of the WPS", "referent CRS", "NA flag value", "lisence URL", "project home")
    
    if(show.env){  return(GSIF.opts)  }
 
}

# load GSIF.opts with some basic information
GSIF.env(show.env = FALSE)

# end of script;
