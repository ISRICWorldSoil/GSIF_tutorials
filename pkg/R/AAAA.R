# Purpose        : Initial settings;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : for more info see [http://cran.r-project.org/doc/manuals/R-exts.html];


################## NEW GSIF CLASSES ##############

### GSIF soil property maps class:
#setClass("GlobalSoilMap", representation (property = 'character', raster = 'RasterBrick', model = 'list', validationPoints = 'SpatialPointsDataFrame', spMetadata = 'SpatialMetadata'), validity <- function(obj) {
#   if(obj@bounds)
#      return('vector with (upper and lower limits) required')
#})

## georecord class:
setClass("geosamples", representation (registry = 'character', methods = 'data.frame', data = 'data.frame'), validity <- function(obj) {
   cnames <- c("sampleid", "producerid", "longitude", "latitude", "TimeSpan.begin", "TimeSpan.end", "altitude", "altitudeMode", "dimension", "volume", "value", "methodid", "measurementError")
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



################## generic functions ##############

if (!isGeneric("getID")){
  setGeneric("getID", function(obj, ...){standardGeneric("getID")})
}

if (!isGeneric("mpspline")){
  setGeneric("mpspline", function(obj, ...){standardGeneric("mpspline")})
}

if (!isGeneric("as.geosamples")){
  setGeneric("as.geosamples", function(obj, ...){standardGeneric("as.geosamples")})
}

## internal methods:


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
    
#    require(plotKML)
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
