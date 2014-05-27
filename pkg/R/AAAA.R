# Purpose        : Initial settings;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Dylan Beaudette (dylan.beaudette@gmail.com); 
# Dev Status     : Pre-Alpha
# Note           : Aqp classes described here -> [http://r-forge.r-project.org/projects/aqp/]; for more info see [http://cran.r-project.org/doc/manuals/R-exts.html];


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
    REST.server = 'http://rest.soilgrids.org/',
    attributes = c("ORCDRC","PHIHOX","SNDPPT","SLTPPT","CLYPPT","CFRVOL","CEC","BLD","TAXGWRB","TAXOUSDA"),
    TimeSpan = list(begin=as.POSIXct("1950-01-01"), end=as.POSIXct("2005-12-30")),
    show.env = TRUE
    ){
    
    md.lst <- list(wps.server=wps.server, ref_CRS=ref_CRS, NAflag=NAflag, license_url=license_url, project_url=project_url, stdepths=stdepths, stsize=stsize, cellsize=cellsize, REST.server=REST.server, attributes=attributes, TimeSpan=TimeSpan)
    
    x <- lapply(names(md.lst), function(x){ assign(x, md.lst[[x]], envir=GSIF.opts) })  
    if(show.env){  
      return(md.lst)  
    } 
}

# load GSIF.opts with some basic information
GSIF.env(show.env = FALSE)

################## NEW GSIF CLASSES ##############

## Copy of the 'SoilProfileCollection' class basically (see [http://aqp.r-forge.r-project.org/aqp-html-manual/]):	
setClass(Class="FAO.SoilProfileCollection", 
  representation=representation(
    idcol='character', # column name containing IDs
    depthcols='character', # 2 element vector with column names for hz top, bottom
    metadata='data.frame', # single-row dataframe with key-value mapping
    horizons='data.frame', # all horizons sorted by ID, top
    site='data.frame', # data about the sampling sites
    sp='SpatialPoints', # (optional) spatial data stored here
    diagnostic='data.frame' # (optional) diagnostic horizons are stored here
  ),
  prototype=prototype(
    idcol='SOURCEID',
    depthcols=c('UHDICM','LHDICM'),
    metadata=data.frame(stringsAsFactors=FALSE), # default units are unkown
    horizons=data.frame(SOURCEID=NA, UDICM=0, LHDICM=200, stringsAsFactors=FALSE),
    site=data.frame(SOURCEID=NA, SPDFAO=1, SOURCEDB=NA, stringsAsFactors=FALSE),
    sp=new('SpatialPoints'),
    diagnostic=data.frame(stringsAsFactors=FALSE)
  ),
  validity=function(object) {

	  ## check horizon logic:
  	h <- object@horizons
  	top <- object@depthcols[1]
    bottom <- object@depthcols[2]	
    if(any(c(is.na(h[[top]]), is.na(h[[bottom]])))) {
  		return("Horizon top and bottom values cannot contain NA values")
 	  }
    test.h <- !h[[top]] < h[[bottom]] 
    if(any(test.h)){
      return("Invalid horizon bottom values found at row:", paste(which(test.h), collapse=", "))
    }
    ## check column names:
    if(any(!names(object@site) %in% soil.vars$varname)|any(!names(object@horizons) %in% soil.vars$varname)){
      test.nm <- !(names(object@site) %in% soil.vars$varname)
      return(paste("Invalid variable name used:", paste(names(object@site)[test.nm], collapse=", ", sep="")))
      test.nm <- !(names(object@horizons) %in% soil.vars$varname)
      return(paste("Invalid variable name used:", paste(names(object@horizons)[test.nm], collapse=", ", sep="")))
    }
    ## check that all required columns are available:
    required <- paste(soil.vars[soil.vars$priority=="required" & (soil.vars$spcslot=="sites"|soil.vars$spcslot=="horizons"),"varname"])
    present <- c(names(object@site), names(object@horizons))
    missing <- !required %in% present
    if(sum(missing)>0){
      return(paste("Missing variable names:", paste(required[missing], collapse=", ", sep="")))
    }
    message("Checking domains...")
    ## munsell colour codes:
    if(any(names(object@horizons) %in% "DCOMNS")){
      if(any(!levels(as.factor(object@horizons$DCOMNS)) %in% levels(munsell$Munsell))){
        message("Removing Munsell colour codes not available in the domain table")
        x <- merge(object@horizons["DCOMNS"], munsell, by.x="DCOMNS", by.y="Munsell", all.x=TRUE, sort=FALSE)
        object@horizons$DCOMNS <- ifelse(is.na(x$R), NA, x$DCOMNS)
      }      
    }
    if(any(names(object@horizons) %in% "MCOMNS")){
      if(any(!levels(as.factor(object@horizons$MCOMNS)) %in% levels(munsell$Munsell))){
        message("Removing Munsell colour codes not available in the domain table")
        x <- merge(object@horizons["MCOMNS"], munsell, by.x="MCOMNS", by.y="Munsell", all.x=TRUE, sort=FALSE)
        object@horizons$MCOMNS <- ifelse(is.na(x$R), NA, x$MCOMNS)
      }      
    }
    ## check domains in the site table:
    for(j in 1:ncol(object@site)){
      vtype <- soil.vars[soil.vars$varname==names(object@site)[j],"type"]
      if(vtype=="factor"){
        DomainId <- soil.vars[soil.vars$varname==names(object@site)[j],"DomainId"]
        if(!is.na(DomainId)){
          levs <- paste(unlist(soil.dom[soil.dom$DomainId == DomainId,"Value"]))
          if(any(!levels(object@site[,j]) %in% levs)){
            return(paste("Invalid domain used for variable:", names(object@site)[j]))
          }
        }
      } else {
        ## remove all values outside the natural range:
        if(!(names(object@site)[j]=="TIMESTRT"|names(object@site)[j]=="TIMEENDR")){
          minval <- soil.vars[soil.vars$varname==names(object@site)[j],"minval"]
          maxval <- soil.vars[soil.vars$varname==names(object@site)[j],"maxval"]
          object@site[,j] <- ifelse(object@site[,j] < minval, NA, ifelse(object@site[,j] > maxval, NA, object@site[,j]))
        }
      }
    }
    ## check metadata slot    
})

## A new class for models fitted in gstat:
setClass("gstatModel", representation(regModel = "ANY", vgmModel = "data.frame", sp = "SpatialPoints"), validity = function(object) {
    ml = c("lm", "glm", "rpart", "randomForest", "lme", "gls", "zeroinfl")
    if(!any(class(object@regModel) %in% ml))
      return(paste("Only models of type", paste(ml, collapse=", "), "are accepted"))
    cn = c("model", "psill", "range", "kappa", "ang1", "ang2", "ang3", "anis1", "anis2")
    if(any(!(names(object@vgmModel) %in% cn)))
      return(paste("Expecting only column names:", paste(cn, collapse=", ")))
    if(!all(cn %in% names(object@vgmModel))){
      x <- cn[!(cn %in% names(object@vgmModel))]
      return(paste("Missing column names:", paste(x, collapse=", "))) 
      }
})

### GSIF soil property maps class:
setClass("SoilGrids", representation (varname = 'character', TimeSpan = 'list', sd1 = 'SpatialPixelsDataFrame', sd2 = 'SpatialPixelsDataFrame', sd3 = 'SpatialPixelsDataFrame', sd4 = 'SpatialPixelsDataFrame', sd5 = 'SpatialPixelsDataFrame', sd6 = 'SpatialPixelsDataFrame'), 
   prototype = list(varname = "NA", TimeSpan = list(begin=Sys.time(), end=Sys.time()), sd1 = NULL, sd2 = NULL, sd3 = NULL, sd4 = NULL, sd5 = NULL, sd6 = NULL), ## will not pass the validity check!
   validity = function(object){
   load(soil.vars)
   if(object@varname %in% soil.vars$varname){
      return(paste("Property", object@varname, "not specified in the Soil Reference Library. See data(soil.vars) for more details."))
   }
   if(!all(sapply(object@TimeSpan, function(x){class(x)[1]})=="POSIXct") & object@TimeSpan[["begin"]] > object@TimeSpan[["end"]]){
      return("'TimeSpan' must indicate 'begin' and 'end' times to which the predictions refer to.") 
   }
   if(ncol(object@sd1)<2|ncol(object@sd2)<2|ncol(object@sd3)<2|ncol(object@sd4)<2|ncol(object@sd5)<2|ncol(object@sd6)<2){
      return("Object in slot 'sd' with at least two realizations (or predictions and variances) required")
   }
   ## check the projection system:
   if(!all(check_projection(object@sd1)|check_projection(object@sd2)|check_projection(object@sd3)|check_projection(object@sd4)|check_projection(object@sd5)|check_projection(object@sd6))){
      ref_CRS = get("ref_CRS", envir = GSIF.opts)
      return(paste("Grids projected in the \"", ref_CRS, "\" projection required.", sep=""))
   }
   ## check the target resolution:
   grd.lst <- get("cellsize", envir = GSIF.opts)
   if(!any(object@sd1@grid@cellsize %in% grd.lst)|!any(object@sd2@grid@cellsize %in% grd.lst)|!any(object@sd3@grid@cellsize %in% grd.lst)|!any(object@sd4@grid@cellsize %in% grd.lst)|!any(object@sd5@grid@cellsize %in% grd.lst)|!any(object@sd6@grid@cellsize %in% grd.lst)){
      return(paste("Grid cell size does not correspond to one of the following:", paste(signif(grd.lst, 4), collapse=", "))) 
   }
   ## check the bounding boxes:
   if(!(any(object@sd1@bbox %in% as.list(object@sd2@bbox, object@sd3@bbox, object@sd4@bbox, object@sd5@bbox, object@sd6@bbox)))){
      return("The bounding box of all 'sd' slots is not standard") 
   }
})


### GlobalSoilMap class (must be 100 m):
setClass("GlobalSoilMap", representation (varname = 'character', TimeSpan = 'list', sd1 = 'SpatialPixelsDataFrame', sd2 = 'SpatialPixelsDataFrame', sd3 = 'SpatialPixelsDataFrame', sd4 = 'SpatialPixelsDataFrame', sd5 = 'SpatialPixelsDataFrame', sd6 = 'SpatialPixelsDataFrame'), 
   prototype = list(varname = "NA", TimeSpan = list(begin=Sys.time(), end=Sys.time()), sd1 = NULL, sd2 = NULL, sd3 = NULL, sd4 = NULL, sd5 = NULL, sd6 = NULL), validity = function(object){
   if(object@TimeSpan.begin > object@TimeSpan.end){
      return("'TimeSpan.begin' must indicate time before or equal to 'TimeSpan.end'") 
   }
   ## check the target resolution:
   grd.lst <- get("cellsize", envir = GSIF.opts)
   if(!all(object@sd1@grid@cellsize == grd.lst[2])|!all(object@sd2@grid@cellsize == grd.lst[2])|!all(object@sd3@grid@cellsize == grd.lst[2])|!all(object@sd4@grid@cellsize == grd.lst[2])|!all(object@sd5@grid@cellsize == grd.lst[2])|!all(object@sd6@grid@cellsize == grd.lst[2])){
      return(paste("Grid cell size does not correspond the prescribed resolution:", paste(signif(grd.lst[2], 4), collapse=", "))) 
   }
   ## check the bounding boxes:
   if(!(any(object@sd1@bbox %in% as.list(object@sd2@bbox, object@sd3@bbox, object@sd4@bbox, object@sd5@bbox, object@sd6@bbox)))){
      return("The bounding box of all 'sd' slots is not standard") 
   }
})


## geosamples class:
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
   ## check if URI exists:
   uri = paste(paste(object@server$URI, "?", sep=""), object@server$version, object@server$service, "request=GetCapabilities", sep="&") 
   require(RCurl)
   try(z <- getURI(uri, .opts=curlOptions(header=TRUE, nobody=TRUE, transfertext=TRUE, failonerror=FALSE)))
   if(!length(x <- grep(z, pattern="404 Not Found"))==0)
      return("Server error: 404 Not Found")
})

## REST class
setClass("REST.SoilGrids", representation (server = 'character', query = 'list', stream = 'list'),
   prototype = list(server=get("REST.server", envir = GSIF.opts), query=list(attributes=get("attributes", envir = GSIF.opts), confidence=c("U","M","L"), depths=c("sd1","sd2","sd3","sd4","sd5","sd6")), stream=list(clipList=NA, param=NA)), ## TH: Might change in future!
   validity = function(object) {
   ## check if URI exists:
   require(RCurl)
   try(z <- getURI(object@server, .opts=curlOptions(header=TRUE, nobody=TRUE, transfertext=TRUE, failonerror=FALSE)))
   if(!length(x <- grep(z, pattern="404 Not Found"))==0){
      return("Server error: 404 Not Found")
   } 
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
   if(!any(names(object@mu) %in% levels(object@predicted@data[,1])))
      return("Class names in the 'predicted' and 'mu' slots do not match")
   # check if the row names in the class.sd, class.c match:
   if(!all(row.names(object@class.c) %in% levels(object@predicted@data[,1])))
      return("Row names in the 'class.c' slot and 'predicted' slots do not match")
   if(!all(row.names(object@class.sd) %in% levels(object@predicted@data[,1])))
      return("Row names in the 'class.sd' slot and 'predicted' slots do not match")
   if(ncol(object@mu@data)<2)
      return("A minimum of two membership maps required")   
   # check if all mu's sum to 1 (plus minus 1%):
   if(!all(rowSums(object@mu@data, na.rm=TRUE)>.99&rowSums(object@mu@data, na.rm=TRUE)<1.01))
      return("Some rows in the 'mu' slot do not sum up to 1")
   # check if the confusion matrix has kappa > 0
   if(length(object@confusion)==0|attr(object@confusion, "error")==0)
      return("Not possible to derive confusion table or no significant match detected")
})



################## generic functions ##############

if(!isClass("ppp")){
  setClass("ppp")
}

if(!isGeneric("getID")){
  setGeneric("getID", function(obj, ...){standardGeneric("getID")})
}

if(!isGeneric("as.data.frame")){
  setGeneric("as.data.frame", function(x, row.names = NULL, optional = FALSE, ...){standardGeneric("as.data.frame")})
}

if(!isGeneric("predict")){
  setGeneric("predict", function(object, ...){standardGeneric("predict")})
}

if(!isGeneric("over")){
  setGeneric("over", function(x, y, ...){standardGeneric("over")})
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

if(!isGeneric("getSpatialTiles")){
  setGeneric("getSpatialTiles", function(obj, ...){standardGeneric("getSpatialTiles")})
}

if(!isGeneric("tile")){
  setGeneric("tile", function(x, ...){standardGeneric("tile")})
}

if(!isGeneric("describe")){
  setGeneric("describe", function(x, ...){standardGeneric("describe")})
}

if(!isGeneric("summary")){
  setGeneric("summary", function(object, ...){standardGeneric("summary")})
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

if (!isGeneric("test.gstatModel")){
  setGeneric("test.gstatModel", function(observations, formulaString, covariates, ...){standardGeneric("test.gstatModel")})
}

if (!isGeneric("fit.regModel")){
  setGeneric("fit.regModel", function(formulaString, rmatrix, predictionDomain, method, ...){standardGeneric("fit.regModel")})
}

if (!isGeneric("fit.vgmModel")){
  setGeneric("fit.vgmModel", function(formulaString, rmatrix, predictionDomain, ...){standardGeneric("fit.vgmModel")})
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

if (!isGeneric("warp")){
  setGeneric("warp", function(obj, ...){standardGeneric("warp")})
}

if (!isGeneric("MaxEnt")){
  setGeneric("MaxEnt", function(occurrences, covariates, ...){standardGeneric("MaxEnt")})
}

# end of script;