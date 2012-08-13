# Purpose        : Converts a SoilProfileCollection to loose records (KML placemarks);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : Hannes I. Reuter;
# Status         : pre-alpha
# Note           : see also "as.data.frame" operation;

## coerce SoilProfileCollection to "geosamples":
setMethod("as.geosamples", signature(obj = "SoilProfileCollection"), 
  function(obj, registry = as.character(NA), sample.area = 1, mxd = 2, dtime = 3600) 
  {

  # reproject if necessary:
  require(aqp)
  require(plotKML)
  if(!check_projection(obj@sp)){
     obj@sp <- reproject(obj@sp)
  }

  # check for duplicates:
  dp <- duplicated(obj@site[,obj@idcol])
  if(sum(dp)>0){
    warning(paste("Duplicated IDs detected in the 'site' slot and will be removed:", paste(obj@site[dp, obj@idcol], collapse=", ", sep="")))
    obj@site <- obj@site[!dp,]
    obj@sp <- obj@sp[!dp,]
  } 

  # estimate thickness in m and depths:
  sampleThickness <- abs(obj@horizons[,obj@depthcols[2]] - obj@horizons[,obj@depthcols[1]])/100 
  depths <- - (obj@horizons[,obj@depthcols[1]] + (obj@horizons[,obj@depthcols[2]] - obj@horizons[,obj@depthcols[1]])/2)/100

  # add the time coordinate if missing:
  if(ncol(obj@sp@coords)==2){
    XYT <- data.frame(cbind(obj@sp@coords, time=rep(NA, nrow(obj@sp@coords))))
  } else {
    XYT <- data.frame(obj@sp@coords)
  }
  names(XYT)[1:2] <- c("x", "y")
  XYT$ID <- profile_id(obj)
  
  # convert site data to geosamples:
  x <- NULL
  site <- obj@site
  # remove columns that are not of interest:
  site[,obj@idcol] <- NULL
  # add the location error:
  locationError = attr(obj@sp@coords, "locationError")
  if(is.null(locationError)) { locationError = rep(as.character(NA), nrow(obj@sp@coords)) }
  XYT$locationError <- locationError
  
  # for each soil variable
  for(j in 1:length(names(site))){
    ll <- length(site[,names(site)[j]])
    observationid = attr(site[,names(site)[j]], "IGSN")
    if(is.null(observationid)) { observationid = rep(as.character(NA), ll) } 
    measurementError = attr(site[,names(site)[j]], "measurementError")
    if(is.null(measurementError)) { measurementError = rep(as.character(NA), ll) }
    sampleArea = attr(site[,names(site)[j]], "sampleArea")
    if(is.null(sampleArea)) { sampleArea = rep(sample.area, ll) }    
    x[[j]] <- data.frame(observationid = as.character(observationid), sampleid = profile_id(obj), longitude = XYT[,1], latitude = XYT[,2], locationError = as.numeric(locationError), TimeSpan.begin = as.POSIXct(XYT[,3]-dtime/2, origin="1970-01-01"), TimeSpan.end = as.POSIXct(XYT[,3]+dtime/2, origin="1970-01-01"), altitude = as.numeric(rep(0, ll)), altitudeMode = rep("relativeToGround", ll), sampleArea = sampleArea, sampleThickness = rep(mxd*sample.area, ll), observedValue = as.character(site[,names(site)[j]]), methodid = rep(names(site)[j], ll), measurementError = as.numeric(measurementError), stringsAsFactors = FALSE) 
  }
  rx <- do.call(rbind, x)
    
  # convert horizon data to geosamples:
  y <- NULL
  hors <- obj@horizons
  # remove columns that are not of interest:
  hors[,obj@idcol] <- NULL
  hors[,obj@depthcols[1]] <- NULL
  hors[,obj@depthcols[2]] <- NULL
  # add coordinates:
  XYTh <- merge(data.frame(ID=obj@horizons[,obj@idcol], dtime=dtime), XYT, by="ID", all.x=TRUE)
  
  # for each soil variable
  for(j in 1:length(names(hors))){
    ll <- length(hors[,names(hors)[j]])
    observationid = attr(hors[,names(hors)[j]], "IGSN")
    if(is.null(observationid)) { observationid = rep(as.character(NA), ll) } 
    measurementError = attr(hors[,names(hors)[j]], "measurementError")
    if(is.null(measurementError)) { measurementError = rep(as.character(NA), ll) }
    sampleArea = attr(hors[,names(hors)[j]], "sampleArea")
    if(is.null(sampleArea)) { sampleArea = rep(sample.area, ll) }
    y[[j]] <- data.frame(observationid = as.character(observationid), sampleid = XYTh$ID, longitude = XYTh$x, latitude = XYTh$y, locationError = as.numeric(XYTh$locationError), TimeSpan.begin = as.POSIXct(XYTh$time-XYTh$dtime/2, origin="1970-01-01"), TimeSpan.end = as.POSIXct(XYTh$time+XYTh$dtime/2, origin="1970-01-01"), altitude = as.numeric(depths), altitudeMode = rep("relativeToGround", ll), sampleArea = sampleArea, sampleThickness = sampleThickness, observedValue = as.character(hors[,names(hors)[j]]), methodid = rep(names(hors)[j], ll), measurementError = as.numeric(measurementError), stringsAsFactors = FALSE) 
  }
  ry <- do.call(rbind, y)
  
  # merge the sites and horizons tables:
  tb <- rbind(rx, ry)
  tb$methodid <- as.factor(tb$methodid)

  # check if the metadata comply with the geosamples standard:
  mnames = c("methodid", "description", "units", "detectionLimit")
  if(any(!(names(obj@metadata) %in% mnames))){ 
    # warning(paste("Missing column names in the 'metadata' slot:", paste(mnames, collapse=", "))) 
    tnames <- levels(tb$methodid)
    metadata <- data.frame(tnames, rep(NA, length(tnames)), rep(NA, length(tnames)), rep(NA, length(tnames)))
    names(metadata) <- mnames
  } else {
    metadata = obj@metadata
  }
 
  # make geosamples:
  gs <- new("geosamples", registry = registry, methods = metadata, data = tb)
    
  return(gs)
})


## subsetting geosamples:
setMethod("subset", signature(x = "geosamples"), function(x, method){
  ret <- x@data[x@data$methodid==method,]
  if(nrow(ret)==0){ warning("Empty object. Methodid possibly not available") }
  attr(ret$methodid, "description") <- x@methods[x@methods$methodid==method,"description"]
  attr(ret$methodid, "units") <- x@methods[x@methods$methodid==method,"units"]
  attr(ret$methodid, "detectionLimit") <- x@methods[x@methods$methodid==method,"detectionLimit"]
  return(ret)  
})


## summary values:
setMethod("show", signature(object = "geosamples"), 
  function(object) 
  {
  cat("  Registry            :", object@registry, "\n")
  cat("  Variables           :", paste(levels(object@data$methodid), collapse=", "), "\n")
  cat("  Total samples       :", nrow(object@data), "\n")
  # create :
  sp <- object@data
  coordinates(sp) <- ~longitude+latitude
  proj4string(sp) <- get("ref_CRS", envir = GSIF.opts)
  cat("  Unique locations    :", length(unique(sp@coords))/2, "\n")
  cat("  Mean location error :", mean(object@data$locationError, na.rm=TRUE), "\n")
  cat("  Min longitude       :", range(object@data$longitude)[1], "\n")  
  cat("  Max longitude       :", range(object@data$longitude)[2], "\n")  
  cat("  Min latitude        :", range(object@data$latitude)[1], "\n")  
  cat("  Max latitude        :", range(object@data$latitude)[2], "\n")  
  # estimate the total area covered by the samples:
  sp.gc <- spTransform(sp, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"))
  Tarea <- signif(diff(sp.gc@bbox[1,])*diff(sp.gc@bbox[2,])/1e6, 4)
  cat("  Total area          :", paste(Tarea), "(square-km)", "\n")
})


## Extract regression matrix:
setMethod("overlay", signature(x = "SpatialPixelsDataFrame", y = "geosamples"), function(x, y, methodid, var.type = "numeric"){
  require(raster)
  require(plotKML)
  
  if(!any(y@data$altitudeMode == "relativeToGround")){
    warning("AltitudeMode accepts only 'relativeToGround' values")
  }
  
  pnts = subset.geosamples(y, method=methodid)
  # reformat observed values:
  if(var.type=="numeric"){
    pnts$observedValue = as.numeric(pnts$observedValue)
  } else { 
      if(var.type=="factor"){
      pnts$observedValue = as.factor(pnts$observedValue)
      }
  }
  
  coordinates(pnts) <- ~longitude+latitude
  proj4string(pnts) <- get("ref_CRS", envir = GSIF.opts) 
  
  index <- overlay(x, spTransform(pnts, x@proj4string))                    
  sel <- !is.na(index)
  out <- cbind(data.frame(pnts[sel,]), x[index[sel],])
  if(nrow(out)==0){ 
    warning("Overlay resulted in an empty table") 
  }
  
  return(out)
  
})

# end of script;