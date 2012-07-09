# Purpose        : Wrapper functions for WPS worldgrids;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : Hannes Reuter (hannes.reuter@wur.nl);
# Status         : tested
# Note           : Location of the server and parameter names might change!;

## get capacities:
setMethod("show", signature(object = "WPS"), function(object){
  uri = paste(paste(object@server$URI, "?", sep=""), paste(object@server$service, object@server$version, "request=GetCapabilities", sep="&"), sep="")
  ret <- xmlTreeParse(uri, useInternalNodes = TRUE)
  ret <- unlist(xmlToList(ret, addAttributes=FALSE))
  # convert to a table:
  ret <- data.frame(field=gsub("\\.", "_", attr(ret, "names")), value=paste(ret), stringsAsFactors = FALSE)
    cat("Raster layers", ":", object@inRastername, "\n")
  for(i in 1:nrow(ret)){
    cat(ret$field[i], ":", ret$value[i], "\n")
  }
})

## get processes:
setMethod("getProcess", signature(x = "WPS"), function(x){
  uri = paste(paste(x@server$URI, "?", sep=""), paste(x@server$service, x@server$version, "request=GetCapabilities", sep="&"), sep="")
  ret <- xmlTreeParse(uri, useInternalNodes = TRUE)
  ret <- unlist(xmlToList(ret, addAttributes=FALSE))
  # convert to a table:
  ret <- data.frame(field=gsub("\\.", "_", attr(ret, "names")), value=paste(ret), stringsAsFactors = FALSE)
  proc <- ret[grep(ret$field, pattern="Process1"),"value"]
  attr(proc, "names") <- ret[grep(ret$field, pattern="Process2"),"value"]
  return(proc)  
})

## get arguments:
setMethod("describe", signature(x = "WPS"), function(x, request = "describeprocess", identifier){
  uri = paste(paste(x@server$URI, "?", sep=""), paste(x@server$service, x@server$version, paste("request=", request, sep=""), paste("identifier=", identifier, sep=""), sep="&"), sep="")
  ret <- xmlTreeParse(uri, useInternalNodes = TRUE)
  ret <- unlist(xmlToList(ret, addAttributes=FALSE))
  # convert to a table:
  ret <- data.frame(field=gsub("\\.", "_", attr(ret, "names")), value=paste(ret), stringsAsFactors = FALSE)
  return(ret)  
})

## overlay points (single point overlay)
setMethod("overlay", signature(x = "WPS", y = "SpatialPoints"), 
  function(x, y) 
  {

  require(RCurl)
  require(sp)
  # point by point
  out <- NULL
  for(i in 1:length(nrow(y@coords))){
    ret <- paste("[x=", y@coords[i,1], ";y=", y@coords[i,2], ";inRastername=", x@inRastername,"]", sep="")
    uri = paste(paste(x@server$URI, "?", sep=""), paste(x@server$service, x@server$version, "request=execute", "identifier=sampler_local1pt_nogml", paste("datainputs=", ret, sep=""), sep="&"), sep="")
    ret <- xmlTreeParse(getURL(uri), useInternalNodes = TRUE)
    nx <- unlist(xmlToList(ret, addAttributes=FALSE))
      if(any(nx %in% "OutData")){
      out[i] <- xmlValue(ret[["//wps:LiteralData"]])
      }
      else {
      out[i] <- NA
      }
  }

  return(out)
})

## subset grid
setMethod("subset", signature(x = "WPS"), 
  function(x, bbox, import = TRUE) 
  {

  require(raster)
  require(RCurl)
  require(rgdal)
  # check that bbox is fine:
  if(nrow(bbox)==2&ncol(bbox)==2&bbox[1,2]<180&bbox[2,2]<90&bbox[1,2]>bbox[1,1]&bbox[1,1]>-180&bbox[2,1]>-90&bbox[2,2]>bbox[2,1]){
    ret <- paste('[bbox=', bbox[1,1], ',', bbox[2,1], ',', bbox[1,2], ',', bbox[2,2], ';inRastername=', x@inRastername, ']', sep="")
    uri = paste(paste(x@server$URI, "?", sep=""), paste(x@server$service, x@server$version, "request=execute", "identifier=subset", paste("datainputs=", ret, sep=""), "responsedocument=OutData=@asreference=true", sep="&"), sep="")
    ret <- xmlTreeParse(getURL(uri), useInternalNodes = TRUE)  
    nx <- unlist(xmlToList(ret, addAttributes=FALSE))
    if(any(nx %in% "Output Subset data")){
       objectname <- xmlAttrs(ret[["//wps:Reference"]], FALSE, FALSE)
       objectname <- objectname[names(objectname)=="href"]
       out <- set.file.extension(paste(x@inRastername, paste(as.vector(bbox), collapse="_"), sep="_"), ".tif")
       try(download.file(objectname, destfile=out, mode="wb"))
      }
      else { stop("Layer not avaialable. Visit www.worldgrids.org for more info.") }  
  }
  else stop("Bounding box required as matrix with LatMin, LatMax, LonMin, LonMax")

  if(import == TRUE){
   out <- flipVertical(readGDAL(out))
   names(out) <- x@inRastername
   return(out)
  }
  else {
    GDALinfo(out)    
  }
})

# end of the script;