# Purpose        : Functions for accessing SoilGrids;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Jorge S. Mendes de Jesus (Jorge.mendesdejesus@wur.nl);
# Status         : tested
# Note           : Location of the server and parameter names might change!


REST.SoilGrids <- function(attributes, depths=paste("sd",1:6,sep=""), confidence=c("L","M","U"), validate=FALSE){
  if(requireNamespace("rjson", quietly = TRUE)){
    if(validate==TRUE){
    ## get the most recent description:
      try( ret <- rjson::fromJSON(file=paste(get("REST.server", envir = GSIF.opts), "query/describe", sep="")), silent = TRUE )
      if(!class(.Last.value)[1]=="try-error" & !length(ret$query)==0){
        if(any(!attributes %in% ret$query$attributes)){
          stop(paste("Requested 'attributes' not present. See '", get("REST.server", envir = GSIF.opts),"' for more info.", sep=""))
        }
        if(any(!depths %in% ret$query$depths)){
          stop(paste("Requested 'depths' not present. See '", get("REST.server", envir = GSIF.opts), "' for more info.", sep=""))
        }
        if(any(!confidence %in% ret$query$confidence)){
          stop(paste("Requested 'confidence' not present. See '", get("REST.server", envir = GSIF.opts),"' for more info.", sep=""))
        }        
      }
    }
    out <- new("REST.SoilGrids", server=get("REST.server", envir = GSIF.opts), query=list(attributes=attributes, confidence=confidence, depths=depths), stream=list(clipList=NA, param=NA))
    return(out)
  }
}


.REST.uri <- function(x, lon, lat){
  ## location and attributes:
  xy <- paste("lon=", lon, "&lat=", lat, sep="")
  at <- NULL
  dpt <- NULL
  cnf <- NULL
  if(!all(get("attributes", envir = GSIF.opts) %in% x@query$attributes)){
    at <- paste("&attributes=", paste(x@query$attributes, collapse=",", sep=""), sep="")
  }
  if(!length(x@query$depths == 6)){
    dpt <- paste("&depths=", paste(x@query$depths, collapse=",", sep=""), sep="")
  }
  if(!all(c("L","M","U") %in% x@query$confidence)){
    cnf <- paste("&confidence=", paste(x@query$confidence, collapse=",", sep=""), sep="")
  }
  uri <- paste(x@server, "query?", xy, at, dpt, cnf, sep="")
}

## overlay method:
setMethod("over", signature(x = "REST.SoilGrids", y = "SpatialPoints"), 
  function(x, y) 
  {

  if(requireNamespace("rjson", quietly = TRUE)){
    ## run point by point
    out <- NULL
    for(i in 1:nrow(y@coords)){
      uri <- .REST.uri(x, lon=y@coords[i,1], lat=y@coords[i,2])
      try( ret <- rjson::fromJSON(file=uri), silent = TRUE )
      if(!class(.Last.value)[1]=="try-error" & !length(ret$properties)==0){
        out[[i]] <- data.frame(ret$properties)
        out[[i]]$lon <- y@coords[i,1]
        out[[i]]$lat <- y@coords[i,2]
      } else {    
        out[[i]] <- data.frame(lon=y@coords[i,1], lat=y@coords[i,2])
      }
    }
    ## bind all elements together:
    out <- plyr::rbind.fill(out)
    return(out)
  }
})

## end of script;