# Purpose        : Derive Spatial Predictive Components for a list of grids;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : Not recommended for large grids;


setMethod("spc", signature(obj = "SpatialPixelsDataFrame", formulaString = "formula"), function(obj, formulaString, scale.=TRUE, silent = FALSE, ...){
 
  require(raster)

  # formula string:
  if(missing(formulaString)) {
     formulaString <- as.formula(paste("~", paste(out@layernames, collapse="+")))
  }
  vars = all.vars(formulaString)
  if(length(vars)< 2){
    stop("At least two covarites required to run Principal Component Analysis")
  }

  # print warning:
  if(nrow(obj)>10e6){
    warning('Operation not recommended for large grids', immediate. = TRUE)
  }
  
  # convert every factor to indicators:
  for(j in 1:length(vars)){
    if(is.factor(obj@data[,vars[j]])){
      ln <- levels(obj@data[,vars[j]])
      for(k in 1:length(ln)){
        vn <- paste(vars[j], k, sep="_")
        obj@data[,vn] <- ifelse(obj@data[,vars[j]]==ln[k], 1, 0)
      }
    message(paste("Converting", vars[j], "to indicators..."))
    }
  } 
  varsn = names(obj)[which(!sapply(obj@data, is.factor))]
  
  out <- brick(obj[varsn])
  # filter the missing values:
  x <- scale(getValues(out)) 
  x[is.na(x)] <- 0 
  
  pcs <- prcomp(formula=formulaString, scale=TRUE, as.data.frame(x))

  # copy values:  
  out@data@values <- pcs$x
  out@layernames <- attr(pcs$x, "dimnames")[2][[1]]
  out <- as(out, "SpatialPixelsDataFrame")
  proj4string(out) <- obj@proj4string
  if(silent==FALSE){
    message(paste("Converting covariates to principal components..."))
    summary(pcs)
  }
 
  pcs <- new("SpatialComponents", predicted = out, pca = pcs[-which(names(pcs)=="x")])
  return(pcs)

}) 

## Extract regression matrix:
setMethod("extract", signature(x = "SpatialComponents", y = "geosamples"), function(x, y, methodid, var.type = "numeric", ...){
  require(raster)
  require(plotKML)
  
  if(!any(y@data$altitudeMode == "relativeToGround")){
    warning("Only 'relativeToGround' values for AltitudeMode expected")
  }
  
  pnts = subset(y, methodid)
  coordinates(pnts) <- ~longitude+latitude
  proj4string(pnts) <- get("ref_CRS", envir = plotKML.opts) 
  pnts <- spTransform(pnts, x@sp@proj4string)
  ov <- extract(brick(x@predicted), pnts)
  out <- cbind(data.frame(pnts), data.frame(ov))
  
  # reformat observed values:
  if(var.type=="numeric"){
    out$observedValue = as.numeric(out$observedValue)
  } 
  else { 
      if(var.type=="factor"){
      out$observedValue = as.factor(out$observedValue)
      }
  }
  
  return(out)
  
})

# end of script;