# Purpose        : Combines a list of spatial pixels (multi-source data merging);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : This is still a rather heuristic approach (a better option would be to use some more robust data assimilation method); weights can also be passed manually;

## merge multiple spatial predictions:
setMethod("merge", signature(x = "SpatialPredictions", y = "SpatialPredictions"), function(x, y, ..., RMSE.l = NULL, silent = TRUE){

  # check inputs:
  r <- list(x, y, ...)
	if (length(r) < 2) {
		  stop("merge needs at least 2 'SpatialPredictions' objects")
	}
	
  # target variable
  variables = sapply(r, FUN=function(x){x@variable})
	# check if the names are consistent:
  if(!length(unique(variables))==1){
	   stop("Merging of objects of class 'SpatialPredictions' requires idential 'variable' slot names")
  }
    
  # variable names:
  cname <- c(unique(variables), "var1.var")
	sel.t <- paste(cname[1], 1:length(r), sep="_")
  sel.s <- paste(cname[2], 1:length(r), sep="_")
	
  # estimate the weights:
	if(is.null(RMSE.l)){
     RMSE.l <- sapply(r, FUN=function(x){var(x@validation$residual, na.rm=T)})
     names(RMSE.l) <- sel.t 
	}
	
  # copy target columns:
  r <- sapply(r, FUN=function(x){slot(x, "predicted")[cname]})
 
  # rename methods for consistency:
  for(j in 1:length(r)){ names(r[[j]]) <- paste(cname, j, sep="_") }
  # copy grid properties:
  cd <- data.frame(t(sapply(r, FUN=function(x){x@grid@cellsize})))
  cs <- data.frame(t(sapply(r, FUN=function(x){x@grid@cells.dim})))

  # resample all grids to the finest resolution:  
  if(all(names(cd)==c("longitude", "latitude", "altitude"))){
    # 3D data:
    stdepth <- cd[1,"altitude"]
    stsize <- cs[1,"altitude"]
    out <- sp3D(r, stdepths = stdepth, stsize = stsize)[[1]]
  
  } else {
    # 2D data:
    if(length(names(cd))==2){ 
    # pick up the most detailed scale:
    cellsize.l <- sapply(r, FUN=function(x){x@grid@cellsize[1]})
    tc <- which(cellsize.l == min(cellsize.l))
    out <- r[[tc[1]]]
    fullgrid(out) <- TRUE
    r[[tc[1]]] <- NULL
    ret <- list(NULL)
    for(j in 1:length(r)){
      # resample all the grids to the finest resolution:
      if(cellsize.l[j] > min(cellsize.l)|!identical(out@bbox, obj[[j]]@bbox)){
         ret[[j]] <- warp(r[[j]], proj4s = proj4string(out), pixsize = min(cellsize.l), GridTopology = out@grid, resampling_method = "cubicspline")
      }
      else {
        ret[[j]] <- r[[j]]
      } 
    }
    ret <- lapply(ret, FUN=function(x){slot(x, "data")})
    out@data <- cbind(out@data, do.call(cbind, ret))
    out <- as(out, "SpatialPixelsDataFrame")
    
    } else {
      stop("2 or 3 dimensional object of class 'SpatialPixelsDataFrame' expected")    
    }    
  } 
	
	# scale the kriging variances using the results of cross-validation:
	out@data[,sel.s] <- RMSE.l / colMeans(out@data[,sel.s]) * out@data[,sel.s]
	
  # derive a weighted mean - we assume that the correlation between the variances (errors) is = 0;
	## GH: See Deutsch, 1965, section 5.3; Searle, 1971, page 89;
	out@data[,cname[1]] <- rowSums(out@data[,sel.t]/out@data[,sel.s]) / rowSums(1/out@data[,sel.s])
	
	if(silent==FALSE){ 
    message(paste("Cross-validation RMSE (type = link):"))
    print(signif(unlist(RMSE.l), 4))
  }
  
	return(out[cname[1]])
})



## merge realizations at multiple scales:
# setMethod("merge", signature(x = "RasterBrickSimulations", y = "RasterBrickSimulations"), function(x, y, ..., RMSE.l = NULL){
#  # weights need to be provided by the user!
#
#  # check inputs:
#  r <- list(x, y, ...)
#	if (length(r) < 2) {
#		  stop("merge needs at least 2 'RasterBrickSimulations' objects")
#	}
#
#  # resample all grids to the finest resolution:  
#
# })

# end of script;