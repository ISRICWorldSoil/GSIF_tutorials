# Purpose        : Run MaxEnt and produce outputs;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : Not recommended for large grids;


## Wrapper function for MaxEnt:
setMethod("MaxEnt", signature(occurrences = "ppp", covariates = "SpatialPixelsDataFrame"), function(occurrences, covariates, nfold = 5, Npoints = 1000, sciname = as.character(NA), period = c(Sys.Date()-1, Sys.Date()),  ...){
  
  require(dismo)
  require(raster)
  # require(plotKML)
  # require(maptools)
  
  # only run if the maxent.jar file is available, in the right folder
  jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
  if (file.exists(jar)) {
    prj <- covariates@proj4string
    sel <- names(covariates)[sapply(covariates@data, is.factor)]
    covariates <- stack(covariates)
    # prepare the occurrence-only records:
    xy <- data.frame(occurrences)[,1:2]
    # fit a MaxEnt model (can take few minutes!):
    if(length(sel)==0){
      me <- maxent(covariates, xy, ...) 
    } else {
      me <- maxent(covariates, xy, factors=sel, ...)     
    }

    # Load rJava
  	if (is.null(getOption('dismo_rJavaLoaded'))) {
		if ( require(rJava) ) {
			.jpackage('dismo')
			options(dismo_rJavaLoaded=TRUE)
		  } else {
			stop('rJava cannot be loaded')
		  }
	  }

    # predict distribution:
    pr <- predict(object=me, x=covariates, ...)  # predict.MaxEnt
    # run cross-validation
    fold <- kfold(xy, k=nfold)
    # randomly take 20% of observations:
    xy.test <- xy[fold == 1,]
    bgp <- randomPoints(covariates, Npoints)  
    ev <- evaluate(me, p=xy.test, a=bgp, x=covariates)
    # this allows estimation of the threshold probability:
    threshold <- ev@t[which.max(ev@TPR + ev@TNR)]
    # prepare data for plotKML:
    hr <- as(calc(pr, fun=function(x){ifelse(x>threshold, 1, NA)}), "SpatialPixelsDataFrame")
    xy <- as(occurrences, "SpatialPoints")
    proj4string(xy) = prj
    # create an object of type "SpatialMaxEntOutput":      
    out <- new("SpatialMaxEntOutput", sciname = sciname, occurrences = xy, TimeSpan.begin = as.POSIXct(period[1]), TimeSpan.end = as.POSIXct(period[2]), maxent = me, sp.domain = hr, predicted = pr)
    return(out)
    
  } else {
    paste("Maxent software could not be located. See 'dismo::maxent' for more info.")
  }
  
})



# end of script;