# Purpose        : Run MaxEnt and produce outputs;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : Not recommended for large grids;


setMethod("MaxEnt", signature(occurrences = "ppp", covariates = "SpatialPixelsDataFrame"), function(occurrences, covariates, nfold = 5, Npoints = 1000, sciname = as.character(NA), period = c(Sys.Date()-1, Sys.Date()),  ...){
  
  require(dismo)
  require(raster)
  require(plotKML)
  
  # only run if the maxent.jar file is available, in the right folder
  jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
  if (file.exists(jar)) {
    prj <- covariates@proj4string
    sel <- names(covariates)[sapply(covariates@data, is.factor)]
    covariates <- stack(covariates)
    # prepare the occurrence-only records:
    xy <- data.frame(occurrences)
    # fit a MaxEnt model (can take few minutes!):
    me <- maxent(covariates, xy, factors=sel) # takes 1-2 mins;
    # predict distribution of Big Foot:
    pr <- predict(me, covariates, ...)
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
  } else {
    paste("Maxent software could not be located. See 'dismo::maxent' for more info.")
  }

  # create an object of type "SpatialMaxEntOutput":      
  out <- new("SpatialMaxEntOutput", sciname = sciname, occurrences = xy, TimeSpan.begin = as.POSIXct(period[1]), TimeSpan.end = as.POSIXct(period[2]), maxent = me, sp.domain = hr, predicted = pr)
  return(out)
  
})

# end of script;