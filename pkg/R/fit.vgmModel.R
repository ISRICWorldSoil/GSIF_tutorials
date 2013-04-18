# Purpose        : Fitting 2D or 3D variograms;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Gerard Heuvelink (gerard.heuvelink@wur.nl); 
# Dev Status     : Pre-Alpha
# Note           : The variogram fitting in geoR is probably more robust, but also more time consuming;


## fit variogram to a 2D or 3D point object:
setMethod("fit.vgmModel", signature(formulaString = "formula", rmatrix = "data.frame", predictionDomain = "SpatialPixelsDataFrame"), function(formulaString, rmatrix, predictionDomain, vgmFun = "Exp", dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]], anis = NULL, subsample = nrow(rmatrix), ivgm, ...){

  ## check input object:
  if(is.na(proj4string(predictionDomain))){ stop("proj4 string required for argument 'predictionDomain'") }

  ## target variable name:
  if(!any(names(rmatrix) %in% all.vars(formulaString))){
    stop("Variables in the 'formulaString' not found in the 'rmatrix' object.")
  }
  
  ## remove missing observations:
  tv = all.vars(formulaString)[1]
  rmatrix <- rmatrix[complete.cases(lapply(all.vars(formulaString), function(x){rmatrix[,x]})),]

  ## spatial coordinates (column names):
  xyn = attr(predictionDomain@bbox, "dimnames")[[1]]
  if(!any(names(rmatrix) %in% xyn)){
       stop(paste("Column names:", paste(xyn[which(!(xyn %in% names(rmatrix)))], collapse=", "), "could not be located in the regression matrix"))
  }
  
  ## add 3D dimension if missing:
  if(dimensions=="3D" & length(xyn)==2){ xyn = c(xyn, "altitude") }  
  
  ## create spatial points:
  coordinates(rmatrix) <- as.formula(paste("~", paste(xyn, collapse = "+"), sep=""))
  proj4string(rmatrix) = predictionDomain@proj4string
  observations = as(rmatrix, "SpatialPoints")

  ## model does not have to be fitted?
  if(vgmFun == "Nug"){
    rvgm <- vgm(nugget=var(rmatrix@data[,tv]), model=vgmFun, range=0, psill=var(rmatrix@data[,tv]))
  } else {

  ## subset if necessary to speed up the computing:
  if(subsample < nrow(rmatrix)){
    pcnt <- subsample/nrow(rmatrix)
    message(paste("Subsetting observations to", signif(pcnt*100, 3), "percent"))
    rmatrix <- rmatrix[runif(nrow(rmatrix))<pcnt,]
  }
   
  ## guess the dimensions:
  if(missing(dimensions)){
      xyn = attr(rmatrix@bbox, "dimnames")[[1]]
      if(length(xyn)==2) {
         dimensions = "2D" 
      } else {
         dimensions = "3D"
      }
  }
  
  if(dimensions == "3D"){
    ## estimate area extent:
    Range = sqrt(areaSpatialGrid(predictionDomain))/3
  
    ## estimate anisotropy:
    if(is.null(anis)){ 
      ## estimate initial range in the vertical direction:
      dr <- abs(diff(range(rmatrix@coords[,3], na.rm=TRUE)))/3
      a2 = 2*dr/Range
      ## estimate anisotropy parameters:
      anis = c(0, 0, 0, 1, a2)
    }
  }
  
  if(dimensions == "2D"){
    ## check if it is projected object:
    if(!is.na(proj4string(predictionDomain))){
      if(!is.projected(predictionDomain)){
        require(fossil)  # Haversine Formula for Great Circle distance
        p.1 <- matrix(c(predictionDomain@bbox[1,1], predictionDomain@bbox[1,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
        p.2 <- matrix(c(predictionDomain@bbox[2,1], predictionDomain@bbox[2,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
        Range = fossil::deg.dist(lat1=p.1[,2], long1=p.1[,1], lat2=p.2[,2], long2=p.2[,1])/2
      } else {
        Range = sqrt(areaSpatialGrid(predictionDomain))/2      
      }
    } else {
      Range = sqrt(areaSpatialGrid(predictionDomain))/2
    }
    ## estimate anisotropy parameters:
    anis = c(0, 1)
  }
  
  ## initial variogram:    
  if(missing(ivgm)){
    if(dimensions == "2D"|dimensions == "3D"){
      ivgm <- vgm(nugget=0, model=vgmFun, range=Range, psill=var(rmatrix@data[,tv]), anis = anis)
    }
  }
  ## TH: 2D+T and 3D+T variogram fitting will be added;

  ## try to fit a variogram:
  try(rvgm <- gstat::fit.variogram(variogram(formulaString, rmatrix), model=ivgm, ...))   
  if(class(.Last.value)[1]=="try-error"){ 
    stop("Variogram model could not be fitted.") 
  }
    
  if(any(!(names(rvgm) %in% c("range", "psill"))) & diff(rvgm$range)==0|diff(rvgm$psill)==0){
    warning("Variogram shows no spatial dependence")     
  }
  
  }
  
  return(list(vgm=rvgm, observations=observations))
  
})

# end of script;