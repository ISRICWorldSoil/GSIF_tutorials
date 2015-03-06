# Purpose        : Fitting 2D or 3D variograms;
# Maintainer     : Bas Kempen (bas.kempen@wur.nl);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl) and Gerard Heuvelink (gerard.heuvelink@wur.nl); 
# Dev Status     : Pre-Alpha
# Note           : The variogram fitting in geoR is probably more robust, but also more time consuming;


## fit variogram to a 2D or 3D point object:
setMethod("fit.vgmModel", signature(formulaString = "formula", rmatrix = "data.frame", predictionDomain = "SpatialPixelsDataFrame"), function(formulaString, rmatrix, predictionDomain, vgmFun = "Exp", dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]], anis = NULL, subsample = nrow(rmatrix), ivgm, cutoff, width = cutoff/15, cressie = FALSE, ...){

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
  proj4string(rmatrix) <- predictionDomain@proj4string
  observations <- rmatrix
  
  ## subset to speed up the computing:
  if(subsample < nrow(rmatrix)){
    pcnt <- subsample/nrow(rmatrix)
    message(paste("Subsetting observations to", signif(pcnt*100, 3), "percent"))
    rmatrix <- rmatrix[runif(nrow(rmatrix))<pcnt,]
  }

  ## model does not have to be fitted?
  if(vgmFun == "Nug"){
    rvgm <- gstat::vgm(nugget=var(rmatrix@data[,tv]), model=vgmFun, range=0, psill=var(rmatrix@data[,tv]))
    svgm <- gstat::variogram(formulaString, rmatrix)
  } else {
     
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
      if(missing(cutoff)) { cutoff = Range }
    }
    
    if(dimensions == "2D"){
      ## check if it is projected object:
      if(!is.na(proj4string(predictionDomain))){
        if(!is.projected(predictionDomain)){
          if(requireNamespace("fossil", quietly = TRUE)){
            ## Haversine Formula for Great Circle distance
            p.1 <- matrix(c(predictionDomain@bbox[1,1], predictionDomain@bbox[1,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
            p.2 <- matrix(c(predictionDomain@bbox[2,1], predictionDomain@bbox[2,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
            Range = fossil::deg.dist(lat1=p.1[,2], long1=p.1[,1], lat2=p.2[,2], long2=p.2[,1])/2
          }
        } else {
          
          ## BK: To avoid problems in variogram fitting (singularity warnings) the default range value should be chosen in proportion to the cut-off value.
          ## BK: To determine the cut-off value of the sample variogram use the default gstat value.
          ## BK: Edzer -> "gstat uses one third of the diagonal of the rectangular (or block for 3D) that spans the data locations."  
          ## BK: The initial value for the range can then be chosen as one-third of the cut-off value
          
          ## compute one third of the diagonal of the bbox of the point locations using Pythagoras theorem
          dbbox <- round((sqrt((rmatrix@bbox[1,2]-rmatrix@bbox[1,1])**2+(rmatrix@bbox[2,2]-rmatrix@bbox[2,1])**2))/3,0)
          Range = dbbox/3
          #Range = sqrt(areaSpatialGrid(predictionDomain))/2      
        }
      } else {
        dbbox <- round((sqrt((rmatrix@bbox[1,2]-rmatrix@bbox[1,1])**2+(rmatrix@bbox[2,2]-rmatrix@bbox[2,1])**2))/3,0)
        Range = dbbox/3
        #Range = sqrt(areaSpatialGrid(predictionDomain))/2
      }
      
      ## estimate anisotropy parameters:
      anis = c(0, 1)
      
      ## BK: if not given, determine sample variogram cutoff value based on gstat default
      if(missing(cutoff)) {
        dbbox <- round((sqrt((rmatrix@bbox[1,2]-rmatrix@bbox[1,1])**2+(rmatrix@bbox[2,2]-rmatrix@bbox[2,1])**2))/3,0)
        cutoff = dbbox
      }
      
    }
    
    ## initial variogram:    
    if(missing(ivgm)){
      if(dimensions == "2D"|dimensions == "3D"){
        ## BK: since most variograms have a nugget effect, set the initial nugget equal to the initial partial sill:
        ivgm <- gstat::vgm(nugget=var(rmatrix@data[,tv])/5, model=vgmFun, range=Range, psill=var(rmatrix@data[,tv])*4/5, anis = anis)
        #ivgm <- vgm(nugget=0, model=vgmFun, range=Range, psill=var(rmatrix@data[,tv]), anis = anis)
      }
    }
    ## TH: 2D+T and 3D+T variogram fitting will be added here;
    
    ## fit sample variogram 
    try( svgm <- gstat::variogram(formulaString, rmatrix, cutoff=cutoff, width=width, cressie=cressie) )
      
    ## try to fit a variogram using default settings:
    try( rvgm <- gstat::fit.variogram(svgm, model=ivgm, ...) )     
    ## BK: the code below does not work for 'singular model' warnings, only when the function does not fit at all
    if(class(.Last.value)[1]=="try-error"){
      try( rvgm <- gstat::fit.variogram(gstat::variogram(formulaString, rmatrix, cutoff=cutoff), model=ivgm, fit.ranges = FALSE, ...) )
      if(class(.Last.value)[1]=="try-error"){    
        warning("Variogram model could not be fitted.") 
      } 
    }
      
    if(any(!(names(rvgm) %in% c("range", "psill"))) & (diff(rvgm$range)==0|diff(rvgm$psill)==0)){
      warning("Variogram shows no spatial dependence")
      rvgm <- NULL     
    }
    
  }
  
  ## BK: sample variogram is now stored as well
  return(list(vgm=rvgm, observations=observations, svgm=svgm))
  
})

# end of script;