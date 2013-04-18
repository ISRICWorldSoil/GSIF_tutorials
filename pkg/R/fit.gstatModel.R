# Purpose        : Fit a 2D or 3D regression-kriging model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : linear models with normally distributed residuals;


## Fit a 'simple' 2D RK model:
setMethod("fit.gstatModel", signature(observations = "SpatialPointsDataFrame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "HB")[[1]], dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]], family = gaussian, stepwise = TRUE, vgmFun = "Exp", subsample = 5000, ...){
 
  ## TH: the function only works with 2D maps at the moment:
  if(length(attr(coordinates(observations), "dimnames")[[2]])>2){
    warning("This method uses only 2D coordinates of the points. For 3D data consider using the 'geosamples-class'.")
  }

  if(!missing(family) & !method=="GLM"){ warning("'family' argument will be ignored. Use 'method=\"GLM\"' instead.")  }

  ## overlay:
  ov <- over(observations, covariates)
  ## all variables of interest:
  tv = all.vars(formulaString)[1]
  seln = names(covariates) %in% all.vars(formulaString)[-1]
  ## check if all covariates are available: 
  if(length(seln)==0){
      stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  ov <- cbind(data.frame(observations[,tv]), ov)  
  
  ## copy coordinate column names for consistency:
  xyn <- which(names(ov) %in% attr(observations@bbox, "dimnames")[[1]])
  if(length(xyn)==2) { 
    names(ov)[xyn] <- attr(covariates@bbox, "dimnames")[[1]][1:2] 
    } else {
    names(ov)[xyn] <- attr(covariates@bbox, "dimnames")[[1]]     
  }

  # check the size of the output:
  if(nrow(ov)==0|is.null(ov[,tv])) {
    stop("The 'over' operations resulted in an empty set.")
  }
  
  # fit/filter the regression model:
  m <- fit.regModel(formulaString = formulaString, rmatrix = ov, predictionDomain = covariates[seln], method = method, dimensions = "2D", family = family, stepwise = stepwise, vgmFun = vgmFun, subsample = subsample, ...)  
    
  return(m)
  
})


## Fit a RK model using geosamples class:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "HB")[[1]], dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]], family = gaussian, stepwise = TRUE, vgmFun = "Exp", subsample = 5000, ...){
  
  ## all columns of interest:
  methodid = all.vars(formulaString)[1]
  seln = names(covariates) %in% all.vars(formulaString)[-1]
  xyn = attr(covariates@bbox, "dimnames")[[1]]
  
  ## prepare regression matrix:
  ov <- over(x=covariates, y=observations, method=methodid, var.type = "numeric")
  if(nrow(ov)==0|is.null(ov$observedValue)) {
    warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
  }
  ## geostats only possible with numeric variables:
  ov[,methodid] = as.numeric(ov$observedValue)
  ov$observedValue = NULL
   
  ## subset to columns of interest:  
  if(length(xyn)==2){ ov <- ov[, c(methodid, names(covariates)[seln], xyn, "altitude")] }
  if(length(xyn)==3){ ov <- ov[, c(methodid, names(covariates)[seln], xyn)] }
  
  ## fit/filter the regression model:
  m <- fit.regModel(formulaString = formulaString, rmatrix = ov, predictionDomain = covariates[seln], method = method, dimensions = "3D", family = family, stepwise = stepwise, vgmFun = vgmFun, subsample = subsample, ...)
  
  return(m)  

})


## Fit a RK model to a list of covariates / formulas:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "list", covariates = "list"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "HB")[[1]], dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]], family = gaussian, stepwise = TRUE, vgmFun = "Exp", subsample = 5000, ...){
    if(!length(formulaString)==length(covariates)){
      stop("'formulaString' and 'covariates' lists of same size expected")
    }
    
    rkm.l <- list(NULL)
    for(l in 1:length(covariates)){   
      rkm.l[[l]] <- fit.gstatModel(observations, formulaString[[l]], covariates[[l]], methodid = methodid, family = family, stepwise = stepwise, subsample = subsample, ...)
    }
    return(rkm.l)  

})


## Fit a RK model and return an object of class "gstatModel" for a list of multiscale grids:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "list"), function(observations, formulaString, covariates, method = list("GLM", "rpart", "randomForest", "HB")[[1]], dimensions = list("3D", "2D", "2D+T", "3D+T")[[1]], family = gaussian, stepwise = TRUE, vgmFun = "Exp", subsample = 5000, ...){

  if(!any(sapply(covariates, class)=="SpatialPixelsDataFrame")){
    stop("List of covariates of class 'SpatialPixelsDataFrame' expected")
  }

  if(!missing(family) & !method=="GLM"){ 
    warning("'family' argument will be ignored. Use 'method=\"GLM\"' instead.")  
  }

  # covariate names:
  covs = unlist(sapply(covariates, FUN=function(x){names(x)}))
  if(!length(unique(covs))==length(covs)){ stop("'Covariates' column names must be unique") }
  
  # target variable:
  methodid = all.vars(formulaString)[1]
     
  # prepare regression matrix:
  ov <- list(NULL)
  for(j in 1:length(covariates)){
    ov[[j]] <- over(x=covariates[[j]], y=observations, method=methodid, var.type="numeric")
      if(nrow(ov[[j]])==0|is.null(ov[[j]]$observedValue)) {
      warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
      }
  }
  xyn = attr(covariates[[1]]@bbox, "dimnames")[[1]]
  ## coordinates in the local space:
  xy <- ov[[1]][,xyn]
  
  gn <- names(observations@data)[!(names(observations@data) %in% c("observationid", xyn))]
  # merge all regression matrices:
  ov <- Reduce(function(x,y) {merge(x[!(names(x) %in% gn)],y[!(names(y) %in% gn)], by="observationid")}, ov)
  ov <- merge(observations@data, ov[,c("observationid", covs)], by="observationid")
  ov <- cbind(ov, xy) 

  # geostats only possible with numeric variables:  
  ov[,methodid] <- as.numeric(ov$observedValue)
  ov$observedValue <- NULL
  
  # check the size of the output object:
  if(nrow(ov)==0|is.null(ov[,methodid])) {
    warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
  }
  
  # fit/filter the regression model:
  m <- fit.regModel(formulaString = formulaString, rmatrix = ov, predictionDomain = covariates[[1]], method = method, dimensions = "3D", family = family, stepwise = stepwise, vgmFun = vgmFun, subsample = subsample, ...)
  
  # save the fitted model:
  return(m) 

})


# end of script;