# Purpose        : Fit a 2D or 3D regression model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : This method will slowly be extended to spatial GLMs, cubist, and Hieararchical Bayes methods;


## Fit a GLM to spatial data:
setMethod("fit.regModel", signature(formulaString = "formula", rmatrix = "data.frame", predictionDomain = "SpatialPixelsDataFrame", method = "character"), function(formulaString, rmatrix, predictionDomain, method = list("GLM", "cubist", "HB")[[1]], dimensions = NULL, family=gaussian, stepwise=TRUE, rvgm=NULL, vgmFun="Exp", subsample = 5000, ...){

  if(method == "GLM"){
  ## target variable name:
  tv = all.vars(formulaString)[1]  
  if(!any(names(rmatrix) %in% tv)){
    stop("Target variable not found in the 'rmatrix' object.")
  }
  
  ## spatial coordinates (column names):
  xyn = attr(predictionDomain@bbox, "dimnames")[[1]]
  ## try to guess the dimensions:
  if(is.null(dimensions)){
    if(length(xyn)==2){ dimensions = "2D" }
    if(length(xyn)==3){ dimensions = "3D" }    
  }  
  if(!any(names(rmatrix) %in% xyn)){
       stop(paste("Column names:", paste(xyn[which(!(xyn %in% names(rmatrix)))], collapse=", "), "could not be located in the regression matrix"))
  }
  
  ## fit/filter the regression model:
  rgm <- glm(formulaString, data=rmatrix, family=family, ...)
  if(stepwise == TRUE){
    rgm <- step(rgm, trace = 0)
  }
   
  ## mask out the missing values:
  if(any(names(rgm) == "na.action")){  rmatrix <- rmatrix[-rgm$na.action,] }
  ## extract the residuals: 
  rmatrix$residual <- residuals(rgm, type="response")
  
  ## test the normality of residuals:
  require(stats)
  st = shapiro.test(rmatrix$residual)
  if(st$p.value < 0.05){
    ## try second test:
    require(nortest)
      at = ad.test(rmatrix$residual)
      if(at$p.value < 0.05){
        warning(paste(st$method, "and", at$method, "report probability of < .05 indicating lack of normal distribution for residuals"), call. = FALSE, immediate. = TRUE)
    }
  }

  ## Fit variogram 2D or 3D:
  if(is.null(rvgm)){
    if(dimensions == "2D"){ rvgm <- fit.vgmModel(residual ~ 1, rmatrix = rmatrix, predictionDomain = predictionDomain, dimensions = "2D", subsample = subsample) }
    if(dimensions == "3D"){ rvgm <- fit.vgmModel(residual ~ 1, rmatrix = rmatrix, predictionDomain = predictionDomain, dimensions = "3D", subsample = subsample) }
  }
  
  ## TH: refit the GLM using the GLS weights? apparently this is possible via the "nlme" package.
  
  ## save the output:
  rkm <- new("gstatModel", regModel = rgm, vgmModel = as.data.frame(rvgm[[1]]), sp = rvgm[[2]])
  return(rkm)
  
  } else {
    stop(paste(method, "method not available at the moment."))
  }

})

## end of script;