# Purpose        : validate a gstatModel;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : works only with linear models with normally distributed residuals;

## cross-validate a "gstatModel" object:
setMethod("validate", signature(obj = "gstatModel"), function(obj, nfold = 5, predictionDomain = NULL, save.gstatModels = FALSE, ...){
   ## nfold has to be =>2
   if(nfold < 2){ stop("'nfold' argument > 2 expected") }

   ## full cross-validation only possible with 'glm' class objects
   if(any(class(obj@regModel)=="glm")){          
    ## get the formString:
    formulaString <- formula(obj@regModel)
    mfamily <- obj@regModel$family
    ## get the regression matrix:
    if(!is.null(obj@regModel$na.action)){
      if(length(obj@regModel$na.action)>0){
       ov <- obj@regModel$data[-obj@regModel$na.action,]
     }
    } else {
       ov <- obj@regModel$data
    }
   
    if(nfold > nrow(ov)){ 
      stop("'nfold' argument must not exceed total number of points") 
    }
    ## get the covariates:
    seln <- all.vars(formulaString)[-1]
    tv <- all.vars(formulaString)[1]
    tm <- obj@regModel$terms[[2]]
     
    ## get the variogram:
    vgmmodel <- obj@vgmModel
    class(vgmmodel) <- c("variogramModel", "data.frame")
    ## predictionDomain:
    if(is.null(predictionDomain)){
      obj2D = data.frame(obj@sp@coords[,1:2])
      coordinates(obj2D) <- names(obj2D)
      message("Estimating the predictionDomain...")
      predictionDomain <- vect2rast(remove.duplicates(obj2D))
      proj4string(predictionDomain) <- obj@sp@proj4string
      predictionDomain <- as(predictionDomain, "SpatialPixelsDataFrame")
    } 

    ## re-fit the data in loops:
    m.l <- list(NULL)
    cv.l <- list(NULL)
    sel <- dismo::kfold(ov, k=nfold)
    message(paste("Running ", nfold, "-fold cross validation with model re-fitting...", sep=""))
    for(j in 1:nfold){
      rmatrix <- ov[!sel==j,]
      nlocs <- ov[sel==j,]
      nlocs <- SpatialPointsDataFrame(obj@sp[sel==j,], data=nlocs)
      if(ncol(nlocs@coords)==2){ 
         dimensions = "2D"
      } else {
         dimensions = "3D"      
      }
      try( m.l[[j]] <- fit.regModel(formulaString=formulaString, rmatrix=rmatrix, predictionDomain=predictionDomain, method="GLM", fit.family=mfamily, dimensions=dimensions, stepwise=TRUE, vgmFun=vgmmodel$model[2]) )
      if(!is.null(m.l[[j]])){
        cv.l[[j]] <- predict.gstatModel(object=m.l[[j]], predictionLocations=nlocs, nfold=0, block=rep(0, ncol(obj@sp@coords)), mask.extra = FALSE, ...)$predicted
        cv.l[[j]]$observed <- eval(tm, nlocs@data)
        cv.l[[j]]$residual <- cv.l[[j]]$observed - cv.l[[j]]$var1.pred
        cv.l[[j]]$zscore <- cv.l[[j]]$residual/sqrt(cv.l[[j]]$var1.var)
        cv.l[[j]]$fold <- rep(j, length(cv.l[[j]]$residual))
        ## clean up:
        cv.l[[j]]@data <- cv.l[[j]]@data[,c("var1.pred", "var1.var", "observed", "residual", "zscore", "fold")]
      }
    }
   
    if(save.gstatModels==TRUE){ 
      cv <- list(do.call(rbind, cv.l), m.l)
      names(cv) <- c("validation", "gstatModels")
    } else {
      cv <- list(do.call(rbind, cv.l))
      names(cv) <- "validation"
    }
    
    return(cv)   
   
   } else { 
    stop(paste("Validate method for a 'regModel' of class", class(obj@regModel), "not available.")) 
   }

})

# end of script;