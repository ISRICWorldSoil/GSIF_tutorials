# Purpose        : evaluate/test gstatModel for model diagnostics;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : ;

.test.gstat.Model <- function(observations, formulaString, covariates, Ns, predictionLocations, save.predictions = TRUE, debug.level = 0, nfold = 5, ...){
  
  ## get method:                       
  methodid <- all.vars(formulaString)[1]
  
  ## derive Ns if not available: 
  if(missing(Ns)){
    if(class(observations)=="SpatialPointsDataFrame"){ 
      Nmax = round(nrow(observations))
      Nmin = round(20 + length(all.vars(formulaString))*10)
    }
    if(class(observations)=="geosamples"){
      Nmax = round(nrow(observations@data[observations@data$methodid == methodid,]))
      Nmin = round(20 + length(all.vars(formulaString))*10)
    }
    ss = round((runif(10)*sqrt(Nmax-Nmin))^2+Nmin)
    Ns <- sort(c(Nmin, ss, Nmax))
  } else {
    Nmax = Ns[length(Ns)]
  }
  ## check numbers:
  if(any(Ns>Nmax)){
    stop("'Ns' argument contains number larger than total number of observations")
  }
  
  ## prepare empty lists:
  m.l <- list(NULL)  
  tvar.l <- as.list(rep(NA, length(Ns)))
  s.l <- as.list(rep(NA, length(Ns))) 
  p.l <- list(NULL)
  ftime <- as.list(rep(NA, length(Ns)))
  ctime <- as.list(rep(NA, length(Ns)))
  cv <- NULL
  
  ## prepare the prediction locations (if missing):
  if(missing(predictionLocations)){
    if(class(observations)=="SpatialPointsDataFrame"){
       predictionLocations = covariates
    }
    if(class(observations)=="geosamples"){
      predictionLocations <- sp3D(covariates)
    }
  }
  
  ## run model testing...
  message(paste("Running model fitting, cross-validation and predictions for", length(Ns), "sampling intensities using N-fold cross-validation..."))
  pb <- txtProgressBar(min=0, max=length(Ns), style=3)
  for(j in 1:length(Ns)){
    
    ## fit models:
    if(class(observations)=="SpatialPointsDataFrame"){
      observations.s <- observations[sample(1:nrow(observations), Ns[j], replace=FALSE),]
      suppressWarnings(suppressMessages( try( ftime[[j]] <- system.time( m.l[[j]] <- fit.gstatModel(observations = observations.s, formulaString = formulaString, covariates = covariates, ...))[[1]])))
    }
    if(class(observations)=="geosamples"){
      ## subset only methods of interest:
      observations.s <- observations@data[observations@data$methodid == methodid,][sample(1:nrow(observations@data[observations@data$methodid == methodid,]), Ns[j], replace=FALSE),]
      metadata <- observations@methods[observations@methods$methodid == methodid,]
      observations.s <- new("geosamples", registry = observations@registry, methods = metadata, data = observations.s)
      suppressWarnings(suppressMessages( try( ftime[[j]] <- system.time( m.l[[j]] <- fit.gstatModel(observations = observations.s, formulaString = formulaString, covariates = covariates, ...))[[1]])))
    }
    
    ## validate models:
    if(any(class(m.l[[j]]@regModel)=="glm")){
      ## cross-validation on a GLM-kriging model:
      suppressWarnings(suppressMessages( try(cv <- validate(m.l[[j]], debug.level = debug.level, nfold = nfold))))
    } else { 
      if(any(class(m.l[[j]]@regModel)=="rpart")|any(class(m.l[[j]]@regModel)=="randomForest")){
      ## TH: validate function does not run on 'rpart' and 'randomForest' type models because they do not include the input regression matrix!
      ## get the regression matrix:
      ov <- over(observations.s, covariates)
      tv = all.vars(formulaString)[1]
      ov <- cbind(data.frame(observations.s[,tv]), ov)
      sel <- kfold(ov, k=nfold)
      ## re-fit using nfold:
      cv.l <- as.list(rep(NA, length(nfold)))
      for(k in 1:nfold){
        rmatrix <- ov[!sel==k,]
        nlocs <- ov[sel==k,]
        if(ncol(m.l[[j]]@sp@coords)==2){ 
           dimensions = "2D"
        } else {
           dimensions = "3D"      
        }
        coordinates(nlocs) <- as.formula(paste("~", paste(attr(m.l[[j]]@sp@coords, "dimnames")[[2]], collapse ="+")))
        proj4string(nlocs) <- m.l[[j]]@sp@proj4string
        if(any(class(m.l[[j]]@regModel)=="rpart")) { 
          suppressWarnings(suppressMessages(try( mm <- fit.regModel(formulaString=formulaString, rmatrix=rmatrix, predictionDomain=covariates, method="rpart", dimensions=dimensions, stepwise=TRUE, vgmFun=m.l[[j]]@vgmModel$model[2]) )))
        }
        if(any(class(m.l[[j]]@regModel)=="randomForest")) {
          suppressWarnings(suppressMessages(try( mm <- fit.regModel(formulaString=formulaString, rmatrix=rmatrix, predictionDomain=covariates, method="randomForest", dimensions=dimensions, stepwise=TRUE, vgmFun=m.l[[j]]@vgmModel$model[2]) )))
        }
        if(!is.null(mm)){ 
          suppressWarnings(suppressMessages( cv.l[[k]] <- predict.gstatModel(object=mm, predictionLocations=nlocs, nfold=0, block=rep(0, ncol(m.l[[j]]@sp@coords)), mask.extra = FALSE, debug.level=debug.level)$predicted ))
          cv.l[[k]]$observed <- eval(formulaString[[2]], nlocs@data)
          cv.l[[k]]$residual <- cv.l[[k]]$observed - cv.l[[k]]$var1.pred
          cv.l[[k]]$zscore <- cv.l[[k]]$residual/sqrt(cv.l[[k]]$var1.var)
          cv.l[[k]]$fold <- rep(j, length(cv.l[[k]]$residual))
          ## clean up:
          cv.l[[k]]@data <- cv.l[[k]]@data[,c("var1.pred", "var1.var", "observed", "residual", "zscore", "fold")]
        }
      }
      ## derive the cross-validation error
      cv <- list(do.call(rbind, cv.l))
      names(cv) <- "validation"
    
    }}
    
    if(!is.null(cv)&is.list(cv)){
      ## RMSE:
      try(tvar.l[[j]] <- sqrt(mean((cv[[1]]$var1.pred-cv[[1]]$observed)^2, na.rm=TRUE)))
      ## failures:
      try(s.l[[j]] <- sum(cv[[1]]$zscore^2 > 1.5 | abs(cv[[1]]$residual) > 3*sd(cv[[1]]$observed, na.rm = TRUE), na.rm = TRUE))
    }
          
    ## test predictions:   
    if(is.list(predictionLocations)){
      suppressWarnings(suppressMessages(try( ctime[[j]] <- system.time( p.l[[j]] <- lapply(predictionLocations, function(x){ predict(m.l[[j]], predictionLocations = x, mask.extra=TRUE, nfold = 0, debug.level = debug.level)}))[[1]]) ))      
    } else {
      suppressWarnings(suppressMessages(try( ctime[[j]] <- system.time( p.l[[j]] <- predict(m.l[[j]], predictionLocations = predictionLocations, mask.extra=TRUE, nfold = 0, debug.level = debug.level))[[1]]) ))
    }

    setTxtProgressBar(pb, j) 
  }
  close(pb)
  cat(j, "\r")
  flush.console() 
  
  out <- data.frame(samples = Ns, RMSE = unlist(tvar.l), pred.sec = unlist(ftime) + unlist(ctime), failures = unlist(s.l))
  
  if(save.predictions == TRUE){
    names(p.l) <- paste("Ns =", as.character(Ns))
    out <- list(performance = out, predictions = p.l)
  } else{
    out <- list(performance = out, predictions = NULL)  
  }
    
  return(out)

}

setMethod("test.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), .test.gstat.Model)
setMethod("test.gstatModel", signature(observations = "SpatialPointsDataFrame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), .test.gstat.Model)

# end of script;
