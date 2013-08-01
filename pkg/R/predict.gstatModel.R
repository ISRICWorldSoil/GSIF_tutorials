# Purpose        : Predict using 2D or 3D regression-kriging model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl) and Gerard Heuvelink (gerard.heuvelink@wur.nl)
# Contributions  : Bas Kempen; 
# Dev Status     : Pre-Alpha
# Note           : works only with unbiased models with normally distributed residuals;


################## prediction #########################
## predict values using a RK model:
predict.gstatModel <- function(object, predictionLocations, nmin = 10, nmax = 30, debug.level = -1, predict.method = c("RK", "KED")[1], nfold = 5, verbose = FALSE, nsim = 0, mask.extra = TRUE, block = predictionLocations@grid@cellsize, zmin = -Inf, zmax = Inf, subsample = length(object@sp), coarsening.factor = 1, vgmmodel = object@vgmModel, subset.observations = !is.na(object@sp@coords[,1]), betas = c(0,1), ...){

  if(nsim<0|!is.numeric(nsim)){
   stop("To invoke conditional simulations set 'nsim' argument to a positive integer number")
  }
  if(!any(coarsening.factor %in% 1:5)){
   stop("'coarsening.factor' must be an integer in the range 1:5")
  }
  
  ## check the projection system:
  if(is.na(proj4string(predictionLocations))){
    stop("proj4string for the 'predictionLocations' object required")
  }
  if(is.na(proj4string(object@sp))){
    stop("proj4string for the 'gstatModel' object required")
  }
  
  ## force SpatialPixels:
  if(!class(predictionLocations)=="SpatialPointsDataFrame"){
    predictionLocations <- as(predictionLocations, "SpatialPixelsDataFrame")
  }
 
  ## target variable name: 
  if(any(class(object@regModel)=="glm"|class(object@regModel)=="lme"|class(object@regModel)=="gls")){
    variable = all.vars(formula(object@regModel))[1]
  }
  
  ## predict regression model (output is a list):
  if(any(class(object@regModel)=="glm")){

    ## filter the missing classes
    ## TH: this is a simplified solution!
    if(any(x <- sapply(object@regModel$model, is.factor))){
      factors <- names(object@regModel$model)[x]
      for(k in 1:length(factors)){
        dom.class <- summary(object@regModel$model[,factors[k]])
        dom.class <- attr(sort(dom.class, decreasing = TRUE)[1], "names")
        fix.c <- levels(predictionLocations@data[,factors[k]])[!(levels(predictionLocations@data[,factors[k]]) %in% levels(object@regModel$model[,factors[k]]))]
        for(j in fix.c){
          predictionLocations@data[,factors[k]][predictionLocations@data[,factors[k]] == j] <- dom.class
        }
      }
    }
    rp <- stats::predict.glm(object@regModel, newdata=predictionLocations, type="response", se.fit = TRUE, na.action = na.pass)
  }
  ## predict outputs from the nlme package:
  if(any(class(object@regModel)=="lme")){
    require(AICcmodavg)
    rp <- AICcmodavg::predictSE.lme(object@regModel, predictionLocations)
  }
  if(any(class(object@regModel)=="gls")){  
    rp <- list(predict(object@regModel, predictionLocations, na.action = na.pass))
  }
  
  if(any(class(object@regModel)=="rpart")){
    rp <- list(predict(object@regModel, predictionLocations))
    variable = all.vars(attr(object@regModel$terms, "variables"))[1]
  }
  if(any(class(object@regModel)=="quantregForest")){
    covs = attr(object@regModel$forest$ncat, "names")
    rp <- list(predict(object@regModel, predictionLocations@data[,covs], quantile=.5)) 
    variable = attr(object@regModel$y, "name")[1]
  }
  if(any(class(object@regModel)=="randomForest")&!any(class(object@regModel)=="quantregForest")){
    rp <- list(predict(object@regModel, predictionLocations, type="response"))
    variable = all.vars(attr(object@regModel$terms, "variables"))[1]
  }  
  ## rename the target variable:   
  if(any(class(object@regModel) %in% c("randomForest", "rpart", "gls"))){
    names(rp)[1] = "fit"
  }
       
  ## vgm for residuals:
  if(missing(vgmmodel)){
    vgmmodel = object@vgmModel
    class(vgmmodel) <- c("variogramModel", "data.frame")
    ## if the variogram is not significant, then use a NULL model:
    if(diff(vgmmodel$range)==0|diff(vgmmodel$psill)==0){ vgmmodel = NULL }
  }

  ## subset observations to the bounding box set by the prediction domain +- 10%:
  if(missing(subset.observations)){
    if(class(predictionLocations)=="SpatialPixelsDataFrame"){
      R = sqrt(areaSpatialGrid(predictionLocations))/3
    } else {
      R = sqrt(diff(predictionLocations@bbox[1,])*diff(predictionLocations@bbox[2,]))/3
    }
    ## 2D:
    if(length(attr(predictionLocations@bbox, "dimnames")[[1]])==2){
      subset.observations = object@sp@coords[,1] > predictionLocations@bbox[1,1]-.1*R & object@sp@coords[,1] < predictionLocations@bbox[1,2]+.1*R & object@sp@coords[,2] > predictionLocations@bbox[2,1]-.1*R & object@sp@coords[,2] < predictionLocations@bbox[2,2]+.1*R
    } else {
    ## 3D:
      Rv = diff(predictionLocations@bbox[3,])
      subset.observations = object@sp@coords[,1] > predictionLocations@bbox[1,1]-.1*R & object@sp@coords[,1] < predictionLocations@bbox[1,2]+.1*R & object@sp@coords[,2] > predictionLocations@bbox[2,1]-.1*R & object@sp@coords[,2] < predictionLocations@bbox[2,2]+.1*R & object@sp@coords[,3] > predictionLocations@bbox[3,1]-1.5*Rv & object@sp@coords[,3] < predictionLocations@bbox[3,2]+1.5*Rv
    }
  }

  ## observed values:
  if(any(class(object@regModel)=="glm")){
    observed <- SpatialPointsDataFrame(object@sp[subset.observations,], data=object@regModel$model[subset.observations,]) 
  }
  if(any(class(object@regModel) %in% c("randomForest", "rpart"))){
    observed <- SpatialPointsDataFrame(object@sp[subset.observations,], data=data.frame(object@regModel$y[subset.observations])) 
  }
  if(any(class(object@regModel) %in% c("gls", "lme"))){
    observed <- SpatialPointsDataFrame(object@sp[subset.observations,], data=data.frame(fitted.values(object@regModel)[subset.observations] + resid(object@regModel)[subset.observations]))
  }
  
  ## TH: Rename the column? is this necessary?
  ## (this assumes that the first variable on the list is always the target var)
  names(observed@data)[1] = variable

  ## check that the proj4 strings match:
  if(!proj4string(observed)==proj4string(predictionLocations)){
    if(!check_projection(observed, ref_CRS=proj4string(predictionLocations))){
      stop("proj4string at observed and predictionLocations don't match")
    } else { ## force the two proj strings to be exactly the same otherwise gstat has problems!
    suppressWarnings(proj4string(observed) <- proj4string(predictionLocations))
    }
  }
  
  ## try to guess physical limits:
  if(any(class(object@regModel)=="glm")){
    if(missing(zmin) & missing(zmax)){
      if(object@regModel$family$family == "gaussian" & object@regModel$family$link == "log"){ zmin = 0; zmax = Inf }
      if(object@regModel$family$link == "log"){ zmin = 0; zmax = Inf }    
      if(object@regModel$family$family == "binomial"){ zmin = 0; zmax = 1 }
      if(object@regModel$family$family == "quasibinomial"){ zmin = 0; zmax = 1 }
      if(object@regModel$family$family == "poisson"){ zmin = 0; zmax = Inf }  
      if(object@regModel$family$family == "Gamma"){ zmin = 0; zmax = Inf }
    }
  }
  if(any(class(object@regModel) %in% c("glm", "lme", "gls"))){
  ## get fitted valus and residuals:
    observed@data[,paste(variable, "modelFit", sep=".")] <- fitted.values(object@regModel)[subset.observations]
    observed@data[,paste(variable, "residual", sep=".")] <- resid(object@regModel)[subset.observations]
  }
  
  if(any(class(object@regModel)=="rpart")){
    observed@data[,paste(variable, "modelFit", sep=".")] <- predict(object@regModel)[subset.observations]
  }
  if(any(class(object@regModel)=="randomForest")){
    observed@data[,paste(variable, "modelFit", sep=".")] <- object@regModel$predicted[subset.observations]
  }
  if(any(class(object@regModel) %in% c("rpart", "randomForest"))&!is.null(object@regModel$y)){
    observed@data[,paste(variable, "residual", sep=".")] <- (object@regModel$y[subset.observations] - observed@data[,paste(variable, "modelFit", sep=".")])
    rp[["residual.scale"]] <- sqrt(mean((observed@data[,paste(variable, "residual", sep=".")])^2, na.rm=TRUE))
    if(is.null(rp[["residual.scale"]])){ rp[["residual.scale"]] = NA }    
  }
  
  ## remove duplicates as they can lead to singular matrix problems:
  if(length(zerodist(observed))>0){
    observed <- remove.duplicates(observed)
  }
  
  ## skip cross-validation if nfold = 0
  if(nfold==0 | !any(class(object@regModel)=="glm")){ 
    cv <- observed[1]
    names(cv) <- "observed"
    cv$var1.pred <- rep(NA, length(cv$observed))
    cv$var1.var <- rep(NA, length(cv$observed))
    cv$residual <- rep(NA, length(cv$observed))
    cv$zscore <- rep(NA, length(cv$observed))
    cv$fold <- rep(1, length(cv$observed))
    cv <- cv[,c("var1.pred", "var1.var", "observed", "residual", "zscore", "fold")]
  }
  
  ## model summary:
  sum.glm = summary(object@regModel)
  
  ## copy predictions to the predictionLocations
  predictionLocations@data[, paste(variable, "modelFit", sep=".")] <- rp[["fit"]]
  
  ## Kriging with External Drift or Universal kriging:
  ## TH: a simplified implementation
  if(predict.method == "KED"){
      ## Calibration of the trend (single predictor)...
      formString <- as.formula(paste(variable, "~", paste(variable, "modelFit", sep="."), sep=""))
      
      if(nsim==0){
        message("Generating predictions using the trend model (KED method)...")
        rk <- gstat::krige(formula=formString, locations=observed, newdata=predictionLocations, model = vgmmodel, beta = betas, nmin = nmin, nmax = nmax, debug.level = debug.level, block = block, ...)
        
        ## mask extrapolation areas:
        rk@data[,paste(variable, "svar", sep=".")] <- rk$var1.var / var(observed@data[,variable], na.rm=TRUE)
        if(mask.extra==TRUE){
          rk$var1.pred <- ifelse(rk@data[,paste(variable, "svar", sep=".")] > 1, NA, rk$var1.pred) 
        }
        ## mask out values outside physical limits:
        rk@data[,variable] <- ifelse(rk$var1.pred > zmax, zmax, ifelse(rk$var1.pred < zmin, zmin, rk$var1.pred))
        if(any(class(object@regModel)=="glm")){
          if(object@regModel$family$family == "poisson"){
            rk@data[,variable] <- as.integer(round(rk@data[,variable], 0))
          }
        }
        
        ## cross-validation (default implementation by gstat):
        if(nfold>0){      
          message(paste("Running ", nfold, "-fold cross validation...", sep=""))
            ## subset if necessary to speed up the computing:
            if(subsample < length(observed)){
              pcnt <- subsample/length(observed)
              message(paste("Subsetting observations to", signif(pcnt*100, 3), "percent"))
                observed.s <- observed[runif(length(observed))<pcnt,]
                cv <- gstat::krige.cv(formString, locations=observed.s, model=vgmmodel, nfold=nfold, verbose=verbose)
            } else {
                cv <- gstat::krige.cv(formString, locations=observed, model=vgmmodel, nfold=nfold, verbose=verbose)            
            }
          proj4string(cv) = observed@proj4string
        }
      }
      else {
        message(paste("Generating", nsim, "conditional simulations using the trend model (KED method)..."))
        rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = nsim, block = block, ...)
        ## mask out values outside physical limits:
        for(i in 1:nsim){
          rk@data[,i] <- ifelse(rk@data[,i] > zmax, zmax, ifelse(rk@data[,i] < zmin, zmin, rk@data[,i]))
          if(any(class(object@regModel)=="glm")){
            if(object@regModel$family$family == "poisson"){
              rk@data[,i] <- as.integer(round(rk@data[,i], 0))
            }
          }
        }    
      }
  
  ## Regression-kriging (the default approach):
  } else {  
   if(predict.method == "RK"){
      if(any(class(object@regModel)=="glm")|any(class(object@regModel)=="lme")){
        ## GH: prediction variance (regression model only)
        ## [http://en.wikipedia.org/wiki/Confidence_and_prediction_bands]
        predictionLocations@data[,"fit.var"] <- rp[["se.fit"]]^2
      } else {
          if(any(class(object@regModel)=="quantregForest")){     
            ## TH: Prediction error for randomForest
            message("Prediction error for 'randomForest' model estimated using the 'quantreg' package.")
            var.rf <- predict(object@regModel, predictionLocations@data[,covs], quantiles=c((1-.682)/2, 1-(1-.682)/2))
            ## TH: this assumes Normal distribution! [https://en.wikipedia.org/wiki/File:Standard_deviation_diagram.svg]
            predictionLocations@data[,"fit.var"] <- ((var.rf[,1] - var.rf[,2])/2)^2
          } else {
            predictionLocations@data[,"fit.var"] <- 0
          }
      }

      ## predict the residuals:
      formString <- as.formula(paste(paste(variable, "residual", sep="."), "~", 1, sep=""))
      if(nsim==0){
        message("Generating predictions using the trend model (RK method)...")        
        ## TH: if the vgmmodel is null, should we use inverse distance interpolation?
        if(is.null(vgmmodel)){
          ## generate empty grid:
          if(any(class(object@regModel)=="glm")){
            rk = predictionLocations["fit.var"] + rp[["residual.scale"]]^2
          }
          if(any(class(object@regModel)=="lme")){
            rk = predictionLocations["fit.var"]      
          }
          names(rk) = "var1.var"
          rk@data[,variable] <- predictionLocations@data[,paste(variable, "modelFit", sep=".")] 
          rk@data[,paste(variable, "svar", sep=".")] <- predictionLocations$var1.var / var(observed@data[,variable], na.rm=TRUE)
          
          ## cross-validation using GLM:
          if(nfold>0 & any(class(object@regModel)=="glm")){
            cv <- observed[variable]
            names(cv) <- "observed"
            cv$var1.pred <- fitted.values(object@regModel)[subset.observations]
            rp.cv <- stats::predict.glm(object@regModel, newdata=object@regModel$data, type="response", se.fit = TRUE)
            cv$var1.var <- (rp.cv[["se.fit"]][-object@regModel$na.action][subset.observations])^2 + (rp.cv[["residual.scale"]])^2
            message("Running GLM cross-validation without any extra model-fitting...")
            try( cv.err <- glm.diag(object@regModel), silent=TRUE)
            if(class(.Last.value)[1]=="try-error") { 
              cv.err <- data.frame(res = rep(NA, length(cv$observed)), rd = rep(NA, length(cv$observed))) 
            }
            cv$residual <- cv.err$res[subset.observations]
            cv$zscore <- cv.err$rd[subset.observations]
            cv$fold <- rep(1, length(cv$observed))
            # print warning:
            message("Reporting cross-validation results using 'boot::glm.diag'")
            cv <- cv[,c("var1.pred","var1.var","observed","residual","zscore","fold")]
          }
        
        } else {
        ## Predict at coarser grid to speed up the processing (only for 2D maps!)
          if(coarsening.factor>1){
            warning("Downscaling ordinary kriging predictions can lead to artifacts", call.=FALSE, immediate.=TRUE)
            predictionLocations.S = spsample(SpatialPixels(SpatialPoints(predictionLocations["fit.var"]@coords[,1:2], predictionLocations@proj4string)), type="regular", cellsize=coarsening.factor*predictionLocations@grid@cellsize[1])
            gridded(predictionLocations.S) <- TRUE
            rk.S <- gstat::krige(formString, locations=observed, newdata=predictionLocations.S, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, ...)
            rk <- gdalwarp(rk.S, GridTopology = predictionLocations@grid, pixsize=predictionLocations@grid@cellsize[1], resampling_method="cubicspline", tmp.file=TRUE)        
            ## convert to SpatialPixels (necessary!):
            rk <- SpatialPixelsDataFrame(predictionLocations@coords[,1:3], data=rk@data[predictionLocations@grid.index,], proj4string=predictionLocations@proj4string, grid=predictionLocations@grid)
          } else {
            rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, block = block, ...)
          }
          ## sum regression and kriging:
          rk@data[,"var1.pred"] <- predictionLocations@data[,paste(variable, "modelFit", sep=".")] + rk@data[,"var1.pred"]
          rk@data[,"var1.var"] <- (predictionLocations@data[,"fit.var"] + rk@data[,"var1.var"])
          ## TH: This formula assumes that the trend and residuals are independent; which is probably not true
          rk@data[,paste(variable, "svar", sep=".")] <- rk@data[,"var1.var"] / var(observed@data[,variable], na.rm=TRUE)  
          ## mask out values outside the physical limits:
          rk@data[,variable] <- ifelse(rk$var1.pred > zmax, zmax, ifelse(rk$var1.pred < zmin, zmin, rk$var1.pred))
          if(any(class(object@regModel)=="glm")){
            if(object@regModel$family$family == "poisson"){
              rk@data[,variable] <- as.integer(round(rk@data[,variable], 0))
            }
          }

          if(nfold>0){
            ## TH: cross-validation for RK model is not implemented in gstat, so we use the KED model for this purpose:
            formString <- as.formula(paste(variable, "~", paste(variable, "modelFit", sep="."), sep=""))
            message(paste("Running ", nfold, "-fold cross validation using 'krige.cv'...", sep=""))
              ## subset if necessary to speed up the computing:
              if(subsample < length(observed)){
                pcnt <- subsample/length(observed)
                message(paste("Subsetting observations to", signif(pcnt*100, 3), "percent"))
                observed.s <- observed[runif(length(observed))<pcnt,]
                cv <- gstat::krige.cv(formString, locations=observed.s, model=vgmmodel, nfold=nfold, verbose=verbose)
              } else {
                cv <- gstat::krige.cv(formString, locations=observed, model=vgmmodel, nfold=nfold, verbose=verbose)
              }
            proj4string(cv) = observed@proj4string
          }
        }
      } ## Simulations using RK model: 
      else {
        if(is.null(vgmmodel)){
        message(paste("Generating", nsim, "'rnorm' simulations using the trend model (RK method)..."))
        # simple rnorm simulations:
        rk = predictionLocations["fit.var"]
        for(i in 1:nsim){ 
          xsim <- rnorm(length(predictionLocations$fit.var), mean=predictionLocations@data[,paste(variable, "modelFit", sep=".")], sd=sqrt(predictionLocations$fit.var + rp[["residual.scale"]]^2))
          rk@data[,i] <- ifelse(xsim > zmax, zmax, ifelse(xsim < zmin, zmin, xsim))
        }
        names(rk) <- paste("sim", 1:nsim, sep="")
        } else {
        message(paste("Generating", nsim, "conditional simulations using the trend model (RK method)..."))
        rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = nsim, block = block, ...)
        for(i in 1:nsim){
          # sum simulated values:
          xsim = rnorm(length(predictionLocations$fit.var), mean=predictionLocations@data[,paste(variable, "modelFit", sep=".")], sd=sqrt(predictionLocations$fit.var))
          rk@data[,i] <- xsim + rk@data[,i]
          ## TH: this does not costs so much time to compute, but it assumes that the trend and residuals are independent!        
          rk@data[,i] <- ifelse(rk@data[,i] > zmax, zmax, ifelse(rk@data[,i] < zmin, zmin, rk@data[,i]))
          if(any(class(object@regModel)=="glm")){
           if(object@regModel$family$family == "poisson"){
             rk@data[,i] <- as.integer(round(rk@data[,i], 0))
           }
          }
        }
        }
      }
  }
  }
   
  # save the output file:
  if(nsim == 0){
    if(class(predictionLocations)=="SpatialPointsDataFrame"){
    rkp <- list(variable = variable, observed = observed, regModel.summary = sum.glm, vgmModel = object@vgmModel, predicted = rk, validation = cv)
    } else {
      rkp <- new("SpatialPredictions", variable = variable, observed = observed, regModel.summary = sum.glm, vgmModel = object@vgmModel, predicted = rk, validation = cv)
    }
  } 
  else {
    if(class(predictionLocations)=="SpatialPointsDataFrame"){
    rkp <- list(variable = variable, observed = observed, regModel.summary = sum.glm, vgmModel = object@vgmModel, predicted = rk, validation = cv)
    } else {
      t1 <- Line(matrix(c(rk@bbox[1,1],rk@bbox[1,2],mean(rk@bbox[2,]),mean(rk@bbox[2,])), ncol=2))
      transect <- SpatialLines(list(Lines(list(t1), ID="t")), observed@proj4string)
      rkp <- new("RasterBrickSimulations", variable = variable, sampled = transect, realizations = brick(rk))
    }
  }
   
  return(rkp)
  
}

setMethod("predict", signature(object = "gstatModel"), predict.gstatModel)


## predict multiple models independently:
predict.gstatModelList <- function(object, predictionLocations, nmin = 10, nmax = 30, debug.level = -1, predict.method = c("RK", "KED")[1], nfold = 5, verbose = FALSE, nsim = 0, mask.extra = TRUE, block = predictionLocations@grid@cellsize, zmin = -Inf, zmax = Inf, subsample = length(object@sp), ...){

    if(is.list(predictionLocations)&!length(object)==length(predictionLocations)){
      stop("'object' and 'predictionLocations' lists of same size expected")
    }
       
    rkp.l <- list(NULL)
    for(l in 1:length(object)){
      rkp.l[[l]] <- predict(object[[l]], predictionLocations[[l]], nmin = nmin, nmax = nmax, debug.level = debug.level, predict.method = predict.method, nfold = nfold, verbose = verbose, nsim = nsim, mask.extra = mask.extra, block = block, zmin = zmin, zmax = zmax, subsample = subsample, ...)
    }
    
    return(rkp.l)
}

setMethod("predict", signature(object = "list"), predict.gstatModelList)

# end of script;