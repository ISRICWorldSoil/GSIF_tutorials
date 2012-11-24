# Purpose        : Predict using 2D or 3D regression-kriging model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : works only with linear models with normally distributed residuals;


################## prediction #########################
## predict values using a RK model:
predict.gstatModel <- function(object, predictionLocations, nmin = 10, nmax = 30, debug.level = -1, predict.method = c("KED", "RK")[1], nfold = 5, verbose = FALSE, nsim = 0, mask.extra = TRUE, block = predictionLocations@grid@cellsize, method = list("GLM", "cubist", "HB")[[1]], zmin = -Inf, zmax = Inf, ...){

  if(nsim<0|!is.numeric(nsim)){
   stop("To invoke conditional simulations set 'nsim' argument to a positive integer number")
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

  if(method == "GLM"){
   
  ## predict regression model:
  rp <- stats::predict.glm(object@regModel, newdata=predictionLocations, type="response", se.fit = TRUE, na.action = na.pass)
  ## target variable name: 
  variable = all.vars(object@regModel$formula)[1]
     
  ## vgm for residuals:
  vgmmodel = object@vgmModel
  class(vgmmodel) <- c("variogramModel", "data.frame")
  ## if the variogram is not significant, then use a NULL model:
  if(attr(vgmmodel, "singular")==TRUE|diff(vgmmodel$range)==0|diff(vgmmodel$psill)==0){ vgmmodel = NULL }

  ## observed values:
  observed <- SpatialPointsDataFrame(object@sp, data=object@regModel$model)
  ## rename the column if necessary
  ## TH: this assumes that the first variable on the list is the target var:
  names(observed@data)[1] = variable

  ## check that the proj4 strings match:
  if(!proj4string(observed)==proj4string(predictionLocations)){
    stop("proj4string at observed and predictionLocations don't match")
  }
  
  ## physical limits:
  if(missing(zmin) & missing(zmax)){
    if(object@regModel$family$family == "gaussian" & object@regModel$family$link == "log"){ zmin = 0; zmax = Inf }
    if(object@regModel$family$link == "log"){ zmin = 0; zmax = Inf }    
    if(object@regModel$family$family == "binomial"){ zmin = 0; zmax = 1 }
    if(object@regModel$family$family == "poisson"){ zmin = 0; zmax = Inf }  
    if(object@regModel$family$family == "gamma"){ zmin = 0; zmax = Inf }    
  }
  
  # values after transformation:
  observed@data[,paste(variable, "glmfit", sep=".")] <- fitted.values(object@regModel)
  observed@data[,paste(variable, "residual", sep=".")] <- residuals(object@regModel)
  
  ## remove duplicates as they can lead to singular matrix problems:
  if(length(zerodist(observed))>0){
    observed <- remove.duplicates(observed)
  }
  
  ## skip cross-validation in nfold = 0
  if(nfold==0){ 
    cv <- observed[1]
    names(cv) <- "observed"
    cv$var1.pred <- rep(NA, length(cv$observed))
    cv$var1.var <- rep(NA, length(cv$observed))
    cv$residual <- rep(NA, length(cv$observed))
    cv$zscore <- rep(NA, length(cv$observed))
    cv$fold <- rep(1, length(cv$observed))
    cv <- cv[,c("var1.pred", "var1.var", "observed", "residual", "zscore", "fold")]
  }
  
  # model summary:
  sum.glm = summary(object@regModel)
  class(sum.glm) = "list"
  
  # copy GLM predictions to the predictionLocations
  predictionLocations@data[, paste(variable, "glmfit", sep=".")] <- rp$fit
  
  if(predict.method == "KED"){
      formString <- as.formula(paste(variable, "~", paste(variable, "glmfit", sep="."), sep=""))
      if(nsim==0){
        message("Generating predictions using the trend model (KED method)...")
        rk <- gstat::krige(formula=formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, block = block, ...)
        
        ## mask extrapolation areas:
        rk@data[,paste(variable, "svar", sep=".")] <- rk$var1.var / var(observed@data[,variable], na.rm=TRUE)
        if(mask.extra==TRUE){
          rk$var1.pred <- ifelse(rk@data[,paste(variable, "svar", sep=".")] > 1, NA, rk$var1.pred) 
        }
        ## mask out values outside physical limits:
        rk@data[,variable] <- ifelse(rk$var1.pred > zmax, zmax, ifelse(rk$var1.pred < zmin, zmin, rk$var1.pred))
        
        ## cross-validation:
        if(nfold>0){      
          message(paste("Running ", nfold, "-fold cross validation...", sep=""))
          cv <- gstat::krige.cv(formString, locations=observed, model=vgmmodel, nfold=nfold, verbose=verbose)
          proj4string(cv) = observed@proj4string
        }
      }
      else {
        message(paste("Generating", nsim, "conditional simulations using the trend model (KED method)..."))
        rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = nsim, block = block, ...)
        ## mask out values outside physical limits:
        for(i in 1:nsim){
          rk@data[,i] <- ifelse(rk@data[,i] > zmax, zmax, ifelse(rk@data[,i] < zmin, zmin, rk@data[,i]))
        }    
      }
  } else {  
   if(predict.method == "RK"){
      predictionLocations@data[, "se.fit"] <- rp$se.fit
      # predict the residuals:
      formString <- as.formula(paste(paste(variable, "residual", sep="."), "~", 1, sep=""))
      if(nsim==0){
        message("Generating predictions using the trend model (RK method)...")
        # TH: if the vgmmodel is null leave the inverse distance interpolation out?
        if(is.null(vgmmodel)){
          # generate empty grid:
          rk = predictionLocations["se.fit"]
          rk@data[,variable] <- predictionLocations@data[,paste(variable, "glmfit", sep=".")] 
          rk@data[,paste(variable, "svar", sep=".")] <- rp$se.fit/var(observed@data[,variable], na.rm=TRUE)
          rk@data[,1] <- NULL
          # cross-validation using GLM:
          if(nfold>0){
            require(boot)
            try(cv.err <- glm.diag(object@regModel), silent=TRUE)
            # cv.err <- boot::cv.glm(data=object@regModel$data, glmfit=object@regModel, K=nfold)
            ## TH: boot::cv.glm fails for unknow reason? 
            if(class(.Last.value)[1]=="try-error") { cv.err <- data.frame(res = rep(NA, length(cv$observed)), rd = rep(NA, length(cv$observed))) }
            cv <- observed[variable]
            names(cv) <- "observed"
            cv$var1.pred <- fitted.values(object@regModel)
            cv$var1.var <- stats::predict.glm(object@regModel, newdata=object@regModel$data, type="response", se.fit = TRUE)$se.fit[-object@regModel$na.action]
            cv$residual <- cv.err$res
            cv$zscore <- cv.err$rd
            cv$fold <- rep(1, length(cv$observed))
            # print warning:
            message("Reporting cross-validation results using 'boot::glm.diag'")
            cv <- cv[,c("var1.pred","var1.var","observed","residual","zscore","fold")]
          }
        } else {
          rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, block = block, ...) 
          ## sum regression and kriging:
          rk@data[,"var1.pred"] <- predictionLocations@data[,paste(variable, "glmfit", sep=".")] + rk@data[,"var1.pred"]
          rk@data[,paste(variable, "svar", sep=".")] <- (predictionLocations$se.fit + rk@data[,"var1.var"]) / var(observed@data[,variable], na.rm=TRUE)  
          ## TH: This formula assumes that the trend and residuals are independent; which is probably not true
          ## mask out values outside the physical limits:
          rk@data[,variable] <- ifelse(rk$var1.pred > zmax, zmax, ifelse(rk$var1.pred < zmin, zmin, rk$var1.pred))

          if(nfold>0){
            formString <- as.formula(paste(variable, "~", paste(variable, "glmfit", sep="."), sep=""))
            message(paste("Running ", nfold, "-fold cross validation...", sep=""))
            cv <- gstat::krige.cv(formString, locations=observed, model=vgmmodel, nfold=nfold, verbose=verbose)
            proj4string(cv) = observed@proj4string
          }
        }
      } 
      else {
        if(is.null(vgmmodel)){
        message(paste("Generating", nsim, "'rnorm' simulations using the trend model (RK method)..."))
        # simple rnorm simulations:
        rk = predictionLocations["se.fit"]
        for(i in 1:nsim){ 
          xsim <- rnorm(length(predictionLocations$se.fit), mean=predictionLocations@data[,paste(variable, "glmfit", sep=".")], sd=predictionLocations$se.fit)
          rk@data[,i] <- ifelse(xsim > zmax, zmax, ifelse(xsim < zmin, zmin, xsim))
        }
        names(rk) <- paste("sim", 1:nsim, sep="")
        } else {
        message(paste("Generating", nsim, "conditional simulations using the trend model (RK method)..."))
        for(i in 1:nsim){
          rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = 1, block = block, ...)
          # sum simulated values:
          xsim = rnorm(length(predictionLocations$se.fit), mean=predictionLocations@data[,paste(variable, "glmfit", sep=".")], sd=predictionLocations$se.fit)
          rk@data[,i] <- xsim + rk@data[,i]
          ## TH: this is inexpensive but it assumes that the trend and residuals are independent!        
          rk@data[,i] <- ifelse(rk@data[,i] > zmax, zmax, ifelse(rk@data[,i] < zmin, zmin, rk@data[,i]))
        }
        names(rk) <- paste("sim", 1:nsim, sep="")

        }
      }
  }
  }
   
  # save the output file:
  if(nsim == 0){
    if(class(predictionLocations)=="SpatialPointsDataFrame"){
    rkp <- list(variable = variable, observed = observed, glm = sum.glm, vgmModel = object@vgmModel, predicted = rk, validation = cv)
    } else {
      rkp <- new("SpatialPredictions", variable = variable, observed = observed, glm = sum.glm, vgmModel = object@vgmModel, predicted = rk, validation = cv)
    }
  } 
  else {
    if(class(predictionLocations)=="SpatialPointsDataFrame"){
    rkp <- list(variable = variable, observed = observed, glm = sum.glm, vgmModel = object@vgmModel, predicted = rk, validation = cv)
    } else {
      t1 <- Line(matrix(c(rk@bbox[1,1],rk@bbox[1,2],mean(rk@bbox[2,]),mean(rk@bbox[2,])), ncol=2))
      transect <- SpatialLines(list(Lines(list(t1), ID="t")), observed@proj4string)
      rkp <- new("RasterBrickSimulations", variable = variable, sampled = transect, realizations = brick(rk))
    }
  }
   
  return(rkp)
  
  } else {
    stop(paste(method, "method not available at the moment."))
  }

}

setMethod("predict", signature(object = "gstatModel"), predict.gstatModel)


## predict multiple models independently:
predict.gstatModelList <- function(object, predictionLocations, nmin = 10, nmax = 30, debug.level = -1, predict.method = c("KED", "RK")[1], nfold = 5, verbose = FALSE, nsim = 0, mask.extra = TRUE, ...){
    if(is.list(predictionLocations)&!length(object)==length(predictionLocations)){
      stop("'object' and 'predictionLocations' lists of same size expected")
    }
       
    rkp.l <- list(NULL)
    for(l in 1:length(object)){
      rkp.l[[l]] <- predict(object[[l]], predictionLocations[[l]], nmin = nmin, nmax = nmax, debug.level = debug.level, predict.method = predict.method, nfold = nfold, verbose = verbose, nsim = nsim, mask.extra = mask.extra, ...)
    }
    
    return(rkp.l)
}

setMethod("predict", signature(object = "list"), predict.gstatModelList)




# end of script;