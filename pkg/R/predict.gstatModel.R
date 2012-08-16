# Purpose        : Fit/predict a 2D or 3D regression-kriging model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : works only with linear models with normally distributed residuals;


################## model fitting #########################
## Fit a GLM to spatial data:
glm.sp <- function(formulaString, rmatrix, predictionDomain, family=gaussian, stepwise=TRUE, rvgm=NULL, vgmFun="Exp", type, ...){
  
  require(gstat)
  require(splines)
  require(plotKML)

  # fit/filter the regression model:
  rgm <- glm(formulaString, data=rmatrix, family=family, ...)
  if(stepwise == TRUE){
    rgm <- step(rgm, trace = 0)
  }
  
  tv = all.vars(formulaString)[1]  
  if(!any(names(rmatrix) %in% tv)){
    stop("Target variable not found in the 'rmatrix' object.")
  } 
  
  # mask out the missing values:
  if(any(names(rgm) == "na.action")){  rmatrix <- rmatrix[-rgm$na.action,] }
  # extract the residuals in the transformed space:
  linkfun <- rgm$family$linkfun
  rmatrix$residual <- linkfun(rmatrix[,tv]) - rgm$linear.predictors
  
  # try to guess the type of rmatrix
  if(missing(type)){
  if(sum(names(rmatrix) %in% c("longitude", "latitude", "altitude"))==3){
    type = "geosamples"
  } else {
    type = "SpatialPointsDataFrame"
  }
  }
  
  # create a 2D / 3D object:
  if(type == "geosamples"){
    coordinates(rmatrix) <- ~ longitude + latitude + altitude
    proj4string(rmatrix) = get("ref_CRS", envir = GSIF.opts)
    suppressWarnings(rmatrix <- spTransform(rmatrix, predictionDomain@proj4string))
  } else {
    if(type == "SpatialPointsDataFrame"){
      xyn = attr(predictionDomain@bbox, "dimnames")[[1]]
      coordinates(rmatrix) <- as.formula(paste("~", paste(xyn, collapse = "+"), sep=""))
      proj4string(rmatrix) = predictionDomain@proj4string
    }
  }
  
  # if the residual variogram is unknown:
  if(is.null(rvgm)){
    if(type == "geosamples"){
      # estimate area extent:
      Range = sqrt(areaSpatialGrid(predictionDomain))/3
      # estimate initial range in the vertical direction:
      dr <- abs(diff(range(rmatrix@coords[,3], na.rm=TRUE)))/3
      a2 = 2*dr/Range

      # estimate anisotropy parameters:
      anis = c(0, 0, 0, 1, a2)    
      ivgm <- vgm(nugget=0, model=vgmFun, range=Range, psill=var(rmatrix$residual), anis = anis)
      # fit the 3D variogram:
      try(rvgm <- gstat::fit.variogram(variogram(residual ~ 1, rmatrix), ivgm, ...))    
    ## TH: This is the most simple implementation - fitting of 3D variograms needs to be improved!
    } else {
    if(type == "SpatialPointsDataFrame"){
      # fit the variogram:
      if(!is.na(proj4string(predictionDomain))){
      if(!is.projected(predictionDomain)){
        require(fossil)  # Haversine Formula for Great Circle distance
        p.1 <- matrix(c(predictionDomain@bbox[1,1], predictionDomain@bbox[1,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
        p.2 <- matrix(c(predictionDomain@bbox[2,1], predictionDomain@bbox[2,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
        Range = fossil::deg.dist(lat1=p.1[,2], long1=p.1[,1], lat2=p.2[,2], long2=p.2[,1])/2
      } else {
        Range = sqrt(areaSpatialGrid(predictionDomain))/2      
      }} else{
        Range = sqrt(areaSpatialGrid(predictionDomain))/2
      }
        ivgm <- vgm(nugget=0, model=vgmFun, range=Range, psill=var(rmatrix$residual))
        try(rvgm <- gstat::fit.variogram(variogram(residual ~ 1, rmatrix), ivgm, ...))
        if(diff(rvgm$range)==0|diff(rvgm$psill)==0){
          warning("Variogram shows no spatial dependence")    
      }
    }
  }}
  
  if(diff(rvgm$range)==0|diff(rvgm$psill)==0){
    warning("Variogram shows no spatial dependence")
  }
  
  out = list(rgm, rvgm, as(rmatrix, "SpatialPoints"))
  names(out) = c("regModel", "vgmModel", "sp")
  
  return(out)

}



## Fit a RK model and return an object of class "gstatModel":
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, methodid, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, ...){
   
  # prepare regression matrix:
  ov <- overlay(x=covariates, y=observations, method=methodid, var.type = "numeric")
  if(nrow(ov)==0|is.null(ov$observedValue)) {
    warning("The overlay operations resulted in an empty set. Check 'methodid' column.")
  }
  # geostats only possible with numeric variables:
  ov$observedValue = as.numeric(ov$observedValue)
  
  # fit/filter the regression model:
  m <- glm.sp(formulaString=formulaString, rmatrix=ov, predictionDomain=covariates, family=family, stepwise=stepwise, rvgm=rvgm, vgmFun=vgmFun, type="geosamples", ...)
  
  # save the fitted model:
  rkm <- new("gstatModel", regModel = m[[1]], vgmModel = as.data.frame(m[[2]]), sp = m[[3]])
  return(rkm)  

})


## Fit a RK model to a list of covariates / formulas:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "list", covariates = "list"), function(observations, formulaString, covariates, methodid, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, ...){
    if(!length(formulaString)==length(covariates)){
      stop("'formulaString' and 'covariates' lists of same size expected")
    }
    
    rkm.l <- list(NULL)
    for(l in 1:length(covariates)){   
      rkm.l[[l]] <- fit.gstatModel(observations, formulaString[[l]], covariates[[l]], methodid = methodid, family = family, stepwise = stepwise, rvgm = rvgm, ...)
    }
    return(rkm.l)  

})


## Fit a RK model and return an object of class "gstatModel" for a list of multiscale grids:
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "list"), function(observations, formulaString, covariates, methodid, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, ...){

  if(!any(sapply(covariates, class)=="SpatialPixelsDataFrame")){
    stop("List of covariates of class 'SpatialPixelsDataFrame' expected")
  }

  # covariate names:
  covs = unlist(sapply(covariates, FUN=function(x){names(x)}))
  if(!length(unique(covs))==length(covs)){ stop("'Covariates' column names must be unique") }
     
  # prepare regression matrix:
  ov <- list(NULL)
  for(j in 1:length(covariates)){
    ov[[j]] <- overlay(x=covariates[[j]], y=observations, method=methodid, var.type = "numeric")
      if(nrow(ov[[j]])==0|is.null(ov[[j]]$observedValue)) {
      warning("The overlay operations resulted in an empty set. Check 'methodid' column.")
      }
  }
  gn <- names(observations@data)[!(names(observations@data) %in% "observationid")]
  # merge all regression matrices:
  ov = Reduce(function(x,y) {merge(x[!(names(x) %in% gn)],y[!(names(y) %in% gn)], by="observationid")}, ov)
  ov = merge(observations@data, ov[,c("observationid", covs)], by="observationid")

  # geostats only possible with numeric variables:  
  ov$observedValue = as.numeric(ov$observedValue)
  
  # check if the size of the output object:
  if(nrow(ov)==0|is.null(ov$observedValue)) {
    warning("The overlay operations resulted in an empty set. Check 'methodid' column.")
  }
  
  # fit/filter the regression model:
  m <- glm.sp(formulaString=formulaString, rmatrix=ov, predictionDomain=covariates[covs], family=family, stepwise=stepwise, rvgm=rvgm, vgmFun=vgmFun, type="geosamples", ...)
  
  # save the fitted model:
  rkm <- new("gstatModel", regModel = m[[1]], vgmModel = as.data.frame(m[[2]]), sp = m[[3]])
  return(rkm) 

})


## Fit a 'simple' 2D RK model:
setMethod("fit.gstatModel", signature(observations = "SpatialPointsDataFrame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, ...){
 
  # the function only works with 2D maps:
  if(length(attr(coordinates(observations), "dimnames")[[2]])>2){
    warning("This method uses only 2D coordinates of the points. For 3D data consider using the 'geosamples-class'.")
  }

  # overlay points:
  index <- overlay(covariates, observations)
  sel <- !is.na(index)
  # all variables of interest:
  tv = all.vars(formulaString)[1]
  seln = names(covariates) %in% all.vars(formulaString)[-1]
  if(length(seln)==0){
      stop("None of the covariates in the 'formulaString' do not match the names in the'covariates' object")
  }
  ov <- cbind(data.frame(observations[sel,tv]), covariates[index[sel],])  

  # check the size of the output:
  if(nrow(ov)==0|is.null(ov[,tv])) {
    warning("The overlay operations resulted in an empty set. Check 'methodid' column.")
  }
  
  # fit/filter the regression model:
  m <- glm.sp(formulaString=formulaString, rmatrix=ov, predictionDomain=covariates[seln], family=family, stepwise=stepwise, rvgm=rvgm, vgmFun=vgmFun, type="SpatialPointsDataFrame", ...)  
    
  # save the fitted model:
  rkm <- new("gstatModel", regModel = m[[1]], vgmModel = as.data.frame(m[[2]]), sp = m[[3]])
  return(rkm)
})


################## cross-validation #########################
## cross-validate a "gstatModel" object:
setMethod("validate", signature(obj = "gstatModel"), function(obj, nfold = 5, predictionDomain = NULL, save.gstatModels = FALSE){

   require(dismo)
   require(gstat)
   require(plotKML)
   if(nfold < 2){ stop("'nfold' argument > 2 expected") }

   # get the formString:
   formulaString <- obj@regModel$formula
   mfamily <- obj@regModel$family
   linkfun <- obj@regModel$family$linkfun
   # get the regression matrix:
   ov <- obj@regModel$data[-obj@regModel$na.action,]
   if(nfold > nrow(ov)){ stop("'nfold' argument must not exceed total number of points") }
   # get the covariates:
   seln = all.vars(formulaString)[-1]
   tv = all.vars(formulaString)[1]
   # get the variogram:
   vgmmodel = obj@vgmModel
   class(vgmmodel) <- c("variogramModel", "data.frame")
   # predictionDomain:
   if(is.null(predictionDomain)){
     obj2D = data.frame(obj@sp@coords[,1:2])
     coordinates(obj2D) <- names(obj2D)
     message("Estimating the predictionDomain...")
     predictionDomain <- vect2rast(remove.duplicates(obj2D))
     proj4string(predictionDomain) <- obj@sp@proj4string
   }

   # re-fit the data in loops:
   m.l <- list(NULL)
   cv.l <- list(NULL)
   sel <- kfold(ov, k=nfold)
   message(paste("Running ", nfold, "-fold cross validation...", sep=""))
   for(j in 1:nfold){
      rmatrix <- ov[!sel==j,]
      nlocs <- ov[sel==j,]
      nlocs <- SpatialPointsDataFrame(obj@sp[sel==j,], data=nlocs)
      m <- glm.sp(formulaString=formulaString, rmatrix=rmatrix, predictionDomain=predictionDomain, family=mfamily, stepwise=TRUE, vgmFun=vgmmodel$model[2])
      m.l[[j]] <- new("gstatModel", regModel = m[[1]], vgmModel = as.data.frame(m[[2]]), sp = m[[3]])
      cv.l[[j]] <- predict.gstatModel(object=m.l[[j]], predictionLocations=nlocs, nfold=0, block=rep(0, ncol(obj@sp@coords)), mask.extra = FALSE)$predicted
      cv.l[[j]]$observed <- linkfun(nlocs@data[,tv])
      cv.l[[j]]$residual <- cv.l[[j]]$observed - cv.l[[j]]$var1.pred
      cv.l[[j]]$zscore <- cv.l[[j]]$residual/sqrt(cv.l[[j]]$var1.var)
      cv.l[[j]]$fold <- rep(j, length(cv.l[[j]]$residual))
      # clean up:
      cv.l[[j]]@data <- cv.l[[j]]@data[,c("var1.pred","var1.var","observed","residual","zscore","fold")]
   }
   
   if(save.gstatModels==TRUE){ 
    cv <- list(do.call(rbind, cv.l), m.l)
    names(cv) <- c("validation", "gstatModels")
   } else {
    cv <- list(do.call(rbind, cv.l))
    names(cv) <- "validation"
   }
    
   return(cv)   

})




################## prediction #########################
## predict values using a RK model:
predict.gstatModel <- function(object, predictionLocations, nmin = 10, nmax = 30, debug.level = -1, method = c("KED", "RK")[1], nfold = 5, verbose = FALSE, nsim = 0, mask.extra = TRUE, block = predictionLocations@grid@cellsize, ...){

  if(nsim<0|!is.numeric(nsim)){
   stop("To invoke conditional simulations set 'nsim' argument to a positive integer number")
  }

  require(gstat)
  require(raster)
  
  # check the projection system:
  if(is.na(proj4string(predictionLocations))){
    stop("proj4string for the 'predictionLocations' object required")
  }
  if(is.na(proj4string(object@sp))){
    stop("proj4string for the 'gstatModel' object required")
  }
  
  # force SpatialPixels:
  if(!class(predictionLocations)=="SpatialPointsDataFrame"){
    predictionLocations <- as(predictionLocations, "SpatialPixelsDataFrame")
  }

  # predict regression model:
  rp <- stats::predict.glm(object@regModel, newdata=predictionLocations, type="link", se.fit = TRUE, na.action = na.pass)
  variable = all.vars(object@regModel$formula)[1] # target variable
  # the vgm for residuals:
  vgmmodel = object@vgmModel
  class(vgmmodel) <- c("variogramModel", "data.frame")
  # if the variogram is not significant, then use a NULL model:
  if(attr(vgmmodel, "singular")==TRUE|diff(vgmmodel$range)==0|diff(vgmmodel$psill)==0){ vgmmodel = NULL }
  # observed values:
  observed <- SpatialPointsDataFrame(object@sp, data=object@regModel$model)
  if(!proj4string(observed)==proj4string(predictionLocations)){
    stop("proj4string at observed and predictionLocations don't match")
  }

  # back-transform function:
  linkfun <- object@regModel$family$linkfun
  invfun <- object@regModel$family$linkinv
  
  # values after transformation:
  observed@data[,paste(variable, "glmfit", sep=".")] <- object@regModel$linear.predictors
  observed@data[,paste(variable, "link", sep=".")] <- linkfun(observed@data[,variable])
  observed@data[,paste(variable, "residual", sep=".")] <- linkfun(observed@data[,variable]) - object@regModel$linear.predictors 
  
  # remove duplicates as they can lead to singular matrix problems:
  if(length(zerodist(observed))>0){
    observed <- remove.duplicates(observed)
  }
  
  # skip cross-validation in nfold = 0
  if(nfold==0){ 
    cv <- observed[variable]
    names(cv) <- "observed"
    cv$var1.pred <- rep(NA, length(cv$observed))
    cv$var1.var <- rep(NA, length(cv$observed))
    cv$residual <- rep(NA, length(cv$observed))
    cv$zscore <- rep(NA, length(cv$observed))
    cv$fold <- rep(1, length(cv$observed))
    cv <- cv[,c("var1.pred","var1.var","observed","residual","zscore","fold")]
  }
  
  # model summary:
  sum.glm = summary(object@regModel)
  class(sum.glm) = "list"
  
  # copy GLM predictions to the predictionLocations
  predictionLocations@data[, paste(variable, "glmfit", sep=".")] <- rp$fit
  
  if(method == "KED"){
      formString <- as.formula(paste(paste(variable, "link", sep="."), "~", paste(variable, "glmfit", sep="."), sep=""))
      if(nsim==0){
        message("Generating predictions using the trend model (KED method)...")
        rk <- gstat::krige(formula=formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, block = block, ...)
        # mask extrapolation areas:
        rk@data[,paste(variable, "svar", sep=".")] <- rk$var1.var / var(observed@data[,paste(variable, "link", sep=".")], na.rm=TRUE)
        if(mask.extra==TRUE){
          rk$var1.pred <- ifelse(rk@data[,paste(variable, "svar", sep=".")] > 1, NA, rk$var1.pred) 
        }
        # back-transform the values:
        rk@data[,variable] = invfun(rk$var1.pred)
        # cross-validation:
        if(nfold>0){      
          message(paste("Running ", nfold, "-fold cross validation...", sep=""))
          cv <- gstat::krige.cv(formString, locations=observed, model=vgmmodel, nfold=nfold, verbose=verbose)
          proj4string(cv) = observed@proj4string
        }
      }
      else {
        message(paste("Generating", nsim, "conditional simulations using the trend model (KED method)..."))
        rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = nsim, block = block, ...)
        # back-transform the values:
        for(i in 1:nsim){
          rk@data[,i] = invfun(rk@data[,i])
        }    
      }
  }

  else {  
   if(method == "RK"){
      predictionLocations@data[, "se.fit"] <- rp$se.fit
      # predict the residuals:
      formString <- as.formula(paste(paste(variable, "residual", sep="."), "~", 1, sep=""))
      if(nsim==0){
        message("Generating predictions using the trend model (RK method)...")
        # TH: if the vgmmodel is null leave the inverse distance interpolation out?
        if(is.null(vgmmodel)){
          # generate empty grid:
          rk = predictionLocations["se.fit"]
          rk@data[,variable] <- invfun(predictionLocations@data[,paste(variable, "glmfit", sep=".")]) 
          rk@data[,paste(variable, "svar", sep=".")] <- rp$se.fit/var(observed@data[,paste(variable, "link", sep=".")], na.rm=TRUE)
          rk@data[,1] <- NULL
          # cross-validation using GLM:
          if(nfold>0){
            require(boot)
            cv.err <- glm.diag(object@regModel)
            # cv.err <- boot::cv.glm(data=object@regModel$data, glmfit=object@regModel, K=nfold)
            ## TH: boot::cv.glm fails for unknow reason? 
            cv <- observed[variable]
            names(cv) <- "observed"
            cv$var1.pred <- object@regModel$linear.predictors
            cv$var1.var <- stats::predict.glm(object@regModel, newdata=object@regModel$data, type="link", se.fit = TRUE)$se.fit[-object@regModel$na.action]
            cv$residual <- cv.err$res
            cv$zscore <- cv.err$rd
            cv$fold <- rep(1, length(cv$observed))
            # print warning:
            warning("Reporting cross-validation results using 'boot::glm.diag'")
            cv <- cv[,c("var1.pred","var1.var","observed","residual","zscore","fold")]
          }
        } else {
          rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, block = block, ...) 
          # sum regression and kriging, and back-transform the values:
          rk@data[,variable] <- invfun(predictionLocations@data[,paste(variable, "glmfit", sep=".")] + rk@data[,"var1.pred"])
          rk@data[,paste(variable, "svar", sep=".")] <- ( predictionLocations$se.fit + rk@data[,"var1.var"] ) / var(observed@data[,paste(variable, "link", sep=".")], na.rm=TRUE)  ## This formula assumes that the trend and residuals are independent; which is probably not true
          if(nfold>0){
            formString <- as.formula(paste(paste(variable, "link", sep="."), "~", paste(variable, "glmfit", sep="."), sep=""))
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
          rk@data[,i] = invfun(xsim)
          }
          names(rk) <- paste("sim", 1:nsim, sep="")
        } else {
        message(paste("Generating", nsim, "conditional simulations using the trend model (RK method)..."))
        for(i in 1:nsim){
          rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = nsim, block = block, ...)
          # back-transform the values:
          rk@data[,i] = invfun(xsim + rk@data[,i])  #TH: this is inexpensive but it assumes that the trend and residuals are independent!        
        }
        }
      }
      }
  }
   
  # save the output file:
  require(plotKML)
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

}

setMethod("predict", signature(object = "gstatModel"), predict.gstatModel)


## predict multiple models independently:
predict.gstatModelList <- function(object, predictionLocations, nmin = 10, nmax = 30, debug.level = -1, method = c("KED", "RK")[1], nfold = 5, verbose = FALSE, nsim = 0, mask.extra = TRUE, ...){
    if(is.list(predictionLocations)&!length(object)==length(predictionLocations)){
      stop("'object' and 'predictionLocations' lists of same size expected")
    }
       
    rkp.l <- list(NULL)
    for(l in 1:length(object)){
      rkp.l[[l]] <- predict(object[[l]], predictionLocations[[l]], nmin = nmin, nmax = nmax, debug.level = debug.level, method = method, nfold = nfold, verbose = verbose, nsim = nsim, mask.extra = mask.extra, ...)
    }
    
    return(rkp.l)
}

setMethod("predict", signature(object = "list"), predict.gstatModelList)


## Get a summary as a data frame:
setMethod("summary", signature(object = "SpatialPredictions"), function(object){
   require(plotKML)
   
   z <- NULL
   z$variable = object@variable
   z$minium = range(object@observed@data[,object@variable])[1]
   z$maximum = range(object@observed@data[,object@variable])[2]
   z$npoints = length(object@observed@data[,object@variable])
   z$area = paste(length(object@predicted[,object@variable]) * object@predicted@grid@cellsize[1] * object@predicted@grid@cellsize[2])
   prj = plotKML::parse_proj4(p4s = proj4string(object@observed), params = as.list("\\+proj="))
      if(prj=="longlat")  {
          areaunits = "square-arcdegrees"
      } 
      else {
          areaunits = "square-m"
      }
   z$area.units = areaunits
#  z$cell.size = object@predicted@grid@cellsize
   z$covariates = all.vars(object@glm$terms)[-1]
   z$family = object@glm$family$family
   z$link = object@glm$family$link
   RMSE <- sqrt(mean((object@validation$var1.pred-object@validation$observed)^2))
   z$RMSE = signif(RMSE, 4)
   tvar <- 1-var(object@validation$residual, na.rm=T)/var(object@validation$observed, na.rm=T)
   z$tvar = signif(tvar*100, 3)
   asint <- as.integer(na.omit(round(object@predicted$var1.pred/(RMSE*.5))))
   tmp <- tempfile()
   save(asint, file=tmp, compress="gzip")
   z$npixels = length(object@predicted[,object@variable])
   # breaks:
   linkfun <- object@glm$family$linkfun
   invfun <- object@glm$family$linkinv
   xz <- range(linkfun(object@predicted@data[,object@variable]), na.rm = TRUE, finite = TRUE)
   xc <- cut(linkfun(object@predicted@data[,object@variable]), breaks = seq(xz[1], xz[2], by=RMSE/2), include.lowest = TRUE)
   z$breaks = invfun(seq(xz[1], xz[2], by=RMSE/2))
   z$bonds = summary(xc)
   z$Bytes = file.info(tmp)$size
   z$compress = "gzip"
   return(z)
})


## Summary for an object of type SpatialPredictions:
setMethod("show", signature(object = "SpatialPredictions"), function(object){
  require(plotKML)
  
  cat("  Variable           :", object@variable, "\n")
  cat("  Minium value       :", range(object@observed@data[,object@variable])[1], "\n")
  cat("  Maximum value      :", range(object@observed@data[,object@variable])[2], "\n")  
  cat("  Size               :", length(object@observed@data[,object@variable]), "\n")  
  # check the projection system:
  Tarea <- length(object@predicted[,object@variable]) * object@predicted@grid@cellsize[1] * object@predicted@grid@cellsize[2]
  prj = plotKML::parse_proj4(p4s = proj4string(object@observed), params = as.list("\\+proj="))
      if(prj=="longlat")  {
          areaunits = "square-arcdegrees"
          lengthunits = "arcdegrees" 
      } 
      else {
          areaunits = "square-m"
          lengthunits = "m"      
      }
  cat("  Total area         :", Tarea, "\n")
  cat("  Total area (units) :", areaunits, "\n")
  cat("  Resolution (x)     :", object@predicted@grid@cellsize[1], "\n")
  cat("  Resolution (y)     :", object@predicted@grid@cellsize[2], "\n")
  cat("  Resolution (units) :", lengthunits, "\n")
  cat("  GLM call formula   :", deparse(object@glm$call$formula), "\n")
  cat("  Family             :", object@glm$family$family, "\n")  
  cat("  Link function      :", object@glm$family$link, "\n")    
  cat("  Vgm model          :", paste(object@vgmModel$model[2]), "\n")
  cat("  Nugget (residual)  :", signif(object@vgmModel$psill[1], 3), "\n")
  cat("  Sill (residual)    :", signif(object@vgmModel$psill[2], 3), "\n")
  cat("  Range (residual)   :", signif(object@vgmModel$range[2], 3), "\n")
  # RMSE at validation points:
  RMSE <- sqrt(mean((object@validation$var1.pred-object@validation$observed)^2))
  cat("  RMSE (validation)  :", signif(RMSE, 4), "\n")
  tvar <- 1-var(object@validation$residual, na.rm=T)/var(object@validation$observed, na.rm=T)
  cat(paste("  Var explained      : ", signif(tvar*100, 3), "% \n", sep=""))
  # Effective bytes:
  asint <- as.integer(na.omit(round(object@predicted$var1.pred/(RMSE*.5))))
  tmp <- tempfile()
  save(asint, file=tmp, compress="gzip");
  cat("  Effective bytes    :", file.info(tmp)$size, "\n")
  cat("  Compression method :", "gzip", "\n")
})

# end of script;