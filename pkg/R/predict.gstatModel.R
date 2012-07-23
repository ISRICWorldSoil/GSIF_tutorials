# Purpose        : Fit/predict a 2D or 3D regression-kriging model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : works only with linear models with normally distributed residuals;



# Fit a RK model and return an object of class "gsm":
setMethod("fit.gstatModel", signature(observations = "geosamples", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, methodid, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, ...){

  require(gstat)
  require(splines)
  require(plotKML)

  # generate formula if missing:
  if(missing(formulaString)) {  
    formulaString <- as.formula(paste(names(observations)[1], "~", paste(names(covariates), collapse="+"), sep=""))
  }
   
  # prepare regression matrix:
  ov <- overlay(x=covariates, y=observations, method=methodid, var.type = "numeric")
  if(nrow(ov)==0|is.null(ov$observedValue)) {
    warning("The overlay operations resulted in an empty set. Check 'methodid' column.")
  }
  # geostats only possible with numeric variables:
  ov$observedValue = as.numeric(ov$observedValue)
  
  # fit/filter the regression model:
  rgm <- glm(formulaString, data=ov, family=family, ...)
  if(stepwise == TRUE){
    rgm <- step(rgm, trace = 0)
  }
  
  # mask out the missing values:
  if(any(names(rgm) == "na.action")){  ov <- ov[-rgm$na.action,] }
  # extract the residuals in the transformed space:
  linkfun <- rgm$family$linkfun
  ov$residual <- linkfun(ov$observedValue) - rgm$linear.predictors
  # create a 3D object:
  coordinates(ov) <- ~ longitude + latitude + altitude
  proj4string(ov) = get("ref_CRS", envir = GSIF.opts)
  suppressWarnings(ov <- spTransform(ov, covariates@proj4string))
  
  # if the residual variogram is unknown:
  if(is.null(rvgm)){
    # estimate area extent:
    Range = sqrt(areaSpatialGrid(covariates))/3
    # estimate initial range in the vertical direction:
    dr <- abs(diff(range(ov@coords[,3], na.rm=TRUE)))/3
    a2 = 2*dr/Range

    # estimate anisotropy parameters:
    anis = c(0, 0, 0, 1, a2)    
    ivgm <- vgm(nugget=0, model=vgmFun, range=Range, psill=var(ov$residual), anis = anis)
    # fit the 3D variogram:
    try(rvgm <- gstat::fit.variogram(variogram(residual ~ 1, ov), ivgm, ...))
    if(diff(rvgm$range)==0|diff(rvgm$psill)==0){
      warning("Variogram shows no spatial dependence")    
    }
    ## TH: This is the most simple implementation - fitting of 3D variograms needs to be improved!

  }
  
  # save the fitted model:
  rkm <- new("gstatModel", regModel = rgm, vgmModel = as.data.frame(rvgm), sp = as(ov, "SpatialPoints"))
  return(rkm)  

})


## Fit a 'simple' RK model (2D points):
setMethod("fit.gstatModel", signature(observations = "SpatialPointsDataFrame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, family = gaussian, stepwise = TRUE, vgmFun = "Exp", rvgm = NULL, ...){

  require(gstat)
  
  # the function only works with 2D maps:
  if(length(attr(coordinates(observations), "dimnames")[[2]])>2){
    stop("Works only with 2D Spatial objects")
  }

  # generate formula if missing:
  if(missing(formulaString)) {  
    formulaString <- as.formula(paste(names(observations)[1], "~", paste(names(covariates), collapse="+"), sep=""))
  }
  
  index <- overlay(covariates, observations)
  sel <- !is.na(index)
  # all variables of interest:
  tv = all.vars(formulaString)[1]
  xyn = attr(coordinates(observations), "dimnames")[[2]]
  seln = names(covariates) %in% all.vars(formulaString)[-1]
  if(length(seln)==0){
      stop("None of the covariates in the 'formulaString' do not match the names in the'covariates' object")
  }
  x <- cbind(data.frame(observations[sel,tv]), covariates[index[sel],])
  # fit the regression model:
  rgm <- glm(formulaString, data=x, family=family, ...)
  if(stepwise == TRUE){
    rgm <- step(rgm, trace = 0)
  }
  # mask out the missing values:
  if(any(names(rgm) == "na.action")){   x <- x[-rgm$na.action,]  }
  # extract the residuals in the transformed space:
  linkfun <- rgm$family$linkfun
  x$residual <- linkfun(x[,tv]) - rgm$linear.predictors
  coordinates(x) <- as.formula(paste("~", paste(xyn, collapse = "+"), sep=""))
  proj4string(x) = observations@proj4string
  
  if(is.null(rvgm)){
    # fit the variogram:
    if(!is.na(proj4string(covariates))){
    if(!is.projected(covariates)){
      require(fossil)  # Haversine Formula for Great Circle distance
      p.1 <- matrix(c(covariates@bbox[1,1], covariates@bbox[1,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
      p.2 <- matrix(c(covariates@bbox[2,1], covariates@bbox[2,2]), ncol=2, dimnames=list(1,c("lon","lat")))  
      Range = deg.dist(lat1=p.1[,2], long1=p.1[,1], lat2=p.2[,2], long2=p.2[,1])/2
    } else{
      Range = sqrt(areaSpatialGrid(covariates))/2
    }}
    else{
      Range = sqrt(areaSpatialGrid(covariates))/2
    }

    ivgm <- vgm(nugget=0, model=vgmFun, range=Range, psill=var(x$residual))
    try(rvgm <- gstat::fit.variogram(variogram(residual ~ 1, x), ivgm, ...))
    if(diff(rvgm$range)==0|diff(rvgm$psill)==0){
      warning("Variogram shows no spatial dependence")    
    }
  }
  
  # save the fitted model:
  rkm <- new("gstatModel", regModel = rgm, sp = as(x, "SpatialPoints"), vgmModel = as.data.frame(rvgm))
  return(rkm)
})


# predict values using a RK model:
predict.gstatModel <- function(object, predictionLocations, nmin = 10, nmax = 30, debug.level = -1, method = c("KED", "RK")[1], nfold = 5, verbose = FALSE, nsim = 0, mask.extra = TRUE, ...){

  if(nsim<0|!is.numeric(nsim)){
   stop("To invoke conditional simulations set 'nsim' argument to a positive integer number")
  }
  nsim <- as.integer(nsim)

  require(gstat)
  require(raster)
  
  # check the projection system:
  if(is.na(proj4string(predictionLocations))){
    stop("proj4string for the 'predictionLocations' required")
  }
  if(is.na(proj4string(object@sp))){
    stop("proj4string for the 'gstatModel' object required")
  }
  # force SpatialPixels:
  predictionLocations <- as(predictionLocations, "SpatialPixelsDataFrame")

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
  
  # skip cross-validation in nfold = 0
  if(nfold==0){ 
    cv <- remove.duplicates(observed[variable])
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
        rk <- gstat::krige(formula=formString, locations=remove.duplicates(observed), newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, ...)
        # mask extrapolation areas:
        rk@data[,paste(variable, "svar", sep=".")] <- rk$var1.var / var(observed@data[,paste(variable, "link", sep=".")], na.rm=TRUE)
        if(mask.extra==TRUE){
          rk$var1.pred <- ifelse(rk@data[,paste(variable, "svar", sep=".")] > 1, NA, rk$var1.pred) 
        }
        # back-transform the values:
        rk@data[,variable] = invfun(rk$var1.pred)
        # skip cross-validation:
        if(nfold>0){      
          message(paste("Running ", nfold, "-fold cross validation...", sep=""))
          cv <- gstat::krige.cv(formString, locations=remove.duplicates(observed), model = vgmmodel, nfold = nfold, verbose = verbose)
          proj4string(cv) = observed@proj4string
        }
      }
      else {
        message(paste("Generating", nsim, "conditional simulations using the trend model (KED method)..."))
        rk <- gstat::krige(formString, locations=remove.duplicates(observed), newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = nsim, ...)
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
            cv.err <- boot::glm.diag(object@regModel)
            # cv.err <- boot::cv.glm(data=object@regModel$data, glmfit=object@regModel, K=nfold) 
            # TH: This one fails for unknow reason? 
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
          observed <- remove.duplicates(observed)
          rk <- gstat::krige(formString, locations=observed, newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, ...) 
          # sum regression and kriging and back-transform the values:
          rk@data[,variable] <- invfun(predictionLocations@data[,paste(variable, "glmfit", sep=".")] + rk@data[,"var1.pred"])
          rk@data[,paste(variable, "svar", sep=".")] <- ( predictionLocations$se.fit + rk@data[,"var1.var"] ) / var(observed@data[,paste(variable, "link", sep=".")], na.rm=TRUE)  ## This formula assumes that the trend and residuals are independent; which is probably not true
          if(nfold>0){
          formString <- as.formula(paste(paste(variable, "link", sep="."), "~", paste(variable, "glmfit", sep="."), sep=""))
          message(paste("Running ", nfold, "-fold cross validation...", sep=""))
          cv <- gstat::krige.cv(formString, locations=observed, model = vgmmodel, nfold = nfold, verbose = verbose)
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
          rk <- gstat::krige(formString, locations=remove.duplicates(observed), newdata=predictionLocations, model = vgmmodel, nmin = nmin, nmax = nmax, debug.level = debug.level, nsim = nsim, ...)
          # back-transform the values:
          rk@data[,i] = invfun(xsim + rk@data[,i])  # this is inexpensive but it assumes that the trend and residuals are independent!        
        }
        }
      }
      }
  }
   
  # save the output file:
  if(nsim == 0){
    rkp <- new("SpatialPredictions", variable = variable, observed = observed, glm = sum.glm, vgmModel = object@vgmModel, predicted = rk, validation = cv)
  } 
  else {
    t1 <- Line(matrix(c(rk@bbox[1,1],rk@bbox[1,2],mean(rk@bbox[2,]),mean(rk@bbox[2,])), ncol=2))
    transect <- SpatialLines(list(Lines(list(t1), ID="t")), observed@proj4string)
    rkp <- new("RasterBrickSimulations", variable = variable, sampled = transect, realizations = brick(rk))
  }
   
  return(rkp)

}

setMethod("predict", signature(object = "gstatModel"), predict.gstatModel)


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