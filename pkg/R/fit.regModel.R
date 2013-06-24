# Purpose        : Fit a 2D or 3D regression model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Bas Kempen (bas.kempen@wur.nl) and Gerard B.M. Heuvelink (gerard.heuvelink@wur.nl); 
# Dev Status     : Pre-Alpha
# Note           : Regression families considered spatial GLMs, CART, and Hieararchical Bayes methods;


## Fit a GLM to spatial data:
setMethod("fit.regModel", signature(formulaString = "formula", rmatrix = "data.frame", predictionDomain = "SpatialPixelsDataFrame", method = "character"), function(formulaString, rmatrix, predictionDomain, method = list("GLM", "rpart", "randomForest", "quantregForest")[[1]], dimensions = NULL, family=gaussian, stepwise=TRUE, rvgm, ...){

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

  ## check if the method exists:
  if(!any(method %in% list("GLM", "rpart", "randomForest", "quantregForest"))){ stop(paste(method, "method not available.")) }
    
  if(method == "GLM"){  
    ## fit/filter the regression model:
    message("Fitting a GLM...")
    rgm <- glm(formulaString, data=rmatrix, family=family)
    if(stepwise == TRUE){
      rgm <- step(rgm, trace = 0)
    }
   
    ## mask out the missing values:
    if(any(names(rgm) == "na.action")){  rmatrix <- rmatrix[-rgm$na.action,] }
    ## extract the response residuals: [http://stackoverflow.com/questions/2531489/understanding-glmresiduals-and-residglm] 
    rmatrix$residual <- resid(rgm, type="response")
  }
  
  if(method == "rpart"){
    ## fit/filter the regression model:
    message("Fitting a regression tree model...")
    rgm <- rpart(formulaString, data=rmatrix)
    if(stepwise == TRUE){
      ## TH: "A good choice of cp for pruning is often the leftmost value for which the mean lies below the horizontal line"
      ## BK: determine row in complexity table with smallest xerror:
      minerror <- min(seq_along(rgm$cptable[,4L])[rgm$cptable[,4L] == min(rgm$cptable[,4L])])
      ## BK: select starting value for evaluation of xerror:
      xerr <- rgm$cptable[1L,4L]
      ## BK: compute 1-SE value:
      dum <- (rgm$cptable[,4L] + rgm$cptable[,5L])[minerror]
      ## BK determine row in complexity table for which xerror is smaller than 1-SE:
      i <- 0
      while (xerr > dum && i <= nrow(rgm$cptable)) {
        i <- i+1L  
        xerr <- rgm$cptable[i,4L]
      }
      # BK: obtain cp parameter and number of splits for selected row:
      cpar <- rgm$cptable[i,1L]
      nsplit <- rgm$cptable[i,2L]
      message(paste("Estimated Complexity Parameter (for prunning):", signif(cpar, 4)))
      rgm <- prune(rgm, cp=cpar)
    }  
    ## extract the residuals:
    if(any(names(rgm) == "na.action")){  rmatrix <- rmatrix[-rgm$na.action,] } 
    rmatrix$residual <- resid(rgm)  
  }
  
  if(method == "randomForest"|method == "quantregForest"){
    ## fit/filter the regression model:
    message("Fitting a randomForest model...")
    ## NA not permitted in response:
    rmatrix <- rmatrix[!is.na(rmatrix[,tv]),]
    if(method == "randomForest"){
      rgm <- randomForest(formulaString, data=rmatrix, na.action=na.pass)
    } else {
      ## TH: the quantreg package developed by Nicolai Meinshausen <meinshausen@stats.ox.ac.uk> is slower but more flexible      
      rgm <- quantregForest(y=eval(formulaString[[2]], rmatrix), x=rmatrix[,all.vars(formulaString)[-1]])
      attr(rgm$y, "name") <- tv  
    }
    ## extract the residuals:
    rmatrix$residual <- rgm$predicted - rgm$y
  }
  
  ## TH: here we will add more regression models...
  
  ## test the normality of residuals:
  require(stats)
  if(length(rmatrix$residual)>4999){
    # subset residuals if necessary...
    x = rmatrix$residual[runif(length(rmatrix$residual))<4000/length(rmatrix$residual)]
  } else {
    x = rmatrix$residual
  }
  st = shapiro.test(x)
  if(st$p.value < 0.05|is.na(st$p.value)){
    ## try second test:
    require(nortest)
      at = ad.test(x)
      if(at$p.value < 0.05|is.na(at$p.value)){
        warning(paste(st$method, "and", at$method, "report probability of < .05 indicating lack of normal distribution for residuals"), call. = FALSE, immediate. = TRUE)
    }
  }

  ## Fit variogram 2D or 3D:
  if(missing(rvgm)){
    if(dimensions == "2D"){ 
      message("Fitting a 2D variogram...")
      rvgm <- fit.vgmModel(residual ~ 1, rmatrix = rmatrix, predictionDomain = predictionDomain, dimensions = "2D", ...) 
    }
    if(dimensions == "3D"){ 
      message("Fitting a 3D variogram...")
      rvgm <- fit.vgmModel(residual ~ 1, rmatrix = rmatrix, predictionDomain = predictionDomain, dimensions = "3D", ...) 
    }
    } else {
       if(is.null(rvgm)){
         rvgm <- fit.vgmModel(residual ~ 1, rmatrix = rmatrix, predictionDomain = predictionDomain, dimensions = dimensions, vgmFun = "Nug", ...)
       } else {
         ## othewise copy the variogram submitted by the user:
         xyn = attr(predictionDomain@bbox, "dimnames")[[1]]
         ## create spatial points:
         coordinates(rmatrix) <- as.formula(paste("~", paste(xyn, collapse = "+"), sep=""))
         proj4string(rmatrix) = predictionDomain@proj4string
         observations = as(rmatrix, "SpatialPoints")
         rvgm <- list(vgm=rvgm, observations=observations)
       }
  }
  
  ## TH: refit the GLM using the GLS weights? apparently this is possible via the "nlme" package.
  
  ## save the output:
  message("Saving an object of class 'gstatModel'...")  
  rkm <- new("gstatModel", regModel = rgm, vgmModel = as.data.frame(rvgm[[1]]), sp = rvgm[[2]])
  return(rkm)

})

## end of script;