# Purpose        : Fit/predict distribution of soil types (memberships);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Bas Kempen (bas.kempen@wur.nl); 
# Dev Status     : Pre-Alpha
# Note           : if the regression model is difficult to fit, it might lead to artifacts;


# Fit a supervised fuzzy kmeans model and predict memberships:
setMethod("spfkm", signature(formulaString = "formula", observations = "SpatialPointsDataFrame", covariates = "SpatialPixelsDataFrame"), function(formulaString, observations, covariates, class.c = NULL, class.sd = NULL, fuzzy.e = 1.2){
  
  ## generate formula if missing:
  if(missing(formulaString)) {  
    formulaString <- as.formula(paste(names(observations)[1], "~", paste(names(covariates), collapse="+"), sep=""))
  }
  ## check the formula string:
  if(!is.formula(formulaString)){
      stop("'formulaString' object of class 'formula' required")
  }
  
  ## selected variables:
  tv = all.vars(formulaString)[1]
  sel = names(covariates) %in% all.vars(formulaString)[-1]
  if(all(sel==FALSE)|length(sel)==0){
      stop("None of the covariates in the 'formulaString' matches the column names in the 'covariates' object")
  }
 
  ## if available, use class centres:
  check_tc <- !is.null(class.c)&!is.null(class.sd)
  if(check_tc){
    if(!class(class.c)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    if(!class(class.sd)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    mout = list(NULL)
  }
  ## otherwise, estimate class centres using the multinomial logistic regression:
  else {
    message("Trying to estimate the class centres using the 'multinom' method...")
    ## multinomial logistic regression:
    rout <- spmultinom(formulaString=formulaString, observations, covariates, class.stats=TRUE, predict.probs=FALSE)
    mout = rout$model
    if(length(unique(rout$fit))<2){ stop("Predictions resulted in <2 classes. See ?multinom for more info") }
    class.c = rout$class.c
    class.sd = rout$class.sd
  }
  
  cl <- as.list(row.names(class.c))
  dsf <- NULL
  ## derive distances in feature space:
  for(c in unlist(cl)){
      dsf[[c]] <- data.frame(lapply(names(covariates)[sel], FUN=function(x){rep(NA, length(covariates@data[,1]))}))
      names(dsf[[c]]) <- names(covariates)[sel]
      for(j in names(covariates)[sel]){
         dsf[[c]][,j] <- ((covariates@data[,j]-class.c[c,j])/class.sd[c,j])^2
      }
  }
  ## sum up distances per class:
  ds <- NULL
  ds <- lapply(dsf, FUN=function(x){sqrt(rowSums(x, na.rm=TRUE, dims=1))})
  names(ds) <- unlist(cl)
  ds <- data.frame(ds)
  ## total sum:
  tt <- rowSums(ds^(-2/(fuzzy.e-1)), na.rm=TRUE, dims=1)
  ## derive the fuzzy membership:
  mm <- covariates[1]
  for(c in unlist(cl)){
    mm@data[,c] <- (ds[,c]^(-2/(fuzzy.e-1))/tt)
  }
  mm@data[,names(covariates)[1]] <- NULL
  
  ## Derive the dominant class:
  maxm <- sapply(data.frame(t(as.matrix(mm@data))), FUN=function(x){max(x, na.rm=TRUE)})
  ## class having the highest membership
  cout <- NULL
  for(c in unlist(cl)){
       cout[which(mm@data[,c] == maxm)] <- c
  }
  cout <- as.factor(cout)
  
  ## construct a map:
  pm <- covariates[1]
  pm@data[,tv] <- cout
  pm@data[,names(covariates)[1]] <- NULL

  ## overlay observations and covariates:
  ov <- over(observations, pm)
  sel.c <- !is.na(ov[,tv]) & !is.na(observations@data[,tv])

  ## kappa statistics:
  require(mda)
  require(psych)
  cf <- confusion(ov[sel.c,tv], as.character(observations@data[sel.c,tv]))
  ## remove missing classes:
  a = attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
  b = attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
  c.kappa = cohen.kappa(cf[a,b])
  message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))  
  
  ## create the output object:
  out <- new("SpatialMemberships", predicted = pm, model = mout, mu = mm, class.c = class.c, class.sd = class.sd, confusion = cf)
  return(out)

})

# end of script;