# Purpose        : Fit/predict distribution of soil types (memberships);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : if the regression model is difficult to fit, it might lead to artifacts;


# Fit a supervised fuzzy kmeans model and predict memberships:
fkmeans <- function(formulaString, observations, covariates, class.c  = NULL, class.sd = NULL, fuzzy.e = 1.2){
   
  # check input classes first:
  if(!class(observations)=="SpatialPointsDataFrame"){
      stop("Class of type 'SpatialPointsDataFrame' required for 'observations'")
  }
  covariates <- as(covariates, "SpatialPixelsDataFrame")
  
  # generate formula if missing:
  if(missing(formulaString)) {  
    formulaString <- as.formula(paste(names(observations)[1], "~", paste(names(covariates), collapse="+"), sep=""))
  }
  # check the formula string:
  if(!is.formula(formulaString)){
      stop("'formulaString' object of class 'formula' required")
  }
  
  # selected variables:
  tv = all.vars(formulaString)[1]
  sel = names(covariates) %in% all.vars(formulaString)[-1]
  if(all(sel==FALSE)|length(sel)==0){
      stop("None of the covariates in the 'formulaString' matches the column names in the'covariates' object")
  }

  # overlay observations and covariates:
  index <- overlay(covariates[sel], observations)  
  ov <- cbind(data.frame(observations[tv]), covariates@data[index,sel])
   
  # if available, use class centres:
  if(!is.null(class.c)&!is.null(class.sd)){
    if(!class(class.c)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    if(!class(class.sd)=="matrix"){ stop("Object of type 'matrix' with column names for covariates and row names correspodning to the class names required") }
    mout = list(NULL)
  }
  # otherwise, estimate class centres using the multinomial logistic regression:
  else {
    # multinomial logistic regression:
    require(nnet)
    message("Fitting a multinomial logistic regression model...")
    mout <- multinom(formulaString, ov)
    cout <- as.factor(predict(mout, newdata=covariates))
    class(mout) = "list"
    # estimate class centres using the results of multinom:
    ca <- aggregate(covariates@data[,sel], by=list(cout), FUN="mean")
    class.c <- as.matrix(ca[-1]); attr(class.c, "dimnames")[[1]] <- ca[,1]
    ca <- aggregate(covariates@data[,sel], by=list(cout), FUN="sd")
    class.sd <- as.matrix(ca[-1]); attr(class.sd, "dimnames")[[1]] <- ca[,1]
    # mask out classes that result in NA:
    for(j in length(nrow(class.c))){ if(any(is.na(class.c[j,]))|any(is.na(class.sd[j,]))){ class.c <- class.c[-j,]; class.sd <- class.sd[-j,] } }
  }
  
  cl <- as.list(row.names(class.c))
  dsf <- NULL
  # derive distances in feature space:
  for(c in unlist(cl)){
      dsf[[c]] <- data.frame(lapply(names(covariates)[sel], FUN=function(x){rep(NA, length(cout))}))
      names(dsf[[c]]) <- names(covariates)[sel]
      for(j in names(covariates)[sel]){
         dsf[[c]][,j] <- ((covariates@data[,j]-class.c[c,j])/class.sd[c,j])^2
      }
  }
  # sum up distances per class:
  ds <- NULL
  ds <- lapply(dsf, FUN=function(x){sqrt(rowSums(x, na.rm=TRUE, dims=1))})
  names(ds) <- unlist(cl)
  ds <- data.frame(ds)
  # total sum:
  tt <- rowSums(ds^(-2/(fuzzy.e-1)), na.rm=TRUE, dims=1)
  # derive the fuzzy membership:
  mm <- covariates[1]
  for(c in unlist(cl)){
    mm@data[,c] <- (ds[,c]^(-2/(fuzzy.e-1))/tt)
  }
  mm@data[,names(covariates)[1]] <- NULL
  
  # if required, derive the dominant class:
  if(!is.null(class.c)&!is.null(class.sd)){
    # highest membership:
    maxm <- sapply(data.frame(t(as.matrix(mm@data))), FUN=function(x){max(x, na.rm=TRUE)})
    # class having the highest membership
    cout <- NULL
    for(c in unlist(cl)){
       cout[which(mm@data[,c] == maxm)] <- c
    }
  }
  
  # kappa statistics:
  require(mda)
  cf <- confusion(as.character(ov[which(!is.na(index)),tv]), cout[index][which(!is.na(index))])
  message(paste("Estimated prediction error:", signif(attr(cf, "error"), 4)))
  
  # construct a map:
  pm <- covariates[1]
  pm@data[,tv] <- cout
  pm@data[,names(covariates)[1]] <- NULL
  
  # create the output object:
  out <- new("SpatialMemberships", predicted = pm, model = mout, mu = mm, class.c = class.c, class.sd = class.sd, confusion = cf)
  return(out)

}


# end of script;