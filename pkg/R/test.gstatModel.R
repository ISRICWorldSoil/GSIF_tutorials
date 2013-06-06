# Purpose        : evaluate/test gstatModel for model diagnostics;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : ;

setMethod("test.gstatModel", signature(observations = "SpatialPointsDataFrame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, Ns, predictionLocations = covariates, ...){
                            
  ## subset:  
  Nmin = round(20 + length(formulaString)*10)
  Nmax = round(nrow(observations))
  if(missing(Ns)){
    ss = round((runif(10)*sqrt(Nmax-Nmin))^2+Nmin)
    Ns <- sort(c(Nmin, ss, Nmax))
  }
  if(any(Ns>Nmax)){
    stop("'Ns' argument contains number larger than size of data")
  }
  
  ## target variable:
  tv = all.vars(formulaString)[1]
  svar = var(observations@data[,tv], na.rm=T)
  
  m.l <- list(NULL)  
  tvar.l <- as.list(rep(NA, length(Ns)))
  s.l <- as.list(rep(NA, length(Ns))) 
  p.l <- list(NULL)
  ctime <- as.list(rep(NA, length(Ns)))
  
  pb <- txtProgressBar(min=0, max=length(Ns), style=3)
  for(j in 1:length(Ns)){
    ## fit models:
    suppressWarnings(suppressMessages( try(m.l[[j]] <- fit.gstatModel(observations = observations[sample(1:nrow(observations), Ns[j], replace=FALSE),], formulaString = formulaString, covariates = covariates, ...))))
    ## validate models:
    suppressWarnings(suppressMessages(try(cv <- validate(m.l[[j]]))))
    ## variance explained:
    try(tvar.l[[j]] <- 1-var(cv[[1]]$residual, na.rm=TRUE)/svar)
    ## failures:
    try(s.l[[j]] <- sum(cv[[1]]$zscore > 3 | cv[[1]]$zscore < -3 ) )
    ## test predictions:
    suppressWarnings(suppressMessages(try( ctime[[j]] <- system.time( p.l[[j]] <- predict(m.l[[j]], predictionLocations = predictionLocations, mask.extra=TRUE, nfold = 0, debug.level = 0))[[1]]) ))
    setTxtProgressBar(pb, j) 
  }
  close(pb)
  
  out <- data.frame(samples = Ns, var.explained = unlist(tvar.l), pred.sec = unlist(ctime), failures = unlist(s.l))
    
  return(out)

})

# end of script;
