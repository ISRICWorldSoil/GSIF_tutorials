# Purpose        : evaluate/test gstatModel for model diagnostics;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : ;

setMethod("test.gstatModel", signature(observations = "SpatialPointsDataFrame", formulaString = "formula", covariates = "SpatialPixelsDataFrame"), function(observations, formulaString, covariates, Ns, predictionLocations = covariates, save.predictions = FALSE, ...){
                            
  ## derive Ns if not available:  
  Nmax = round(nrow(observations))
  if(missing(Ns)){
    Nmin = round(20 + length(formulaString)*10)
    ss = round((runif(10)*sqrt(Nmax-Nmin))^2+Nmin)
    Ns <- sort(c(Nmin, ss, Nmax))
  }
  if(any(Ns>Nmax)){
    stop("'Ns' argument contains number larger than size of data")
  }
  
  ## empty lists:
  m.l <- list(NULL)  
  tvar.l <- as.list(rep(NA, length(Ns)))
  s.l <- as.list(rep(NA, length(Ns))) 
  p.l <- list(NULL)
  ctime <- as.list(rep(NA, length(Ns)))
  
  message(paste("Running model fitting, cross-validation and predictions for", length(Ns), "sampling intensities using N-fold cross-validation..."))
  pb <- txtProgressBar(min=0, max=length(Ns), style=3)
  for(j in 1:length(Ns)){
    ## fit models:
    suppressWarnings(suppressMessages( try(m.l[[j]] <- fit.gstatModel(observations = observations[sample(1:nrow(observations), Ns[j], replace=FALSE),], formulaString = formulaString, covariates = covariates, ...))))
    ## validate models:
    suppressWarnings(suppressMessages(try(cv <- validate(m.l[[j]]))))
    if(is.list(cv)){
      ## variance explained:
      svar = var(cv[[1]]$observed, na.rm=T)
      try(tvar.l[[j]] <- 1-var(cv[[1]]$residual, na.rm=TRUE)/svar)
      ## failures:
      try(s.l[[j]] <- sum(cv[[1]]$zscore^2 > 1.5 | abs(cv[[1]]$residual) > 3*sqrt(svar), na.rm = TRUE))
    }
    
    ## test predictions:
    suppressWarnings(suppressMessages(try( ctime[[j]] <- system.time( p.l[[j]] <- predict(m.l[[j]], predictionLocations = predictionLocations, mask.extra=TRUE, nfold = 0, debug.level = 0))[[1]]) ))
    setTxtProgressBar(pb, j) 
  }
  close(pb)
  cat(j, "\r")
  flush.console() 
  
  out <- data.frame(samples = Ns, var.explained = unlist(tvar.l), pred.sec = unlist(ctime), failures = unlist(s.l))
  
  if(save.predictions == TRUE){
    names(p.l) <- paste("Ns =", as.character(Ns))
    out <- list(performance = out, predictions = p.l)
  } else{
    out <- list(performance = out, predictions = NULL)  
  }
    
  return(out)

})

# end of script;
