# Purpose        : Derive Spatial Predictive Components for a list of grids;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : Not recommended for large grids;


setMethod("spc", signature(obj = "SpatialPixelsDataFrame", formulaString = "formula"), function(obj, formulaString, scale. = TRUE, silent = FALSE, ...){

  # formula string:
  if(missing(formulaString)) {
     formulaString <- as.formula(paste("~", paste(out@layernames, collapse="+")))
  }
  vars = all.vars(formulaString)
  if(length(vars)< 2){
    stop("At least two covarites required to run Principal Component Analysis")
  }
  obj@data <- obj@data[,vars]

  # print warning:
  if(silent==FALSE){
  if(nrow(obj)>10e6){
    warning('Operation not recommended for large grids', immediate. = TRUE)
  }}
  
  # convert every factor to indicators:
  for(j in 1:length(vars)){
    if(is.factor(obj@data[,vars[j]])){
      # remove classes without pixels:
      obj@data[,vars[j]] <- as.factor(paste(obj@data[,vars[j]]))
      ln <- levels(obj@data[,vars[j]])
      for(k in 1:length(ln)){
        vn <- paste(vars[j], k, sep="_")
        obj@data[,vn] <- ifelse(obj@data[,vars[j]]==ln[k], 1, 0)
      }
    message(paste("Converting", vars[j], "to indicators..."))
    }
  } 
  varsn = names(obj)[which(!sapply(obj@data, is.factor))]
  obj@data <- obj@data[,varsn]

  ## filter the missing values:
  if(scale. == TRUE){
    x <- scale(obj@data) 
    x[is.na(x)] <- 0
    x <- as.data.frame(x)
    sd.l <- lapply(x, sd)
    x0 <- sd.l==0
    if(any(x0)){
      message(paste("Columns with zero variance removed:", names(x)[which(x0)]), immediate. = TRUE)
      formulaString.f = as.formula(paste("~", paste(varsn[-which(x0)], collapse="+")))
      ## principal component analysis:
      pcs <- prcomp(formula=formulaString.f, x)
    } else {
      pcs <- prcomp(formula=formulaString, x)
    }
  } else {
    pcs <- prcomp(formula=formulaString, obj@data) 
  }

  ## copy values: 
  obj@data <- as.data.frame(pcs$x)
  proj4string(obj) <- obj@proj4string
  if(silent==FALSE){
    message(paste("Converting covariates to principal components..."))
    summary(pcs)
  }
 
  pcs <- new("SpatialComponents", predicted = obj, pca = pcs[-which(names(pcs)=="x")])
  return(pcs)

}) 


setMethod("spc", signature(obj = "list", formulaString = "list"), function(obj, formulaString, scale. = TRUE, silent = FALSE, ...){
    if(!length(obj)==length(formulaString)){
      stop("'obj' and 'formulaString' lists of same size expected")
    }
    
    pcs.l <- list(NULL)
    for(l in 1:length(obj)){
      pcs.l[[l]] <- spc(obj = obj[[l]], formulaString = formulaString[[l]], ...)
    }
    
    return(pcs.l)
})


# end of script;