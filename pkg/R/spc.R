# Purpose        : Derive Spatial Predictive Components for a list of grids;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : Not recommended for large grids;


setMethod("spc", signature(obj = "SpatialPixelsDataFrame", formulaString = "formula"), function(obj, formulaString, scale. = TRUE, silent = FALSE, ...){
 
  require(raster)

  # formula string:
  if(missing(formulaString)) {
     formulaString <- as.formula(paste("~", paste(out@layernames, collapse="+")))
  }
  vars = all.vars(formulaString)
  if(length(vars)< 2){
    stop("At least two covarites required to run Principal Component Analysis")
  }

  # print warning:
  if(nrow(obj)>10e6){
    warning('Operation not recommended for large grids', immediate. = TRUE)
  }
  
  # convert every factor to indicators:
  for(j in 1:length(vars)){
    if(is.factor(obj@data[,vars[j]])){
      ln <- levels(obj@data[,vars[j]])
      for(k in 1:length(ln)){
        vn <- paste(vars[j], k, sep="_")
        obj@data[,vn] <- ifelse(obj@data[,vars[j]]==ln[k], 1, 0)
      }
    message(paste("Converting", vars[j], "to indicators..."))
    }
  } 
  varsn = names(obj)[which(!sapply(obj@data, is.factor))]
  
  out <- obj[varsn]

  # filter the missing values:
  x <- scale(out@data) 
  x[is.na(x)] <- 0 
  
  pcs <- prcomp(formula=formulaString, scale=TRUE, as.data.frame(x))

  # copy values: 
  out@data <- as.data.frame(pcs$x)
  proj4string(out) <- obj@proj4string
  if(silent==FALSE){
    message(paste("Converting covariates to principal components..."))
    summary(pcs)
  }
 
  pcs <- new("SpatialComponents", predicted = out, pca = pcs[-which(names(pcs)=="x")])
  return(pcs)

}) 


# end of script;