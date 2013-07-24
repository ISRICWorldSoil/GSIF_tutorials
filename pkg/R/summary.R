# Purpose        : Get a summary of a SpatialPredictions class;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : ;


## Get a summary as a data frame:
setMethod("summary", signature(object = "SpatialPredictions"), function(object){
   
   z <- NULL
   z$variable = object@variable
   z$minium = range(object@observed@data[,object@variable])[1]
   z$maximum = range(object@observed@data[,object@variable])[2]
   z$npoints = length(object@observed@data[,object@variable])
   z$area = paste(length(object@predicted[,object@variable]) * object@predicted@grid@cellsize[1] * object@predicted@grid@cellsize[2])
   try(prj <- plotKML::parse_proj4(p4s = rgdal::CRSargs(CRS(proj4string(object@observed))), params = as.list("\\+proj=")))
      if(prj=="longlat")  {
          areaunits = "square-arcdegrees"
      } else {  
          areaunits = "square-m"
      }
   z$area.units = areaunits
#  z$cell.size = object@predicted@grid@cellsize
   if(any(class(object@regModel.summary)=="summary.glm")){
     z$covariates = all.vars(object@regModel.summary$terms)[-1]
     z$family = object@regModel.summary$family$family
   } else {
     z$covariates = NA
     z$family = NA
   }
   RMSE <- sqrt(mean((object@validation$var1.pred-object@validation$observed)^2))
   z$RMSE = signif(RMSE, 4)
   ## if cross-validation results are available:
   if(is.na(RMSE)|is.null(RMSE)){
    z$tvar = NA
    z$breaks = NA
    z$bonds = NA
    z$Bytes = NA
    z$compress = NA
   } else {
    tvar <- 1-var(object@validation$residual, na.rm=T)/var(object@validation$observed, na.rm=T)
    z$tvar = signif(tvar*100, 3)
    asint <- as.integer(na.omit(round(object@predicted@data[,object@variable]/(RMSE*.5), 0)))
    tmp <- tempfile()
    save(asint, file=tmp, compress="gzip")
    # breaks:
    xz <- range(object@predicted@data[,object@variable], na.rm = TRUE, finite = TRUE)
    xc <- cut(object@predicted@data[,object@variable], breaks = seq(xz[1], xz[2], by=RMSE/2), include.lowest = TRUE)
    z$breaks = seq(xz[1], xz[2], by=RMSE/2)
    z$bonds = summary(xc)
    z$Bytes = file.info(tmp)$size
    z$compress = "gzip"
   }
   z$npixels = length(object@predicted[,object@variable])
   
   return(z)
})


## Summary for an object of type SpatialPredictions:
setMethod("show", signature(object = "SpatialPredictions"), function(object){
  
  cat("  Variable           :", object@variable, "\n")
  cat("  Minium value       :", range(object@observed@data[,object@variable])[1], "\n")
  cat("  Maximum value      :", range(object@observed@data[,object@variable])[2], "\n")  
  cat("  Size               :", length(object@observed@data[,object@variable]), "\n")  
  # check the projection system:
  Tarea <- length(object@predicted[,object@variable]) * object@predicted@grid@cellsize[1] * object@predicted@grid@cellsize[2]
  try(prj <- plotKML::parse_proj4(p4s = rgdal::CRSargs(CRS(proj4string(object@observed))), params = as.list("\\+proj=")))
      if(prj=="longlat")  {
          areaunits = "square-arcdegrees"
          lengthunits = "arcdegrees" 
      } else {
          areaunits = "square-m"
          lengthunits = "m"      
      }
  cat("  Total area         :", Tarea, "\n")
  cat("  Total area (units) :", areaunits, "\n")
  cat("  Resolution (x)     :", object@predicted@grid@cellsize[1], "\n")
  cat("  Resolution (y)     :", object@predicted@grid@cellsize[2], "\n")
  cat("  Resolution (units) :", lengthunits, "\n")
  if(any(class(object@regModel.summary)=="summary.glm")){
    cat("  GLM call formula   :", deparse(object@regModel.summary$call$formula), "\n")
    cat("  Family             :", object@regModel.summary$family$family, "\n")  
    cat("  Link function      :", object@regModel.summary$family$link, "\n")    
  }
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
  asint <- as.integer(na.omit(round(object@predicted@data[,object@variable]/(RMSE*.5), 0)))
  tmp <- tempfile()
  save(asint, file=tmp, compress="gzip");
  cat("  Effective bytes    :", file.info(tmp)$size, "\n")
  cat("  Compression method :", "gzip", "\n")
})

# end of script;