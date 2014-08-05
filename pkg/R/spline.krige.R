# Purpose        : Speading up ordinary kriging (e.g. of the regression residuals);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : this function is ONLY useful for highly clustered point data sets;

spline.krige <- function(formula, locations, newdata, newlocs=NULL, model, te=as.vector(newdata@bbox), file.name, silent=FALSE, t_cellsize=newdata@grid@cellsize[1], optN=20, quant.nndist=.5, nmax=30, ...){
  if(!class(locations)=="SpatialPointsDataFrame"){
    stop("Object 'locations' of class 'SpatialPointsDataFrame' expected")
  }
  if(!class(newdata)=="SpatialPixelsDataFrame"){
    stop("Object 'newdata' of class 'SpatialPixelsDataFrame' expected")
  }
  if(is.null(newlocs)){ 
     newlocs <- resample.grid(locations, newdata, silent=silent, t_cellsize=t_cellsize, quant.nndist=quant.nndist)$newlocs
  }
  if(silent==FALSE){
    message("Predicting at variable grid...")
  }
  s_te <- as.vector(newdata@bbox)
  if(missing(formula)){
    formula <- as.formula(paste(names(locations)[1], 1, sep="~"))
  }
  class(model) <- c("variogramModel", "data.frame")
  ok <- krige(formula, locations=locations[!is.na(locations@data[,1]),], newdata=newlocs, model=model, nmax=nmax, debug.level=-1, ...)
  ## write points to a shape file:
  tmp <- list(NULL)
  tmp.out <- list(NULL)
  for(k in 1:ncol(ok@data)){
    tmp[[k]] <- set.file.extension(tempfile(), ".shp")
    writeOGR(ok[k], names(ok)[k], dsn=tmp[[k]], "ESRI Shapefile")
    if(missing(file.name)){
      tmp.out[[k]] <- tempfile()
    } else {
      tmp.out[[k]] <- paste(file.name, k, sep="_")
    }
    ## point to grid (spline interpolation):
    rsaga.geoprocessor(lib="grid_spline", module=6, param=list(SHAPES=tmp[[k]], FIELD=0, TARGET=0, NPMIN=3, NPMAX=nmax, USER_XMIN=te[1]+t_cellsize/2, USER_XMAX=te[3]+t_cellsize/2, USER_YMIN=te[2]-t_cellsize/2, USER_YMAX=te[4]-t_cellsize/2, USER_SIZE=t_cellsize, USER_GRID=set.file.extension(tmp.out[[k]], ".sgrd")), show.output.on.console = FALSE)
    ## fill gaps:
    rsaga.geoprocessor(lib="grid_tools", module=7, param=list(INPUT=set.file.extension(tmp.out[[k]], ".sgrd"), RESULT=set.file.extension(tmp.out[[k]], ".sgrd")), show.output.on.console = FALSE)
    if(!all(te==s_te)|t_cellsize<newdata@grid@cellsize[1]){
      if(silent==FALSE){ message(paste("Resampling band", k, "to the target resolution and extent...")) }
      if(t_cellsize<newdata@grid@cellsize[1]){
        rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT=set.file.extension(tmp.out[[k]], ".sgrd"), TARGET=0, SCALE_DOWN_METHOD=4, USER_XMIN=s_te[1]+t_cellsize/2, USER_XMAX=s_te[3]+t_cellsize/2, USER_YMIN=s_te[2]-t_cellsize/2, USER_YMAX=s_te[4]-t_cellsize/2, USER_SIZE=t_cellsize, USER_GRID=set.file.extension(tmp.out[[k]], ".sgrd")), show.output.on.console=FALSE)
      } else {
        ## resample:
        rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT=set.file.extension(tmp.out[[k]], ".sgrd"), TARGET=0, SCALE_DOWN_METHOD=0, SCALE_UP_METHOD=0, USER_XMIN=s_te[1]+t_cellsize/2, USER_XMAX=s_te[3]+t_cellsize/2, USER_YMIN=s_te[2]-t_cellsize/2, USER_YMAX=s_te[4]-t_cellsize/2, USER_SIZE=t_cellsize, USER_GRID=set.file.extension(tmp.out[[k]], ".sgrd")), show.output.on.console=FALSE)
      }
    }
    if(missing(file.name)){
      if(k==1){
        out <- readGDAL(set.file.extension(tmp.out[[k]], ".sdat"), silent=TRUE)
        proj4string(out) <- newdata@proj4string
        names(out) <- names(ok)[k]
      } else {
        out@data[,names(ok)[k]] <- readGDAL(set.file.extension(tmp.out[[k]], ".sdat"), silent=TRUE)$band1
      }
      unlink(paste0(tmp.out[[k]],".*"))
    } else {
      if(silent==FALSE){ message(paste("Created output file:", tmp.out[[k]])) }
    }
  }
  return(out)
}

## resample using variable sampling intensity:
resample.grid <- function(locations, newdata, silent=FALSE, t_cellsize, optN, quant.nndist, nstrata=4){
    if(silent==FALSE){
      message("Deriving density map...")
    }
    require(spatstat)
    require(maptools)
    ## derive density map:
    W <- as.matrix(newdata[1])
    W <- ifelse(is.na(W), FALSE, TRUE)
    suppressWarnings( locs.ppp <- ppp(x=locations@coords[,1], y=locations@coords[,2], xrange=newdata@bbox[1,], yrange=newdata@bbox[2,], mask=t(W)[ncol(W):1,]) )
    dist.locs <- nndist(locs.ppp)
    dmap <- maptools::as.SpatialGridDataFrame.im(density(locs.ppp, sigma=quantile(dist.locs, quant.nndist)))
    dmap.max <- max(dmap@data[,1], na.rm=TRUE)
    dmap@data[,1] <- signif(dmap@data[,1]/dmap.max, 3)
    ## TH: not sure if here is better to use quantiles or a regular split?
    breaks.d <- seq(0, 1, by=1/(nstrata))
    #breaks.d <- quantile(dmap@data[,1], seq(0, 1, by=1/(nstrata+1)), na.rm=TRUE)
    if(sd(dmap@data[,1])==0){ stop("Density map shows no variance. See '?resample.grid' for more information.") }
    dmap$strata <- cut(x=dmap@data[,1], breaks=breaks.d, include.lowest=TRUE, labels=paste0("L", 1:nstrata))
    proj4string(dmap) = locations@proj4string
    ## regular sampling proportional to the sampling density (rule of thumb: one sampling point can be used to predict 'optN' grids):
    newlocs <- list(NULL)
    for(i in 1:length(levels(dmap$strata))){
      im <- dmap[dmap$strata==paste0("L",i),"strata"]
      im <- as(im, "SpatialPixelsDataFrame")
      if(i==length(levels(dmap$strata))){ 
        Ns <- round(sum(!is.na(im$strata)) * newdata@grid@cellsize[1]/t_cellsize)
      } else {
        ov <- over(locations, im)
        if(sum(!is.na(ov$strata))==0){ 
          Ns <- optN
        } else {
          Ns <- round(sum(!is.na(ov$strata)) * optN)
        }
      }
      newlocs[[i]] <- sp::spsample(im, type="regular", n=Ns)
    }
    newlocs <- do.call(rbind, newlocs)
    if(silent==FALSE){
      message(paste("Generated:", length(newlocs), "prediction locations."))
    }
    return(list(newlocs=newlocs, density=dmap))
}

## end of script;