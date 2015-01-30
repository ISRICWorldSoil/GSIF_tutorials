# Purpose        : Speading up ordinary kriging (e.g. of the regression residuals);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : this function is ONLY useful for highly clustered point data sets;

spline.krige <- function(formula, locations, newdata, newlocs=NULL, model, te=as.vector(newdata@bbox), file.name, silent=FALSE, t_cellsize=newdata@grid@cellsize[1], optN=20, quant.nndist=.5, nmax=30, predictOnly=FALSE, resample=TRUE, saga.env, saga.lib=c("grid_spline","grid_tools"), saga.module=c(4,0), ...){
  if(!class(locations)=="SpatialPointsDataFrame"){
    stop("Object 'locations' of class 'SpatialPointsDataFrame' expected")
  }
  if(!class(newdata)=="SpatialPixelsDataFrame"){
    stop("Object 'newdata' of class 'SpatialPixelsDataFrame' expected")
  }
  if(is.null(newlocs)){ 
    newlocs <- resample.grid(locations, newdata, silent=silent, t_cellsize=t_cellsize, quant.nndist=quant.nndist)$newlocs
  }
  if(missing(saga.env)){
    saga.env <- rsaga.env() 
  }
  s_te <- as.vector(newdata@bbox)
  if(silent==FALSE){
    message("Predicting at variable grid...")
  }
  if(missing(formula)){
    formula <- as.formula(paste(names(locations)[1], 1, sep="~"))
  }
  class(model) <- c("variogramModel", "data.frame")
  tvar <- all.vars(formula)[1]
  ok <- krige(formula, locations=locations[!is.na(locations@data[,tvar]),], newdata=newlocs, model=model, nmax=nmax, debug.level=-1, ...)
  ## write points to a shape file:
  tmp <- list(NULL)
  tmp.out <- list(NULL)
  if(predictOnly==TRUE){
    ok <- ok["var1.pred"]
  }
  for(k in 1:ncol(ok@data)){
    tmp[[k]] <- set.file.extension(tempfile(), ".shp")
    writeOGR(ok[k], names(ok)[k], dsn=tmp[[k]], "ESRI Shapefile")
    if(missing(file.name)){
      tmp.out[[k]] <- tempfile()
    } else {
      tmp.out[[k]] <- paste(file.name, k, sep="_")
    }
    ## point to grid (spline interpolation):
    suppressWarnings( rsaga.geoprocessor(lib=saga.lib[1], module=saga.module[1], param=list(SHAPES=tmp[[k]], FIELD=0, TARGET=0, METHOD=1, LEVEL_MAX=14, USER_XMIN=te[1]+t_cellsize/2, USER_XMAX=te[3]-t_cellsize/2, USER_YMIN=te[2]+t_cellsize/2, USER_YMAX=te[4]-t_cellsize/2, USER_SIZE=t_cellsize, USER_GRID=set.file.extension(tmp.out[[k]], ".sgrd")), show.output.on.console = FALSE, env=saga.env) )
    if(resample==TRUE){
      if(!all(te==s_te)|t_cellsize<newdata@grid@cellsize[1]){
        if(silent==FALSE){ message(paste("Resampling band", k, "to the target resolution and extent...")) }
        if(t_cellsize<newdata@grid@cellsize[1]){
          suppressWarnings( rsaga.geoprocessor(lib=saga.lib[2], module=saga.module[2], param=list(INPUT=set.file.extension(tmp.out[[k]], ".sgrd"), TARGET=0, SCALE_DOWN_METHOD=4, USER_XMIN=s_te[1]+t_cellsize/2, USER_XMAX=s_te[3]-t_cellsize/2, USER_YMIN=s_te[2]+t_cellsize/2, USER_YMAX=s_te[4]-t_cellsize/2, USER_SIZE=t_cellsize, USER_GRID=set.file.extension(tmp.out[[k]], ".sgrd")), show.output.on.console=FALSE, env=saga.env) )
        } else {
          ## upscale:
          suppressWarnings( rsaga.geoprocessor(lib=saga.lib[2], module=saga.module[2], param=list(INPUT=set.file.extension(tmp.out[[k]], ".sgrd"), TARGET=0, SCALE_DOWN_METHOD=0, SCALE_UP_METHOD=0, USER_XMIN=s_te[1]+t_cellsize/2, USER_XMAX=s_te[3]-t_cellsize/2, USER_YMIN=s_te[2]+t_cellsize/2, USER_YMAX=s_te[4]-t_cellsize/2, USER_SIZE=t_cellsize, USER_GRID=set.file.extension(tmp.out[[k]], ".sgrd")), show.output.on.console=FALSE, env=saga.env) )
        }
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
      if(silent==FALSE){ message(paste0("Created output SAGA GIS grid: ", tmp.out[[k]], ".sdat")) }
    }
  }
  if(missing(file.name)){
    return(out)
  }
}

## resample using variable sampling intensity:
resample.grid <- function(locations, newdata, silent=FALSE, n.sigma, t_cellsize, optN, quant.nndist=.5, breaks.d=NULL){
  if(silent==FALSE){
    message("Deriving density map...")
  }
  if(requireNamespace("spatstat", quietly = TRUE)&requireNamespace("maptools", quietly = TRUE)){
    ## derive density map:
    W <- as.matrix(newdata[1])
    W <- ifelse(is.na(W), FALSE, TRUE)
    suppressWarnings( locs.ppp <- spatstat::ppp(x=locations@coords[,1], y=locations@coords[,2], xrange=newdata@bbox[1,], yrange=newdata@bbox[2,], mask=t(W)[ncol(W):1,]) )
    if(missing(n.sigma)){
      dist.locs <- spatstat::nndist(locs.ppp)
      n.sigma <- quantile(dist.locs, quant.nndist)
    }
    if(n.sigma < 2*t_cellsize){ 
      warning(paste0("Estimated 'Sigma' too small. Using 2 * newdata cellsize."))
      n.sigma = 2*newdata@grid@cellsize[1]
    }
    if(n.sigma > sqrt(length(newdata)*newdata@grid@cellsize[1]*newdata@grid@cellsize[2]/length(locations))){ 
      warning(paste0("'Sigma' set at ", signif(n.sigma, 3), ". This is possibly an unclustered point sample. See '?resample.grid' for more information.")) 
    }
    dmap <- maptools::as.SpatialGridDataFrame.im(density(locs.ppp, sigm=n.sigma))
    dmap.max <- max(dmap@data[,1], na.rm=TRUE)
    dmap@data[,1] <- signif(dmap@data[,1]/dmap.max, 3)
    if(is.null(breaks.d)){
      ## TH: not sure if here is better to use quantiles or a regular split?
      breaks.d <- expm1(seq(0, 3, by=3/10))/expm1(3)
      #breaks.d <- quantile(dmap@data[,1], seq(0, 1, by=1/5)), na.rm=TRUE)
    }
    if(sd(dmap@data[,1], na.rm=TRUE)==0){ stop("Density map shows no variance. See '?resample.grid' for more information.") }
    dmap$strata <- cut(x=dmap@data[,1], breaks=breaks.d, include.lowest=TRUE, labels=paste0("L", 1:(length(breaks.d)-1)))
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
}

## end of script;