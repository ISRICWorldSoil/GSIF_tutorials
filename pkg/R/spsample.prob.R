# Purpose        : estimate occurrence probabilities;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ;
# Status         : pre-alpha
# Note           : Not recommended for large grids;


setMethod("spsample.prob", signature(observations = "SpatialPoints", covariates = "SpatialPixelsDataFrame"), function(observations, covariates, quant.nndist=.95, n.sigma, ...){
   
  ## mask out missing combinations:
  covariates <- covariates[stats::complete.cases(covariates@data),]
  ov <- over(observations, covariates)
  observations <- observations[stats::complete.cases(ov),]
   
  if(requireNamespace("spatstat", quietly = TRUE)&requireNamespace("maxlike", quietly = TRUE)){
    mg_owin <- spatstat::as.owin(data.frame(x = covariates@coords[,1], y = covariates@coords[,2], window = TRUE))
    suppressWarnings( locs.ppp <- spatstat::ppp(x=coordinates(observations)[,1], y=coordinates(observations)[,2], window=mg_owin) )
    dist.locs <- spatstat::nndist(locs.ppp)                    
    ## inlcusion probabilities geographical space:
    if(missing(n.sigma)){
      n.sigma <- quantile(dist.locs, quant.nndist)
    }
    if(n.sigma < 0.5*sqrt(length(covariates)*covariates@grid@cellsize[1]*covariates@grid@cellsize[2]/length(observations))){ 
        warning(paste0("'Sigma' set at ", signif(n.sigma, 3), ". Consider increasing the value.")) 
    }
    message(paste("Deriving kernel density map using sigma", signif(n.sigma, 3), "..."))
    dmap <- maptools::as.SpatialGridDataFrame.im(density(locs.ppp, sigm=n.sigma, ...))
    ## Pixel values are estimated intensity values, expressed in 'points per unit area' (hence multiply by area).
    dmap.max <- max(dmap@data[,1], na.rm=TRUE)
    dmap@data[,1] <- signif(dmap@data[,1]/dmap.max, 3)
    
    ## occurrence probabilities in feature space:
    message("Deriving inclusion probabilities using MaxLike analysis...")
    fm <- as.formula(paste("~", paste(names(covariates), collapse="+")))
    ml <- maxlike::maxlike(formula=fm, rasters=stack(covariates), points=observations@coords, method="BFGS", savedata=TRUE)
    ## bug in "maxlike" (https://github.com/rbchan/maxlike/issues/1); need to replace this 'by hand':
    ml$call$formula <- fm
    ## TH: this operation can be time consuming and is not recommended for large grids!
    ml.p <- as(predict(ml), "SpatialPixelsDataFrame")
    ## sum two occurrence probabilities (masks for the two maps need to be exactly the same):
    covariates$iprob <- signif((ml.p@data[,1] + dmap@data[,1])/2, 3)
     
    out <- list(prob=covariates["iprob"], observations=as(observations, "SpatialPoints"), density=dmap, maxlike=ml.p, maxlikeFit=ml[-which(names(ml)=="rasters")])
    return(out) 
  }
})