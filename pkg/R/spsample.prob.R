# Purpose        : estimate inclusion probabilities;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ;
# Status         : pre-alpha
# Note           : Not recommended for large grids;


setMethod("spsample.prob", signature(observations = "SpatialPoints", covariates = "SpatialPixelsDataFrame"), function(observations, covariates, test.SRS = FALSE, quant.nndist=.95, n.sigma, ...){
   
  ## mask out missing combinations:
  covariates <- covariates[stats::complete.cases(covariates@data),]
  ov <- over(observations, covariates)
  observations <- observations[stats::complete.cases(ov),]
   
  require(spatstat)
  mg_owin <- spatstat::as.owin(data.frame(x = data.frame(covariates)[,"x"], y = data.frame(covariates)[,"y"], window = TRUE))
  suppressWarnings( locs.ppp <- spatstat::ppp(x=coordinates(observations)[,1], y=coordinates(observations)[,2], window=mg_owin) )
  dist.locs <- spatstat::nndist(locs.ppp)                    
  ## test Complete Spatial Randomness:
  if(test.SRS == TRUE){
    message("Testing spatial randomness of points...")
    env <- spatstat::envelope(locs.ppp, fun=Gest)
  } else {
    env <- NULL
  }
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
  
  ## inclusion probabilities feature space:
  require(dismo)
  message("Deriving inclusion probabilities using MaxEnt analysis...")
  me <- MaxEnt(occurrences=locs.ppp, covariates=covariates)
  ## TH: this operation can be time consuming and is not recommended for large grids!
  ## sum two inclusion probabilities (this assumes that masks are exactly the same):
  covariates$iprob <- (as(me@predicted, "SpatialPixelsDataFrame")@data[,1] + dmap@data[,1])/2
   
  out <- list(prob=covariates["iprob"], observations=as(observations, "SpatialPoints"), density=dmap, envelope=env, maxent=me)
  return(out) 
})