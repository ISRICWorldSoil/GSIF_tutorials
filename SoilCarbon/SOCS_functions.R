
## Get global 30 m resolution soil covariates
get_30m_covariates <- function(te, tr, p4s, input.vrt=c("/mnt/cartman/GlobalForestChange2000-2014/first.vrt", "/mnt/cartman/GlobalForestChange2000-2014/treecover2000.vrt", "/mnt/cartman/GlobalSurfaceWater/occurrence.vrt", "/mnt/cartman/GlobalSurfaceWater/extent.vrt", "/mnt/cartman/SRTMGL1/SRTMGL1.2.tif", "/mnt/cartman/Landsat/bare2010.vrt"), l.bands=c("REDL00", "NIRL00", "SW1L00", "SW2L00")){
  ## paste0("/mnt/cartman/GlobCover30/2010/glc", seq(10,100,by=10),".vrt"))
  for(j in 1:length(input.vrt)){
    if(basename(input.vrt[j])=="first.vrt"){
      tmp <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".tif")
      system(paste0('gdalwarp ', input.vrt[j], ' ', tmp, ' -co \"COMPRESS=DEFLATE\" -t_srs \"', p4s, '\" -tr ', tr, ' ', tr, ' -te ', te))
      for(k in 1:4){
        outname = paste0(gsub(".vrt", paste0("_", l.bands[k]), gsub("/", "_", gsub("/mnt/cartman/", "", (input.vrt[j])))), ".tif")
        if(!file.exists(outname)){ system(paste0('gdal_translate ', tmp, ' ', outname, ' -co \"COMPRESS=DEFLATE\" -b ', k)) }
      }
      unlink(tmp)
    } else {
      outname = RSAGA::set.file.extension(gsub("/", "_", gsub("/mnt/cartman/", "", (input.vrt[j]))), ".tif")
      if(!file.exists(outname)){ system(paste0('gdalwarp ', input.vrt[j], ' ', outname, ' -co \"COMPRESS=DEFLATE\" -t_srs \"', p4s, '\" -tr ', tr, ' ', tr, ' -te ', te)) }
    }
  }
}

saga_DEM_derivatives <- function(INPUT, MASK=NULL, sel=c("SLP","TWI","CRV","VBF","VDP","OPN","DVM"), saga_cmd="saga_cmd"){
  if(!tools::file_ext(INPUT)=="sgrd"){ 
    system(paste0('gdal_translate ', INPUT, ' -of \"SAGA\" ', RSAGA::set.file.extension(INPUT, ".sdat")))
    INPUT = RSAGA::set.file.extension(INPUT, ".sgrd")
  }
  if(!is.null(MASK)){
    ## Fill in missing DEM pixels:
    suppressWarnings( system(paste0(saga_cmd, ' grid_tools 25 -GRID=\"', INPUT, '\" -MASK=\"', MASK, '\" -CLOSED=\"', INPUT, '\"')) )
  }
  ## Slope:
  if(any(sel %in% "SLP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 0 -ELEVATION=\"', INPUT, '\" -SLOPE=\"', gsub(".sgrd", "_slope.sgrd", INPUT), '\" -C_PROF=\"', gsub(".sgrd", "_cprof.sgrd", INPUT), '\"') ) ) )
  }
  ## TWI:
  if(any(sel %in% "TWI")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_hydrology 15 -DEM=\"', INPUT, '\" -TWI=\"', gsub(".sgrd", "_twi.sgrd", INPUT), '\"') ) ) )
  }
  ## MrVBF:
  if(any(sel %in% "VBF")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 8 -DEM=\"', INPUT, '\" -MRVBF=\"', gsub(".sgrd", "_vbf.sgrd", INPUT), '\" -T_SLOPE=10 -P_SLOPE=3') ) ) )
  }
  ## Valley depth:
  if(any(sel %in% "VDP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_channels 7 -ELEVATION=\"', INPUT, '\" -VALLEY_DEPTH=\"', gsub(".sgrd", "_vdepth.sgrd", INPUT), '\"') ) ) )
  }
  ## Openess:
  if(any(sel %in% "OPN")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_lighting 5 -DEM=\"', INPUT, '\" -POS=\"', gsub(".sgrd", "_openp.sgrd", INPUT), '\" -NEG=\"', gsub(".sgrd", "_openn.sgrd", INPUT), '\" -METHOD=0' ) ) ) )
  }
  ## Deviation from Mean Value:
  if(any(sel %in% "DVM")){
    suppressWarnings( system(paste0(saga_cmd, ' statistics_grid 1 -GRID=\"', INPUT, '\" -DEVMEAN=\"', gsub(".sgrd", "_devmean.sgrd", INPUT), '\" -RADIUS=11' ) ) )
  }
}

## Spatiotemporal overlay
ov_st <- function(x, y, profs, vnames, col.names=c("SOURCEID","YEAR","OCDENS","DEPTH.f","LONWGS84","LATWGS84")){
  sel = which(profs$YEAR_c==x)
  pnts = profs[sel,]
  if(nrow(pnts)>0){
    sfInit(parallel=TRUE, cpus=56)
    sfExport("pnts", "y")
    sfLibrary(rgdal)
    sfLibrary(raster)
    out <- data.frame(sfClusterApplyLB(y, function(i){try( raster::extract(y=pnts, x=raster(i)) )}))
    sfStop()
    names(out) = vnames
    out <- cbind(as.data.frame(pnts)[,col.names], out)
    return(out)
  }
}

## Convert soil horizon data to x,y,d regression matrix for 3D modeling:
hor2xyd = function(x, U="UHDICM", L="LHDICM", treshold.T=15){
  x$DEPTH <- x[,U] + (x[,L] - x[,U])/2
  x$THICK <- x[,L] - x[,U]
  sel = x$THICK < treshold.T
  ## begin and end of the horizon:
  x1 = x[!sel,]; x1$DEPTH = x1[,L]
  x2 = x[!sel,]; x2$DEPTH = x1[,U]
  y = do.call(rbind, list(x, x1, x2))
  return(y)
}


