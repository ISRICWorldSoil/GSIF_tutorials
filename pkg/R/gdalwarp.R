# Purpose        : resampling with FWTools;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : FWTools is NOT an R package!;

.gdalwarp.SpatialPixels <- function(obj, proj4s = proj4string(obj), GridTopology = NULL, pixsize, resampling_method = "bilinear", NAflag = get("NAflag", envir = GSIF.opts), tmp.file = FALSE, show.output.on.console = FALSE, program){
  
  if(missing(program)){
    if(.Platform$OS.type == "windows") {
      fw.dir = .FWTools.path()        
      program = shQuote(shortPathName(normalizePath(file.path(fw.dir, "bin/gdalwarp.exe"))))
    } else {
      program = "gdalwarp"
  }}
  
  if(!nchar(program)==0){
    require(stringr)
    message(paste('Resampling', length(names(obj)), 'layers to CRS(\"', stringr::str_trim(substr(proj4s, 1, 20)), ' ... ', '\") with grid cell size:', pixsize, '...'))
    
    if(any(class(obj)=="SpatialPixelsDataFrame")){
      size = ncol(obj) 
    } else { if(any(class(obj)=="RasterLayer")){ 
      size = 1
    }}

    pb <- txtProgressBar(min=0, max=size, style=3)

    for(i in 1:size){
  
    ## name the temp file:
    if(any(class(obj)=="RasterLayer")){
        if(raster::inMemory(obj)==TRUE & tmp.file==TRUE){
           tf <- tempfile()
           extension <- ".tif"
        } else {
           ## file extension:
           require(tools)
           extension <- paste(".", file_ext(raster::filename(obj)), sep="")
           tf <- strsplit(raster::filename(obj), extension)[[1]]
        }
    }
    if(any(class(obj)=="SpatialPixelsDataFrame")){
        extension <- ".tif"
        if(tmp.file==TRUE){
           tf <- tempfile()
        } else {
           tf <- paste(normalizeFilename(deparse(substitute(obj, env = parent.frame()))), names(obj)[i], sep="_")
        }
    }
    
    ## check if it is factor or numeric:
    if(any(class(obj)=="RasterLayer")){
      isfactor <- is.factor(obj)
    } else {
      if(any(class(obj)=="SpatialPixelsDataFrame")){
        isfactor <- is.factor(obj@data[,i])
      }
    }
    
    ## write to a file if necessary:
    if(any(class(obj)=="SpatialPixelsDataFrame")){
       if(isfactor){
        x <- writeRaster(raster(obj[i]), filename=paste(tf, extension, sep=""), format="GTiff", overwrite=TRUE)
        if(raster::maxValue(x)==0) { warning(paste("Layer", names(obj[i]), "is of type 'factor' but contains no levels")) }
       } else {
         writeGDAL(obj[i], paste(tf, extension, sep=""), "GTiff", mvFlag = NAflag)
        }
    } else {
      if(any(class(obj)=="RasterLayer")){
        if(raster::inMemory(obj)==TRUE){
          x <- writeRaster(obj, filename=paste(tf, extension, sep=""), format="GTiff", overwrite=TRUE)
        }
        if(isfactor){ 
          if(raster::maxValue(x)==0) { warning(paste("Layer", names(obj), "is of type 'factor' but contains no levels")) }
          } 
    }}
      
    ## resample to WGS84 system:
    if(isfactor){
      if(is.null(GridTopology)){
        if(is.na(proj4string(obj))|proj4string(obj)=="NA"){
          system(paste(program, ' ', tf, extension, ' ', tf, '_ll.tif -dstnodata \"255\" -r near -ot \"Byte\" -tr ', pixsize, ' ', pixsize, sep=""), show.output.on.console = show.output.on.console)
        } else {
          system(paste(program, ' ', tf, '.tif', ' -t_srs \"', proj4s, '\" ', tf, '_ll.tif -dstnodata \"255\" -r near -ot \"Byte\" -tr ', pixsize, ' ', pixsize, sep=""), show.output.on.console = show.output.on.console)
        } 
      } else {
        if(class(GridTopology)=="GridTopology"){
          bbox = bbox(SpatialGrid(GridTopology))
            if(is.na(proj4string(obj))|proj4string(obj)=="NA"){
              system(paste(program, ' ', tf, extension, ' ', tf, '_ll.tif -dstnodata \"255\" -ot \"Byte\" -r near', ' -te ', bbox[1,1], ' ', bbox[2,1], ' ', bbox[1,2], ' ', bbox[2,2], ' -ts ', GridTopology@cells.dim[1], ' ', GridTopology@cells.dim[2], sep=""), show.output.on.console = show.output.on.console)
            } else {
              system(paste(program, ' ', tf, '.tif', ' -t_srs \"', proj4s, '\" ', tf, '_ll.tif -dstnodata \"255\" -ot \"Byte\" -r near', ' -te ', bbox[1,1], ' ', bbox[2,1], ' ', bbox[1,2], ' ', bbox[2,2], ' -ts ', GridTopology@cells.dim[1], ' ', GridTopology@cells.dim[2], sep=""), show.output.on.console = show.output.on.console)
            }
        } else {
          stop("'GridTopology-class' object required")
        }
      }
    }
       else {
          if(is.null(GridTopology)){ 
            if(is.na(proj4string(obj))|proj4string(obj)=="NA"){
              system(paste(program, ' ', tf, extension, ' ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r ', resampling_method, ' -tr ', pixsize, ' ', pixsize, sep=""), show.output.on.console = show.output.on.console) 
            } else {
              system(paste(program, ' ', tf, '.tif', ' -t_srs \"', proj4s, '\" ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r ', resampling_method, ' -tr ', pixsize, ' ', pixsize, sep=""), show.output.on.console = show.output.on.console)            
            }
          } else {
              if(class(GridTopology)=="GridTopology"){
                bbox = bbox(SpatialGrid(GridTopology))
                if(is.na(proj4string(obj))|proj4string(obj)=="NA"){
                  system(paste(program, ' ', tf, extension, ' ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r ', resampling_method, ' -te ', bbox[1,1], ' ', bbox[2,1], ' ', bbox[1,2], ' ', bbox[2,2], ' -ts ', GridTopology@cells.dim[1], ' ', GridTopology@cells.dim[2], sep=""), show.output.on.console = show.output.on.console)
                } else {
                  system(paste(program, ' ', tf, '.tif', ' -t_srs \"', proj4s, '\" ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r ', resampling_method, ' -te ', bbox[1,1], ' ', bbox[2,1], ' ', bbox[1,2], ' ', bbox[2,2], ' -ts ', GridTopology@cells.dim[1], ' ', GridTopology@cells.dim[2], sep=""), show.output.on.console = show.output.on.console)                
                }
              } else {
                stop("'GridTopology-class' object required")
              }
            }
        }
        
        ## read images back to R:
        if(i==1){
          if(any(class(obj)=="RasterLayer")){ 
            if(raster::inMemory(obj)==TRUE){
              res <- readGDAL(paste(tf, "_ll.tif", sep=""), silent = FALSE)
              names(res) <- names(obj)[i]
            } else {
              res <- raster(paste(tf, "_ll.tif", sep=""))
            }
          } else {
              if(any(class(obj)=="SpatialPixelsDataFrame")) {
                res <- readGDAL(paste(tf, "_ll.tif", sep=""), silent = FALSE)
                names(res) <- names(obj)[i]          
              }
          }
        } else{
          if(any(class(obj)=="SpatialPixelsDataFrame")) { 
            res@data[,names(obj)[i]] <- readGDAL(paste(tf, "_ll.tif", sep=""), silent = FALSE)$band1
          }
        }
        
        # reformat to the original factors:
        if(isfactor & any(class(obj)=="SpatialPixelsDataFrame")){
          res@data[,i] <- as.factor(res@data[,i])
          levels(res@data[,i]) = levels(obj@data[,i])
        }
        ## clean up:
        if(any(class(obj)=="RasterLayer")){ 
          if(raster::inMemory(obj)==FALSE){
            message(paste("\n", paste(tf, "_ll.tif", sep="")))
          }
        } else {
            unlink(paste(tf, "_ll.tif", sep=""))
            unlink(paste(tf, extension, sep=""))
        }        
    setTxtProgressBar(pb, i)          
    }
  close(pb)
  cat(i, "\r")
  flush.console()
  
  } else { 
    stop("Could not locate FWTools. First install and test FWTools. See 'plotKML.env()' for more info.") 
  }
  
  return(res)

}

setMethod("gdalwarp", signature(obj = "SpatialPixelsDataFrame"), .gdalwarp.SpatialPixels)
setMethod("gdalwarp", signature(obj = "RasterLayer"), .gdalwarp.SpatialPixels)

## end of script;
