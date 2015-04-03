# Purpose        : Overlay points/lines/polygons over any list of files e.g. Landsat scenes;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : in the case multiple files follow the same pattern the values are aggregated;

extract.list <- function(x, y, path=".", ID="SOURCEID", method="simple", is.pattern=FALSE, force.projection = TRUE, NAflag = "", show.progress=TRUE, isFactor=FALSE, ...){
  if(requireNamespace("reshape", quietly = TRUE)){
    if(path=="."){ path <- getwd() }
    if(!file.exists(path)|!file.info(path)[1,"isdir"]) { stop(paste("Directory:", path, "does not exists")) }
    if(class(x)=="SpatialPoints"){
      x <- SpatialPointsDataFrame(x, data.frame(ID=as.factor(as.character(1:nrow(coordinates(x))))))
      names(x) <- paste(ID)
    } else {
      if(is.null(x@data[,ID])) { stop("'ID' column not found") }
    }
    ov <- list(NULL)
    ## simple case -> file list
    if(is.pattern==FALSE){
        message(paste("Extracting values for", length(x), "points from", length(y), "rasters..."))
        if (show.progress) { pb <- txtProgressBar(min=0, max=length(y), style=3) }
        for(i in 1:length(y)){
          if(.Platform$OS.type == "windows"){
            fname <- normalizePath(paste(path, y[i], sep="\\"), winslash="/")
          } else {
            fname <- normalizePath(paste(path, y[i], sep="/"))
          }
          try( r <- raster(fname) )
          if(!class(.Last.value)[1]=="try-error"){
            if(force.projection==TRUE){ 
              raster::projection(r) <- proj4string(x)
            }
            ov[[i]] <- extract(r, x, method=method, ...)
            names(ov)[i] <- basename(y[i])
          }
          if (show.progress) { setTxtProgressBar(pb, i) }
        }
        if (show.progress) { close(pb) }
  
        ## format to a data.frame:
        sel <- sapply(ov, is.null)
        ov.df <- cbind(x@data[ID], as.data.frame(ov[!sel]))
  
    ## look for pattern in files (e.g. Landsat scenes):
    } else {
        message(paste("Extracting values for", length(x), "points using pattern matching..."))
        if (show.progress) { pb <- txtProgressBar(min=0, max=length(y), style=3) }
        for(i in 1:length(y)){
          ## normalize var name:
          vname <- make.names(gsub("[[:punct:]]", "", y[i]))
          ## list all files (each can have a different coordinate system):
          lst <- list.files(path=path, pattern=glob2rx(y[i]), recursive=TRUE)
          if(length(lst)>0){
            tmp <- list(NULL)
            for(k in 1:length(lst)){
              ## overlay per file:
              if(.Platform$OS.type == "windows"){
                r <- raster(normalizePath(paste(path, lst[k], sep="\\"), winslash="/"))
              } else{
                r <- raster(normalizePath(paste(path, lst[k], sep="/")))
              }
              x.t <- spTransform(x, CRS(projection(r)))
              ## select only points within the bounding box to speed up:
              subs <- x.t@coords[,1] > extent(r)@xmin & x.t@coords[,1] < extent(r)@xmax & x.t@coords[,2] > extent(r)@ymin & x.t@coords[,2] < extent(r)@ymax
              if(sum(subs)>0){
                tmp[[k]] <- data.frame(x.t@data[,ID][subs], extract(r, x.t[subs,], method=method, ...))
                names(tmp[[k]]) <- paste(c(ID, vname))
                ## remove missing values:
                if(!NAflag==""){ tmp[[k]] <- tmp[[k]][!tmp[[k]][,2]==NAflag,] }
              }
            }
            ## fill-in missing pixels:
            tmp <- do.call(rbind, tmp)
            if(isFactor==FALSE){
              frm <- as.formula(paste(vname, "~", ID))
              message("Aggregating values using the ID column...")
              ov[[i]] <- aggregate(frm, data=tmp, mean, na.rm=TRUE)
            } else {
              ov[[i]] <- tmp
            }
            names(ov)[i] <- vname
          }
          if (show.progress) { setTxtProgressBar(pb, i) }
        }
        if (show.progress) { close(pb) }
        if(length(y)>1){
          ov.df <- reshape::merge_recurse(ov)
        } else {
          ov.df <- ov[[1]]
        }
    }
    return(ov.df)
  }
}

setMethod("extract", signature(x = "SpatialPoints", y = "character"), extract.list)
setMethod("extract", signature(x = "SpatialPointsDataFrame", y = "character"), extract.list)

## end of script;