# Purpose        : Overlay points/lines/polygons over any list of files e.g. Landsat scenes;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Status         : pre-alpha
# Note           : in the case multiple files follow the same pattern the values are aggregated;

extract.list <- function(x, y, path=".", ID="SOURCEID", method="simple", is.pattern=FALSE, force.projection = TRUE, NAflag = "", ...){
  require(reshape)
  if(path=="."){ path <- getwd() }
  if(!file.info(path)[1,"isdir"]) { stop(paste("Directory:", path, "does not exists")) }
  if(length(y)<2){ stop("'list' argument must contain more than one element") }
  if(class(x)=="SpatialPoints"){
    x <- SpatialPointsDataFrame(x, data.frame(ID=as.factor(as.character(1:nrow(coordinates(x))))))
    names(x) <- paste(ID)
  }
  ov <- list(NULL)
  ## simple case -> file list
  if(is.pattern==FALSE){
      message(paste("Extraction values from", length(x), "points for", length(y), "rasters..."))
      pb <- txtProgressBar(min=0, max=length(y), style=3)
      for(i in 1:length(y)){
        fname <- normalizePath(paste(path, y[i], sep="\\"), winslash="/")
        try( r <- raster(fname) )
        if(!class(.Last.value)[1]=="try-error"){
          if(force.projection==TRUE){ 
            raster::projection(r) <- proj4string(x)
          }
          ov[[i]] <- extract(r, x, method=method, ...)
          names(ov)[i] <- y[i]
        }
        setTxtProgressBar(pb, i)
      }
      close(pb)

      ## format to a data.frame:
      sel <- sapply(ov, is.null)
      ov.df <- cbind(x@data[ID], as.data.frame(ov[!sel]))

  ## look for pattern in files (e.g. Landsat scenes):
  } else {
      message(paste("Extraction values from", length(x), "points using pattern matching..."))
      pb <- txtProgressBar(min=0, max=length(y), style=3)
      for(i in 1:length(y)){
        ## normalize var name:
        vname <- gsub("[[:punct:]]", "", y[i])
        ## list all files (each can have a different coordinate system!):
        lst <- list.files(path=path, pattern=glob2rx(y[i]), recursive=TRUE)
        if(length(lst)>0){
          tmp <- list(NULL)
          for(k in 1:length(lst)){
            ## overlay per file:
            r <- raster(normalizePath(paste(path, lst[k], sep="\\"), winslash="/"))
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
          ## fill-in all missing pixels:
          tmp <- do.call(rbind, tmp)
          frm <- as.formula(paste(vname, "~", ID))
          ov[[i]] <- aggregate(frm, data=tmp, mean, na.rm=TRUE)
          names(ov)[i] <- vname
        }
        setTxtProgressBar(pb, i)
      } 
      close(pb)
      ov.df <- reshape::merge_recurse(ov)
  }
  return(ov.df)
}

setMethod("extract", signature(x = "SpatialPoints", y = "character"), extract.list)
setMethod("extract", signature(x = "SpatialPointsDataFrame", y = "character"), extract.list)

## end of script;