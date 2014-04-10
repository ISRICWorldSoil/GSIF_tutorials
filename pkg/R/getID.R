# Purpose        : Locate the cells for a given bounding box;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : Pierre Roudier;
# Status         : tested
# Note           : Land mask will be continously updated;


setMethod("getID", signature(obj = "SpatialPolygons"), function(obj, pixsize = 3/3600, empty.tif = FALSE, compress = FALSE, zipname = set.file.extension(tempfile(tmpdir = getwd()), "zip")){
 
  if(is.projected(obj)){
    stop('Object of type "SpatialPolygons" projected in the geographical coordinates (WGS84) required')
  }

  # over polygons and landmask to get a list of IDs:
  load(system.file("data/landmask.rda", package="GSIF"))
  gridded(landmask) <- ~x+y
  proj4string(landmask) <- "+proj=longlat +datum=WGS84"
  list.id <- over(landmask, obj)
  names.id <- NULL
  
  for(j in which(!is.na(list.id))){ 
    X0 <- floor(landmask[j,]@coords[,1])
    Y0 <- floor(landmask[j,]@coords[,2])
    fname <- paste(ifelse(X0<0, "W", "E"), abs(X0), "_", ifelse(Y0<0, "S", "N"), abs(Y0), sep="")

  if(empty.tif == TRUE) {
    if(is.na(file.info(paste(getwd(), "/", fname, ".tif", sep=""))$size)){
      grid.1deg <- expand.grid(Lon=seq(X0+pixsize/2, X0+1-pixsize/2, by=pixsize), Lat=seq(Y0+pixsize/2, Y0+1-pixsize/2, by=pixsize), KEEP.OUT.ATTRS=FALSE)
      grid.1deg$blank <- rep(1, length(grid.1deg[,1]))
      coordinates(grid.1deg) <- ~Lon+Lat
      gridded(grid.1deg) <- TRUE
      proj4string(grid.1deg) <- CRS("+proj=longlat +datum=WGS84")
      writeGDAL(grid.1deg["mask"], paste(fname, ".tif", sep=""), "GTiff", type="Byte", mvFlag=0)
      }
      
      if(compress == TRUE) {
      zip(zipfile=zipname, files=paste(fname, ".tif", sep=""))
      unlink(paste(fname, ".tif", sep=""))  
      } 
  } 

  names.id[j] <- fname
  
  } 

  names.lst <- levels(as.factor(names.id))
  return(names.lst)

}) 

# end of script;