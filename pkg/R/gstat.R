# Purpose        : Uses windows binaries for gstat.exe (ordinary / universal kriging predictions and simulations);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : This code was developed on Windows 7 OS. Important note: Gstat.exe is not used mainted by the author of gstat any more (read more at http://www.gstat.org/). Possibly a better option is to use Python interface to SGeMS (http://sgems.sourceforge.net/) and/or HPGL library (http://hpgl.aoizora.org);


## Writes geosamples to a GEO-EAS text file
write.data.geosamples <- function(obj, covariates = NULL, methodid, outfile, outformat = "GeoEAS", digits = 4) {
  
  # prepare regression matrix:
  if(is.null(covariates)){
    # subset to XYZV:
    ov <- obj@data[,c("longitude", "latitude", "altitude", "observedValue")]
  } else {
    ov <- over(x=covariates, y=obj, methodid=methodid, var.type = "numeric")
    ov <- ov[,c(attr(covariates@coords, "dimnames")[[2]], "altitude", "observedValue", names(covariates))]
    names(ov)[which(names(ov) %in% attr(covariates@coords, "dimnames")[[2]])] <- c("longitude", "latitude")
  }

  # geostats only possible with numeric variables:
  ov$observedValue = as.numeric(ov$observedValue)
  ov <- ov[!is.na(ov$observedValue),]
  # round up numbers:
  ov$longitude <- signif(ov$longitude, digits+2)
  ov$latitude <- signif(ov$latitude, digits+2)
  ov$altitude <- signif(ov$altitude, digits)
  ov$observedValue <- signif(ov$observedValue, digits)
  for(j in names(covariates)){ if(is.numeric(ov[,j])) { ov[,j] <- signif(ov[,j], digits) } }

  if(outformat == "GeoEAS"){
    # Write the header
    cat('GeoEAS file created in R\n', file=outfile)
    cat(length(names(ov)), file=outfile, append=TRUE)
    cat('\n', file=outfile, append=TRUE)
    write(cbind(names(ov)), file=outfile, append=TRUE)
    write.table(ov, file=outfile, append=TRUE, sep='\t', col.names=FALSE, row.names=FALSE)
    
  } else {
    stop("Format currently not supported")
  }
 
}

setMethod("write.data", signature="geosamples", write.data.geosamples)

## write spatial points:
write.data.SpatialPoints <- function(obj, covariates = NULL, methodid, outfile, outformat = "GeoEAS", digits = 4) {
  
  # prepare regression matrix:
  if(is.null(covariates)){
    # subset to XYZV:
    ov <- data.frame(obj[methodid])
    # re-order:
    ov <- ov[,c(attr(obj@coords, "dimnames")[[2]], methodid)]
  } else {
    ov <- sp::over(x=obj, y=covariates)
    ov <- cbind(data.frame(obj[methodid]), ov)
    # re-order:
    ov <- ov[,c(attr(obj@coords, "dimnames")[[2]], methodid, names(covariates))]
  }

  # geostats only possible with numeric variables:
  ov[,methodid] = as.numeric(ov[,methodid])
  ov <- ov[!is.na(ov[,methodid]),]
  # round up numbers:
  ov[,methodid] <- signif(ov[,methodid], digits)
  for(j in names(covariates)){ if(is.numeric(ov[,j])) { ov[,j] <- signif(ov[,j], digits) } }

  if(outformat == "GeoEAS"){
    # Write the header
    cat('GeoEAS file created in R\n', file=outfile)
    cat(length(names(ov)), file=outfile, append=TRUE)
    cat('\n', file=outfile, append=TRUE)
    write(cbind(names(ov)), file=outfile, append=TRUE)
    write.table(ov, file=outfile, append=TRUE, sep='\t', col.names=FALSE, row.names=FALSE)
    
  } else {
    stop("Format currently not supported")
  }
 
}

setMethod("write.data", signature="SpatialPointsDataFrame", write.data.SpatialPoints)


## generate a script (this will work only with formulas where the predictors are actual gridded maps):
makeGstatCmd <- function(formString, vgmModel, outfile, easfile, nsim = 0, nmin = 20, nmax = 40, radius, zmap = 0, predictions = "var1.pred.hdr", variances = "var1.svar.hdr", xcol = 1, ycol = 2, zcol = 3, vcol = 4, Xcols){

  methodid <- all.vars(formString)[1]
  masks <- all.vars(formString)[-1] 
  out = file(outfile, "w")

  # 2D or 3D? see page 39 and 40 in the gstat manual [http://www.gstat.org/gstat.pdf]
  if(vgmModel$ang2[2]==0&vgmModel$ang3[2]==0&vgmModel$anis2[2]==1){ 
    dims = 2 
  } else {
    dims = 3
  }
  # covariate columns:
  if(missing(Xcols)){ 
    if(dims==3){ 
      Xcols = paste(5:(length(masks)+4), collapse="&", sep="") 
    } else {
      Xcols = paste(4:(length(masks)+3), collapse="&", sep="")     
    }
  }
  # estimate the radius:
  if(missing(radius)){
     radius = vgmModel$range[2]*1.6
  }

  if(dims == 3){
    write(paste("data(", methodid, "): '", easfile, "', x=", xcol, ", y=", ycol, ", z=", zcol, ", v=", vcol, ", min=", nmin, ", max=", nmax, ", radius=", radius, ", X=", Xcols, ", average;", sep=""), file=out, append=TRUE)
    write(paste("variogram(", methodid, "): ", vgmModel$psill[1], " Nug(0) + ", vgmModel$psill[2], " ", vgmModel$model[2], "(", vgmModel$range[2], ",", vgmModel$ang1[2], ",", vgmModel$ang2[2], ",", vgmModel$ang3[2], ",", vgmModel$anis1[2], ",", vgmModel$anis2[2], ");", sep=""), file=out, append=TRUE)  
    write(paste("set zmap=", zmap, ";", sep=""), file=out, append=TRUE)
    } else {
    write(paste("data(", methodid, "): '", easfile, "', x=", xcol, ", y=", ycol, ", v=", vcol, ", min=", nmin, ", max=", nmax, ", radius=", radius, ", X=", Xcols, ", average;", sep=""), file=out, append=TRUE)
    write(paste("variogram(", methodid, "): ", vgmModel$psill[1], " Nug(0) + ", vgmModel$psill[2], " ", vgmModel$model[2], "(", vgmModel$range[2], ",", vgmModel$ang1[2], ",", 0, ",", 0, ",", vgmModel$anis1[2], ",", 1, ");", sep=""), file=out, append=TRUE)      
    }
    write(paste("mask: ", paste("'", masks, "'", collapse=",", sep=""), ";", sep=""), file=out, append=TRUE) 
    if(nsim==0){
      write(paste("predictions(", methodid, "): '", predictions, "';", sep=""), file=out, append=TRUE) 
      write(paste("variances(", methodid, "): '", variances, "';", sep=""), file=out, append=TRUE) 
    } else { 
    # simulations:
    if(nsim>0){
      write("method: gs;", out, append=TRUE) 
      write(paste("set nsim=", nsim, ";", sep=""), file=out, append=TRUE)     
      write(paste("predictions(", methodid, "): '", predictions, "';", sep=""), file=out, append=TRUE) 
    }
  }
     
  close(out)
  message("You can now run the gstat command script using the 'system' command")
  
}

# end of script;

