# Purpose        : Geotif validator
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ; 
# Dev Status     : Beta
# Note           : For more details see [https://docs.google.com/document/d/1APpIyOFdXlJhXGd_huejJjX3EnCTUj2q_K8SfljL62E/pub]

SoilGrid.validator <- function(obj, domain, ground.truth, N.sample=2000, xml.file, z.lim, md.type="INSPIRE", test.URL=FALSE){

  if(!class(obj)=="GDALobj"|!class(domain)=="GDALobj"){
    stop("Objects of class \"GDALobj\" required.")
  }
  if(obj[["bands"]]>1){ stop("File contains multiple layers. Layer with a single band required.") }
  if(!missing(ground.truth)){
    if(!class(ground.truth)=="SpatialPointsDataFrame"){
      stop("'ground.truth' object of class \"SpatialPointsDataFrame\" required.")
    }
  }
  
  if(requireNamespace("RCurl", quietly = TRUE)&requireNamespace("Hmisc", quietly = TRUE)&requireNamespace("XML", quietly = TRUE)){
         
    GDALobj.file <- attr(obj, "file")
    domain.file <- attr(domain, "file")
    prj4_string <- attr(obj, "projection")
    rows <- obj[["rows"]]
    columns <- obj[["columns"]]
    GDType <- attr(obj, "df")$GDType[1]
    res.x <- obj[["res.x"]]
    res.y <- obj[["res.y"]]
    ll.x <- obj[["ll.x"]]
    ll.y <- obj[["ll.y"]]
  
    ## C_MTD | Metadata file is available (in the same directory)?
    if(missing(xml.file)){ xml.file <- paste0(GDALobj.file, ".xml") }
    if(file.exists(xml.file)){
      message(paste("Reading the metadata file: ", xml.file, sep=""))
      ret <- XML::xmlTreeParse(xml.file)  
      ml <- XML::xmlRoot(ret)
      a <- XML::xmlTreeParse(xml.file, useInternalNodes = TRUE)
      if(md.type=="INSPIRE"){
        C_MTD <- XML::xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["title"]][[1]])
      } else {
        C_MTD <- NA
      }
    } else {
      stop("Metadata file missing.")
    }
    
    ## C_MTS | Metadata file passes schema validation (see INSPIRE metadata validator)?
    if(md.type=="INSPIRE"){ 
      if(test.URL==TRUE){
        message("Connecting to the INSPIRE schema...")
        try(xsd <- XML::xmlTreeParse(get("inspire_xsd", envir=plotKML.opts), isSchema=TRUE, useInternalNodes = TRUE))
        try(xsd.val <- XML::xmlSchemaValidate(xsd, a))
        ## print the results of validation:
        x <- NULL
        try(x <- sapply(xsd.val$errors, function(x){x["msg"][[1]]}))
        if(nzchar(x)){
          C_MTS <- paste("Warning:", length(xsd.val$errors), "errors detected")
        } else {
          C_MTS <- "No problems detected"
        }
      } else {
        C_MTS <- "Not tested"
      }
    }
  
    ## C_URL | All layers in the documentation are also available on the server (the URL in the metadata file exists)?
    try( C_URL <- XML::xmlValue(ml[["distributionInfo"]][["MD_Distribution"]][["transferOptions"]][["MD_DigitalTransferOptions"]][["onLine"]][["CI_OnlineResource"]][["linkage"]][["URL"]]) )
    if(!is.na(C_URL)&exists("C_URL")){
      if(RCurl::url.exists(C_URL)){
        ## test downloading the file:
        if(test.URL==TRUE){
          message("Testing download speed...")
          x <- system.time( download.file(C_URL, destfile=basename(C_URL)) ) 
          download.time <- x[["elapsed"]]
        } else {
          download.time <- NA
        }
      } else {
        C_URL <- "URL does not exist"
        download.time <- NA
      }
    } else {
      C_URL <- "Missing URL in metadata file"
      download.time <- NA
    }
   
    ## C_URD | URLs to find full data description have been provided in the Metadata file?
    try( C_URD <- XML::xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["pointOfContact"]][["CI_ResponsibleParty"]][["contactInfo"]][["CI_Contact"]][["onlineResource"]][["CI_OnlineResource"]][["linkage"]][[1]]) )
    if(!is.na(C_URD)&exists("C_URD")){
      if(RCurl::url.exists(C_URD)){
        C_URD <- "URL validated"
      } else {
        C_URD <- "URL does not exist"
      }
    } else {
      C_URD <- "Missing URL in metadata file"
    }
  
    ## C_SMK | Spatial predictions are available for at least 95% of the soil mask?
    ri <- raster(GDALobj.file)
    if(fromDisk(ri)){
      rpnts <- spsample(SpatialGrid(GridTopology(c(ll.x,ll.y), c(res.x,res.y), c(columns,rows)), proj4string = CRS(prj4_string)), n=N.sample, type="random")
      message(paste("Overlaying", N.sample, "random points..."))
      ov <- data.frame(z=extract(x=ri, y=rpnts), d=extract(x=raster(domain.file), y=rpnts))
    } else {
      ## smaller files can be imported to memory:
      try( s <- stack(c(GDALobj.file, domain.file)) )
      if(!class(.Last.value)[1]=="try-error"){
        ov <- data.frame(z=getValues(s[[1]]), d=getValues(s[[2]]))
      } else {
        stop("Overlay resulted in an empty set") 
      }
    }
    NA.mask <- !is.na(ov$d)
    C_SMK <- round(sum(!is.na(ov$z[NA.mask]))/sum(NA.mask)*100, 1)
    mask.size <- round(sum(NA.mask)/length(NA.mask)*100, 1)
    
    ## N_NRC | All rasters conform to the same grid (raster stack) i.e. all have exactly the same number of rows and columns, and the same coordinate system? 
    if(any(!rows==domain[["rows"]]|!columns==domain[["columns"]]|!res.x==domain[["res.x"]]|!ll.x==domain[["ll.x"]]|!ll.y==domain[["ll.y"]]|!prj4_string==attr(domain, "projection"))){
      N_NRC <- "Inconsistent grid"
    } else {
      N_NRC <- "Consistent"
    }
    
    ## N_PRJ | Proj4 string is provided (embedded in the GeoTIFF)? 
    N_PRJ <- ifelse(prj4_string=="NA"|nchar(prj4_string)==0, "Missing", prj4_string)
    ## N_PRS | Provided Proj4 string is standard - i.e. it is registered under some coordinate system registry e.g. http://www.prj2epsg.org/ (and available in different sofware formats)?
    wkt <- showWKT(prj4_string)
    if(requireNamespace("rjson", quietly = TRUE)){
      if(test.URL==TRUE){
      try( ret <- rjson::fromJSON(file=paste0("http://www.prj2epsg.org/search.json?mode=wkt&terms=", RCurl::curlEscape(wkt))), silent = TRUE )
        if(!class(.Last.value)[1]=="try-error" & !length(ret)==0){
          N_PRS <- ifelse(ret$exact, ret$codes[[1]]$code, "Could not determine EPSG code")
        } else {
          N_PRS <- "Not existent in the http://www.prj2epsg.org/"
        }
      }
    } else {
      N_PRS <- "Not checked"
    }
  
    ## N_LIM | For numeric variables - No spatial predictions outside of physical limits (upper and lower limits) e.g. pH < 2 or pH > 11? / For factor-type variables correspond to the original domain e.g. no predicted values outside the original legend?
    message("Deriving value range...")
    z.rn <- range(ov$z[NA.mask], na.rm=TRUE)
    ## R_RME | RMSE i.e. Attribute Measurement Resolution or numeric resolution (RMSE/2) has been defined in the metadata file?
    try( numeric.resolution <- as.numeric(XML::xmlValue(ml[["eainfo"]][["detailed"]][["attr"]][["attrdomv"]][["rdom"]][["attrmres"]])) )
    if(!is.na(numeric.resolution)&exists("numeric.resolution")){
      R_RME <- "RMSE defined in metadata"      
    } else {
      numeric.resolution <- (z.rn[2]-z.rn[1])/nclass.Sturges(ov$z[NA.mask])
      R_RME <- NA
    }
    if(!missing(z.lim)){
      if(is.numeric(z.lim)&length(z.lim)==2){
        N_LIM <- ifelse(z.rn[1]<z.lim[1]|z.rn[2]>z.lim[2], "Predictions outside natural limits", "Predictions within limits")
      } else {
        N_LIM <- NA
      }
    } else {
      N_LIM <- NA
    }
     
    ## R_RMS | Spatial overlays using independent point data sets (at least 3 study areas / subsets representing the domain of interest) do not show RMSE > 2 times higher than the original RMSE i.e. four times the numeric resolution specified (extra analysis required)?
    if(!missing(ground.truth)){
      message("Testing RMSE using ground truth data (first column) ...")
      if(!proj4string(ground.truth)==prj4_string){
         ground.truth <- spTransform(ground.truth, CRS(prj4_string))
      }
      ov.cv <- data.frame(meas=ground.truth$z, pred=extract(x=ri, y=ground.truth))
      RMSE <- signif( sqrt(mean((ov.cv$meas - ov.cv$pred)^2, na.rm=TRUE)) , 3)
      R_RMS <- ifelse(RMSE>4*numeric.resolution, "RMSE exceeds 2 x declared RMSE", "RMSE within limits")
      ## N_HST | Histogram of the predicted variable corresponds to the sampled histogram (Chi-square test)?
      hbb <- Hmisc::histbackback(ov.cv$meas, ov.cv$pred, prob=TRUE)
      ks <- ks.test(hbb$left, hbb$right)
      N_HST <- ifelse(ks$p.value < 0.05, "Two-sample Kolmogorov-Smirnov test failed", "No significant difference")
    } else {
      RMSE <- NA
      R_RMS <- NA
      N_HST <- NA
      ov.cv <- NA
    }
    
    out <- list(summaries=list(file.name=basename(GDALobj.file), prj4_string=N_PRJ, rows=rows, columns=columns, GDType=GDType, ll.x=ll.x, ll.y=ll.y, Z_MIN=signif(min(ov$z[NA.mask], na.rm=TRUE), 3), Z_MAX=signif(max(ov$z[NA.mask], na.rm=TRUE), 3), Z_STD=signif(sd(ov$z[NA.mask], na.rm=TRUE), 3), numeric.resolution=signif(numeric.resolution, 3), N_LIM=N_LIM, mask.size=mask.size, C_SMK=C_SMK, C_MTD=C_MTD, C_MTS=C_MTS, C_URL=C_URL, C_URD=C_URD, R_RME=R_RME, RMSE=RMSE, R_RMS=R_RMS, N_HST=N_HST, download.time=download.time), cross.val=ov.cv)
    message("Finished!")
    return(out)
  }

}