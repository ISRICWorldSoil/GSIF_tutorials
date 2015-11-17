# title         : exportSoilGrids1km.R
# purpose       : exports maps and table data;
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL.
# inputs        : Predicted soil properties and classes as geotifs;
# outputs       : txt, KML files;
# remarks 1     : mosaicking takes about 1-2 hours but compression takes about 4-5 hours!!

library(rgdal)
library(GSIF)
library(RSAGA)
library(plotKML)
library(XML)

fw.path <- utils::readRegistry("SOFTWARE\\WOW6432NODE\\FWTools")$Install_Dir
gdalwarp <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gtifout.lst <- sapply(list.files(pattern=glob2rx("*.tif.gz$"), recursive=TRUE, full.names=FALSE), function(x){strsplit(x, ".tif.gz")[[1]][1]})
load("m_TAXGWRB.rda"); load("m_TAXOUSDA.rda")
load("world_country.rda")
kml.outdir <- "G:/SoilGrids1km/KML"

## make legends for SAGA GIS:
load("soil.legends.rda")
#data(soil.legends)
sel.leg <- which(names(soil.legends) %in% c("ORCDRC","PHIHOX","PHIKCL","BLD","CEC","SNDPPT","SLTPPT","CLYPPT","CRFVOL"))
for(j in sel.leg){
    Group = as.factor(paste(signif(soil.legends[[j]]$MIN, 3), signif(soil.legends[[j]]
    $MAX, 3), sep=" - "))
    unlink(set.file.extension(paste(names(soil.legends)[j]), ".txt"))
    makeSAGAlegend(x=Group, col_pal=soil.legends[[j]]$COLOR, filename=set.file.extension(paste(names(soil.legends)[j]), ".txt"), MINIMUM=soil.legends[[j]]$MIN, MAXIMUM=soil.legends[[j]]$MAX)
}
## categories separately because some classes can get 'kicked-out' during modelling:
USDAf <- data.frame(Group=m_TAXOUSDA$lev, levs=1:length(m_TAXOUSDA$lev))
USDA.col <-  merge(USDAf, soil.legends[["TAXOUSDA"]])
unlink("TAXOUSDA.txt")
makeSAGAlegend(x=as.factor(paste(USDA.col$Group)), MINIMUM=USDA.col$levs, MAXIMUM=USDA.col$levs+1, col_pal=USDA.col$COLOR, filename="TAXOUSDA.txt")
WRBf <- data.frame(Group=m_TAXGWRB$lev, levs=1:length(m_TAXGWRB$lev))
WRB.col <- merge(WRBf, soil.legends[["TAXGWRB"]])
unlink("TAXGWRB.txt")
makeSAGAlegend(x=as.factor(paste(WRB.col$Group)), MINIMUM=WRB.col$levs, MAXIMUM=WRB.col$levs+1, col_pal=WRB.col$COLOR, filename="TAXGWRB.txt")

legends.sdl <- soil.legends[sel.leg]
USDA.col$MIN <- USDA.col$levs; USDA.col$MAX <- USDA.col$levs+1
WRB.col$MIN <- WRB.col$levs; WRB.col$MAX <- WRB.col$levs+1
legends.sdl[["TAXOUSDA"]] <- USDA.col
legends.sdl[["TAXGWRB"]] <- WRB.col

## SLD files:
library(XML)
for(j in 1:length(legends.sdl)){
  sld.file <- file(set.file.extension(names(legends.sdl)[j], ".sld"), "w", blocking=FALSE)
  l1 = newXMLNode("StyledLayerDescriptor", attrs=c("xsi:schemaLocation" = "http://www.opengis.net/sld StyledLayerDescriptor.xsd", version="1.0.0"), namespaceDefinitions=c("http://www.opengis.net/sld", "xsi" = "http://www.w3.org/2001/XMLSchema-instance", "ogc" = "http://www.opengis.net/ogc", "gml" = "http://www.opengis.net/gml"))
  l2 <- newXMLNode("NamedLayer", parent = l1)
  l3 <- newXMLNode("Name", "SoilGrids1km", parent = l2)
  l3b <- newXMLNode("UserStyle", parent = l2)
  l4 <- newXMLNode("Title", paste("SoilGrids1km", names(legends.sdl)[j], sep="_"), parent = l3b)
  l4b <- newXMLNode("FeatureTypeStyle", parent = l3b)
  l5 <- newXMLNode("Rule", parent = l4b)
  l6 <- newXMLNode("RasterSymbolizer", parent = l5)
  if(names(legends.sdl)[j] %in% c("TAXOUSDA","TAXGWRB")){
    Group = paste(legends.sdl[[j]]$Group)
  } else {
    Group = paste(signif(legends.sdl[[j]]$MIN, 3), signif(legends.sdl[[j]]$MAX, 3), sep=" - ") ## \226
  }
  levs = legends.sdl[[j]]$MAX   
  l7 <- newXMLNode("ColorMap", attrs=c(type="intervals"), parent = l6)
  if(names(legends.sdl)[j] %in% c("TAXOUSDA","TAXGWRB")){
    txt <- sprintf('<ColorMapEntry color="%s" quantity="%.0f" label="%s" opacity="%.1f"/>', c("#FFFFFF", strtrim(legends.sdl[[j]]$COLOR, 7)), c(1, levs), c("NODATA", Group), c(0.0, rep(.7, length(legends.sdl[[j]]$COLOR))))  
  } else {
    if(names(legends.sdl)[j] %in% c("PHIHOX","SNDPPT","SLTPPT","CLYPPT")){
       txt <- sprintf('<ColorMapEntry color="%s" quantity="%.0f" label="%s" opacity="%.1f"/>', c("#FFFFFF", strtrim(legends.sdl[[j]]$COLOR, 7)), ceiling(c(0, levs)), c("NODATA", Group), c(0.0, rep(.7, length(legends.sdl[[j]]$COLOR))))
    } else {
       txt <- sprintf('<ColorMapEntry color="%s" quantity="%.1f" label="%s" opacity="%.1f"/>', c("#FFFFFF", strtrim(legends.sdl[[j]]$COLOR, 7)), c(0, levs), c("NODATA", Group), c(0.0, rep(.7, length(legends.sdl[[j]]$COLOR))))
    }
  }
  parseXMLAndAdd(txt, l7)
  saveXML(l1, sld.file)
  close(sld.file)
}  

load("tiles.pol.rda")
## GeoJeson of the tiles:
sel = tiles.pol$TNAMES %in% sapply(rda.lst, function(x){strsplit(x, "/")[[1]][[2]]})
tiles.geo <- tiles.pol[sel,]
writeOGR(tiles.geo["TNAMES"], "tiles.geojson", "tiles", driver="GeoJSON")

## list of layers (METADATA):
s.list <- data.frame(Filename=gtifout.lst, Description=NA, stringsAsFactors=FALSE)
for(j in names(soil.legends)){
  if(j %in% c("TAXOUSDA","TAXGWRB")){
    for(k in grep(s.list$Filename, pattern=j)){
      tname <- strsplit(s.list$Filename[k], "_")[[1]][2]
      if(is.na(tname)){ 
        Description <- paste("Predicted classes (as integers)")
      }else{
        Description <- paste("Predicted probability of occurence (0-100%) for class ", tname, sep="")
      }
      s.list[k,"Description"] <- Description
    }
  } else {
    for(k in grep(s.list$Filename, pattern=j)){
      typei <- strsplit(s.list$Filename[k], "_")[[1]][3]
      depthi <- as.integer(strsplit(strsplit(s.list$Filename[k], "_")[[1]][2], "sd")[[1]][2])
      if(j == "ORCDRC"){
        Description <- paste("Soil organic carbon content (fine earth fraction) in permilles (", ifelse(typei=="M", "mean estimate", ifelse(typei=="L", "90% lower limit", "90% upper limit")), ") for ", get("stdepths", envir = GSIF.opts)[depthi]*-100, " cm depth", sep="")  
      }
      if(j == "PHIHOX"){
        Description <- paste("Soil pH x 10 in H2O (", ifelse(typei=="M", "mean estimate", ifelse(typei=="L", "90% lower limit", "90% upper limit")), ") for ", get("stdepths", envir = GSIF.opts)[depthi]*-100, " cm depth", sep="")  
      }
      if(j == "CRFVOL"){
        Description <- paste("Coarse fragments volumetric in percent (", ifelse(typei=="M", "mean estimate", ifelse(typei=="L", "90% lower limit", "90% upper limit")), ") for ", get("stdepths", envir = GSIF.opts)[depthi]*-100, " cm depth", sep="")  
      }
      if(j %in% c("SNDPPT", "SLTPPT", "CLYPPT")){
        Description <- paste("Soil texture fraction ", ifelse(j=="SNDPPT", "sand", ifelse(j=="SLTPPT", "silt", "clay")), " in percent (", ifelse(typei=="M", "mean estimate", ifelse(typei=="L", "90% lower limit", "90% upper limit")), ") for ", get("stdepths", envir = GSIF.opts)[depthi]*-100, " cm depth", sep="")  
      }
      if(j == "BLD"){
        Description <- paste("Bulk density in kg / cubic-meter (", ifelse(typei=="M", "mean estimate", ifelse(typei=="L", "90% lower limit", "90% upper limit")), ") for ", get("stdepths", envir = GSIF.opts)[depthi]*-100, " cm depth", sep="")  
      }
      if(j == "CEC"){
        Description <- paste("Cation exchange capacity in cmolc/kg (", ifelse(typei=="M", "mean estimate", ifelse(typei=="L", "90% lower limit", "90% upper limit")), ") for ", get("stdepths", envir = GSIF.opts)[depthi]*-100, " cm depth", sep="")  
      }
    s.list[k,"Description"] <- Description
    }
  }  
}
write.csv(s.list, "SoilGrids1km_list.csv")

## Create PNGs whole world at 5 km resolution:
makeWorldPNG <- function(tif, tifout=tempfile(), png.width=7200, png.height=png.width/2, cuts, class.labels, pal, resample="near", bbox=matrix(c(-180,-90,180,90), nrow=2), mvFlag, over.lines=world_country){
  ## resample if necessary:
  system(paste(gdalwarp, tif, tifout, '-r', resample, '-te', bbox[1,1], bbox[2,1],bbox[1,2], bbox[2,2], '-tr', 360/png.width, 180/png.height))
  img <- readGDAL(tifout)
  if(!missing(mvFlag)){
    img$band1 <- ifelse(img$band1==mvFlag, NA, img$band1)
  }
  ## classify to the legend:
  img$out <- cut(img$band1, breaks = cuts, labels = class.labels)
  ## plot the image
  require(RSAGA)
  raster_name = set.file.extension(tif, ".png")
  png(filename = raster_name, bg = "transparent", type="cairo-png", width=png.width, height=png.height)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  raster::image(img["out"], col = pal, frame.plot = FALSE, main="")
  lines(over.lines, col="black")
  legend("left", legend=rev(class.labels), fill=rev(pal), cex=3)
  dev.off()
}

sel.lst <- c("ORCDRC","PHIHOX","CLYPPT","BLD")
depths.lst <- 1:6
daten = "25_nov_2013"
cl.lst <- c("L","M","U")
for(k in sel.lst){
  cuts = c(soil.legends[[k]]$MIN, soil.legends[[k]]$MAX[nrow(soil.legends[[k]])])
  class.labels = paste(soil.legends[[k]]$MIN, " - ", soil.legends[[k]]$MAX)
  pal = soil.legends[[k]]$COLOR
  for(j in depths.lst){
    for(cl in cl.lst){
      fname = paste(k, "_sd", j, "_", cl, "_", daten, ".tif", sep="")
      if(!file.exists(set.file.extension(fname, ".png"))){
        makeWorldPNG(tif=fname, cuts=cuts, class.labels=class.labels, pal=pal, mvFlag=0)
      }
    }
  }
}


## KML plots:
tifout.lst <- read.csv("tifout.lst.csv")
library(raster)
library(RSAGA)

for(k in 1:nrow(tifout.lst)){
  ## list all prediction locations (tiles):
  spred = strsplit(paste(tifout.lst[k,2]), "_")[[1]][3]=="M"
  sdep = any(strsplit(paste(tifout.lst[k,2]), "_")[[1]][2] %in% c("sd1","sd3","sd5"))
  svar = any(strsplit(paste(tifout.lst[k,2]), "_")[[1]][1] %in% c("ORCDRC","PHIHOX","CLYPPT","BLD","CRFVOL"))
  if(spred&sdep&svar){
    tvar.name <- strsplit(paste(tifout.lst[k,2]), "_")[[1]][1:2]
    tif3D.lst <- normalizePath(list.files(pattern=glob2rx(paste(tvar.name[1], "_", tvar.name[2], "_M_1km_*.tif", sep="")), recursive=TRUE, full.names=TRUE))
    c.dir <- getwd()
    ## start writing the KML:
    setwd(kml.outdir)
    kml.file <- paste(tvar.name[1], "_", tvar.name[2], "_1km.kml", sep="")
    kml_open(kml.file)
    ## add description:
    desc <- "<iframe width =\"1200\" height=\"740\" frameborder=\"0\" src=\"http://www.isric.org/content/soilgrids/\">www.SoilGrids.org</iframe>"
    netw <- paste("<name>soilgrids1km:", paste(tvar.name[1], tvar.name[2], sep="_"), "</name>
  		<open>1</open>
	 	  <Url>
			 <href>http://85.214.37.88:80/geoserver/soilgrids1km/wms?height=1024&amp;width=1024&amp;layers=soilgrids1km:", paste(tvar.name[1], tvar.name[2], sep="_"), "&amp;request=GetMap&amp;service=wms&amp;styles=", tvar.name[1], "&amp;format_options=SUPEROVERLAY:false;KMPLACEMARK:false;KMSCORE:40;KMATTR:true;&amp;srs=EPSG:4326&amp;format=application/vnd.google-earth.kmz&amp;transparent=false&amp;version=1.1.1&amp;</href>
			 <viewRefreshMode>onStop</viewRefreshMode>
			 <viewRefreshTime>1</viewRefreshTime>
		  </Url>", sep="")
    kml.out <- get("kml.out", envir=plotKML.fileIO)
    description_txt <- sprintf('<description><![CDATA[%s]]></description>', desc)
    parseXMLAndAdd(description_txt, parent=kml.out[["Document"]])  
    network_txt <- sprintf('<NetworkLink>%s</NetworkLink>', netw)
    parseXMLAndAdd(network_txt, parent=kml.out[["Document"]])
    assign('kml.out', kml.out, envir=plotKML.fileIO)
    ## add logo and a legend:
    logo = "http://www.isric.org/sites/default/files/SoilGrids_banner.png"
    legenda = paste("http://worldgrids.org/maps/", tvar.name[1], "_legend.png", sep="")
    kml_screen(image.file = logo, position = "UR", sname = "ISRIC logo")
    kml_screen(image.file = legenda, position = "UL", sname = paste(tvar.name[1], "legend"))
    #pal = soil.legends[[tvar.name[1]]]$COLOR
    #for(i in 1:length(tif3D.lst)){
#      raster_name = set.file.extension(basename(tif3D.lst[[i]]), ".png")
#      if(!file.exists(raster_name)){
#        tif <- readGDAL(tif3D.lst[[i]])
#        tif$value <- cut(tif@data[,1], breaks=c(soil.legends[[tvar.name[1]]]$MIN[1], soil.legends[[tvar.name[1]]]$MAX), include.lowest=TRUE)
#        pw = tif@grid@cells.dim[1]*3
#        ph = tif@grid@cells.dim[2]*3
#        kml_layer(tif, colour=value, colour_scale=pal, html.table=desc, raster_name=raster_name, plot.legend=FALSE, png.width=pw, png.height=ph)
#      }
#    }
    kml_close(kml.file)
    setwd(c.dir)
  }
}

## Sinusoidal projection:
sin.csy <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
pol <- Polygons(list(Polygon(data.frame(x=c(-179.9,179.9,rep(179.9,200),179.9,-179.9,rep(-179.9,200),-179.9),y=c(83.5,83.5,seq(83.5,-58,length.out=200),-58,-58,seq(-58,83.5,length.out=200),83.5)))), "id1")
ll.bbox <- SpatialPolygonsDataFrame(SpatialPolygons(list(pol), proj4string=CRS("+proj=longlat +datum=WGS84")), data=data.frame(mask=1, row.names="id1"))
sin.bbox <- spTransform(ll.bbox, CRS(sin.csy))
writeOGR(sin.bbox, "sin_bbox.shp", ".", "ESRI Shapefile")

#folder_name = { x <- strsplit(strsplit(raster_name, ".png")[[1]][1], "_")[[1]]; x[length(x)] }

## end of script;