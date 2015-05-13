## Depth to bedrock / depth to saprolite (>80% rocks):

library(aqp)
library(GSIF)
library(plotKML)
library(sp)
library(plyr)
library(raster)
#global
dir <- "E:\\data\\soildata\\depth\\points\\profs\\"
setwd(dir)
##############################################
## SOIL PROFILE DATA
##############################################
## Mexican soil profile observations:
#LIM_ROCA:!!!
#An X appears if the physical limitation of the soil is the bedrock.
#LIM_REGO:
#An X appears if the physical limitation of the soil is the regolith.
#LIM_CEME:????
#An X appears if the physical limitation of the soil is a cementation.
#LIM_NIVF:!!!!
#An X appears if the physical limitation of the soil is the water table.
inegi <- read.csv(".\\soil\\INEGIv1_1km.csv")
inegi.xy <- inegi[, c("IDENTIFI", "X_COORD", "Y_COORD")]
coordinates(inegi.xy) <- ~ X_COORD + Y_COORD
mx.csy <- paste0("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=", 
    "-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs")
proj4string(inegi.xy) <- CRS(mx.csy)
inegi.ll <- as.data.frame(reproject(inegi.xy))
inegi.ll[!(inegi$X_COORD > 0&inegi$Y_COORD > 0), 2:3] <- NA
inegi$LONWGS84 <- inegi.ll$x
inegi$LATWGS84 <- inegi.ll$y
## Depth to bedrock i.e. 'R' horizon:
sel.r1 <- grep(pattern="x", inegi$LIM_ROCA, ignore.case=FALSE) 
sel.r2 <- grep(pattern="R", inegi$HSIMBOLO, ignore.case=FALSE)
#inegi$HSIMBOLO[sel.r2]
sel.R  <- unique(sel.r1, sel.r2)
inegi$BDRICM <- NA
inegi$BDRICM[sel.R] <- inegi$PROFUNDI[sel.R]
#View(inegi[, c("IDENTIFI", "PROFUNDI", "HSIMBOLO", "LIM_ROCA", "BDRICM")])
bdr.d <- aggregate(inegi$BDRICM, list(inegi$IDENTIFI), max, na.rm=TRUE)
names(bdr.d) <- c("IDENTIFI", "BDRICM")
bdr.d$BDRICM <- ifelse(is.infinite(bdr.d$BDRICM), NA, bdr.d$BDRICM)
MX.spdb <- plyr::join(bdr.d, inegi[, c("IDENTIFI", "LONWGS84", "LATWGS84")], 
            type="left", match ="first")
MX.spdb$SOURCEID <- paste0("MX_", MX.spdb$IDENTIFI)
MX.spdb <- MX.spdb[complete.cases(MX.spdb), ]
#ploting
MX.ll <- MX.spdb[,2:4]
coordinates(MX.ll) <- ~ LONWGS84+LATWGS84
proj4string(MX.ll) <- CRS("+proj=longlat +datum=WGS84")
#plotKML(MX.ll)
#hist(MX.spdb$BDRICM,100)

## USA data:
load(".\\soil\\NCSS_all.rda")
sel.r  <- grep("R", NCSS_all$horizons$hzn_desgn, ignore.case = FALSE)
sel.rn <- grep("IR", NCSS_all$horizons$hzn_desgn, ignore.case = FALSE) 
sel.r  <- sel.r[!(sel.r %in% sel.rn)]   
length(sel.r)
#View(sort(NCSS_all$horizons[sel.r, c("hzn_desgn")]
#NCSS_all$horizons <- NCSS_all$horizons[sel.r, ]
NCSS_all$horizons$BDRICM <- NA
## fix typos:
max.d <- aggregate(NCSS_all$horizons$hzn_bot, 
        list(NCSS_all$horizons$site_key), max, na.rm=TRUE)
names(max.d)[1] <- "site_key"
NCSS_all$horizons$hzn_topF <- ifelse(NCSS_all$horizons$layer_sequence>1&
            NCSS_all$horizons$hzn_top==0, 
            plyr::join(NCSS_all$horizons["site_key"], max.d)$x, 
            NCSS_all$horizons$hzn_top) 
NCSS_all$horizons$BDRICM[sel.r] <- NCSS_all$horizons$hzn_topF[sel.r]
bdr2.d <- aggregate(NCSS_all$horizons$BDRICM[sel.r], 
            list(NCSS_all$horizons$site_key[sel.r]), min, na.rm=TRUE)
names(bdr2.d) <- c("site_key", "BDRICM")
bdr2.d$BDRICM <- ifelse(is.infinite(bdr2.d$BDRICM), NA, bdr2.d$BDRICM)
US.spdb <- plyr::join(bdr2.d, NCSS_all$sites[, c("site_key", "LON", "LAT")], 
            type="left", match ="first")
names(US.spdb)[3:4] <- c("LONWGS84", "LATWGS84")
US.spdb$SOURCEID <- paste0("US_", US.spdb$site_key)
US.spdb <- US.spdb[complete.cases(US.spdb), ]
#BDRICM == 1189 is checked
US.ll <- US.spdb[,2:4]
coordinates(US.ll) <- ~ LONWGS84+LATWGS84
proj4string(US.ll) <- CRS("+proj=longlat +datum=WGS84")
#plotKML(US.ll)
#hist(US.spdb$BDRICM,200)

## Canada:
load(".\\SOIL\\CanSIS.rda")
#H2O_TABLE water table depth
#PFROST_DEPTH Permafrost - Depth to
#FAMILY_DEPTH class: Family Criteria - Soil Depth 
names(CanSIS) <- c("sites", "horizons")
summary(CanSIS$horizons$PS_VCSAND)  ## Coarse fragments
## depth to bedrock (R horizon):
summary(as.factor(CanSIS$sites$FAMILY_DEPTH))
View(CanSIS$sites[!CanSIS$sites$FAMILY_DEPTH ==" ", 
        c("SOURCEID","FAMILY_DEPTH")])
#bulk density > 2.67
#sel.b <- which(!is.na(CanSIS$horizons$BD) &CanSIS$horizons$BD>2.67)
sel.r  <- grep("R", CanSIS$horizons$HORIZON, ignore.case=FALSE)
sel.rn <- c(grep ("DRIFT", CanSIS$horizons$HORIZON, ignore.case=FALSE),
           grep ("GIR", CanSIS$horizons$HORIZON, ignore.case=FALSE)) 
sel.r  <- sel.r[!(sel.r %in% sel.rn)]    
CanSIS$horizons[!is.na(CanSIS$horizons$BD) & CanSIS$horizons$BD>2.67, ]     
#CanSIS$horizons[sel.r,"HORIZON"]
#not nure, but use the deepest first
sel.n <- which(CanSIS$horizons$SOURCEID %in% CanSIS$sites$SOURCEID[grep("lit", 
        CanSIS$sites$FAMILY_DEPTH, ignore.case=TRUE)])         
View(CanSIS$horizons[sel.r,c("SOURCEID","UDEPTH.x","LDEPTH.x","HORIZON")])
CanSIS$horizons$BDRICM <- NA
CanSIS$horizons$BDRICM[sel.n] <- CanSIS$horizons$LDEPTH.x[sel.n]
CanSIS$horizons$BDRICM[sel.r] <- CanSIS$horizons$UDEPTH.x[sel.r]
bdr3.d <- aggregate(CanSIS$horizons$BDRICM, list(CanSIS$horizons$SOURCEID), max, na.rm=TRUE)
names(bdr3.d) <- c("SOURCEID", "BDRICM")
bdr3.d$BDRICM <- ifelse(is.infinite(bdr3.d$BDRICM), NA, bdr3.d$BDRICM)
CA.spdb <- plyr::join(bdr3.d, CanSIS$sites[,c("SOURCEID", "LONWGS84", "LATWGS84")], type="left", match ="first")
CA.spdb <- CA.spdb[complete.cases(CA.spdb), ]
CA.ll <- CA.spdb[,2:4]
coordinates(CA.ll) <- ~ LONWGS84+LATWGS84
proj4string(CA.ll) <- CRS("+proj=longlat +datum=WGS84")
#plotKML(CA.ll)
#hist(CA.spdb$BDRICM,40)


##############################################
## WELL DATA (DRILLINGS)
##############################################

## Wells data compiled by Wei:
dir <- "E:\\data\\soildata\\depth\\points\\profs\\"
setwd(dir)
us <- read.csv(".\\well\\wells_us.txt", sep="\t")
ca <- read.csv(".\\well\\wells_ca.txt", sep="\t")
us$BDRICM <- us$D_BR*100
names(us)[2:3] <- c("LONWGS84", "LATWGS84")
us$SOURCEID <- paste0("USWELL_", us$Source)
#hist(us$BDRICM, breaks=6000, col="grey", xlim = c(0,2000)) 
#hist(log1p(us$BDRICM), breaks=40, col="grey") 
ca$BDRICM <- ca$D_BR*100
names(ca)[2:3] <- c("LONWGS84", "LATWGS84")
ca$SOURCEID <- paste0("CAWELL_", ca$Source)
#hist(ca$BDRICM, breaks=6000, col="grey", xlim = c(0,2000)) 
#hist(log1p(ca$BDRICM), breaks=40, col="grey")

#us.ll <- us[log1p(us$BDRICM)<4,c("BDRICM","LONWGS84", "LATWGS84","SOURCEID")]
#coordinates(us.ll) <- ~ LONWGS84+LATWGS84
#proj4string(us.ll) <- CRS("+proj=longlat +datum=WGS84")
#plotKML(us.ll[runif(length(us.ll))<0.05,])

## bind all together:
na.bbox <- matrix(c(-175,14,-58,70), nrow=2)
NA.pnts <- do.call(rbind, list( 
        ca[,c("SOURCEID","LONWGS84","LATWGS84","BDRICM")], 
        us[,c("SOURCEID","LONWGS84","LATWGS84","BDRICM")], 
        US.spdb[,c("SOURCEID","LONWGS84","LATWGS84","BDRICM")], 
        MX.spdb[,c("SOURCEID","LONWGS84","LATWGS84","BDRICM")], 
        CA.spdb[,c("SOURCEID","LONWGS84","LATWGS84","BDRICM")]))
#-------------------------------check values        
#for( i in 2:4)
#{
#    print(colnames(NA.pnts[i]))
#    print(max(NA.pnts[i], na.rm = TRUE))
#    print(min(NA.pnts[i], na.rm = TRUE))
#}
## subset to NorthAmerica continent:
NA.pnts <- NA.pnts[complete.cases(NA.pnts), ]
NA.pnts <- NA.pnts[NA.pnts$LATWGS84 > na.bbox[2,1] & 
        NA.pnts$LATWGS84 < na.bbox[2,2] & NA.pnts$LONWGS84 < na.bbox[1,2] &
        NA.pnts$LONWGS84 > na.bbox[1,1], ]
plot(NA.pnts[, c("LONWGS84", "LATWGS84")])
str(NA.pnts)
hist(log1p(NA.pnts$BDRICM), breaks=40, col="grey")
## Concentrated values around 0, 1 feet and 2 feet!! Are these all ok. 
save(NA.pnts, file="NA.pnts.rda")
NA.sp <- NA.pnts[, c("BDRICM","LONWGS84", "LATWGS84","SOURCEID")]
coordinates(NA.sp) <- ~ LONWGS84+LATWGS84
proj4string(NA.sp) <- CRS("+proj=longlat +datum=WGS84")
save(NA.sp, file="NA.sp.rda")

grd <- vect2rast(NA.sp["BDRICM"], cell.size=.1, bbox=na.bbox)
plot(log1p(raster(grd)), col=SAGA_pal[[1]])
grd.pol <- grid2poly(as(grd, "SpatialPixelsDataFrame"))
kml(grd.pol, colour=log1p(BDRICM), colour_scale=SAGA_pal[[1]])
save.image(paste("prepare.RData", sep = ""))
## end of script;

