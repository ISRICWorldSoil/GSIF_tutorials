# title         : rw_AFSP.R
# purpose       : Reading and writing of Africa Soil Profile data (ca 12,000 profiles);
# reference     : Africa Soil Profiles Database [http://africasoils.net/data/legacyprofile]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL.
# inputs        : "http://www.isric.org/sites/default/files/datasets/AfSP012Qry_ISRIC.7z" DBF tables
# outputs       : R data frame for SoilProfileCollection;
# remarks 1     : afsp.rda is only a subset of the complete database;

library(aqp)
library(plyr)
library(foreign)
library(GSIF)
library(plotKML)

## Download the database:
download.file("http://www.isric.org/sites/all/modules/pubdlcnt/pubdlcnt.php?file=http://www.isric.org/sites/default/files/datasets/AfSP012Qry_ISRIC.7z&nid=557", destfile=paste(getwd(),"AfSP012Qry_ISRIC.7z", sep="/"))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
system("7za -e AfSP012Qry_ISRIC.7z")

# import tables:
Profiles <- read.dbf("AfSP012Qry_ISRIC/GIS_Dbf/AfSP012Qry_Profiles.dbf", as.is=TRUE)
Layers <- read.dbf("AfSP012Qry_ISRIC/GIS_Dbf/AfSP012Qry_Layers.dbf", as.is=TRUE)
#View(Profiles[1:200,])
#View(Layers[1:200,])
#summary(Profiles[,c("WRB06", "WRB06rg", "FAO88", "FAO74", "USDA")])
## depth to bedrock compare with "RockDpth":
Profiles[Profiles$RockDpth=="1.1","ProfileID"]
Layers[Layers$ProfileID=="ML 7076_721",c("UpDpth","LowDpth","HorDes")]
Layers[Layers$ProfileID=="MW lrep_MC17",c("UpDpth","LowDpth","HorDes")]
Layers[Layers$ProfileID=="TZ 13578_MP26",c("UpDpth","LowDpth","HorDes")]

# select columns of interest:
s.lst <- c("ProfileID", "Reliab", "X_LonDD", "Y_LatDD", "XYAccur", "T_Year", "WRB06", "WRB06rg", "FAO88", "USDA", "LocalCls", "Location", "Drain", "RockDpth")
h.lst <- c("ProfileID", "LayerID", "LayerNr", "UpDpth", "LowDpth", "HorDes", "ColorM", "ColorD", "FldTxtr", "CfPc", "Sand", "Silt", "Clay", "BlkDens", "PHH2O", "PHKCl", "EC", "OrgC", "TotalN", "Ecec", "Bsat", "CecSoil", "LabTxtr", "VolAWC", "ExCa", "ExNa", "ExMg", "ExK", "ExBases", "ExAl", "ExAcid")

## correct soil site and horizon data:
sites <- Profiles[!duplicated(Profiles$ProfileID),]
sites <- sites[!is.na(sites$X_LonDD)&!is.na(sites$Y_LatDD), s.lst]
sites$XYAccur <- ifelse(sites$XYAccur==-9999, NA, sites$XYAccur)
sites$T_Year <- ifelse(sites$T_Year==-9999, NA, sites$T_Year)
sites$T_Year <- as.Date(paste(sites$T_Year), format="%Y")
sites$WRB06[which(sites$WRB06=="NA")] <- NA
sites$WRB06rg[which(sites$WRB06rg=="NA")] <- NA
sites$FAO88[which(sites$FAO88=="NA")] <- NA
sites$USDA[which(sites$USDA=="NA")] <- NA
sites$LocalCls[which(sites$LocalCls=="NA")] <- NA
sites$Location[which(sites$Location=="NA")] <- NA
sites$Drain[which(sites$Drain=="NA")] <- NA
sites$RockDpth[which(sites$RockDpth=="NA")] <- NA 
RockDpth.sel <- grep(pattern="^>", sites$RockDpth, ignore.case=FALSE, fixed=FALSE)
RockDpth.selP <- as.numeric(sapply(sites$RockDpth[RockDpth.sel], function(x){strsplit(x, ">")[[1]][2]})) >= 1.8
sites$RockDpth.f <- 100 * as.numeric(sites$RockDpth)
sites$RockDpth.f[RockDpth.sel[RockDpth.selP]] <- 200
sites$RockDpth <- NULL
summary(!is.na(sites$WRB06rg))
summary(!is.na(sites$USDA))
summary(!is.na(sites$LocalCls))
summary(!is.na(sites$Drain))
summary(!is.na(sites$RockDpth.f))
str(sites)
## plot "RockDpth"
RockDpth <- sites[!is.na(sites$RockDpth.f),c("RockDpth.f","X_LonDD","Y_LatDD","ProfileID")]
coordinates(RockDpth) <- ~ X_LonDD+Y_LatDD
proj4string(RockDpth) <- "+proj=longlat +datum=WGS84"
plotKML(RockDpth)

horizons <- Layers[,h.lst]
horizons$LayerNr <- as.integer(horizons$LayerNr)
horizons$FldTxtr[which(horizons$FldTxtr=="NA")] <- NA
horizons$HorDes[which(horizons$HorDes=="NA")] <- NA
horizons$ColorM[which(horizons$ColorM=="NA")] <- NA
horizons$ColorD[which(horizons$ColorD=="NA")] <- NA
horizons$ColorM[which(horizons$ColorM=="<Null>")] <- NA
horizons$ColorD[which(horizons$ColorD=="<Null>")] <- NA
horizons$FldTxtr[which(horizons$FldTxtr=="NA")] <- NA
# fix the Munsell color codes?
# replace "-9999":
for(j in c("CfPc", "Sand", "Silt", "Clay", "BlkDens", "PHH2O", "PHKCl", "EC", "OrgC", "TotalN", "Ecec", "Bsat", "CecSoil", "LabTxtr", "VolAWC", "ExCa", "ExNa", "ExMg", "ExK", "ExBases", "ExAl", "ExAcid")){
  if(is.numeric(horizons[,j])){
    horizons[,j] <- signif(ifelse(horizons[,j]==-9999, NA, horizons[,j]), 2)
  } else {
    horizons[,j] <- ifelse(horizons[,j]==-9999, NA, horizons[,j])
  }
}
str(horizons)

# rename columns:
names(sites) <- c("SOURCEID", "RLBISRIC", "LONWGS84", "LATWGS84", "STDXYZ", "TIMESTRR","TAXNWRB", "TAXGWRB", "TAXNFAO", "TAXNUSDA", "TAXN", "LOCNAME", "DRAINFAO", "RockDpth.f")
names(horizons) <- c("SOURCEID", "SAMPLEID", "LSQINT", "UHDICM", "LHDICM", "HZDFAO", "MCOMNS", "DCOMNS", "TEXMHT", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHOX", "PHIKCL", "ECN", "ORCDRC", "NTO", "ECE", "BST", "CEC", "TEXMHL", "AWCIMM", "ECAPH7", "ENAPH7", "EMGPH7", "EXKPH7", "EXBPH7", "EALKCL", "EACKCL")
## Depth to bedrock i.e. 'R' horizon (there is a column 'RockDpth' in the Profiles table but it is empty):
sel.r <- grep(pattern="^R", horizons$HZDFAO, ignore.case=FALSE, fixed=FALSE)
#sel.r2 <- grep(pattern="*/R", horizons$HZDFAO, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- horizons$CRFVOL >= 95 & !is.na(horizons$CRFVOL)
#sel.r3 <- grep(pattern="CR", horizons$HZDFAO, ignore.case=FALSE, fixed=FALSE)
horizons$BDRICM <- ifelse(horizons$LHDICM > 200, 200, NA)
horizons$BDRICM[sel.r] <- horizons$UHDICM[sel.r]
horizons$BDRICM[sel.r2] <- horizons$UHDICM[sel.r2]
#horizons$BDRICM[sel.r3] <- horizons$LHDICM[sel.r3]
bdr.d <- aggregate(horizons$BDRICM, list(horizons$SOURCEID), min, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
bdr.d$BDRICM <- ifelse(bdr.d$BDRICM==Inf, NA, bdr.d$BDRICM)
sites.m <- merge(sites, bdr.d, all.y=FALSE)
sites.m$BDRICM <- ifelse(is.na(sites.m$BDRICM) & !is.na(sites.m$RockDpth.f), sites.m$RockDpth.f, sites.m$BDRICM)
sites.m$BDRICM <- ifelse(sites.m$BDRICM >200, 200, sites.m$BDRICM)
hist(sites.m$BDRICM)

## some points have 0, 0 coordinates:
sel.xy <- sites.m$LONWGS84==0 & sites.m$LATWGS84==0
## Save object:
afsp <- list(sites=sites.m[!sel.xy,c("SOURCEID", "RLBISRIC", "LONWGS84", "LATWGS84", "STDXYZ", "TIMESTRR","TAXNWRB", "TAXNFAO", "TAXGWRB", "TAXNUSDA","BDRICM","DRAINFAO")], horizons=horizons[,c("SOURCEID", "LSQINT", "HZDFAO", "UHDICM", "LHDICM", "MCOMNS", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHOX", "PHIKCL", "ECN", "ORCDRC", "CEC", "EACKCL")])
afsp$sites$TAXNWRB <- iconv(afsp$sites$TAXNWRB, to="ASCII", sub="byte")
afsp$sites$TAXNWRB <- as.factor(afsp$sites$TAXNWRB)
afsp$sites$TAXNUSDA <- iconv(afsp$sites$TAXNUSDA, to="ASCII", sub="byte")
afsp$sites$TAXNUSDA <- as.factor(afsp$sites$TAXNUSDA)
afsp$sites$TAXNFAO <- iconv(afsp$sites$TAXNFAO, to="ASCII", sub="byte")
afsp$sites$TAXNFAO <- as.factor(afsp$sites$TAXNFAO)
afsp$sites$DRAINFAO <- as.factor(afsp$sites$DRAINFAO)
## check on more time if there are any missing coordinates:
sel.Cs <- !(afsp$sites$LONWGS84 == 0 & afsp$sites$LATWGS84 == 0) & !is.na(afsp$sites$LONWGS84) & !is.na(afsp$sites$LATWGS84)
summary(sel.Cs)
afsp$sites <- afsp$sites[sel.Cs,]
afsp$horizons <- afsp$horizons[which(afsp$horizons$SOURCEID %in% afsp$sites$SOURCEID),]
str(afsp)
save(afsp, file="afsp.rda", compress="xz")

## export depth to bedrock:
BDR <- afsp$sites[!is.na(afsp$sites$BDRICM),c("BDRICM","LONWGS84","LATWGS84","SOURCEID")]
coordinates(BDR) <- ~ LONWGS84 + LATWGS84
proj4string(BDR) <- "+proj=longlat +datum=WGS84"
plotKML(BDR)

## save sample data for WorlSoilProfiles
sites.m$SOURCEDB = "AfSPDB"
wsp10 <- list(sites=sites.m[!sel.xy,c("SOURCEID", "LONWGS84", "LATWGS84", "TIMESTRR", "TAXGWRB", "TAXNUSDA", "BDRICM", "SOURCEDB")], horizons=horizons[,c("SOURCEID", "UHDICM", "LHDICM", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHOX", "PHIKCL", "ECN", "ORCDRC", "CEC")])
str(wsp10)
wsp10$sites$TAXGWRB <- as.character(wsp10$sites$TAXGWRB)
wsp10$sites$TAXNUSDA <- as.character(wsp10$sites$TAXNUSDA)
lapply(as.list(wsp10$sites), function(x){sum(!is.na(x))})
lapply(as.list(wsp10$horizons), function(x){sum(!is.na(x))})
save(wsp10, file="D:/SPDB/wsp10.rda", compress="xz")

afsis.leg <- list(sites=sites.m[,c("SOURCEID", "LONWGS84", "LATWGS84", "TIMESTRR", "TAXGWRB", "TAXNUSDA", "BDRICM", "SOURCEDB", "LOCNAME", "DRAINFAO")], horizons=horizons[,c("SOURCEID", "SAMPLEID", "UHDICM", "LHDICM", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHOX", "PHIKCL", "ECN", "ORCDRC", "NTO", "ECE", "BST", "CEC", "TEXMHT", "ECAPH7", "ENAPH7", "EMGPH7", "EXKPH7", "EXBPH7", "EALKCL", "EACKCL")])
afsis.leg$sites$TAXGWRB <- as.character(afsis.leg$sites$TAXGWRB)
afsis.leg$sites$TAXNUSDA <- as.character(afsis.leg$sites$TAXNUSDA)
save(afsis.leg, file="afsis.leg.rda", compress="xz")

## try to plot soil colors:
afsp.spc <- join(afsp$horizons, afsp$sites, type='inner')
depths(afsp.spc) <- SOURCEID ~ UHDICM + LHDICM
site(afsp.spc) <- ~ LONWGS84+LATWGS84+STDXYZ+TIMESTRR+TAXNWRB+TAXGWRB+TAXNUSDA
coordinates(afsp.spc) <- ~ LONWGS84 + LATWGS84
proj4string(afsp.spc) <- "+proj=latlong +datum=WGS84"
str(afsp.spc)
# get colors:
sc <- data.frame(mcolor=afsp.spc@horizons$MCOMNS, index=1:length(afsp.spc@horizons$MCOMNS))
# reformat color codes:
sc$Munsell <- sub(" ", "", sub("/", "_", sc$mcolor))
hue.lst <- expand.grid(c("2.5", "5", "7.5", "10"), c("YR","GY","BG","YE","YN","YY","R","Y","B","G"))
hue.lst$mhue <- paste(hue.lst$Var1, hue.lst$Var2, sep="") 
for(j in hue.lst$mhue[1:24]){ 
  sc$Munsell <- sub(j, paste(j, "_", sep=""), sc$Munsell, fixed=TRUE)
  sc$Munsell <- sub("__", "_", sc$Munsell, fixed=TRUE) 
}
# match MunsellRGB table:
load(file("http://gsif.isric.org/lib/exe/fetch.php?media=munsell_rgb.rdata"))
Munsell.rgb <- merge(sc, munsell.rgb, by="Munsell", all.x=TRUE, all.y=FALSE)
Munsell.rgb$R <- ifelse(is.na(Munsell.rgb$R), 255, Munsell.rgb$R)
Munsell.rgb$G <- ifelse(is.na(Munsell.rgb$G), 255, Munsell.rgb$G)
Munsell.rgb$B <- ifelse(is.na(Munsell.rgb$B), 255, Munsell.rgb$B)
Munsell.rgb <- Munsell.rgb[order(Munsell.rgb$index),]
afsp.spc@horizons$m_color <- rgb(red=Munsell.rgb$R, green=Munsell.rgb$G, blue=Munsell.rgb$B, maxColorValue = 255)

# kml(afsp.spc[1:200], color.name ="m_color", var.name="ORCDRC")
## It would take a lot of time to export all profiles at once!!

afsp.df <- as.data.frame(afsp.spc)
write.csv(afsp.df, file="afsp.csv")
coordinates(afsp.df) <- ~ LONWGS84 + LATWGS84
proj4string(afsp.df) <- "+proj=latlong +datum=WGS84"
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(afsp.df[,"ORCDRC_A"], file = "afsp_ORCDRC_A.kml", colour = ORCDRC_A, shape = shape, labels = ORCDRC_A, size = log1p(ORCDRC_A))
# takes 15 minutes!!

# end of script;