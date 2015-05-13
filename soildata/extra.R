
# load the Munsell color codes table:
load(file("http://globalsoilmap.net/data/munsell_rgb.RData"))
str(munsell.rgb)
# pick a random colour:
munsell.rgb[round(runif(1)*2350, 0),]
# convert to HSI values:
library(colorspace)
as(RGB(R=munsell.rgb[1007,"R"]/255, G=munsell.rgb[1007,"G"]/255, B=munsell.rgb[1007,"B"]/255), "HSV")

# plot the ISRIC WISE soil colors in Google Earth:

library(RODBC)
# Before you can connect, you need to set-up an ODBC connection channel on your PC:
# 1. Install ODBC driver for windows: [http://dev.mysql.com/downloads/connector/odbc/]
# 2. Set-up an ODBC connection; call it "WISE"
#   db: spatiala_WISE
#   username: spatiala_soildb
#   pw: readonly_12
# Now connect to the database:
cWISE <- odbcConnect(dsn="WISE", DBMSencoding="UTF-8")
odbcGetInfo(cWISE)
sqlTables(cWISE)$TABLE_NAME
# get soil properties:
mcolor.tbl <- sqlQuery(cWISE, "SELECT wise3_id, honu, topdep, botdep, mcolor FROM horizon", na.strings=c("<NA>", "NA"))
str(mcolor.tbl)


# levels(mcolor.tbl$mcolor)
summary(mcolor.tbl$mcolor)
# reformat color codes:
mcolor.tbl$Munsell <- sub(" ", "", sub("/", "_", mcolor.tbl$mcolor))
hue.lst <- expand.grid(c("2.5", "5", "7.5", "10"), c("YR","GY","BG","YE","YN","YY","R","Y","B","G"))
hue.lst$mhue <- paste(hue.lst$Var1, hue.lst$Var2, sep="") 
for(j in hue.lst$mhue[1:24]){ 
mcolor.tbl$Munsell <- sub(j, paste(j, "_", sep=""), mcolor.tbl$Munsell, fixed=TRUE) 
}
mcolor.tbl$depth <- mcolor.tbl$topdep + (mcolor.tbl$botdep-mcolor.tbl$topdep)/2 
# merge munsell colors and RGB values:
mcolor.RGB <- merge(mcolor.tbl, munsell.rgb, by="Munsell")
str(mcolor.RGB)
mcolor.RGB$Rc <- round(mcolor.RGB$R/255*100, 0)
mcolor.RGB$Gc <- round(mcolor.RGB$G/255*100, 0)
mcolor.RGB$Bc <- round(mcolor.RGB$B/255*100, 0)
# reformat colors suited for Google Earth:
mcolor.RGB$FBGR <- paste("ff", ifelse(mcolor.RGB$Bc<10, paste("0", mcolor.RGB$Bc, sep=""), ifelse(mcolor.RGB$Bc==100, "ff", mcolor.RGB$Bc)), ifelse(mcolor.RGB$Gc<10, paste("0", mcolor.RGB$Gc, sep=""), ifelse(mcolor.RGB$Gc==100, "ff", mcolor.RGB$Gc)), ifelse(mcolor.RGB$Rc<10, paste("0", mcolor.RGB$Rc, sep=""), ifelse(mcolor.RGB$Rc==100, "ff", mcolor.RGB$Rc)), sep="")

# get coordinates of points:
site.tbl <- sqlQuery(cWISE, "SELECT wise3_id, latit, latdeg, latmin, 	latsec, longi, londeg, lonmin, lonsec FROM site", na.strings=c("<NA>","NA"))
# reformat coordinates:
site.tbl$latsec <- ifelse(is.na(site.tbl$latsec), 0, site.tbl$latsec)
site.tbl$lonsec <- ifelse(is.na(site.tbl$lonsec), 0, site.tbl$lonsec)
site.tbl$latmin <- ifelse(is.na(site.tbl$latmin), 0, site.tbl$latmin)
site.tbl$lonmin <- ifelse(is.na(site.tbl$lonmin), 0, site.tbl$lonmin)
site.tbl <- subset(site.tbl, !is.na(site.tbl$latdeg)|!is.na(site.tbl$londeg))  # no point in using profiles that have no geo-reference!
# concatenate to WGS84 coordinates:
site.tbl$LAT <- cols2dms(site.tbl$latdeg, site.tbl$latmin, site.tbl$latsec, site.tbl$latit)
site.tbl$LON <- cols2dms(site.tbl$londeg, site.tbl$lonmin, site.tbl$lonsec, site.tbl$longi)
str(site.tbl)

# convert to a 3D table:
mcolor.XYd <- merge(site.tbl[,c("wise3_id", "LON", "LAT")], mcolor.RGB[,c("wise3_id", "honu", "depth", "Munsell", "FBGR")], by="wise3_id") 
str(mcolor.XYd)