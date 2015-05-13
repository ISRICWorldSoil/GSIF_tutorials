rm(list = ls(all = TRUE))
library(foreign)
library(maptools)
library(rgdal)
library(stringr)
#library(RODBC)
library(plotKML)
library(RSAGA)
library(XML)
#dir_out <- "E:\\data\\soildata\\depth\\points\\"
#load(paste0(dir_out, "us.RData"))

#define global, change if needed
options(stringsAsFactors = FALSE) # not usful for readShapePoints, 
                                # fix s_n == 10
# change dir* if needed
dir_f <- "E:\\data\\soildata\\depth\\points\\codegsifb\\head\\"
source(paste(dir_f, "functions.r", sep = ""))
source(paste(dir_f, "crawl.r", sep = ""))
dir1  <- "E:\\data\\soildata\\depth\\points\\states\\"
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy") 
setwd(dir1)
dirs   <- shell("dir * /B", intern = T)
states <- cbind(1:length(dirs), dirs)
colnames(states) <- c("s_code", "states")
states
dir_out <- "E:\\data\\soildata\\depth\\points\\"#      s_code states             
# [1,] "1"    "Alaska"      
# [2,] "2"    "Arizona"     
# [3,] "3"    "Indiana"     
# [4,] "4"    "Iowa"        
# [5,] "5"    "Kansas"      
# [6,] "6"    "Kentuky"     
# [7,] "7"    "Maine"       
# [8,] "8"    "Minnesota"   
# [9,] "9"    "Missouri"    
#[10,] "10"   "Nevada"      
#[11,] "11"   "NewHampshire"
#[12,] "12"   "Newyork"     
#[13,] "13"   "Ohio"        
#[14,] "14"   "Other"       
#[15,] "15"   "Pennsylvania"
#[16,] "16"   "Tennessee"   
#[17,] "17"   "Vermont"     

swap2n <- function(a)
{
    t <- a[, 1]
    a[, 1] <- a[, 2]
    a[, 2] <- t
    a
}


rec.lst <- as.list(rep(NA, length(dirs)))
names(rec.lst) <- dirs

#------------------------all in feet
#-------------------------------Alaska
#0 are not kept, 2043 out of 3320 are 0
s_n <- 1 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- read.csv(".\\tmp\\all.csv")
tmp <- tmp[tmp[, "Bedrock.Depth"] > 0 & !is.na(tmp[, "Bedrock.Depth"]), 
        c("Meridian", "Township", "Range", "Section", "Quarter.Sections",
        "Bedrock.Depth", "Static.Water.Depth", "Hole.Depth")]
tmp <- tmp[complete.cases(tmp[, c("Meridian", "Township", "Range", 
        "Section")]), ]
tmp <- tmp[tmp[,"Township"]!= "", ]
####takes too long to read alt_first?????????????????????        
#tmp2 <- readShapePoly(".\\tmp\\alt_AK\\alt_first")
#Copper River Meridian 1905 61¡ã49¡ä04¡åN 145¡ã18¡ä37¡åW Alaska 
#Fairbanks Meridian 1910 64¡ã51¡ä50.048¡åN 147¡ã38¡ä25.94¡åW Alaska 
#Kateel River Meridian 1956 65¡ã26¡ä16.374¡åN 158¡ã45¡ä31.01¡åW Alaska 
#Seward Meridian 1911 60¡ã07¡ä37¡åN 149¡ã21¡ä26¡åW Alaska 
#Umiat Meridian 1956 69¡ã23¡ä29.654¡åN 152¡ã00¡ä04.55¡åW Alaska 
base <- data.frame(n = c(12, 13, 44, 28, 45), 
        Meridian = c("C", "F", "K", "S", "U"))
tmp <- join(tmp, base[, c("Meridian", "n")], type = "left", match = "first")
tmp$Township2 <- str_c(str_sub(tmp$Township, 1, -2), "0", 
                str_sub(tmp$Township, -1))
tmp$Range2 <- str_c(str_sub(tmp$Range, 1, -2), "0", str_sub(tmp$Range, -1))
tmp$lndkey <- str_c("AK", tmp$n, "T", tmp$Township2, tmp$Range2)
tmp$mtrs   <- str_c(tmp$lndkey, ifelse(tmp$Section < 10, "00", "0"), 
            tmp$Section)

# get tht bounding box of a mtrs
shps <- shell("dir .\\tmp\\alt_AK\\split\\*.shp /B", intern = T)
shps <- str_replace(shps, ".shp", "")
tmp5 <- NULL
for(i in 1:length(shps))
{
    print(i)
    #split the file by SAGA 10*4
    tmp2 <- readShapePoly(paste0(".\\tmp\\alt_AK\\split\\", shps[i]))
    #get lat and long according to indkey from shp file
    tmp2$mtrs <- as.character(tmp2$mtrs)
    
    tmp3 <- NULL
    for(j in which(tmp2$mtrs %in% tmp$mtrs))
    {
    
    
        mtrs   <- tmp2$mtrs[j]
        coords <- tmp2@polygons[[j]]@Polygons[[1]]@coords 
        tmp4 <- getcoords(coords)
        tmp4 <- data.frame(xmax = tmp4[1], xmin = tmp4[2],
                ymax = tmp4[3], ymin = tmp4[4])  
        tmp3 <- rbind(tmp3, cbind(mtrs, tmp4))    
    }
    
    tmp5 <- rbind(tmp5, tmp3)
}
tmp2 <- aggregate(tmp5[, 2:5], by = list(tmp5$mtrs), FUN = "max" )
tmp2 <- rename(tmp2, c(Group.1 = "mtrs"))
tmp  <- tmp[tmp$mtrs %in% tmp2$mtrs, ]
tmp  <- join(tmp, tmp2, by = "mtrs", type = "left")
save("tmp.rda",list("tmp"))
# get bounding box of Quarter.Sections
for(i in 1:dim(tmp)[1])
{

    tmp3 <- unlist(str_split(tmp$Quarter.Sections[i], " ")) 
    if(!(length(tmp3 == 1) & tmp3 == "")) 
    {
        for(j in length(tmp3):1)
        {
            tmp[i, c("xmax", "xmin", "ymax", "ymin")] <- 
                qsec_box(tmp[i, c("xmax", "xmin", "ymax", "ymin")], tmp3[j])
        }
    }
}
# get long, lat and acc
tmp$long <- (tmp$xmax + tmp$xmin)/2
tmp$lat  <- (tmp$ymax + tmp$ymin)/2
tmp$xdis <- (tmp$xmax - tmp$xmin) * 40075000 * sin((90-(tmp$ymax + tmp$ymin)/2)*pi/180)/360
tmp$ydis <- (tmp$ymax - tmp$ymin) * 110940
tmp$acc  <- (tmp$xdis + tmp$ydis)/4
tmp$acc  <- ifelse(tmp$acc < 55, 50, ifelse(tmp$acc < 110, 100, ifelse(tmp$acc < 210, 
        200, ifelse(tmp$acc < 420, 400, 800))))
tmp$acc  <- as.character(tmp$acc) 
tmp  <- tmp[, c("long", "lat", "Bedrock.Depth", 
    "Static.Water.Depth", "Hole.Depth", "acc")]
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"0 out of 1249 is 0"
del_unzip()

#to_deg <- function(t)
#{
#  
#    t <- unlist(str_split(t, "¡ã"))
#    t <- unlist(str_split(t, "¡ä"))
#    t <- str_replace(t, pattern = "¡å", replacement = "")
#    t <- as.numeric(t)
#    t <- t[1] + t[2]/60 + t[3]/60/60
#    t
#}

#base <- data.frame(n = c(12, 13, 44, 28, 45), Meridian = c("C", "F", "K", "S", "U"), 
#        long = c("145¡ã18¡ä37¡å", "147¡ã38¡ä25.94¡å", "158¡ã45¡ä31.01¡å",
#                "149¡ã21¡ä26¡å", "152¡ã00¡ä04.55¡å"),
#        lat  = c("61¡ã49¡ä04¡å", "64¡ã51¡ä50.048¡å", "65¡ã26¡ä16.374¡å",
#                "69¡ã23¡ä29.654¡å", "69¡ã23¡ä29.654¡å"))
#base$long <- sapply(base$long, to_deg)
#base$lat  <- sapply(base$lat, to_deg)
#tmp <- join(tmp, base[, c("Meridian", "n")], type = "left", match = "first")
#rsaga.sgrd.to.esri(".\\alt_ak\\alt_AK\\alt_twnsh2.sgrd", 
#    "t", ".\\alt_ak\\alt_AK\\")
#tmp3<- readGDAL(fname=".\\alt_ak\\alt_AK\\t.asc") 
#rsaga.get.modules(libs = "shapes_tools")
#rsaga.get.modules(Save Grid As)
#rsaga.get.usage(lib = "shapes_tools", module = "Split Shapes Layer")
#t<-rsaga.geoprocessor(lib = "shapes_tools", module = "Split Shapes Layer", 
#    param = list(SHAPES = ".\\tmp\\alt_AK\\alt_twnshp.shp",
#        NX = 2, NY =2))
      
#-------------------------------Arizona
#0 are kept
s_n <- 2 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
dir2 <- "\\tmp\\DGM-52(DepthToBedrock)v.1.0\\Shapefiles\\"
tmp  <- readShapePoints(paste(getwd(),dir2,"dtb_Boreholes", sep = ""))
proj4string(tmp) <- 
    "+proj=utm +zone=12 +ellps=clrk66 +datum=NAD27 +units=m +no_defs"
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords,tmp@data$DEPTH, NA, NA, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"9 out of 562 is 0"
del_unzip()

#Projected Coordinate System:    NAD_1927_UTM_Zone_12N
#Projection: Transverse_Mercator
#False_Easting:  500000.00000000
#False_Northing: 0.00000000
#Central_Meridian:   -111.00000000
#Scale_Factor:   0.99960000
#Latitude_Of_Origin: 0.00000000
#Linear Unit:    Meter

#-------------------------------Indiana
#have position accuracy by description
#too many 0 value, 0 are not kept
s_n <- 3
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(paste(getwd(),"\\tmp\\Waterwells_IDNR_IN", 
        sep = ""))
proj4string(tmp) <- 
    "+proj=utm +zone=16 +ellps=GRS80 +units=m +no_defs"
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords,tmp@data$DBLBEDROCK, tmp@data$DBLSTATIC, 
    tmp@data$DBLDEPTH, tmp@data$LOC_TYPE)
tmp <- tmp[tmp[, 3] > 0, ]
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"0 out of 67308 is 0"
del_unzip()
#Grid_Coordinate_System_Name: Universal Transverse Mercator
#Universal_Transverse_Mercator:
#D_North_American_1983
#UTM_Zone_Number: 16
#Transverse_Mercator:
#Scale_Factor_at_Central_Meridian: 0.999600
#Longitude_of_Central_Meridian: -87.000000
#Latitude_of_Projection_Origin: 0.000000
#False_Easting: 500000.000000
#False_Northing: 0.000000

#-------------------------------Iowa
#have "WNUMBER", "Depth(Total)","Depth(Bedrock)","Static Water (ft)"
#some kmz are missing, use the processed data
#too many 0 value, 0 are not kept
# get processed data
s_n <- 4
setwd(paste(dir1, dirs[s_n], sep = ""))
tmp <- read.table("hor1.txt", header = TRUE)
tmp <- tmp[, c(5, 6, 3, 4, 2)]
tmp <- cbind(tmp, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"0 out of 43957 is 0"

#setwd(paste(dir1, dirs[s_n], "\\kmz", sep = ""))
#do_unzip()
#kml_files <- shell("dir .\\tmp\\* /B", intern = T)
#tmp <- data.frame(Long = character(),
#                    Lat = character(), 
#                    Depth = character())
#for(i in 1:length(kml_files))
#{
#
#    anames <- c("WNUMBER", "Depth(Total)", "Depth(Bedrock)", 
#                "Static Water (ft)")  
#    tmp <- paste(getwd(), "/tmp/", kml_files[i], sep = "")
#    if(kml_files[i] == "Iowa.kml")
#    {
#        tmp <- readOGR(dsn = tmp, layer = "Iowa.hp")
#    }else if(kml_files[i] == "Kossuth.kml")
#    {
#        tmp <- readOGR(dsn = tmp, layer = "Kossuth2")
#    }else
#    {
#        tmp <- readOGR(dsn = tmp, layer = str_sub(kml_files[i], end = -5))
#    }
#    n   <- dim(tmp@coords)[1]  
#    tmp2 <- data.frame(Long = character(),
#                     Lat = character(), 
#                     Depth = character())
#    for(j in 1:n)
#    {
#        tmp2 <- rbind(tmp2, getkml(tmp@data$Description[j], tmp@coords[j, ],
#            anames)[c(5, 6, 3)])
#    }
#    colnames(tmp2) <- colnames(rec_IA)
#    rec_IA <- rbind(rec_IA, tmp2)
#}
#i <- sapply(rec_IA, is.character)
#rec_IA[i] <- lapply(rec_IA[i], as.numeric)
#rec_IA <- cbind(rec_IA,3)
#colnames(rec_IA) <- c_names
#print_0(rec.lst[[s_n]][, 4])
##"22575 out of 55485 is 0"
#rec_IA <- rec_IA[rec_IA[, 3] > 0, ]
#row.names(rec_IA) <- NULL
#del_unzip()
#rm(tmp2)

#-------------------------------Kansas
#0 are kept
s_n <- 5
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(paste(getwd(),"\\tmp\\SouthCentral_KS_Bedrock_Wells_dd", 
        sep = ""))
tmp <- cbind(tmp@coords,tmp@data$BED_DEPTH, NA, 
    tmp@data$WELL_DEPTH, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"62 out of 14491 is 0"
del_unzip()

#-------------------------------Kentuky
#too many 0 value, 0 are not kept
#save xlsx as csv first in Excel
s_n <- 6
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- read.csv(paste(getwd(),"\\tmp\\waterwells\\waterwell.csv", 
        sep = ""), header= TRUE)
tmp <- tmp[, c("lon27", "lat27", "Depth_to_Bedrock", 
    "Static_Water_Level", "Total_Depth")]
tmp <- tmp[!is.na(tmp[, 3]), ]
tmp <- cbind(tmp, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"21144 out of 66175 is 0"
rec.lst[[s_n]] <- rec.lst[[s_n]][rec.lst[[s_n]][, 4] > 0, ]
del_unzip()

#-------------------------------Maine
#0 are kept
s_n <- 7
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
dir2 <- "\\tmp\\mapping_applications\\well_database\\_python\\"
tmp <- read.csv(paste(getwd(), dir2, "mgs_wells.csv", 
        sep = ""), header= TRUE)
tmp  <- tmp[, c("LONGITUDE", "LATITUDE", "OVERBURDEN_THICKNESS_FT", 
    "WELL_STATIC_LEVEL_FT", "WELL_DEPTH_FT", "LOCATION_ACCURACY")]
tmp <- tmp[!is.na(tmp[, 3]), ]
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"374 out of 41629 is 0"
del_unzip()

#-------------------------------Minnesota
#0 are  kept
s_n <- 8
setwd(paste(dir1, dirs[s_n], "\\dtb", sep = ""))
do_unzip()
dir2 <- "\\tmp\\OFR06_02\\"
tmp <- readShapePoints(paste(getwd(), dir2, "CWIpt_loc", 
        sep = ""))
proj4string(tmp) <- 
    "+proj=utm +zone=15 +ellps=GRS80 +units=m +no_defs"
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords, tmp@data$DEPTH2BDRK, NA, tmp@data$DEPTH_DRLL, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"613 out of 54112 is 0"
del_unzip()
#Projected Coordinate System:    NAD_1983_UTM_Zone_15N
#Projection: Transverse_Mercator
#False_Easting:  500000.00000000
#False_Northing: 0.00000000
#Central_Meridian:   -93.00000000
#Scale_Factor:   0.99960000
#Latitude_Of_Origin: 0.00000000
#Linear Unit:    Meter

#-------------------------------Missouri
#0 are  kept
s_n <- 9
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(paste(getwd(),  "\\tmp\\MO_2006_Well_Logs_shp", 
        sep = ""))
proj4string(tmp) <- 
    "+proj=utm +zone=15 +ellps=GRS80 +units=m +no_defs"
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords, tmp@data$DEPTHTOBED, NA, tmp@data$DRILLDEPTH, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"0 out of 22679 is 0"
rec.lst[[s_n]] <- rec.lst[[s_n]][rec.lst[[s_n]][, 4] > -0.01, ]
del_unzip()
#Projected Coordinate System:    NAD_1983_UTM_Zone_15N
#Projection: Transverse_Mercator
#False_Easting:  500000.00000000
#False_Northing: 0.00000000
#Central_Meridian:   -93.00000000
#Scale_Factor:   0.99960000
#Latitude_Of_Origin: 0.00000000
#Linear Unit:    Meter


#-------------------------------Nevada
#problems in reading mdb files ?????????????????????????????
#0 are  kept
#export out.dbf from mdb.....using the following SQL
#SELECT dbo_wlog.longitude, dbo_wlog.latitude, dbo_wlog.depth_bedrock
#FROM dbo_wlog
#where dbo_wlog.depth_bedrock>-0.001;
s_n <- 10
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
#mdb_file <- paste(getwd(), "/tmp/wlog.mdb", sep = "")
#tmp <- odbcConnect("E:\\wlog.mdb")
#c_NV <- odbcConnect(dsn="NV")
#sqlTables(cNCSS)$TABLE_NAME
tmp  <- foreign :: read.dbf(paste(getwd(),"\\tmp\\out.dbf", sep = ""))
tmp[, 1]  <- -tmp[, 1] 
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"2 out of 680 is 0"
del_unzip()
#-------------------------------NewHampshire
#0 are  kept
s_n <- 11
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- read.csv(".\\tmp\\all.txt", sep = "\t", colClasses = "character")
tmp <- tmp[tmp$DEPTH_TO_BEDROCK != "", c("WRB_NUMBER", "ADDRESS_1",  
    "TOWN", "DEPTH_TO_BEDROCK", "STATIC_WATER_LEVEL", "TOTAL_DEPTH")]
tmp[, c("DEPTH_TO_BEDROCK", "STATIC_WATER_LEVEL", "TOTAL_DEPTH")] <-  
    sapply(tmp[, c("DEPTH_TO_BEDROCK", "STATIC_WATER_LEVEL", 
    "TOTAL_DEPTH")], pattern = " ft", replacement = "",str_replace)
tmp2 <- read.csv(".\\tmp\\ll.csv") 
tmp2 <- tmp2[tmp2$lon != "lon",]
tmp2 <- tmp2[!is.na(tmp2$lon), ]
tmp <- join(tmp2, tmp, type = "left", by = "WRB_NUMBER")
tmp  <- cbind(tmp[, c("lon", "lat", "DEPTH_TO_BEDROCK", "STATIC_WATER_LEVEL", 
            "TOTAL_DEPTH")], NA)              
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
# "0 out of 19782 is 0"
del_unzip()

#tmp4 <- tmp[!(tmp$WRB_NUMBER %in% tmp2$WRB_NUMBER), ]
##### block 1
#tmp2 <- NULL
#tmp3 <- NULL
#library(XML)
#library(stringr)
#for(i in 1:500)
#{
#    print(i)
#    tmp2 <- rbind(tmp2, get_coords2(tmp4$WRB_NUMBER[i]))
#    if(i %% 1000 == 0)
#    {
#        tmp3 <- rbind(tmp3,tmp2)
#        tmp2 <- NULL 
#        save.image(paste(dir_out, "ll1.RData", sep = ""))   
#    }
#}
#tmp3<- rbind(tmp3,tmp2)
#write.csv(tmp3, "1.csv")

---------------------Newyork
#0 are  kept
s_n <- 12
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(paste(getwd(), "\\tmp\\WaterWellProgram", 
        sep = ""))
proj4string(tmp) <- 
    "+proj=utm +zone=18 +ellps=GRS80 +units=m +no_defs"
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords, as.character(tmp@data$Rock_Depth), 
    as.character(tmp@data$GW_Depth), as.character(tmp@data$Well_Depth), "")
tmp <- tmp[!is.na(tmp[, 3]), ]
tmp <- tmp[tmp[, 2]>40,]  
tmp[tmp[, 3]>70,]  
rec.lst[[s_n]] <- form_rec(tmp, s_n)
rec.lst[[s_n]] <- rec.lst[[s_n]][rec.lst[[s_n]][,4] > -0.01, ]
print_0(rec.lst[[s_n]][, 4])
#"165 out of 14652 is 0"
del_unzip()


#PROJCS["NAD_1983_UTM_Zone_18N",GEOGCS["GCS_North_American_1983",
#DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],
#PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],
#PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],
#PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-75.0],
#PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],
#UNIT["Meter",1.0]]

#-------------------------------Ohio
#0 are  kept
#determine depth ==0 if depth_CALC is near 0
#Attribute_Label: DEPTH_SFX -----useful but not used here
#Attribute_Definition: Indicates a well that does not reach bedrock.
s_n <- 13
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(paste(getwd(), "\\tmp\\bedrock\\BT24K\\btpoints", 
        sep = ""))
proj4string(tmp) <- 
    "+proj=lcc +lat_1=38.73333333333333 +lat_2=40.03333333333333 +lat_0=38 
    +lon_0=-82.5 +x_0=600000.0000000001 +y_0=0 +ellps=GRS80 +datum=NAD83 
    +to_meter=0.3048006096012192 +no_defs "
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords, tmp@data$DEPTH, tmp@data$DEPTH_CALC)
tmp <- tmp[tmp[, 3] > -0.01, ] 
tmp <- tmp[tmp[, 3] > 0 | abs(tmp[, 3]-tmp[, 4]) < 3, ] 
tmp <- cbind(tmp[, c(1,2,3)], NA, NA, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"309 out of 56578 is 0"
del_unzip()
#PROJCS["NAD_1983_StatePlane_Ohio_South_FIPS_3402_Feet",
#GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",
#SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],
#UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],
#PARAMETER["False_Easting",1968500.0],PARAMETER["False_Northing",0.0],
#PARAMETER["Central_Meridian",-82.5],
#PARAMETER["Standard_Parallel_1",38.73333333333333],
#PARAMETER["Standard_Parallel_2",40.03333333333333],
#PARAMETER["Latitude_Of_Origin",38.0],UNIT["Foot_US",0.3048006096012192]]

#-------------------------------other
#0 are  kept
#problems in reading mdb files ?????????????????????????????
#export PTS_DEPT.DBF from mdb
s_n <- 14
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- foreign :: read.dbf(paste(getwd(),"\\tmp\\PTS_DEPT.DBF", sep = ""))
tmp <- cbind(tmp[, c("X_NAD83", "Y_NAD83", "DEPTH")], as.numeric(NA), 
        as.numeric(NA), "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
"38 out of 16950 is 0"
del_unzip()

#-------------------------------Pennsylvania
#0 are  kept
#BedrockNotReached-----useful but not used here
#problems in reading mdb files ?????????????????????????????
#export TBLGENSI.DBF from mdb
s_n <- 15
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- foreign :: read.dbf(paste(getwd(),"\\tmp\\TBLGENSI.DBF", sep = ""))
tmp <- cbind(tmp[, c("LONGITUDED", "LATITUDEDD", "DEPTHTOBED")], as.numeric(NA), 
    tmp[, c("WELLDEPTH", "LATLONGACC")])
tmp <- tmp[!(is.na(tmp[, 1]) | is.na(tmp[, 2]) | is.na(tmp[, 3]) | 
        tmp[, 1] == 0 | tmp[, 2] == 0), ]
# fix some long and lat in order
tmp$LATITUDEDD[tmp$LATITUDEDD > 100] <- tmp$LATITUDEDD[tmp$LATITUDEDD > 100]/100
tmp <- tmp[tmp$LATITUDEDD < 80, ]
tmp2 <- which(tmp$LATITUDEDD < -70)
tmp[tmp2, c("LONGITUDED", "LATITUDEDD")] <- 
    swap2n(tmp[tmp2, c("LONGITUDED", "LATITUDEDD")])
tmp2 <- which(tmp$LATITUDEDD < -30)
tmp$LATITUDEDD[tmp2] <- - tmp$LATITUDEDD[tmp2]
tmp2 <- which(tmp$LONGITUDED < -10000)
tmp$LONGITUDED[tmp2] <- tmp$LONGITUDED[tmp2]/100000
tmp2 <- which(tmp$LONGITUDED < -100)
tmp$LONGITUDED[tmp2] <- tmp$LONGITUDED[tmp2]/100
tmp <- tmp[tmp$LONGITUDED > -81, ]
tmp2 <- which(tmp$LONGITUDED > 70)
tmp$LONGITUDED[tmp2] <- - tmp$LONGITUDED[tmp2]
tmp <- tmp[tmp$LONGITUDED < -60, ]
#range(tmp$LONGITUDED) 
# fix end
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"1146 out of 147115 is 0"
del_unzip()

#-------------------------------Tennessee
#0 are  kept
s_n <- 16
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(paste(getwd(), "\\tmp\\regolith", 
        sep = ""))
proj4string(tmp) <- 
    "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords, tmp@data$DEPTH_TO_B, as.numeric(NA), 
        as.numeric(NA), "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"0 out of 29914 is 0"
del_unzip()
#GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",
#SPHEROID["GRS_1980",6378137.0,298.257222101]],
#PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]

#-------------------------------Vermont
#0 are not kept
s_n <- 17
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(paste(getwd(), "\\tmp\\Water_PVTWELLS_point", 
        sep = ""))
proj4string(tmp) <- 
    "+proj=tmerc +lat_0=42.5 +lon_0=-72.5 +k=0.9999642857142857 
    +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs  "
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords, tmp@data$OverBurden, tmp@data$StaticWate, 
    tmp@data$WellDepth, tmp@data$LocationDe)
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"8212 out of 100208 is 0"
rec.lst[[s_n]] <- rec.lst[[s_n]][rec.lst[[s_n]][, 4] > 0, ]
del_unzip()
#PROJCS["NAD_1983_StatePlane_Vermont_FIPS_4400",
#GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",
#SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],
#UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],
#PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],
#PARAMETER["Central_Meridian",-72.5],
#PARAMETER["Scale_Factor",0.9999642857142858],
#PARAMETER["Latitude_Of_Origin",42.5],UNIT["Meter",1.0]]


#library(foreign)
#library(maptools)
#library(rgdal)
#library(stringr)
##library(RODBC)
#library(plotKML)
#library(RSAGA)
#dir_out <- "E:\\data\\soildata\\depth\\points\\"
#load(paste(dir_out, "us1.RData", sep = ""))

#-------------------------------merge rec.lst
wells <- rec.lst[[1]]
for(i in 2:length(rec.lst))
{
    wells <- rbind(wells,rec.lst[[i]])
      
}
i <- 4:6
wells[i] <- lapply(wells[i], function(x) x <- x*0.3048)

setwd(dir1)
#-------------------------------check values
for(i in 1:length(rec.lst))
{
    png(filename = paste0("..\\pic\\", names(rec.lst)[i], ".png"))      
    hist(rec.lst[[i]][, 4], breaks = 200,
        xlab = names(rec.lst)[i])
    dev.off()
    png(filename = paste0("..\\pic\\", 2,names(rec.lst)[i], ".png"))      
    hist(rec.lst[[i]][, 4], 
        breaks = seq(from = -1, to = max(rec.lst[[i]][, 4])+10, by = 1 ), 
        xlim = c(0,200), xlab = names(rec.lst)[i])
    dev.off()
    print(names(rec.lst)[i])
    print_0(rec.lst[[i]][, 4]) 
}
#tmp <- unique(rec.lst[[1]][, 2:3])
#wellsp <- rec.lst[[15]][rec.lst[[15]][,4]<5, 2:4]
#coordinates(wellsp) <- ~ Long+Lat
#proj4string(wellsp) <- CRS("+proj=longlat +datum=WGS84")
#plotKML(wellsp)
# too much 0 value: Arizona is ok
# too much vuale near 0: kentuky is OK, maine is OK
# too much 1, nevada may be OK? 
#for(i in 1:length(rec.lst))
#{
#    print(states[i,])
#    for( j in 2:6)
#    {
#        print(colnames(rec.lst[[i]][j]))
#        print(max(rec.lst[[i]][j], na.rm = TRUE))
#        print(min(rec.lst[[i]][j], na.rm = TRUE))
#    }
#}

for( i in 2:6)
{
    print(colnames(wells[i]))
    print(max(wells[i], na.rm = TRUE))
    print(min(wells[i], na.rm = TRUE))
}
#wells[wells[,3]<20,]
wells <- wells[wells[, 4] < 4000, ]
print_0(wells[,4])
#"2718 out of 661970 is 0"

#-------------------------------get the list of accuracy of position
#levels(as.factor(wells[, 7]))
#levels(as.factor(wells[wells[, 7] != "", 1]))
#Alaska
s_n = 1
tmp  <- c("50", "100", "200", "400", "800" )
acc <- cbind(s_n, tmp, tmp)
#Indiana
s_n = 3
tmp  <- c("1", "2", "3")
tmp2 <- c("Actual Location", "Estimated Location/Geocoding addresses",
 "Estimated Location/TRS_quarter sections_county")
acc <- rbind(acc, cbind(s_n, tmp, tmp2))
#Maine
s_n = 7
tmp  <- c("WITHIN 100M", "WITHIN 10M", "WITHIN 15M", "WITHIN 30M", 
        "WITHIN 5M")
acc <- rbind(acc, cbind(s_n, tmp, tmp))
#Nevada
s_n = 10
tmp  <- c("1", "5", "F", "G", "H", "M", "T" )
acc <- rbind(acc, cbind(s_n, tmp, tmp))
#s_n = 13
tmp  <- c("1", "F", "M", "S", "T" )
tmp2 <- c("UNKNOWN", "ACCURATE TO +5 SECONDS",
    "ACCURATE TO +1 MINUTE", "ACCURATE TO +1 SECOND",
    "ACCURATE TO +10 SECONDS")
acc <- rbind(acc, cbind(s_n, tmp, tmp2))
#Pennsylvania
s_n =15
tmp  <- as.character(1:7)
tmp2 <- c("1:24,000 USGS Digitized",  
    "calcualted coordinate measured off 1:24,000 USGS Quad",
    "E911 Address", "field located with GPS unit - uncorrected",            
    "GPS location", "screen digitized", "Welldriller/Clarion")
acc <- rbind(acc, cbind(s_n, tmp, tmp2))
colnames(acc) <- c("source", "flag", "descrip")

#-------------------------------output
#rm(rec.lst)
write.table(cbind(states,count(wells$Source)[2]), paste(dir_out, "states_us.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
write.table(acc, paste(dir_out, "acc_us.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
write.table(wells, paste(dir_out, "wells_us.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
write.table(count(wells$Source), paste(dir_out, "wells_us.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
save.image(paste(dir_out, "us.RData", sep = ""))
