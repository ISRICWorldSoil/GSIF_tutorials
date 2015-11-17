rm(list = ls(all = TRUE))
library(rgdal)
library(maptools)
library(plotKML)
library(stringr)

#dir_out <- "E:\\data\\soildata\\depth\\points\\profs\\well\\"
#load(paste(dir_out, "ca.RData", sep = ""))

#define global, change if needed
options(stringsAsFactors = FALSE) # not usful for readShapePoints, 
                                # fix s_n == 1
# change dir* if needed
dir_f <- "E:\\data\\soildata\\depth\\points\\codegsifb\\head\\"
source(paste(dir_f, "functions.r", sep = ""))
dir1 <- "E:\\data\\soildata\\depth\\points\\canada\\"
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy") 
setwd(dir1)
dirs   <- shell("dir * /B", intern = T)
states <- cbind(1:length(dirs), dirs)
colnames(states) <- c("s_code", "states")
states
dir_out <- "E:\\data\\soildata\\depth\\points\\"
#     s_code states           
#[1,] "1"    "Britishcolumbia"
#[2,] "2"    "GINwaterwell"   
#[3,] "3"    "Nova Scotia"    
#[4,] "4"    "ontario"        
#[5,] "5"    "quebec"  

rec.lst <- as.list(rep(NA, length(dirs)))
names(rec.lst) <- dirs
#-------------------------------Britishcolumbia, feet
#0 are kept
s_n <- 1 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp  <- readShapePoints(paste(getwd(), 
    "\\tmp\\GW_WW_WRBC\\GW_WW_WRBC_point", sep = ""))
proj4string(tmp) <- 
    "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 
    +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords[, c(1, 2)], tmp@data$BDRCKDPTH, tmp@data$WTRDPTH, 
    tmp@data$DPTHWLLDRL, as.character(tmp@data$TMCCRCCD))
tmp <- tmp[!is.na(tmp[,3]), ]
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
# "1498 out of 32627 is 0"
del_unzip()

#-------------------------------GINwaterwell
# can be derived but bot reliable to use

#-------------------------------Nova Scotia, feet
#export well.txt from mde file
#0 are kept
s_n <- 3 
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- read.table(paste(getwd(),"\\tmp\\well.txt", sep = ""), 
        sep = ",", head = T)
tmp <- tmp[, c("Easting", "Northing", "DepthToBedrock",  
    "wyStaticLevel", "TotalOrFinishedDepth", "UTMAccuracy")] 
tmp <- tmp[!(is.na(tmp[, 1]) | is.na(tmp[,2])| is.na(tmp[,3])), ]   
coordinates(tmp) <- ~ Easting+Northing
proj4string(tmp) <- 
    "+proj=utm +zone=20 +ellps=WGS84 +units=m +no_defs"
tmp <- reproject(tmp)
tmp <- cbind(tmp@coords, tmp@data[, c("DepthToBedrock",  
    "wyStaticLevel", "TotalOrFinishedDepth", "UTMAccuracy")])

rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
# "359 out of 91071 is 0"
del_unzip()

#-------------------------------ontario, meter
#0 are not kept
#ODH data not used, as it need some process
s_n <- 4 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp  <- readShapePoints(paste(getwd(), "\\tmp\\wwis_out", sep = ""))
tmp <- cbind(tmp@coords, tmp@data$DP_BEDROCK, tmp@data$STATIC_LEV, 
    tmp@data$DEPTH, "")
tmp <- tmp[!is.na(tmp[, 3]), ]
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
# "21501 out of 366207 is 0"
rec.lst[[s_n]] <- rec.lst[[s_n]][rec.lst[[s_n]][, 4] > 0, ]
del_unzip()

#-------------------------------quebec, meter
#0 are kept
s_n <- 5 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
setwd(paste(getwd(), "\\tmp\\qc", sep = ""))
kml_files <- shell("dir *.kml /B", intern = T)
#Well depth, Depth of the rock
anames <- c("PROF_PUITS", "PROF_ROC")
tmp3 <- data.frame(Depth = character(),
                    DTB = character(),
                    Long = character(),
                    Lat = character())
for(j in 1:(length(kml_files)))
#for(j in 1:2)
{
    tmp <- readOGR(dsn = kml_files[j], layer = str_sub(kml_files[j], 
            end = -5))
    p_n <- dim(tmp@coords)[1]
    tmp2 <- NULL   
    for(i in 1:p_n)
    {
        tmp2 <- rbind(tmp2, getkml2(tmp@data$Description[i], 
            tmp@coords[i, ], anames))
    }
    colnames(tmp2) <- colnames(tmp3)
    tmp3 <- rbind(tmp3, tmp2)
    tmp3 <- cbind(tmp3[tmp3[, 2] != "", ])    
}
tmp <- cbind(tmp3[, c(3, 4, 2)], as.numeric(NA), tmp3[, 1], "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"2736 out of 111659 is 0"
setwd(paste(dir1, dirs[s_n], sep = ""))
del_unzip()
rm(tmp)
rm(tmp2)
rm(tmp3)

#-------------------------------merge rec.lst
setwd(dir1)
wells <- NULL
for(i in 1:length(rec.lst))
{
    if(i==1 || i==3)
    {
        tmp <- rec.lst[[i]]
        j <- c(4,5,6)
        tmp[j] <- lapply(rec.lst[[i]][j], function(x) x <- x*0.3048)
        wells <- rbind(wells,tmp)
    }else
    {
        wells <- rbind(wells,rec.lst[[i]])
    }
}
wells <- wells[!is.na(wells[, 4]), ]
#-------------------------------check values
#for(i in c(1,3:length(rec.lst)))
#{
#    png(filename = paste0("..\\pic\\", names(rec.lst)[i], ".png"))      
#    if(i == 4 | i == 5)
#    {
#        hist(rec.lst[[i]][, 4]/0.3048, breaks = 200,
#            xlab = names(rec.lst)[i])
#        dev.off()
#        png(filename = paste0("..\\pic\\", 2,names(rec.lst)[i], ".png"))      
#        hist(rec.lst[[i]][, 4]/0.3048, 
#            breaks = seq(from = -1, to = max(rec.lst[[i]][, 4])/0.3048+10, by = 1 ), 
#            xlim = c(0,200), xlab = names(rec.lst)[i])
#        dev.off()
#    }else
#    {
#        hist(rec.lst[[i]][, 4], breaks = 200,
#            xlab = names(rec.lst)[i])
#        dev.off()
#        png(filename = paste0("..\\pic\\", 2,names(rec.lst)[i], ".png"))      
#        hist(rec.lst[[i]][, 4], 
#            breaks = seq(from = -1, to = max(rec.lst[[i]][, 4])+10, by = 1 ), 
#            xlim = c(0,200), xlab = names(rec.lst)[i])
#        dev.off()
#    }
#    print(names(rec.lst)[i])
#    print_0(rec.lst[[i]][, 4]) 
#}
#tmp <- count(rec.lst[[5]][, 4])
#wellsp <- rec.lst[[4]][rec.lst[[4]][,4]<2, 2:4]
#wellsp <- wellsp[runif(dim(wellsp)[1])<0.1, ]
#coordinates(wellsp) <- ~ Long+Lat
#proj4string(wellsp) <- CRS("+proj=longlat +datum=WGS84")
#plotKML(wellsp)
## too much 0 value: ontario may be deleted?!!!! 
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
print_0(wells[, 4])
#"4593 out of 580063 is 0"
#-------------------------------get the list of accuracy of position
#levels(as.factor(wells[, 7]))
s_n = 1
tmp  <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
acc <- cbind(s_n, tmp, tmp)
s_n  <- 3
tmp  <- c("1", "2", "3", "4")
tmp2 <- c("50", "800", "1500", "8000")
acc <- rbind(acc, cbind(s_n, tmp, tmp2))
colnames(acc) <- c("source", "flag", "descrip")

#-------------------------------output
#rm(rec.lst)
tmp <- count(wells$Source)
names(tmp) <- c("s_code","freq")
write.table(join(as.data.frame(states),tmp), paste(dir_out, "states_ca.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
write.table(acc, paste(dir_out, "acc_ca.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
write.table(wells, paste(dir_out, "wells_ca.txt", sep = ""), 
    row.names = FALSE, sep = "\t")

save.image(paste(dir_out, "ca.RData", sep = ""))

