rm(list = ls(all = TRUE))


#define global, change if needed
dir_f <- "E:\\data\\soildata\\depth\\points\\codegsifb\\head\\"
source(paste(dir_f, "functions.r", sep = ""))
dir1  <- "E:\\data\\soildata\\depth\\points\\russia\\"
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy") 
setwd(dir1)
dirs   <- shell("dir * /B", intern = T)
states <- cbind(1:length(dirs), dirs)
colnames(states) <- c("s_code", "states")
states
dir_out <- "E:\\data\\soildata\\depth\\points\\"#      s_code states        



###
rec.lst <- as.list(rep(NA, length(dirs)))
names(rec.lst) <- dirs



#-------------------------------russia & alaska soil
#0 are  kept,
s_n <- 1 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- read.csv(".\\tmp\\desc.dat", sep ="\t", header = F)

