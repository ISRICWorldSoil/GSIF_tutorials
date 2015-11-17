rm(list = ls(all = TRUE))

library(foreign)
library(plyr)
library(stringr)

#global
options(stringsAsFactors = FALSE)
dir <- "ngis/"
dir_out <- "depth/points/"
## unzip large files:
unzip("lookup/litlist1.zip", exdir = "./lookup")
unzip("lookup/litlist5.zip", exdir = "./lookup")

# files of unique records
list_files <- c("hgulist.txt", "conlist.txt", "keywords.txt",
    "litlist1.txt", "litlist2.txt", "litlist3.txt", "litlist4.txt", "litlist5.txt")
#files to seprate lithology logs
lit_seps <- c("lit_del.txt", "lit_keep.txt","lit_qsure.txt",
     "lit_nsure.txt", "lit_noinf.txt", "lit_nouse.txt", "lit_qdel.txt")
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy")


p_dim <- function(lst)
{
    t <- list(del = dim(lst[["del"]])[1], 
         keep  = dim(lst[["keep"]])[1],
         qsure = dim(lst[["qsure"]])[1], 
         nsure = dim(lst[["nsure"]])[1],
         noinf = dim(lst[["noinf"]])[1], 
         nouse = dim(lst[["nouse"]])[1],
         qdel  = dim(lst[["qdel"]])[1])       
    t <- unlist(t)
    print(t)
    print(sum(t))
} 
key_del<-function(del,rec) # get record index to be deleted in del list
{
    tmp <- grep(pattern = del[1], rec, ignore.case = TRUE)
    for (i in 2:length(del))
    {
        tmp <- c(tmp, grep(pattern = del[i], rec, 
            ignore.case = TRUE))
    }
    unique(tmp)
}
key_del2 <- function(p1, p2,  lit)
{
    tmp  <- grepl(pattern = p1, lit$Descriptio, ignore.case = TRUE)   
    tmp2 <- lit[tmp,]
    for(i in length(p2))
    {
        tmp2[, "Descriptio"] <- gsub (pattern = p2[i], replacement = "", 
            tmp2[, "Descriptio"], ignore.case = TRUE)
    }     
    tmp  <- grepl(pattern = p1, tmp2$Descriptio, ignore.case = TRUE)
    tmp3 <- rownames(lit) %in% rownames(tmp2[tmp, ])
    print(sum(tmp))
    del  <- tmp2[tmp, ]
    lit  <- lit[!tmp3,]
    list(lit = lit, del = del)
}
key_keep <- function(rec_nsure, rec_keep, rec_del)
{
    tmp <- grep(pattern = "consolidate", rec_nsure$Descriptio, ignore.case = TRUE)
    if(exists("rec_keep"))
    {
        rec_keep <- rbind(rec_keep, rec_nsure[tmp, ])
    }else
    {
        rec_keep  <- rec_nsure[tmp, ]
    }
    rec_nsure <- rec_nsure[-tmp, ]
    
    tmp <- grep(pattern = "consolidate", rec_del$Descriptio, ignore.case = TRUE)
    rec_keep <- rbind(rec_keep, rec_del[tmp, ])
    rec_del <- rec_del[-tmp,]
    
    tmp <- grep(pattern = "unconsolidate", rec_keep$Descriptio, ignore.case = TRUE)
    rec_del <- rbind(rec_del, rec_keep[tmp, ])
    rec_keep <- rec_keep[-tmp, ]
    list(rec_nsure, rec_keep, rec_del)
} 
key_nouse <- key_del 
key_noinf <- key_del 
key_qdel  <- key_del
key_qsure <- key_del 
str_sep <- function(s)
{
    s <- str_replace_all(s, "\"", "")
    unlist(str_split(s, "\t"))
}

#get the depth to bedrock
get_dtb <- function(bores, keep, index)
{
    tmp <- keep[index,][ , c("BoreID", "FromDepth")] 
    tmp <- join(tmp, bores[ ,c("BoreID", "Latitude", "Longitude", 
            "BoreDepth")], by = c("BoreID"), type = "left" )
    tmp <- data.frame( tmp["BoreID"], tmp["Longitude"], 
            tmp["Latitude"], tmp["FromDepth"], NA,
            tmp["BoreDepth"],  NA)        
    colnames(tmp) <- c_names
    tmp
}
#----------------------------------------unzip files
setwd(dir)
#zip_files <- system("ls *.zip", intern = T) # linux
zip_files  <- shell("dir *.zip /B", intern = T)
fold_names <- sub(".zip", "", zip_files)
dbfs <- c("NGIS_Bores.dbf", "NGIS_BoreholeLog.dbf", 
        "NGIS_ConstructionLog.dbf","NGIS_LithologyLog.dbf")
#for(i in 1:length(zip_files))
#{
#    unzip(zip_files[i],exdir = fold_names[i])
#}

#-------------------------------------------
#read the dbfs
col.lst <- list(
    bores = c("HydroID", "Latitude", "Longitude", "RefElev", "BoreDepth"),
    hgus  = c("BoreID", "FromDepth", "ToDepth", "HGUID", "Descriptio"),
    cons  = c("BoreID", "FromDepth", "ToDepth", "Material"),
    lits  = c("BoreID", "FromDepth", "ToDepth", "MajorLithC", "MinorLithC", 
    "Descriptio"))
dbf.lst <- list(bores = NULL, hgus= NULL, cons = NULL, lits = NULL)
for(i in 1:length(fold_names))
{
    dir_2 <- paste(dir, "\\", fold_names[i],"\\", fold_names[i], "\\", sep = "")
    for(j in 1:length(dbf.lst))
    {
        dbf.lst[[j]] <- rbind(dbf.lst[[j]], 
            read.dbf(paste(dir_2, dbfs[j], sep = ""))[col.lst[[j]]])
    } 
} 
for( i in 2:4)
{
    j <- sapply(dbf.lst[[i]], is.factor)
    dbf.lst[[i]][j] <- lapply(dbf.lst[[i]][j], as.character)
}

#----------------------------------------remove unzipped
#for (i in 1:length(fold_names))
#{ 
#    cmd <- paste("rmdir /S /Q ", fold_names[i], sep = "")
#    shell(cmd)
#}

#-----------------------------------------------  
#deal with cons
#tmp <- summary(cons[4], maxsum = 200) 
#tmp <- count(cons[4]) 
#write.table(tmp, paste(dir, "\\", "conlist2.txt", sep = ""), sep = "\t", 
#    row.names = FALSE)
# store as conlist.txt and indicate which one to keep mannually:
#0: delete, 1: keep, 2: not sure, 3: not clear what it is
conls <- read.table(paste(dir, "\\", list_files[2], sep = ""),header = TRUE)
conls <- conls[conls$mark == 1, ]
cons_out <- match_df(dbf.lst[[3]], conls, on = "Material")

#-----------------------------------------------
#deal with lits
#list to seprate the lits
ls.lst <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
lit.lst <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
#check  the FromDepth from different area:
#something may be wrong for BoreID > 70000000
#for (i in 1:7)
#{
#    i <- i*10000000
#    tmp3 <- dbf.lst[[4]][(dbf.lst[[4]]$BoreID < i+10000000 & 
#    dbf.lst[[4]]$BoreID > i), "FromDepth"]
#    hist(tmp3, xlim = c(0, 300), breaks = seq(-30, 3700,by = 30), 
#       freq = FALSE, xlab = i)
#    hist(tmp3, xlim = c(0, 50), breaks = seq(-30, 3700,by = 2), 
#       freq = FALSE, xlab = i)
#}

#----------------------------------------corrections
names(dbf.lst[[4]])
#"BoreID"     "FromDepth"  "ToDepth"    "MajorLithC" "MinorLithC" "Descriptio"
for(i in 4:6)
{
    dbf.lst[[4]][is.na(dbf.lst[[4]][, i]), i] <- ""
}
dbf.lst[[4]][dbf.lst[[4]]$MinorLithC == "NA", "MinorLithC"]  <- ""
dbf.lst[[4]][dbf.lst[[4]]$MinorLithC == "Na", "Descriptio"]  <- ""
tmp2 <- c("\\", "\"")
for(i in 1:length(tmp2))
{
    tmp <- str_sub(dbf.lst[[4]]$Descriptio, -1) == tmp2[i] 
    dbf.lst[[4]][tmp, "Descriptio"] <- 
        str_sub(dbf.lst[[4]][tmp, "Descriptio"], 1, -2)
}
#unique combination of MajorLithC, MinorLithC and Descriptio
#tmp <- count(dbf.lst[[4]][4:6]) 
#write.table(tmp, paste(dir, "\\", "litlist12.txt", sep = ""), sep = "$", 
#    row.names = FALSE)
#store as litlist1.txt

#-----------------------------------------------
#classify according to unique MajorLithC
#tmp <- count(dbf.lst[[4]][4])
#write.table(tmp, paste(dir, "\\", "litlist22.txt", sep = ""), sep = "\t",
#     row.names = FALSE)
#store as litlist2.txt and indicate which one to keep mannually, 
#juded according to description:
#0: delete, 1: keep, 2: quite sure, 3: not sure, 4: not clear what it is,
#5: if this layer exit, the whole points is not used, human disturbance
#6: qure sure delete,
#note that NULL and the numbers should be converted in text format in saving
tmp <- read.table(paste(dir, "\\", "litlist2.txt", sep = ""), header = TRUE)
tmp2 <- c(0, 1, 3, 4, 5)
for(i in tmp2)
{
    ls.lst[[i+1]] <- tmp[tmp$token == i, 1:2]
    lit.lst[[i+1]] <- match_df(dbf.lst[[4]], ls.lst[[i+1]], 
                    on = "MajorLithC")
}
str(dbf.lst[[4]])
p_dim(lit.lst)
#   del    keep   nsure   noinf   nouse 
# 747568     497 2435873      62    8116 
#rm(dbf.lst[[4]]) # save storeage

#-----------------------------------------------
#classify according to unique MinorLithC
#tmp <- count(lit.lst[["nsure"]][5])
#write.table(tmp, paste(dir, "\\", "litlist32.txt", sep = ""), sep = "\t", 
#   row.names = FALSE)
tmp <- read.table(paste(dir, "\\", "litlist3.txt", sep = ""), header = TRUE)
ls.lst <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
tmp2 <- c(0, 1, 4, 5)
for(i in tmp2)
{
    ls.lst[[i+1]] <- tmp[tmp$token == i, 1:2]
    tmp3 <- match_df(lit.lst[["nsure"]], ls.lst[[i+1]], on = "MinorLithC")
    lit.lst[[i+1]] <- rbind(lit.lst[[i+1]], tmp3)
    lit.lst[["nsure"]] <- subset(lit.lst[["nsure"]], 
        !(rownames(lit.lst[["nsure"]]) %in% rownames(tmp3)))
}
p_dim(lit.lst)
#    del    keep   nsure   noinf   nouse 
# 757345     544 2426044      65    8118 

#-----------------------------------------------   
#classify according to unique Descriptio

#check the records
#sel <- grep(pattern = "wash", dbf.lst[[4]]$Descriptio, ignore.case = TRUE)
#sel <-  dbf.lst[[4]]$Descriptio == "Slippery back"
#ii  <- order(dbf.lst[[4]][sel, "FromDepth"])
#tmp3 <- dbf.lst[[4]][sel,]
#tmp3 <- tmp3[ii,]
#View(tmp3)
#-----------------------------------------------
# classify Descriptio with more than 100 records
#tmp <- count(lit.lst[["nsure"]][6])
#write.table(tmp, paste(dir, "\\", "litlist42.txt", sep = ""), sep = "\t", 
#   row.names = FALSE)
tmp <- read.table(paste(dir, "\\", "litlist4.txt", sep = ""), header = TRUE, 
   sep = "\t")
ls.lst <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
tmp2 <- c(0, 1, 2, 4, 5)
for(i in tmp2)
{
    ls.lst[[i+1]] <- tmp[tmp$token == i, 1:2]
    tmp3 <- match_df(lit.lst[["nsure"]], ls.lst[[i+1]], on = "Descriptio")
    lit.lst[[i+1]] <- rbind(lit.lst[[i+1]], tmp3)
    lit.lst[["nsure"]] <- subset(lit.lst[["nsure"]], 
        !(rownames(lit.lst[["nsure"]]) %in% rownames(tmp3)))
}
p_dim(lit.lst)   
#    del    keep   qsure   nsure   noinf   nouse 
#1551698   34018  276506 1290023   29617   10254 
#-----------------------------------------------
#classify Descriptio with keywords  
#library(foreign)
#library(plyr)
#library(stringr)
#options(stringsAsFactors = FALSE)
#dir <- "E:\\data\\soildata\\depth\\points\\australia\\ngis"
#load(paste0(dir,"\\as2.RData"))
save.image(paste0(dir_out,"\\as1.RData"))
tmp <- read.table(paste(dir, "\\", "keywords.txt", sep = ""), header = TRUE, 
    sep = "\t") 
i <- sapply(tmp, is.factor)
tmp[i]  <- lapply(tmp[i], as.character)
ls.lst  <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
for(tmp2 in names(ls.lst))
{
    ls.lst[[tmp2]] <- tmp[!(tmp[tmp2] == ""), tmp2] 
}

# check the ls.lst[["del"]] manually to make sure it right
#for (i in 1:length(ls.lst[["del"]]))
#{
#    tmp <- grep(pattern = ls.lst[["del"]][i], lit.lst[["nsure"]]$Descriptio,  
#       ignore.case=TRUE)
#    tmp2 <- lit.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    if (str_sub(ls.lst[["del"]][i], -1) == "/")
#    {
#        tmp3 <- str_c(str_sub(ls.lst[["del"]][i], 1, -2), "1")
#    }else
#    {
#        tmp3 <- ls.lst[["del"]][i]
#    }
#    write.table(tmp2, paste(dir, "\\del\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}

# The order of operation should be: del, keep,(+ mud, clay, silt) nouse, noinfor, 
# qdel, qsure and pkeep
#--------------------------del
tmp <- key_del(ls.lst[["del"]], lit.lst[["nsure"]]$Descriptio)
lit.lst[["del"]] <- rbind(lit.lst[["del"]], lit.lst[["nsure"]][tmp, ])
lit.lst[["nsure"]] <- lit.lst[["nsure"]][-tmp, ]
p_dim(lit.lst)
#--------------------------keep
#only one to keep "consolidate"

tmp <- key_keep(lit.lst[["nsure"]], lit.lst[["keep"]], lit.lst[["del"]])
lit.lst[["nsure"]] <- tmp[[1]]
lit.lst[["keep"]]  <- tmp[[2]]
lit.lst[["del"]]  <- tmp[[3]]
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2393312   50457  276506  431970   29617   10254 

#--------------------------correction 
#Caly to clay, slit to silt,
#clay stone,clay-stone(clastone) to claystone
#silt stone,silt-stone to siltstone, 
#mud stone to mudstone

tmp <- grepl(pattern = "caly", lit.lst[["nsure"]]$Descriptio, 
        ignore.case = TRUE) & !grepl(pattern = "scaly",
        lit.lst[["nsure"]]$Descriptio, ignore.case = TRUE)
lit.lst[["nsure"]][tmp, "Descriptio"] <- gsub (pattern = "caly", 
        replacement = "clay", lit.lst[["nsure"]][tmp, "Descriptio"],
         ignore.case = TRUE)
tmp2 <- c("slit", "silthe", "CALCISIlTITE", "silt stone", "silt-stone", 
            "silstone", "siltsone", "siltone", "clay stone", "clay-stone", 
            "claysone", "mud stone", "mud-stone", "mud rock")
tmp3 <- c("silt", "slithe", "CALCISLITITE", "siltstone", "siltstone", 
            "siltstone", "siltstone", "siltstone", "claystone", "claystone", 
            "claystone", "mudstone", "mudstone", "mudrock")
for(i in length(tmp2))
{
    lit.lst[["nsure"]][, "Descriptio"] <- gsub (pattern = tmp2[i], 
        replacement = tmp3[i], lit.lst[["nsure"]][, "Descriptio"],
        ignore.case = TRUE)
}

#delete clay, but not claystone etc.    
# the following keywords are kept for "clay"
tmp <- c("clayston", "claybound", "clayish", "clayee", "clayed", "clayband")
tmp <- key_del2("clay", tmp, lit.lst[["nsure"]])
lit.lst[["nsure"]] <- tmp[["lit"]]
lit.lst[["del"]] <- rbind(lit.lst[["del"]], tmp[["del"]])
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2534619   50457  276506  290663   29617   10254 


#delete silt but not siltstone etc. 
# the following keywords are kept for "silt"
tmp <- c("siltston", "siltit", "siltrex", "Silted", "Siltband", "siltier")
tmp <- key_del2("silt", tmp, lit.lst[["nsure"]])
lit.lst[["nsure"]] <- tmp[["lit"]]
lit.lst[["del"]] <- rbind(lit.lst[["del"]], tmp[["del"]])
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2559077   50457  276506  266205   29617   10254 
            
#delete mud but not mudstone etc
# the following keywords are kept for "mud"
tmp <- c("mudston","mudso", "mudtone", "mudst", "Mudsdtone", 
        "Mudshale", "Muderong", "mudsstone")
tmp <- key_del2("mud", tmp, lit.lst[["nsure"]])
lit.lst[["nsure"]] <- tmp[["lit"]]
lit.lst[["del"]] <- rbind(lit.lst[["del"]], tmp[["del"]])
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2570362   50457  276506  254920   29617   10254 

# check the ls.lst[["nouse"]] manually to make sure it right
#for (i in 1:length(ls.lst[["nouse"]]))
#{
#    tmp <- grep(pattern = ls.lst[["nouse"]][i], 
#       lit.lst[["nsure"]]$Descriptio, ignore.case=TRUE)
#    tmp2 <- lit.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["nouse"]][i]
#    write.table(tmp2, paste(dir, "\\nouse\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}

#--------------------------nouse
tmp <- key_nouse(ls.lst[["nouse"]], lit.lst[["nsure"]]$Descriptio)
lit.lst[["nouse"]] <- rbind(lit.lst[["nouse"]], lit.lst[["nsure"]][tmp, ])
lit.lst[["nsure"]] <- lit.lst[["nsure"]][-tmp, ]
p_dim(lit.lst)
#   del    keep   qsure   nsure   noinf   nouse 
#2570362   50457  276506  252379   29617   12795 

# check the ls.lst[["noinf"]]o manually to make sure it right
#for (i in 1:length(ls.lst[["noinf"]]))
#{
#    tmp <- grep(pattern = ls.lst[["noinf"]][i], 
#        lit.lst[["nsure"]]$Descriptio, ignore.case=TRUE)
#    tmp2 <- lit.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["noinf"]][i]
#    write.table(tmp2, paste(dir, "\\noinf\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}

#--------------------------noinf

tmp <- key_noinf(ls.lst[["noinf"]], lit.lst[["nsure"]]$Descriptio)
lit.lst[["noinf"]] <- rbind(lit.lst[["noinf"]], lit.lst[["nsure"]][tmp, ])
lit.lst[["nsure"]] <- lit.lst[["nsure"]][-tmp, ]
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2570362   50457  276506  243005   38991   12795 

# check the ls.lst[["qdel"]] manually to make sure it right
#for (i in 1:length(ls.lst[["qdel"]]))
#{
#    tmp <- grep(pattern = ls.lst[["qdel"]][i], lit.lst[["nsure"]]$Descriptio, 
#        ignore.case=TRUE)
#    tmp2 <- lit.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["qdel"]][i]
#    write.table(tmp2, paste(dir, "\\qdel\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}
#--------------------------qdel
tmp <- key_qdel(ls.lst[["qdel"]], lit.lst[["nsure"]]$Descriptio)
lit.lst[["qdel"]]  <- rbind(lit.lst[["qdel"]], lit.lst[["nsure"]][tmp, ])
lit.lst[["nsure"]] <- lit.lst[["nsure"]][-tmp, ]
p_dim(lit.lst)
# del    keep   qsure   nsure   noinf   nouse    qdel 
#2570362   50457  276506  207955   38991   12795   35050 

# check the ls.lst[["qdel"]] manually to make sure it right
#for (i in 1:length(ls.lst[["qsure"]]))
#{
#    tmp <- grep(pattern = ls.lst[["qsure"]][i], 
#        lit.lst[["nsure"]]$Descriptio, ignore.case=TRUE)
#    tmp2 <- lit.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["qsure"]][i]
#    write.table(tmp2, paste(dir, "\\qsure\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}

#--------------------------qsure
tmp <- key_qsure(ls.lst[["qsure"]], lit.lst[["nsure"]]$Descriptio)
lit.lst[["qsure"]] <- rbind(lit.lst[["qsure"]], lit.lst[["nsure"]][tmp, ])
lit.lst[["nsure"]] <- lit.lst[["nsure"]][-tmp, ]
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse    qdel 
#2570362   50457  442078   42383   38991   12795   35050 

tmp <- key_qsure(ls.lst[["qsure"]], lit.lst[["qdel"]]$Descriptio)
lit.lst[["qsure"]] <- rbind(lit.lst[["qsure"]], lit.lst[["qdel"]][tmp, ])
lit.lst[["qdel"]] <- lit.lst[["qdel"]][-tmp, ]
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse    qdel 
#2570362   50457  470368   42383   38991   12795    6760 
#--------------------------nsure
#get rid of meaningless character
i <- 1
while(i < 7)
{
    tmp <- str_sub(lit.lst[["nsure"]]$Descriptio, 1, 1)
    tmp <- (tmp <= "/") | (tmp >= ":" & tmp <= "@") | 
        (tmp >= "[" & tmp <= "`") | (tmp >= "{" & tmp <= "~")
    lit.lst[["nsure"]][tmp, "Descriptio"] <- str_sub(
            lit.lst[["nsure"]][tmp, "Descriptio"], 2)  
    i <- i + 1    
}
                             


#tmp <- count(lit.lst[["nsure"]]$Descriptio)
#write.table(tmp, paste(dir, "\\", "litlist52.txt", sep = ""), sep = "\t", 
#    row.names = FALSE)
#store as litlist5.txt and indicate which one to keep mannually, 
#juded according to description:
#0: delete, 1: keep, 2: quite sure, 3: not sure, 4: not clear what it is,
#5: if this layer exit, the whole points is not used, human disturbance
#6: qure sure delete,
#note that NULL and the numbers should be converted in text format in saving
# only  classified the records with freq bigger than 10  and fist few ordered
tmp <- readLines(paste(dir, "\\", "litlist5.txt", sep = "")) 
tmp <- tmp[2:length(tmp)]
tmp <- data.frame(t(sapply(tmp, str_sep)))
row.names(tmp) <- NULL
tmp[, 3] <- as.numeric(tmp[, 3])
colnames(tmp) <- c("Descriptio", "freq", "token")
ls.lst  <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
tmp2 <- c(0, 2, 4, 5)
for(i in tmp2)
{
    ls.lst[[i+1]] <- tmp[tmp$token == i, 1:2]
    tmp3 <- match_df(lit.lst[["nsure"]], ls.lst[[i+1]], on = "Descriptio")
    lit.lst[[i+1]] <- rbind(lit.lst[[i+1]], tmp3)
    lit.lst[["nsure"]] <- subset(lit.lst[["nsure"]], 
        !(rownames(lit.lst[["nsure"]]) %in% rownames(tmp3)))
}
p_dim(lit.lst)
#    del    keep   qsure   nsure   noinf   nouse    qdel 
#2570633   50457  470628   36713   44134   12844    6707 

#write.table(lit.lst[["del"]], paste(dir, "\\", lit_seps[1], sep = ""), 
#   sep = "$", row.names = FALSE)
#write.table(lit.lst[["keep"]], paste(dir, "\\", lit_seps[2], sep = ""),  
#   sep = "$", row.names = FALSE)
#write.table(lit.lst[["qsure"]], paste(dir, "\\", lit_seps[3], sep = ""), 
#   sep = "$", row.names = FALSE)
#write.table(lit.lst[["nsure"]], paste(dir, "\\", lit_seps[4], sep = ""), 
#   sep = "$", row.names = FALSE)
#write.table(lit.lst[["noinf"]], paste(dir, "\\", lit_seps[5], sep = ""),  
#   sep = "$", row.names = FALSE)
#write.table(lit.lst[["nouse"]], paste(dir, "\\", lit_seps[6], sep = ""),  
#   sep = "$", row.names = FALSE)
#check the records manually
tmp <- grepl(pattern = "Grav", lit.lst[["nsure"]]$Descriptio, ignore.case = TRUE, 
        fixed = FALSE) #& !grepl(pattern = "mudston", lit.lst[["nsure"]]$Descriptio, ignore.case = TRUE, 
        #fixed = FALSE) #& !grepl(pattern = "mudrock", lit.lst[["nsure"]]$Descriptio, ignore.case = TRUE, 
        #fixed = FALSE) #& !grepl(pattern = "silts", lit.lst[["nsure"]]$Descriptio, ignore.case = TRUE, 
        #fixed = FALSE)& !grepl(pattern = "Silted", lit.lst[["nsure"]]$Descriptio, ignore.case = TRUE, 
        #fixed = FALSE)
#tmp <-  dbf.lst[[4]]$Descriptio == "Slippery back"
ii <- order(lit.lst[["nsure"]][tmp, "FromDepth"])
tmp3 <- lit.lst[["nsure"]][tmp, ]
tmp3 <- tmp3[ii, ]
View(tmp3)
#library(foreign)
#library(plyr)
#library(stringr)
#options(stringsAsFactors = FALSE)
#dir <- "E:\\data\\soildata\\depth\\points\\australia\\ngis"
#load(paste0(dir,"\\as2.RData"))
save.image(paste0(dir_out,"\\as2.RData"))

#-----------------------------------------------  
#deal with hgus
#write.table(hgus, paste(dir, "\\", "hgus.txt", sep = ""), sep = "\t", 
#    row.names = FALSE)
#-----------------------------------------------
#classify according to keywords
tmp <- read.table(paste(dir, "\\", "keywords.txt", sep = ""), header = TRUE, 
    sep = "\t") 
ls.lst  <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
hgu.lst <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
for(tmp2 in names(ls.lst))
{
    ls.lst[[tmp2]] <- tmp[!(tmp[tmp2] == ""), tmp2] 
}

# The order of operation should be: del, keep,(+ clay, silt) nouse, noinfor, 
# qdel, qsure and pkeep
#--------------------------del
tmp <- key_del(ls.lst[["del"]], dbf.lst[["hgus"]]$Descriptio)
hgu.lst[["del"]]   <- dbf.lst[["hgus"]][tmp, ]
hgu.lst[["nsure"]] <- dbf.lst[["hgus"]][-tmp, ]
p_dim(hgu.lst)
#   del  nsure 
#170140 340686 
#--------------------------keep
#only one to keep "consolidate"
tmp <- key_keep(hgu.lst[["nsure"]], hgu.lst[["keep"]], hgu.lst[["del"]])
hgu.lst[["nsure"]] <- tmp[[1]]
hgu.lst[["keep"]]  <- tmp[[2]]
hgu.lst[["del"]]   <- tmp[[3]]
p_dim(hgu.lst)
#   del   keep  nsure 
#169829    582 340415 
#--------------------------correction 
#Caly to clay, slit to silt,
#clay stone,clay-stone(clastone) to claystone
#silt stone,silt-stone to siltstone, 
#mud stone to mudstone

tmp <- grepl(pattern = "caly", hgu.lst[["nsure"]]$Descriptio, 
        ignore.case = TRUE) & !grepl(pattern = "scaly",
        hgu.lst[["nsure"]]$Descriptio, ignore.case = TRUE)
hgu.lst[["nsure"]][tmp, "Descriptio"] <- gsub (pattern = "caly", 
        replacement = "clay", hgu.lst[["nsure"]][tmp, "Descriptio"],
         ignore.case = TRUE)
tmp2 <- c("slit", "silthe", "CALCISIlTITE", "silt stone", "silt-stone", 
            "silstone", "siltsone", "siltone", "clay stone", "clay-stone", 
            "claysone", "mud stone", "mud-stone", "mud rock")
tmp3 <- c("silt", "slithe", "CALCISLITITE", "siltstone", "siltstone", 
            "siltstone", "siltstone", "siltstone", "claystone", "claystone", 
            "claystone", "mudstone", "mudstone", "mudrock")
for(i in length(tmp2))
{
    hgu.lst[["nsure"]][, "Descriptio"] <- gsub (pattern = tmp2[i], 
        replacement = tmp3[i], hgu.lst[["nsure"]][, "Descriptio"],
        ignore.case = TRUE)
}

#delete clay, but not claystone etc.    
# the following keywords are kept for "clay"
tmp <- c("clayston", "claybound", "clayish", "clayee", "clayed", "clayband")
tmp <- key_del2("clay", tmp, hgu.lst[["nsure"]])
hgu.lst[["nsure"]] <- tmp[["lit"]]
hgu.lst[["del"]] <- rbind(hgu.lst[["del"]], tmp[["del"]])
p_dim(hgu.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2534619   50457  276506  290663   29617   10254 


#delete silt but not siltstone etc. 
# the following keywords are kept for "silt"
tmp <- c("siltston", "siltit", "siltrex", "Silted", "Siltband", "siltier")
tmp <- key_del2("silt", tmp, hgu.lst[["nsure"]])
hgu.lst[["nsure"]] <- tmp[["lit"]]
hgu.lst[["del"]] <- rbind(hgu.lst[["del"]], tmp[["del"]])
p_dim(hgu.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2559077   50457  276506  266205   29617   10254 
            
#delete mud but not mudstone etc
# the following keywords are kept for "mud"
tmp <- c("mudston","mudso", "mudtone", "mudst", "Mudsdtone", 
        "Mudshale", "Muderong", "mudsstone")
tmp <- key_del2("mud", tmp, hgu.lst[["nsure"]])
hgu.lst[["nsure"]] <- tmp[["lit"]]
hgu.lst[["del"]] <- rbind(hgu.lst[["del"]], tmp[["del"]])
p_dim(hgu.lst)
#    del    keep   qsure   nsure   noinf   nouse 
#2570362   50457  276506  254920   29617   10254 

# check the ls.lst[["nouse"]] manually to make sure it right
#for (i in 1:length(ls.lst[["nouse"]]))
#{
#    tmp <- grep(pattern = ls.lst[["nouse"]][i], 
#       hgu.lst[["nsure"]]$Descriptio, ignore.case=TRUE)
#    tmp2 <- hgu.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["nouse"]][i]
#    write.table(tmp2, paste(dir, "\\nouse\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}

#--------------------------nouse
tmp <- key_nouse(ls.lst[["nouse"]], hgu.lst[["nsure"]]$Descriptio)
hgu.lst[["nouse"]] <- rbind(hgu.lst[["nouse"]], hgu.lst[["nsure"]][tmp, ])
hgu.lst[["nsure"]] <- hgu.lst[["nsure"]][-tmp, ]
p_dim(hgu.lst)
#   del    keep   qsure   nsure   noinf   nouse 
#2570362   50457  276506  252379   29617   12795 

# check the ls.lst[["noinf"]]o manually to make sure it right
#for (i in 1:length(ls.lst[["noinf"]]))
#{
#    tmp <- grep(pattern = ls.lst[["noinf"]][i], 
#        hgu.lst[["nsure"]]$Descriptio, ignore.case=TRUE)
#    tmp2 <- hgu.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["noinf"]][i]
#    write.table(tmp2, paste(dir, "\\noinf\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}

#--------------------------noinf

tmp <- key_noinf(ls.lst[["noinf"]], hgu.lst[["nsure"]]$Descriptio)
hgu.lst[["noinf"]] <- rbind(hgu.lst[["noinf"]], hgu.lst[["nsure"]][tmp, ])
hgu.lst[["nsure"]] <- hgu.lst[["nsure"]][-tmp, ]
p_dim(hgu.lst)
#   del   keep  nsure  noinf  nouse 
#188791    582 298557  20886   2010  

# check the ls.lst[["qdel"]] manually to make sure it right
#for (i in 1:length(ls.lst[["qdel"]]))
#{
#    tmp <- grep(pattern = ls.lst[["qdel"]][i], hgu.lst[["nsure"]]$Descriptio, 
#        ignore.case=TRUE)
#    tmp2 <- hgu.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["qdel"]][i]
#    write.table(tmp2, paste(dir, "\\qdel\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}
#--------------------------qdel
tmp <- key_qdel(ls.lst[["qdel"]], hgu.lst[["nsure"]]$Descriptio)
hgu.lst[["qdel"]]  <- rbind(hgu.lst[["qdel"]], hgu.lst[["nsure"]][tmp, ])
hgu.lst[["nsure"]] <- hgu.lst[["nsure"]][-tmp, ]
p_dim(hgu.lst)
#   del   keep  nsure  noinf  nouse   qdel 
#188791    582 260042  20886   2010  38515 

# check the ls.lst[["qdel"]] manually to make sure it right
#for (i in 1:length(ls.lst[["qsure"]]))
#{
#    tmp <- grep(pattern = ls.lst[["qsure"]][i], 
#        hgu.lst[["nsure"]]$Descriptio, ignore.case=TRUE)
#    tmp2 <- hgu.lst[["nsure"]][tmp, ]
#    ii <- order(tmp2[, "FromDepth"])
#    tmp2 <- tmp2[ii, ]
#    tmp3 <- ls.lst[["qsure"]][i]
#    write.table(tmp2, paste(dir, "\\qsure\\", tmp3, ".txt", sep=""), 
#       row.names = FALSE,sep = "\t")
#}

#--------------------------qsure
tmp <- key_qsure(ls.lst[["qsure"]], hgu.lst[["nsure"]]$Descriptio)
hgu.lst[["qsure"]] <- rbind(hgu.lst[["qsure"]], hgu.lst[["nsure"]][tmp, ])
hgu.lst[["nsure"]] <- hgu.lst[["nsure"]][-tmp, ]
p_dim(hgu.lst)
#   del   keep  qsure  nsure  noinf  nouse   qdel 
#188791    582 208050  51992  20886   2010  38515 

tmp <- key_qsure(ls.lst[["qsure"]], hgu.lst[["qdel"]]$Descriptio)
hgu.lst[["qsure"]] <- rbind(hgu.lst[["qsure"]], hgu.lst[["qdel"]][tmp, ])
hgu.lst[["qdel"]] <- hgu.lst[["qdel"]][-tmp, ]
p_dim(hgu.lst)
#   del   keep  qsure  nsure  noinf  nouse   qdel 
#213211    582 251875  16677  20883   1466   6132 

####################################!!!!!!!!!!!!!!!!!!!!!
#tmp <- count(hgu.lst[["nsure"]]$Descriptio)
#names(tmp) <- c("Descriptio", "freq")
#write.table(tmp, paste(dir, "\\", "hgulist12.txt", sep = ""), sep = "\t", 
#    row.names = FALSE)
#write.table(hgu.lst[["nsure"]], paste(dir, "\\", "hgus.txt", sep = ""), sep = "\t", 
#    row.names = FALSE)
#store as hgulist1.txt and indicate which one to keep mannually, 
#juded according to description:
#0: delete, 1: keep, 2: quite sure, 3: not sure, 4: not clear what it is,
#5: if this layer exit, the whole points is not used, human disturbance
#6: qure sure delete,
#note that NULL and the numbers should be converted in text format in saving
# only  classified the records with freq bigger than 10  and fist few ordered
#according to http://dbforms.ga.gov.au/pls/www/geodx.strat_units.int
tmp <- readLines(paste(dir, "\\", "hgulist1.txt", sep = "")) 
tmp <- tmp[2:length(tmp)]   
tmp <- data.frame(t(sapply(tmp, str_sep)))
row.names(tmp) <- NULL
tmp[, 3] <- as.numeric(tmp[, 3])
colnames(tmp) <- c("Descriptio", "freq", "token")
ls.lst  <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)
tmp2 <- c(0, 2, 4, 5)
for(i in tmp2)
{
    ls.lst[[i+1]] <- tmp[tmp$token == i, 1:2]
    tmp3 <- match_df(hgu.lst[["nsure"]], ls.lst[[i+1]], on = "Descriptio")
    hgu.lst[[i+1]] <- rbind(hgu.lst[[i+1]], tmp3)
    hgu.lst[["nsure"]] <- subset(hgu.lst[["nsure"]], 
        !(rownames(hgu.lst[["nsure"]]) %in% rownames(tmp3)))
}
p_dim(hgu.lst)
#   del   keep  qsure  nsure  noinf  nouse   qdel 
#213416    582 253539   2361  33330   1466   6132 

#library(foreign)
#library(plyr)
#library(stringr)
#options(stringsAsFactors = FALSE)
#dir <- "E:\\data\\soildata\\depth\\points\\australia\\ngis"
#load(paste0(dir,"\\as3.RData"))
save.image(paste0(dir_out,"\\as3.RData"))
str(dbf.lst[["bores"]])
str(cons_out)
str(lit.lst)
str(hgu.lst)
dbf.lst[["bores"]] <- rename(dbf.lst[["bores"]], c("HydroID"="BoreID"))

maxd.lst <- list(del = NULL, keep = NULL, qsure = NULL, nsure = NULL, 
    noinf = NULL, nouse = NULL, qdel =NULL)

tmp <- rbind(cons_out[, c("BoreID", "FromDepth")], 
        lit.lst[[2]][, c("BoreID", "FromDepth")], 
        hgu.lst[[2]][, c("BoreID", "FromDepth")])    

 
for(i in 1:7)
{
    if(i == 2 )
    {
        tmp <- rbind(cons_out[, c("BoreID", "FromDepth", "ToDepth")], 
            lit.lst[[i]][, c("BoreID", "FromDepth", "ToDepth")], 
            hgu.lst[[i]][, c("BoreID", "FromDepth", "ToDepth")])
    }else
    {
        tmp <- rbind(
            lit.lst[[i]][, c("BoreID", "FromDepth", "ToDepth")], 
            hgu.lst[[i]][, c("BoreID", "FromDepth", "ToDepth")])    
    }
    maxd.lst[[i]] <- aggregate(tmp, by = list(tmp$BoreID), max)[ , 2:4]
}


#use keep - qsure  - nsure - nouse - noinf  -qdel (strictest)
tmp  <- maxd.lst[["keep"]]
tmp2 <- tmp$BoreID %in% c(maxd.lst[["qsure"]]$BoreID, 
        maxd.lst[["nsure"]]$BoreID, maxd.lst[["noinf"]]$BoreID, 
        maxd.lst[["nouse"]]$BoreID, maxd.lst[["qdel"]]$BoreID)       
wells1 <- get_dtb(dbf.lst[["bores"]], tmp, !tmp2)

# use keep + qsure - nouse (most relaxed)
tmp  <- rbind(maxd.lst[["keep"]], maxd.lst[["qsure"]])
tmp  <- aggregate(tmp, by = list(tmp$BoreID), max)[ , 2:4]
tmp2 <- tmp$BoreID %in% maxd.lst[["nouse"]]$BoreID 
wells2 <- get_dtb(dbf.lst[["bores"]], tmp, !tmp2)
# use keep + qsure - nsure - nouse - noinf  - qdel 
tmp  <- rbind(maxd.lst[["keep"]], maxd.lst[["qsure"]])
tmp  <- aggregate(tmp, by = list(tmp$BoreID), max)[ , 2:4]
tmp2 <- tmp$BoreID %in% c( 
        maxd.lst[["nsure"]]$BoreID, maxd.lst[["noinf"]]$BoreID, 
        maxd.lst[["nouse"]]$BoreID, maxd.lst[["qdel"]]$BoreID)       
wells3 <- get_dtb(dbf.lst[["bores"]], tmp, !tmp2)

write.table(wells3, paste(dir_out, "wells_as.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
write.table(dim(wells3)[1], paste(dir_out, "states_as.txt", sep = ""), 
    row.names = FALSE, sep = "\t")    
    
write.table(wells1, paste(dir_out, "wells_as2.txt", sep = ""), 
    row.names = FALSE, sep = "\t")
write.table(dim(wells1)[1], paste(dir_out, "states_as2.txt", sep = ""), 
    row.names = FALSE, sep = "\t")  
    

 
#library(foreign)
#library(plyr)
#library(stringr)
#options(stringsAsFactors = FALSE)
#dir <- "E:\\data\\soildata\\depth\\points\\australia\\ngis"
#load(paste0(dir,"\\as3.RData"))
save.image(paste0(dir_out,"\\as.RData"))
