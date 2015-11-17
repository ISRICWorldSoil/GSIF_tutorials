# title         : TAXG_1km_AfSIS.R
# purpose       : Mapping soil types (probabilities) for Africa (1 km resolution);
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.isric.org/]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles "afsp.rda" and WorldGrids maps
# outputs       : Predictions of pH, clay content and organic carbon at 6 depths (point support);
# remarks 1     : This is LARGE data! Not recommended for a PC without at least 16GB of RAM and multicore processor. The script is extra long because at many places we were forced to run things in loops i.e. tile by tile;


#setwd("E:\\gsif\\AFRICA\\1km")
#install.packages("GSIF", repos=c("http://R-Forge.R-project.org"))
#library(GSIF)
library(maptools)
library(rgdal)

#load(file="TAXGWRB_AfSIS.RData")
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
fw.path <- utils::readRegistry("SOFTWARE\\WOW6432NODE\\FWTools")$Install_Dir
gdalwarp <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))

## Load input data
## Points:
data(afsp)
s.lst <- c("SOURCEID", "LONWGS84", "LATWGS84", "TAXGWRB", "TAXNFAO")
afsp.xy <- afsp$sites[,s.lst]
coordinates(afsp.xy) <- ~ LONWGS84 + LATWGS84
proj4string(afsp.xy) <- "+proj=latlong +datum=WGS84"
## Covariates:
if(is.na(file.info("afgrid1km.rda")$size)){
  download.file("http://worldgrids.org/rda/afgrid1km.rda", "afgrid1km.rda")
}
load("afgrid1km.rda")
afsp.xy <- spTransform(afsp.xy, afgrid1km@proj4string)
summary(afsp.xy$TAXNFAO) ## TH: Too many classes!
summary(afsp.xy$TAXGWRB)

## clean up:
afsp.xy$TAXGWRB[afsp.xy$TAXGWRB == "Vr"] <- NA
afsp.xy$TAXGWRB <- as.factor(ifelse(afsp.xy$TAXGWRB == "NA", NA, paste(afsp.xy$TAXGWRB)))
summary(afsp.xy$TAXGWRB)

## Multinomial logistic regression:
formulaString2 = as.formula(paste('TAXGWRB ~ DEMSRE3a + SLPSRT3a + TDMMOD3a + TDSMOD3a + TNMMOD3a + TNSMOD3a + EVMMOD3a + EVSMOD3a + PREGSM1a +', paste("G0", c(1:7,9), "ESA3a", sep="", collapse="+"), '+', paste("G", 10:21, "ESA3a", sep="", collapse="+"), ' + SGEUSG3x + TWISRE3a + L3POBI3a'))
sel2 = names(afgrid1km) %in% all.vars(formulaString2)[-1]
ov <- over(afsp.xy, afgrid1km[sel2])  
ov <- cbind(data.frame(afsp.xy), ov)
gc()
str(ov)
summary(is.na(ov$TDMMOD3a)&!is.na(ov$TAXGWRB))
## 857 points fall outside the soil mask!
    
require(nnet)
mout <- nnet::multinom(formulaString2, ov, MaxNWts = 2300)
## predict dominant class:
L = 10  # number of tiles
sel <- c(seq(1, nrow(afgrid1km), by=round(nrow(afgrid1km)/L)), nrow(afgrid1km))

afgrid1km$TAXGWRB.p <- NULL
gc()
for(j in 1:L){  ## takes about 10 mins and > 8GB!
    gc()
    afgrid1km@data[sel[j]:sel[j+1], "TAXGWRB.p"] <- as.factor(paste(predict(mout, newdata=afgrid1km[sel[j]:sel[j+1],], na.action = na.pass)))
    gc()
}
GWRB.lvs <- levels(afgrid1km$TAXGWRB.p)
afgrid1km$TAXGWRB.p <- as.integer(afgrid1km$TAXGWRB.p)
gc()
## write to GDAL:
writeGDAL(afgrid1km["TAXGWRB.p"], "TAXGWRB3a.tif", "GTiFF", type="Byte", mvFlag=0)
gc()
## write a table with class names:
write.table(data.frame(Codes=1:length(GWRB.lvs), TAXGWRB=GWRB.lvs), file="TAXGWRB_legend.txt")

## predict probabilities:
probs <- predict(mout, newdata=afgrid1km[1:1000,], type="probs", na.action = na.pass)
mm <- afgrid1km[1:length(attr(probs, "dimnames")[[2]])]
names(mm) <- attr(probs, "dimnames")[[2]]
gc()
for(j in 1:L){
    gc()
    gc()
    mm@data[sel[j]:sel[j+1],] <- round(data.frame(predict(mout, newdata=afgrid1km[sel[j]:sel[j+1],], type="probs", na.action = na.pass))*100)
    gc()
}
str(mm@data)
## write probabilities to geotifs:
for(i in 1:ncol(mm)){
  gc()
  writeGDAL(mm[i], paste("TAXGWRB3a_", names(mm)[i], ".tif", sep=""), "GTiFF", type="Byte", mvFlag=255)
  gc()
}

## kappa statistics:
require(mda)
require(psych)
cout.m <- as.factor(paste(predict(mout, newdata=ov, na.action = na.pass)))
cf <- confusion(cout.m, as.character(ov[,"TAXGWRB"]))
## remove missing classes:
a = attr(cf, "dimnames")[[1]] %in% attr(cf, "dimnames")[[2]] 
b = attr(cf, "dimnames")[[2]] %in% attr(cf, "dimnames")[[1]]
c.kappa = cohen.kappa(cf[a,b])
ac <- sum(diag(cf))/sum(cf)*100
## number of subjects: 5592
message(paste("Estimated Cohen Kappa (weighted):", signif(c.kappa$weighted.kappa, 4)))  ##  17%
message(paste("Map purity:", signif(ac, 3)))  ## 21%


## compress all produced maps:
system("7za a TAXGWRB_1km_mnr.tif.7z TAXGWRB3a*.tif")
## save the model:
save(mout, file="TAXG.m.rda", compress="xz")

#save.image(file="TAXGWRB_AfSIS.RData")

# end of script; 