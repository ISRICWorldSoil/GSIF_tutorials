# title         : rw_AFSP.R
# purpose       : Reading and writing of Africa Soil Profile data;
# reference     : Africa Soil Profiles Database [http://africasoils.net/data/legacyprofile]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Jun 2012.
# inputs        : "http://www.isric.org/sites/default/files/private/datasets/AfSP01Qry_ISRIC.zip" DBF tables
# outputs       : R data frame for SoilProfileCollection;
# remarks 1     : afsp.rda is only a subset of the complete database;

library(aqp)

## Download the database:
download.file("http://www.isric.org/sites/default/files/private/datasets/AfSP01Qry_ISRIC.zip", destfile=paste(getwd(),"AfSP01Qry_ISRIC.zip", sep="/"))
unzip("AfSP01Qry_ISRIC.zip")

# import tables:
Profiles <- read.dbf("AfSP01Qry_ISRIC/GIS_Dbf/AfSP01Qry_Profiles.dbf")
Layers <- read.dbf("AfSP01Qry_ISRIC/GIS_Dbf/AfSP01Qry_Layers.dbf")
# View(Profiles$dbf[1:200,])
# View(Layers$dbf[1:200,])

# select columns of interest:
s.lst <- c("ProfileID", "Reliab", "X_LonDD", "Y_LatDD", "XYAccur", "T_Year", "WRB06", "WRB06rg", "FAO88", "USDA", "LocalCls")
h.lst <- c("ProfileID", "LayerNr", "UpDpth", "LowDpth", "HorDes", "ColorM", "ColorD", "FldTxtr", "CfPc", "Sand", "Silt", "Clay", "BlkDens", "PHH2O", "PHKCl", "OrgC", "TotalN", "Ecec", "Bsat", "CecSoil", "LabTxtr", "VolAWC")

# correct soil site and horizon data:
sites <- Profiles$dbf[!duplicated(Profiles$dbf$ProfileID),]
sites <- sites[!is.na(sites$X_LonDD)&!is.na(sites$Y_LatDD), s.lst]
sites$XYAccur <- ifelse(sites$XYAccur==-9999, NA, sites$XYAccur)
sites$T_Year <- as.Date(sites$T_Year)
sites$WRB06[which(sites$WRB06=="NA")] <- NA
sites$WRB06rg[which(sites$WRB06rg=="NA")] <- NA
sites$FAO88[which(sites$FAO88=="NA")] <- NA
sites$USDA[which(sites$USDA=="NA")] <- NA
sites$LocalCls[which(sites$LocalCls=="NA")] <- NA
summary(!is.na(sites$WRB06rg))
summary(!is.na(sites$USDA))
summary(!is.na(sites$LocalCls))
str(sites)

horizons <- Layers$dbf[,h.lst]
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
for(j in c("CfPc", "Sand", "Silt", "Clay", "BlkDens", "PHH2O", "PHKCl", "OrgC", "TotalN", "Ecec", "Bsat", "CecSoil", "LabTxtr", "VolAWC")){
  horizons[,j] <- signif(ifelse(horizons[,j]==-9999, NA, horizons[,j]), 2)
}
str(horizons)

# rename columns:
names(sites) <- c("SOURCEID", "RLBISRIC", "LONWGS84", "LATWGS84", "STDXYZ", "TIMESTRR","TAXNWRB", "TAXGWRB", "TAXNFAO", "TAXNUSDA", "TAXN")
names(horizons) <- c("SOURCEID", "LSQINT", "UHDICM", "LHDICM", "HZDFAO", "MCOMNS", "DCOMNS", "TEXMHT", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHO5", "PHIKCL", "ORCDRC", "NTO", "ECE", "BST", "CEC", "TEXMHL", "AWCIMM")

# Save object:
afsp <- list(sites=sites[,c("SOURCEID", "RLBISRIC", "LONWGS84", "LATWGS84", "STDXYZ", "TIMESTRR","TAXNWRB", "TAXGWRB", "TAXNUSDA")], horizons=horizons[,c("SOURCEID", "LSQINT", "UHDICM", "LHDICM", "MCOMNS", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHO5", "PHIKCL", "ORCDRC")])
str(afsp)
save(afsp, file="afsp.rda", compress="xz")

# end of script;