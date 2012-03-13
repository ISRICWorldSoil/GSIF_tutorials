# title         : rw_NCSS.R
# purpose       : Reading and writing of NSCD data;
# reference     : NCSS Characterization Database [http://ssldata.nrcs.usda.gov/]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Feb 2012.
# inputs        : "Repository2010.mdb" MS Access dbs
# outputs       : R data frames for SoilProfileCollection;
# remarks 1     : LARGE dataset!

## Download the database:
# download.file("http://globalsoilmap.net/data/Repository2010.7z", destfile=paste(getwd(),"Repository2010.7z",sep="/"))

library(RODBC)
library(aqp)
library(sp)
# define a new function to merge the degree, min, sec columns:
cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}
# load("Repository2010.RData")

# ------------------------------------------------------------
# Fetch tables - NCSS
# ------------------------------------------------------------

cNCSS <- odbcConnect(dsn="NCSS")
sqlTables(cNCSS)$TABLE_NAME
# get tables:
site <- sqlFetch(cNCSS, "site", stringsAsFactors=FALSE)
str(site)  # 37,333 profiles!!!
pedon <- sqlFetch(cNCSS, "pedon", stringsAsFactors=FALSE)
str(pedon)
pedon$observation_date <- as.character(as.Date(pedon$observation_date))
# has to remove Date format otherwise it gets problems with NA's
PSDA <- sqlFetch(cNCSS, "PSDA and Rock Fragments", stringsAsFactors=FALSE)   ## gravel content ("wpg2") is in mass percentage and needs to be converted to volume %
str(PSDA)
Organic <- sqlFetch(cNCSS, "Organic", stringsAsFactors=FALSE)
str(Organic)
CEC <- sqlFetch(cNCSS, "CEC and Bases", stringsAsFactors=FALSE)
str(CEC)
Carbon <- sqlFetch(cNCSS, "Carbon and Extractions", stringsAsFactors=FALSE)
str(Carbon)
BulkDens <- sqlFetch(cNCSS, "Bulk Density and Moisture", stringsAsFactors=FALSE)
str(BulkDens)  ## we need "Bulk Density, <2mm Fraction, Ovendry"
pH <- sqlFetch(cNCSS, "pH and Carbonates", stringsAsFactors=FALSE)
str(pH)
layer <- sqlFetch(cNCSS, "layer", stringsAsFactors=FALSE)
str(layer)
tax <- sqlFetch(cNCSS, "Taxonomy_New", stringsAsFactors=FALSE)
tax$last_correlated_date <- as.character(as.Date(tax$last_correlated_date))
Phosphorus <- sqlFetch(cNCSS, "Phosphorus", stringsAsFactors=FALSE)
MajorElements <- sqlFetch(cNCSS, "Major Elements", stringsAsFactors=FALSE)
str(MajorElements)
Salt <- sqlFetch(cNCSS, "Salt", stringsAsFactors=FALSE)
str(Salt)
TaxonomyE <- sqlFetch(cNCSS, "Taxonomy_Error", stringsAsFactors=FALSE)
str(TaxonomyE)
TaxonomyE$last_correlated_date <- as.character(as.Date(TaxonomyE$last_correlated_date))
procs <- sqlFetch(cNCSS, "analysis_procedure", stringsAsFactors=FALSE)
# save(list=c("Organic", "site", "pedon", "PSDA", "CEC", "Carbon", "BulkDens", "pH", "layer", "tax", "Phosphorus", "MajorElements", "Salt", "TaxonomyE"), file="cNCSS.RData")


# ------------------------------------------------------------
# Re-format columns
# ------------------------------------------------------------

# add missing columns:
layer$DEPTH <- layer$hzn_top + (layer$hzn_bot - layer$hzn_top)/2
# mask out profiles with missing coordinates:
site <- site[!is.na(site$longitude_degrees)&!is.na(site$latitude_degrees),]
site$longitude_it <- ifelse(site$longitude_direction=="west", "W", "E")
site$latitude_it <- ifelse(site$latitude_direction=="north"|site$latitude_direction=="North", "N", "S")
site$LAT <- cols2dms(site$latitude_degrees, site$latitude_minutes, site$latitude_seconds, site$latitude_it)
site$LON <- cols2dms(site$longitude_degrees, site$longitude_minutes, site$longitude_seconds, site$longitude_it)
summary(site$LAT) # 253 NA's
summary(site$LON)
summary(site$latitude_seconds)
str(site)

# ------------------------------------------------------------
# create the horizon and site tables (takes time!)
# ------------------------------------------------------------
 
h1 <- merge(PSDA, CEC[,!(names(CEC) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
h2 <- merge(Organic[,!(names(Organic) %in% c("result_source_key", "prep_code", "c_tot", "oc", "n_tot", "c_n_ra"))], Carbon, by=c("natural_key"), all=TRUE)
h3 <- merge(h1, h2[,!(names(h2) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
h4 <- merge(h3, layer, by=c("natural_key"), all=TRUE)
h5 <- merge(h4[,!(names(h4) %in% c("db_od"))], BulkDens[,!(names(BulkDens) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
horizon <- merge(h5[,!(names(h5) %in% c("ph_h2o"))], pH[,!(names(pH) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
# names(horizon)
## fix some columns:
horizon$wpg2 <- ifelse(horizon$wpg2 > 100|horizon$wpg2 <0 , NA, horizon$wpg2)
## merge BulkDens and PSDA and derive GRAVEL content:
mBD <- (horizon$wpg2/100 * 2.6 + (1-horizon$wpg2/100) * horizon$db_od)
horizon$GRAVEL <- horizon$wpg2*mBD/2.6
## check visually:
# plot(y=horizon$GRAVEL, x=horizon$wpg2, xlim=c(0,100), ylab="GRAVEL (vol %)", xlab="GRAVEL (mass %)", main="NCSS (250K measurements)", pch="+", cex=.7)


# subset tables:
hs <- subset(horizon[,c("site_key", "natural_key", "layer_sequence", "hzn_desgn", "hzn_bot", "hzn_top", "hzn_vert_subdvn", "clay_tot_psa", "sand_tot_psa", "silt_tot_psa", "pyr_col", "wpg2", "caco3", "ph_hist", "oc", "c_tot", "base_sum", "cec_sum", "cec_nh4", "ph_h2o", "ph_kcl", "db_od")], !is.na(horizon$DEPTH)&!is.na(horizon$site_key)&horizon$layer_sequence<15)
# strange - there are layer_sequence numbers up to 259?!
str(hs)
# summary(as.factor(horizon$layer_sequence))
hs$site_key <- as.factor(hs$site_key)
# remove duplicate site keys (there are many!!):
hs$layer_ID <- as.factor(paste(hs$site_key, hs$layer_sequence, sep="_"))
hs <- hs[!duplicated(hs$layer_ID),]
length(levels(as.factor(hs$layer_ID)))
mean(hs$oc, na.rm=TRUE) ## Organic Carbon in %!
# estimate OC by correcting for CaCO3:
hs$ORCDRC <- 10*ifelse(!is.na(hs$c_tot), ifelse((hs$ph_h2o > 7)&!is.na(hs$caco3), hs$c_tot - .12 * hs$caco3, hs$c_tot), hs$oc)  
hs$ORCDRC <- ifelse(hs$ORCDRC < 0, 0, hs$ORCDRC) 
# stats:
summary(hs$ORCDRC)

pedon$site_key <- as.factor(pedon$site_key)
# there are also pedons with multiple site IDs!?
s1 <- merge(pedon, tax[,!(names(tax) %in% c("natural_key"))], by=c("pedon_key"), all=TRUE)
# summary(s1$site_key)  
s1 <- s1[!duplicated(s1$site_key),]
str(s1)
# merge the site table:
site$site_key <- as.factor(site$site_key)
site.s <- merge(site, s1, by=c("site_key"), all.x=TRUE)
site.s <- site.s[!is.na(site.s$LAT)&!is.na(site.s$LON), c("site_key","user_site_id","LAT","LON","horizontal_datum_name","observation_date","sampled_taxon_name","correlated_taxon_name","correlated_taxon_kind","correlated_class_name","SSL_class_name")]
str(site.s)

# ------------------------------------------------------------
# SoilProfileCollection
# ------------------------------------------------------------

NCSS_all <- list(sites=site.s, horizons=hs)
str(NCSS_all)
# save(NCSS_all, file="NCSS_all.rda", compress="xz")
# Promote to SoilProfileCollection
NCSS.spc <- join(NCSS_all$horizons, NCSS_all$sites, type='inner')
depths(NCSS.spc) <- site_key ~ hzn_top + hzn_bot
# TAKES 4-5 mins!
# extract site data
site(NCSS.spc) <- ~ LON + LAT + horizontal_datum_name + observation_date + sampled_taxon_name + correlated_taxon_name + correlated_taxon_kind + correlated_class_name + SSL_class_name
# generate SpatialPoints
coordinates(NCSS.spc) <- ~ LON + LAT
# assign CRS data
proj4string(NCSS.spc) <- "+proj=latlong +datum=NAD83"
str(NCSS.spc)
# library(plotKML)
# kml(NCSS.spc, var.name="ORCDRC")

# subset to Indiana state:
sel <- site.s$LAT>37.2&site.s$LAT<42.2&site.s$LON< -84.5&site.s$LON> -87.5
site.s_in <- site.s[sel,]
hs_in <- hs[hs$site_key %in% site.s_in$site_key,]
NCSS <- list(sites=site.s_in, horizons=hs_in)
# round up the numbers:
NCSS$horizons$ORCDRC <- round(NCSS$horizons$ORCDRC, 1)
NCSS$horizons$clay_tot_psa <- round(NCSS$horizons$clay_tot_psa, 0)
NCSS$horizons$sand_tot_psa <- round(NCSS$horizons$sand_tot_psa, 0)
NCSS$horizons$silt_tot_psa <- round(NCSS$horizons$silt_tot_psa, 0)
NCSS$horizons$oc <- round(NCSS$horizons$oc, 2)
NCSS$horizons$c_tot <- round(NCSS$horizons$c_tot, 2)
NCSS$horizons$base_sum <- round(NCSS$horizons$base_sum, 1)
NCSS$horizons$caco3 <- round(NCSS$horizons$caco3, 2)
NCSS$horizons$cec_sum <- round(NCSS$horizons$cec_sum, 1)
NCSS$horizons$cec_nh4 <- round(NCSS$horizons$cec_nh4, 1)
NCSS$horizons$ph_h2o <- round(NCSS$horizons$ph_h2o, 1)
NCSS$horizons$ph_kcl <- round(NCSS$horizons$ph_kcl, 1)
NCSS$horizons$db_od <- round(NCSS$horizons$db_od, 2)
# Save object:
save(NCSS, file="NCSS.rda", compress="xz")
# save(NCSS, file="NCSS.rda", compress="gzip")
save.image("Repository2010.RData")

# end of script;