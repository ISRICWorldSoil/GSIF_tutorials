# title         : samples_AfSIS.R
# purpose       : Preparing soil samples for Africa;
# reference     : Methodology for global soil mapping from GBIF package [http://gsif.r-forge.r-project.org/]
# producer      : Prepared by T. Hengl (tom.hengl@wur.nl)
# address       : In Wageningen, NL.
# inputs        : Africa soil profiles (afsp) from the GSIF package
# outputs       : Rda file with soil profiles formated as "geosamples";
# remarks 1     : ;

library(aqp)
library(GSIF)
library(plyr)

data(afsp)
s.lst <- c("SOURCEID", "LONWGS84", "LATWGS84", "TAXNWRB", "TAXNFAO")
h.lst <- c("SOURCEID", "UHDICM", "LHDICM", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "PHIHO5", "ORCDRC", "CEC")
# create object of type "SoilProfileCollection"
afsp.spc <- join(afsp$horizons[,h.lst], afsp$sites[,s.lst], type='inner')
depths(afsp.spc) <- SOURCEID ~ UHDICM + LHDICM
site(afsp.spc) <- ~ LONWGS84 + LATWGS84 + TAXNWRB + TAXNFAO
# convert to geosamples:
coordinates(afsp.spc) <- ~ LONWGS84 + LATWGS84
proj4string(afsp.spc) <- "+proj=latlong +datum=WGS84"
str(afsp.spc)
afsp.geo <- as.geosamples(afsp.spc)
str(afsp.geo)
## 569,620 observations;

## save file:
save(afsp.geo, file="afsp.geo.rda", compress="xz")

## save as a CSV file:
#write.csv(afsp.geo@data, "afsp.geo.csv")
#system("7za a -tgzip afsp.geo.csv.gz afsp.geo.csv")
#unlink("afsp.geo.csv")


# end of script;