# Purpose        : Derive Soil organic carbon stock in kg / m^2;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Gerard Heuvelink (gerard.heuvelink@wur.nl);
# Dev Status     : Stable
# Note           : Original formula available at e.g. [http://www.eoearth.org/view/article/156087/]


OCSKGM <- function(ORCDRC, BLD=1400, CRFVOL=0, HSIZE, ORCDRC.sd=10, BLD.sd=100, CRFVOL.sd=5, se.prop=TRUE){
    if(any(ORCDRC[!is.na(ORCDRC)]<0)|any(BLD[!is.na(BLD)]<0)|any(CRFVOL[!is.na(CRFVOL)]<0)){
       warning("Negative values for 'ORCDRC', 'BLD', 'CRFVOL' found")
    }
    OCSKG <- ORCDRC/1000 * HSIZE/100 * BLD * (100-CRFVOL)/100
    if(se.prop==TRUE){
       if(any(ORCDRC.sd[!is.na(ORCDRC.sd)]<0)){
          ORCDRC.sd = ifelse(is.na(ORCDRC.sd)|ORCDRC.sd<0, 0, ORCDRC.sd); warning("Replacing negative values for 'ORCDRC.sd'")
       }
       if(any(BLD.sd[!is.na(BLD.sd)]<0)){
          BLD.sd = ifelse(is.na(BLD.sd)|BLD.sd<0, 0, BLD.sd); warning("Replacing negative values for 'BLD.sd'")
       }
       if(any(CRFVOL.sd[!is.na(CRFVOL.sd)]<0)){
          CRFVOL.sd = ifelse(is.na(CRFVOL.sd)|CRFVOL.sd<0, 0, CRFVOL.sd); warning("Replacing negative values for 'CRFVOL.sd'")
       }
       ## Formula derived by Gerard. See also: [http://books.google.nl/books?id=C\_XWjSsboeUC]
       OCSKG.sd <- 1E-7*HSIZE*sqrt(BLD^2*(100-CRFVOL)^2*ORCDRC.sd^2 + ORCDRC^2*(100-CRFVOL)^2*BLD.sd^2 + ORCDRC^2*BLD^2*CRFVOL.sd^2)
       attr(OCSKG, "measurementError") <- signif(OCSKG.sd, 3)
       attr(OCSKG, "units") <- "kilograms per square-meter"
    }
    return(OCSKG)
}

## end of script;