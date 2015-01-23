\name{LRI}
\alias{LRI}
\alias{ERZDICM}
\title{Limiting Rootability}
\description{Derive Limiting Rootability using observed soil properties at at least three depths.}
\usage{LRI(UHDICM, LHDICM, SNDPPT, SLTPPT, CLYPPT, CRFVOL, BLD,
   ORCDRC, ECN, CEC, ENA, EACKCL, EXB, PHIHOX, CRB, GYP, tetaS, 
   fix.values=TRUE, thresholds, print.thresholds=FALSE)
}
\arguments{
  \item{UHDICM}{numeric; upper horizon depth in cm}
  \item{LHDICM}{numeric; lower horizon depth in cm}
  \item{SNDPPT}{numeric; sand content in percent}
  \item{SLTPPT}{numeric; silt content in percent}
  \item{CLYPPT}{numeric; clay content in percent}
  \item{CRFVOL}{numeric; volume percentage of coarse fragments (> 2 mm)}
  \item{BLD}{numeric; bulk density in kg per cubic-meter for the horizon/solum}
  \item{ORCDRC}{numeric; soil organic carbon concentration in permille or g per kg}
  \item{ECN}{numeric; electrical conductivity in dS per m}
  \item{CEC}{numeric; Cation Exchange Capacity in cmol per kilogram}
  \item{ENA}{numeric; exchangable Na in cmol per kilogram}
  \item{EACKCL}{numeric; exchangable acidity in cmol per kilogram}  
  \item{EXB}{numeric; exchangable bases in cmol per kilogram}  
  \item{PHIHOX}{numeric; soil pH in water suspension}
  \item{CRB}{numeric; CaCO3 (carbonates) in g per kg}  
  \item{GYP}{numeric; CaSO4 (gypsum) in g per kg}
  \item{tetaS}{numeric; volumetric percentage (optional; if not provided it will be derived using the \code{AWCPTF} Pedo-Transfer Function)}
  \item{fix.values}{logical; specifies whether to correct values of textures and bulk density to avoid creating nonsensical values}
  \item{thresholds}{data.frame; optional table containing threshold values for \code{"CRFVOL"}, \code{"tetaS"} (volumetric percentage), \code{"BLD.f"} (clay-adjusted BLD), \code{"SNDPPT"}, \code{"CLY.d"} (difference in clay between horizons), \code{"SND.d"} (difference in sand between horizons), \code{"PHIHOX.L"} (lower limits for pH), \code{"PHIHOX.H"} (upper limits for pH), \code{"ECN"}, \code{"ENA.f"} (exchangable saturated Na), \code{"ENA"}, \code{"EACKCL.f"} (exchangable saturated acidity), \code{"CRB"} (carbonates), and \code{"GYP"} (gypsum)}
  \item{print.thresholds}{logical; specifies whether to attach the threshold values to the output object}
}
\value{
Returns a vector with \code{TRUE} / \code{FALSE} values where \code{FALSE} indicates rooting not possible. Threshold values used to derive Limiting Rootability scores are set based on common soil agricultural productivity tresholds (e.g. in this case for maize), and can be adjusted via the \code{thresholds} argument. This functions also accounts for textural changes (sudden changes in sand and clay content) and saturated water content.
}
\note{
Horizons need to be sorted by depth e.g. 0-5, 5-15, 15-30\dots For each soil property at least three depths are needed otherwise the function reports an error. Missing values are automatically replaced using smoothing splines.
}
\author{ Johan Leenaars and Maria Ruiperez Gonzalez  }
\seealso{ \code{\link{AWCPTF}} }
\examples{
## sample profile from Nigeria (ISRIC:NG0017):
UHDICM = c(0, 18, 36, 65, 87, 127)
LHDICM = c(18, 36, 65, 87, 127, 181)
SNDPPT = c(66, 70, 54, 43, 35, 47)
SLTPPT = c(13, 11, 14, 14, 18, 23)
CLYPPT = c(21, 19, 32, 43, 47, 30)
CRFVOL = c(17, 72, 73, 54, 19, 17)
BLD = c(1.57, 1.60, 1.52, 1.50, 1.40, 1.42)*1000
PHIHOX = c(6.5, 6.9, 6.5, 6.2, 6.2, 6.0)
CEC = c(9.3, 4.5, 6.0, 8.0, 9.4, 10.9)
ENA = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.2)
EACKCL = c(0.1, 0.1, 0.1, NA, NA, 0.5)
EXB = c(8.9, 4.0, 5.7, 7.4, 8.9, 10.4)
ORCDRC = c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
x <- LRI(UHDICM=UHDICM, LHDICM=LHDICM, SNDPPT=SNDPPT, 
   SLTPPT=SLTPPT, CLYPPT=CLYPPT, CRFVOL=CRFVOL, 
   BLD=BLD, ORCDRC=ORCDRC, CEC=CEC, ENA=ENA, EACKCL=EACKCL, 
   EXB=EXB, PHIHOX=PHIHOX, print.thresholds=TRUE)
x
## Most limiting: BLD.f and CRFVOL, but nothing < 20

## Effective Rootable Depth:
sel <- x==FALSE
if(!all(sel==FALSE)){ 
  UHDICM[which(sel==TRUE)[1]] 
} else {
  max(LHDICM)
}

xI <- attr(x, "minimum.LRI")
## derive Effective rooting depth:
ERDICM(UHDICM=UHDICM, LHDICM=LHDICM, minimum.LRI=xI, DRAINFAO="M")
}