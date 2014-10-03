\name{AWCPTF}
\encoding{latin1}
\alias{AWCPTF}
\title{Available soil water capacity }
\description{Derive available soil water capacity  (in cubic-meter per cubic-meter) based on a Pedo-Transfer Function developed using the Africa Soil Profile Database (\enc{Wösten}{Wosten} et al. 2013).}
\usage{AWCPTF(SNDPPT, SLTPPT, CLYPPT, ORCDRC, 
    BLD=1682, CEC, PHIHOX, h1=-10, h2=-20, h3=-31.6, 
    pwp=-1585, PTF.coef, fix.values=TRUE )
}
\arguments{
  \item{SNDPPT}{numeric; sand content in percent}
  \item{SLTPPT}{numeric; silt content in percent}
  \item{CLYPPT}{numeric; clay content in percent}
  \item{ORCDRC}{numeric; soil organic carbon concentration in permille or g / kg}
  \item{BLD}{numeric; bulk density in kg / cubic-meter for the horizon/solum}
  \item{CEC}{numeric; Cation Exchange Capacity in cmol per kilogram}
  \item{PHIHOX}{numeric; soil pH in water suspension}  
  \item{h1}{numeric; moisture potential in kPa e.g. -10 (pF 2.0)}
  \item{h2}{numeric; moisture potential in kPa e.g. -20 (pF 2.3)}
  \item{h3}{numeric; moisture potential in kPa e.g. -31.6 (pF 2.5)}
  \item{pwp}{numeric; moisture potential at wilting point in kPa e.g. -1585 (pF 4.2)}
  \item{PTF.coef}{data.frame; optional conversion coefficients (Pedo-Transfer Function)}
  \item{fix.values}{logical; specifies whether to correct values of textures and bulk density to avoid creating non-sensical values}
}
\value{
Returns a data frame with the following columns:
 \itemize{
   \item \code{AWCh1}: available soil water capacity (volumetric fraction) for h1; 
   \item \code{AWCh2}: available soil water capacity (volumetric fraction) for h2;
   \item \code{AWCh3}: available soil water capacity (volumetric fraction) for h3;     
   \item \code{WWP}: available soil water capacity (volumetric fraction) until wilting point;   
}
}
\note{
Pedotransfer coefficients (\code{PTF.coef}) developed by Hodnett and Tomasella (2002). \code{fix.values} will correct sand, silt and clay fractions so they sum up to 100, and will replace bulk density values using global minimum maximum values.
}
\author{ Johan Leenaars and Maria Ruiperez Gonzalez  }
\references{
\itemize{
  \item Hodnett, M. G., & Tomasella, J. (2002). \href{http://dx.doi.org/10.1016/S0016-7061(02)00105-2}{Marked differences between van Genuchten soil water-retention parameters for temperate and tropical soils: a new water-retention pedo-transfer functions developed for tropical soils}. Geoderma, 108(3), 155-180.
  \item \enc{Wösten}{Wosten}, J. H. M., Verzandvoort, S. J. E., Leenaars, J. G. B., Hoogland, T., & Wesseling, J. G. (2013). \href{http://dx.doi.org/10.1016/j.geoderma.2012.11.021}{Soil hydraulic information for river basin studies in semi-arid regions}. Geoderma, 195, 79-86.
}
}
\examples{
SNDPPT = 30 
SLTPPT = 25 
CLYPPT = 48 
ORCDRC = 23 
BLD = 1200 
CEC = 12 
PHIOX = 6.4
x <- AWCPTF(SNDPPT, SLTPPT, CLYPPT, ORCDRC, BLD, CEC, PHIHOX)
x
## predict AWC for the whole AfSP DB:
data(afsp)
names(afsp$horizons) 
BLDf <- ifelse(is.na(afsp$horizons$BLD), 
   mean(afsp$horizons$BLD, na.rm=TRUE), afsp$horizons$BLD)
xd <- AWCPTF(afsp$horizons$SNDPPT, afsp$horizons$SLTPPT, 
  afsp$horizons$CLYPPT, afsp$horizons$ORCDRC, 
  BLD=BLDf*1000, afsp$horizons$CEC, afsp$horizons$PHIHOX)
str(xd)
}