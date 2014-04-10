\name{OCSKGM}
\alias{OCSKGM}
\title{Soil organic carbon stock calculus}
\description{Compute soil organic carbon stock (in kilograms per square-meter) for a given horizon/solum depth and based on soil organic carbon concentration, horizon/solum thickness, bulk density and percentage of coarse fragments.}
\usage{OCSKG(ORCDRC, BLD=1682, CRFVOL=0, HSIZE, 
     ORCDRC.sd=10, BLD.sd=100, CRFVOL.sd=5, se.prop=TRUE)
}
\arguments{
  \item{ORCDRC}{numeric; soil organic carbon concentration in permille or g / kg}
  \item{BLD}{numeric; bulk density in kg / cubic-meter}
  \item{CRFVOL}{numeric; percentage of coarse fragments in the horizon/solum}
  \item{HSIZE}{numeric; thickness of the horizon/solum in cm}
  \item{ORCDRC.sd}{numeric; standard error of estimating \code{ORCDRC} must be positive number}
  \item{BLD.sd}{numeric; standard error of estimating \code{BLD} must be positive number}
  \item{BLD.sd}{numeric; standard error of estimating \code{BLD} must be positive number}
}
\value{
\href{http://www.eoearth.org/view/article/156087}{Soil organic carbon stock} in kilograms per square-meter. To convert to tones per hectar multiply by 10.
}
\note{
Propagated error (attached as an attribute) is estimated using the Taylor Series Method and shows only an approximate estimate. A more robust way to estimate the propagated uncertainty would be to use (geo)statistical simulations. See Heuvelink (1998) for more info.
}
\author{ Tomislav Hengl, Niels Batjes and Gerard Heuvelink }
\references{
\itemize{
  \item Heuvelink, G. B. (1998) Error propagation in environmental modelling with GIS. CRC Press, 150 p.
  \item Nelson, D.W., and L.E. Sommers (1982) Total carbon, organic carbon, and organic matter. p. 539–580. In A.L. Page et al. (ed.) Methods of soil Analysis. Part 2. 2nd ed. Agron. Monogr. 9. ASA and SSSA, Madison, WI.
}
}
\examples{
Area <- 1E4  ## 1 ha
HSIZE <- .3 ## 0--30 cm
ORCDRC <- 50  ## 5%
ORCDRC.sd <- 10  ## +/-1%
BLD <- 1500  ## 1.5 tones per cubic meter
BLD.sd <- 100  ## +/-0.1 tones per cubic meter
CRFVOL <- 10  ## 10%
CRFVOL.sd <- 5  ## +/-5%
x <- OCSKGM(ORCDRC, BLD, CRFVOL, HSIZE, ORCDRC.sd, BLD.sd, CRFVOL.sd)
x  ## 20.25 +/-4.41 kg/m^2
## in tones per ha:
x[[1]] * Area / 1000
}