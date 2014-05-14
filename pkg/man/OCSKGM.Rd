\name{OCSKGM}
\alias{OCSKGM}
\title{Soil organic carbon stock}
\description{Derive soil organic carbon stock / storage (in kilograms per square-meter) and propagated uncertainty for a given horizon/solum depth and based on soil organic carbon concentration, horizon/solum thickness, bulk density and percentage of coarse fragments.}
\usage{OCSKGM(ORCDRC, BLD=1682, CRFVOL=0, HSIZE, 
     ORCDRC.sd=10, BLD.sd=100, CRFVOL.sd=5, se.prop=TRUE)
}
\arguments{
  \item{ORCDRC}{numeric; soil organic carbon concentration in permille or g / kg}
  \item{BLD}{numeric; bulk density in kg / cubic-meter for the horizon/solum}
  \item{CRFVOL}{numeric; percentage of coarse fragments (above 2 mm in diameter) in the sample}
  \item{HSIZE}{numeric; thickness of the horizon/solum in cm}
  \item{ORCDRC.sd}{numeric; standard error of estimating \code{ORCDRC} (must be positive number)}
  \item{BLD.sd}{numeric; standard error of estimating \code{BLD} (must be positive number}
  \item{CRFVOL.sd}{numeric; standard error of estimating \code{CRFVOL} (must be positive number)}
  \item{se.prop}{logical; specifies whether to derive propagated error}
}
\value{
\href{http://www.eoearth.org/view/article/156087}{Soil organic carbon stock} in kilograms per square-meter. To convert to tonnes per hectar multiply by 10.
}
\note{
Propagated error (attached as an attribute) is estimated using the Taylor Series Method and shows only an approximate estimate. A more robust way to estimate the propagated uncertainty would be to use (geo)statistical simulations. See Heuvelink (1998) for more info.
}
\author{ Tomislav Hengl, Niels Batjes and Gerard Heuvelink }
\references{
\itemize{
  \item Heuvelink, G. B. (1998) Error propagation in environmental modelling with GIS. CRC Press, 150 p.
  \item Nelson, D.W., and L.E. Sommers (1982) Total carbon, organic carbon, and organic matter. p. 539-580. In A.L. Page et al. (ed.) Methods of soil Analysis. Part 2. 2nd ed. Agron. Monogr. 9. ASA and SSSA, Madison, WI.
}
}
\examples{
Area <- 1E4  ## 1 ha
HSIZE <- 30 ## 0--30 cm
ORCDRC <- 50  ## 5%
ORCDRC.sd <- 10  ## +/-1%
BLD <- 1500  ## 1.5 tonnes per cubic meter
BLD.sd <- 100  ## +/-0.1 tonnes per cubic meter
CRFVOL <- 10  ## 10%
CRFVOL.sd <- 5  ## +/-5%         
x <- OCSKGM(ORCDRC, BLD, CRFVOL, HSIZE, ORCDRC.sd, BLD.sd, CRFVOL.sd)
x  ## 20.25 +/-4.41 kg/m^2
## in tonnes per ha:
x[[1]] * Area / 1000
}