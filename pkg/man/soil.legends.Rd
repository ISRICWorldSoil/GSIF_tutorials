\name{soil.legends}
\docType{data}
\encoding{latin1}
\alias{soil.legends}
\alias{soil.dom}
\alias{soil.vars}
\alias{makeSAGAlegend}
\title{Standard color palettes for soil properties and classes}
\description{Standard color palettes for soil properties and classes that can be used to display global soil data.}
\usage{data(soil.legends)}
\format{
Contains a list of color palettes (data frames with class names / break points, and cumulative probabilities) for:
\describe{
  \item{\code{ORCDRC}}{numeric; soil organic carbon content in permille}
  \item{\code{PHIHOX}}{numeric; pH index measured in water solution}
  \item{\code{PHIKCL}}{numeric; pH index measured in KCl solution}
  \item{\code{BLD}}{numeric; bulk density in kg per cubic meter}
  \item{\code{CEC}}{numeric; Cation Exchange Capacity}
  \item{\code{SNDPPT}}{numeric; weight percentage of the sand particles (0.05--2 mm)}
  \item{\code{SLTPPT}}{numeric; weight percentage of the silt particles (0.0002--0.05 mm)}
  \item{\code{CLYPPT}}{numeric; weight percentage of the clay particles (<0.0002 mm)}
  \item{\code{CRFVOL}}{numeric; volumetric percentage of coarse fragments (>2 mm)}
  \item{\code{TAXGWRB}}{factor; World Reference base groups}
  \item{\code{TAXOUSDA}}{factor; Keys to Soil Taxonomy suborders}
}
}
\note{Breaks for continuous soil properties were determined using the \code{quantiles} function and by visually inspecting the histograms to maximize the contrast in output maps. Based on a compilation of global soil profile data (\url{http://soilprofiles.org}).}
\author{Tomislav Hengl}
\references{
\itemize{
\item Global Soil Information Facilities (\url{http://gsif.isric.org})
}
}
\examples{
data(soil.legends)
pal <- soil.legends$ORCDRC$COLOR
names(pal) <- signif((soil.legends$ORCDRC$MAX + 
    soil.legends$ORCDRC$MIN)/2, 3)
pal
data(soil.vars)
soil.vars[soil.vars$varname=="ORCDRC",]
## make SAGA GIS palette:
makeSAGAlegend(x=as.factor(names(pal)), col_pal=pal, 
    filename="ORCDRC.txt")
}
\keyword{datasets}