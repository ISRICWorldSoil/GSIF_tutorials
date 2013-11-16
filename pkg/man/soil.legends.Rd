\name{soil.legends}
\docType{data}
\encoding{latin1}
\alias{soil.legends}
\title{Standard color palettes for soil properties and classes}
\description{Standard color palettes for soil properties and classes that can be used to display global soil data.}
\usage{data(soil.legends)}
\format{
Contains a list of color palettes (data frames with class names / break points, and cumulative probabilities) for:
\describe{
  \item{\code{ORCDRC}}{numeric; soil organic carbon content in permille}
  \item{\code{PHIHO5}}{numeric; pH index measured in water solution}
  \item{\code{PHIHO5}}{numeric; pH index measured in KCl solution}
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
\author{Global Soil Information Facilities (\url{http://gsif.isric.org})}
\examples{
data(soil.legends)
pal <- soil.legends$ORCDRC$COLOR
names(pal) <- (soil.legends$ORCDRC$MAX + soil.legends$ORCDRC$MIN)/2
leg.width <- (max(nchar(names(pal)))*20+150)/100
par(mar=c(.0,.0,.0,.5))
plot(x=rep(1, length(pal)), y=1:length(pal), axes=FALSE, 
  xlab='', ylab='', pch=15, cex=1.5, col=pal, 
  xlim=c(0,.6*leg.width), asp=.6)
text(x=rep(1, length(pal)), y=1:length(pal), 
  labels=names(pal), cex=.5, pos=4, offset=1)
}
\keyword{datasets}