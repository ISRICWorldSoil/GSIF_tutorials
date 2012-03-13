\name{WPSworldgrids}
\alias{WPSworldgrids}
\title{Web Processing Service for worldgrids.org}
\description{Overlay, subset and aggregation (wrapper) functions for the WPS \href{http://worldgrids.org/doku.php?id=wiki:functions.r}{worldgrids.org}.}
\usage{
WPSworldgrids(layername, \dots)
}
\arguments{
  \item{layername}{object}
  \item{\dots}{additional arguments}
}
\details{WPS is implemented via the \href{http://pywps.wald.intevation.org}{PyWPS} module. The standard format for the images on the worldgrids.org repository is \code{"GeoTiff"}. Complete grids (global coverage) can be obtained directly from worldgrids.org. We do not advise using the "bbox" argument to obtain grids that cover more than 30 percent of the global coverage.}
\references{
\itemize{
\item text
\item text
}
}
\author{ Tomislav Hengl \& Hannes I. Reuter}
\seealso{ \code{\link{landmask}}
}
\examples{
\dontrun{
WPSworldgrids(layername="RNMGSP1", xcoord=10, ycoord=10)
WPSworldgrids(layername="TD1MOD1", bbox=matrix(c(-20,30,40,60), nrow=2))
# 1 km resolution:
WPSworldgrids(layername="LN1DSM3", bbox=matrix(c(7,51,9,54), nrow=2))
}
}
\keyword{text}
