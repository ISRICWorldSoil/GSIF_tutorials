\name{GlobalSoilMap-class}
\docType{class}
\encoding{latin1}
\alias{GlobalSoilMap-class}
\alias{GlobalSoilMap}
\title{A class for GlobalSoilMap soil property maps}
\description{A class containing predictions of target soil property at six standard depths following the \href{http://globalsoilmap.net/specifications}{GlobalSoilMap.net specifications}: \code{sd1} = 2.5 cm (0--5), \code{sd2} = 10 cm (5--15), \code{sd3} = 22.5 cm (15--30), \code{sd4} = 45 cm (30--60), \code{sd5} = 80 cm (60--100), \code{sd6} = 150 cm (100--200).}
\section{Slots}{
  \describe{    
    \item{\code{varname}:}{object of class \code{"character"}; abbreviated variable name registered in the Global Soil Data registry}
    \item{\code{TimeSpan}:}{object of class \code{"list"}; contains begin and end of the sampling period of class \code{"POSIXct"}}
    \item{\code{sd1}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 2.5 cm (0--5)}
    \item{\code{sd2}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 10 cm (5--15)}
    \item{\code{sd3}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 22.5 cm (15--30)}
    \item{\code{sd4}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 45 cm (30--60)}
    \item{\code{sd5}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 80 cm (60--100)}
    \item{\code{sd6}:}{object of class \code{"SpatialPixelsDataFrame"}; predictions and variances, or number of realizations of the target variable at depth 150 cm (100--200)}
  }
}
\references{
\itemize{
  \item Hartemink, A. E., Hempel, J., Lagacherie, P., McBratney, A., McKenzie, N., MacMillan, R. A., ... & Zhang, G. L. (2010). \href{http://dx.doi.org/10.1007/978-90-481-8863-5_33}{GlobalSoilMap.net --- A New Digital Soil Map of the World}. In Digital Soil Mapping (pp. 423--428). Springer Netherlands.
  \item Sanchez, P. A., S. Ahamed, F. Carre, A. E. Hartemink, J. Hempel, J. Huising, P. Lagacherie, A. B. McBratney, N. J. McKenzie, M L. de \enc{Mendonça}{Mendonca}-Santos, et al., (2009) \href{http://dx.doi.org/10.1126/science.1175084}{Digital Soil Map of the World}. Science, 325(5941): 680--681.
}
}
\seealso{ \code{\link{SoilGrids-class}}, \code{\link{SpatialComponents-class}}, \code{\link{geosamples-class}} 
}
\keyword{classes}