\name{geosamples-class}
\docType{class}
\alias{geosamples-class}
\alias{subset,geosamples-method}
\alias{subset.geosamples}
\alias{show,geosamples-method}
\alias{overlay,geosamples-method}
\alias{write.data}
\alias{write.data,geosamples-method}
\alias{write.data.geosamples}
\alias{write.data,SpatialPointsDataFrame-method}
\alias{write.data.SpatialPoints}
\title{A class for spatially and temporally referenced samples}
\description{A class for spatially and temporally referenced samples with fixed column names (standardized geosamples). Corresponds to the point \code{"Placemark"} in the KML schema.}
\section{Slots}{
  \describe{
    \item{\code{registry}:}{object of class \code{"character"}; URI of the online registry i.e. the URL where the \code{"producerid"} column can be linked to all other connected metadata}
    \item{\code{methods}:}{object of class \code{"data.frame"}; a table with method names (\code{"methodid"}), a one sentence description of each method (\code{"description"}), measurement units or levels (\code{"units"}), and associated detection limits (\code{"detectionLimit"})}
    \item{\code{data}:}{object of class \code{"data.frame"}; a standardized table with fixed column names: \code{"observationid"} (unique observation ID; as specified in the data registry service), \code{"sampleid"} (producer's ID; usually site ID and horizon ID or sequence number), \code{"longitude"} (longitude on the WGS84 ellipsoid), \code{"latitude"} (latitude on the WGS84 ellipsoid), \code{"locationError"} (error radius in meters), \code{"TimeSpan.begin"}, \code{"TimeSpan.end"}, \code{"altitude"}, \code{"altitudeMode"}, \code{"sampleArea"} (spatial support in square meters), \code{"sampleThickness"} (thickness of horizons in meters or vertical support), \code{"observedValue"} (measured value), \code{"methodid"} (method name), \code{"measurementError"} (estimated measurement error for that specific observation)}
  }
The column names in the data slot largely reflect the \href{https://developers.google.com/kml/documentation/kmlreference}{KML schema elements}. It also aims at being compatible with the \href{http://www.opengeospatial.org/standards/om/}{OGC Observations and measurements specifications}.\cr 
Geosamples-class can be used to store and manipulate geological, hydrological, geochemical, biodiversity, soil science and similar field samples. Geological and soil samples can be registered via the \href{http://www.geosamples.org/}{geosamples.org}, in which case the \code{"observationid"} will correspond to the unique sample identifier (IGSN). \code{"sampleid"} column allows linking geosamples to the original ID's. 
}
\section{Methods}{
  \describe{
	\item{show}{\code{signature(obj = "geosamples")}: summarize object by listing methods, total number of observations, total area covered etc. }
	\item{subset}{\code{signature(obj = "geosamples")}: subset to a single variable type; returns a data frame }
  \item{overlay}{\code{signature(covariates = "SpatialPixelsDataFrame", observations = "geosamples")}: overlay geosamples and spatial pixels }
  \item{write.data}{\code{signature(obj = "geosamples")}: write geosamples to an external format e.g. GeoEAS }
  }
}
\author{ Tomislav Hengl and Hannes I. Reuter }
\references{
\itemize{
\item Dyson, E., (2003) \href{http://cdn.oreilly.com/radar/r1/09-03.pdf}{Online Registries: The DNS and Beyond\dots}. Edventure, Vol 21(8).
\item International Geo Sample Number (\url{http://en.wikipedia.org/wiki/International_Geo_Sample_Number})
\item KML Reference (\url{https://developers.google.com/kml/documentation/kmlreference})
\item OGC Observations and Measurements standard (\url{http://www.opengeospatial.org/standards/om/})
\item SESAR, the System for Earth Sample Registration (\url{http://www.geosamples.org})
}
}
\seealso{ \code{\link{as.geosamples}} }
\keyword{classes}
