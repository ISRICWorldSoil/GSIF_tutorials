\name{as.geosamples}
\docType{methods}
\alias{as.geosamples}
\alias{as.geosamples,SpatialProfileCollection-method}
\title{Converts an object of class \code{"SpatialProfileCollection"} to a standardized data frame}
\description{Converts an object of class \code{"SpatialProfileCollection"} to an object of class \code{"geosamples"} with all measurements broken separately as individual records. Geo-records are standardized universal spatially and temporally referenced samples from the Earth's surface.}
\usage{
as.geosamples(obj)
}
\arguments{
  \item{obj}{object of class \code{"SoilProfileCollection"}}
}
\details{Returns an object of type \code{"geosamples"} with the following slots:
  \describe{
   \item{\code{registry}}{character; an unique URI defining the global registry with all metadata and attribute lists (domains)}
   \item{\code{methods}}{data.frame; a table with method names (\code{"methodid"}), a detailed description of each method (\code{"description"}), measurement units or levels (\code{"units"}), and associated detection limits (\code{"detectionLimit"})}
   \item{\code{data}}{data.frame; a standardized table with fixed column names: \code{"sampleid"} (universal ID; as specified in the data registry service), \code{"producerid"} (producer's ID; usually site ID and horizon ID or sequence number), \code{"longitude"}, \code{"latitude"}, \code{"locationError"}, \code{"TimeSpan.begin"}, \code{"TimeSpan.end"}, \code{"altitude"}, \code{"altitudeError"}, \code{"altitudeMode"}, \code{"dimension"}, \code{"volume"}, \code{"value"} (measured value), \code{"methodid"}, \code{"measurementError"}}
  }
The column names in the data slot largely reflect the KML schema elements. For more information about what each column means, refer to the \href{https://developers.google.com/kml/documentation/kmlreference}{KML reference}.\cr 
Many columns required by the \code{"geosamples"} class can be missing and will result in \code{NA} values. These will likely need to be inserted afterwards, within or outside R.}
\note{The basic idea of this class is to outsource the metadata management to some global publicly available service. Georecord class can be used to prepare geological, hydrological, geochemical, biodiversity, soil science and similar samples. Geological and soil samples can be registered via the \href{http://www.geosamples.org/}{geosamples.org}, for example, so that the \code{"sampleid"} corresponds to the IGSN (unique sample identifier). The methods used in the field of soil science could be all implemented via the \href{http://soilinformationstandards.org}{SoilLM} scheme. \code{"producerid"} column is recommended as it allows linking to the original data set.}
\author{ Tomislav Hengl and Hannes I. Reuter }

\seealso{ \code{\link{as.data.frame}}, \code{\link{geosamples-class}} }
\examples{
\dontrun{

}
}
\keyword{methods}
