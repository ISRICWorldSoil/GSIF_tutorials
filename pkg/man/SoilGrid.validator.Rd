\name{SoilGrid.validator}
\alias{SoilGrid.validator}
\title{Validate SoilGrid (spatial predictions)}
\description{Validate SoilGrid (spatial predictions) i.e. soil property maps following the GSIF validation protocol.}
\usage{
SoilGrid.validator(obj, domain, ground.truth, N.sample=2000, 
   xml.file, z.lim, md.type="INSPIRE", test.URL=FALSE)
}
\arguments{
  \item{obj}{\code{"GDALobj"} object i.e. a pointer to a spatial layer of interest (single slice)}
  \item{domain}{\code{"GDALobj"} object i.e. a pointer to a spatial layer contain soil mask}
  \item{ground.truth}{\code{"SpatialPointsDataFrame"}; contains values of the target variable at exactly the same depth / same support size sampled either using Simple Random Sampling or regular sampling on a grid}  
  \item{N.sample}{integer; random sampling size}  
  \item{xml.file}{character; metadata file (should have the same name as \code{obj} file)}    
  \item{z.lim}{numeric; upper and lower physical limits}
  \item{md.type}{character; metadata standard (currently \href{http://inspire.ec.europa.eu/}{INSPIRE})}    
  \item{test.URL}{logical; specifies whether to validate XML schema / test download times and proj4 string}
}
\value{Returns a data frame with validation results. Explanation of codes is available in the SoilGrids.org data validation protocol.}
\note{One SoilGrid layer (2D slice) basically contains predictions on a regular grid for a specific soil depth (at either point or block support). The \code{ground.truth} data must refer to the exactly the same depth and the same support size and should ideally be collected using some probability spatial sampling (see e.g. \code{sp::spsample}). To estimate values of soil properties at standard depths, consider using \code{\link{mpspline}} function.\cr
Numeric resolution is derived as estimated RMSE/2. Numeric resolution can be best specified as \code{Attribute_Measurement_Resolution} (the smallest unit increment to which an attribute value is measured). \cr
Increasing \code{N.sample} can lead to more precise results at the cost of higher computing time.
}
\author{Tomislav Hengl }
\seealso{\code{plotKML::spMetadata} }