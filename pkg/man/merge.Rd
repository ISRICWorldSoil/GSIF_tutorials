\name{merge}
\docType{methods}
\alias{merge}
\alias{merge,SpatialPredictions,SpatialPredictions-method}
\alias{merge,RasterBrickSimulations,RasterBrickSimulations-method}
\title{Merge multiple predictions}
\description{Merges objects of class \code{"SpatialPredictions"} or \code{"RasterBrickSimulations"} and produces average predictions where the two objects overlap spatially. If the predictions are available at different resolutions, then it downscales all other grids to the smallest grid cell size using bicubic splines (for predictions) i.e. nearest neighbor algorithm (for simulations). Weigths can be passed via the \code{RMSE.l} argument, otherwise they will be estimated from \code{validation} slot (in the case of the \code{"SpatialPredictions"} class objects).
}
\usage{
\S4method{merge}{SpatialPredictions,SpatialPredictions}(x, y, \dots, 
          RMSE.l = NULL, silent = TRUE)
}
\arguments{
  \item{x}{object of class \code{"SpatialPredictions"} or \code{"RasterBrickSimulations"}}
  \item{y}{object of class \code{"SpatialPredictions"} or \code{"RasterBrickSimulations"}}
  \item{...}{additional objects of class \code{"SpatialPredictions"} or \code{"RasterBrickSimulations"}}
  \item{RMSE.l}{numeric; list of mean prediction errors for each object (these are used as weights during the averaging)}
  \item{silent}{logical; specifies whether to print out the progress and used RMSE's}
} 
\details{Merging of multiple predictions using weighted averaging is a heuristic approach to mapping. The main problem is that this method assumes that the predictions are completely independent (independent covariates, independent models), but this not might be the case and hence the predictions will be sub-optimal. Merging multiple predictions is however attractive for situations where the predictions do not have the same extent, so that a more global map can be used to fill in the gaps in locally produced predictions.}
\value{
Returns an object of type \code{"SpatialPixelsDataFrame"} or \code{"RasterBrickSimulations"} that contains only the merged values.
}
\references{
\itemize{ 
\item Heuvelink, G.B.M., Bierkens, M.F.P. (1992) \href{http://dx.doi.org/10.1016/0016-7061(92)90002-O}{Combining soil maps with interpolations from point observations to predict quantitative soil properties}. Geoderma 55(1-2): 1-15.
}
}
\author{Tomislav Hengl and Gerard B.M. Heuvelink}
\keyword{methods}