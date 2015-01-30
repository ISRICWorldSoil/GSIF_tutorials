\name{ERDICM}
\alias{ERDICM}
\title{Effective Rooting Zone depth}
\description{Derive Effective Rooting Zone depth i.e. an effective depth suitable for plant growth. Usually minimum depth of soil out of three standard rooting depths: limiting soil properties, depth to water-stagnating layer and depth to bedrock.}
\usage{ERDICM(UHDICM, LHDICM, minimum.LRI, DRAINFAO, BDRICM, 
    threshold.LRI=20, srd=150, drain.depths, smooth.LRI=TRUE)
}
\arguments{
  \item{UHDICM}{numeric; upper horizon depth in cm}
  \item{LHDICM}{numeric; lower horizon depth in cm}
  \item{minimum.LRI}{numeric; minimum Limiting Rootability index}
  \item{DRAINFAO}{factor; FAO drainage class e.g. \code{"V"}, \code{"P"}, \code{"I"}, \code{"M"}, \code{"W"}, \code{"S"}, \code{"E"}}
  \item{BDRICM}{numeric; depth to bedrock in cm}
  \item{threshold.LRI}{numeric; treshold index for LRI}
  \item{srd}{numeric; maximum depth of interest}
  \item{drain.depths}{data.frame; estimate effective rooting depth per drainage class (\code{DRAINFAO})}
  \item{smooth.LRI}{logical; specify whether to smooth LRI values using splines}
}
\value{
Returns a vector of effective rooting depth in cm.
}
\author{ Johan Leenaars and Maria Ruiperez Gonzalez  }
\seealso{ \code{\link{LRI}} }