\name{makeGstatCmd}
\alias{makeGstatCmd}
\encoding{latin1}
\title{Make a gstat command script}
\description{Generates a command script based on the regression model and variogram. This can then be used to run predictions/simulations by using the pre-compiled binary \code{gstat.exe}.}
\usage{
makeGstatCmd(formString, vgmModel, outfile, easfile, 
        nsim = 0, nmin = 20, nmax = 40, radius, zmap = 0, 
        predictions = "var1.pred.hdr", variances = "var1.svar.hdr", 
        xcol = 1, ycol = 2, zcol = 3, vcol = 4, Xcols)
}
\arguments{
  \item{formString}{object of class \code{"formula"} --- regression model} 
  \item{vgmModel}{object of class \code{"vgmmodel"} or \code{"data.frame"}} 
  \item{outfile}{character; output file for the command script}
  \item{easfile}{character; file name for the GeoEAS file with observed values}
  \item{nsim}{integer; number of simulations}
  \item{nmin}{integer; smallest number of points in the search radius (see gstat user's manual)}
  \item{nmax}{integer; largest number of points in the search radius (see gstat user's manual)}
  \item{radius}{numeric; search radius (see gstat user's manual)}
  \item{zmap}{numeric; fixed value for the 3D dimension in the case of 3D kriging}
  \item{predictions}{character; output file name for predictions}
  \item{variances}{character; output file name for kriging variances}
  \item{xcol}{integer; position of the x column in the GeoEAS file}
  \item{ycol}{integer; position of the y column in the GeoEAS file}
  \item{zcol}{integer; position of the z column in the GeoEAS file}
  \item{vcol}{integer; position of the target variable column in the GeoEAS file}
  \item{Xcols}{integer; column numbers for the list of covariates}
}
\details{To run the script under Windows OS you need to obtain the pre-compiled \code{gstat.exe} program from the www.gstat.org website, and put it in some directory e.g. \code{c:/gstat/}. Then add the program to your path (see environmental variable under Windows > Control panel > System > Advanced > Environmental variables), or copy the exe program directly to some windows system directory.}
\note{The advantage of using \code{gstat.exe} is that it loads large grids much faster to memory than if you use gstat in R, hence it is potentially more suited for computing with large grids. The draw back is that you can only pass simple linear regression models to \code{gstat.exe}. The stand-alone gstat is not maintained by the author of gstat any more.}
\author{ Tomislav Hengl }
\references{
\itemize{
\item Bivand, R.S., Pebesma, E.J., and \enc{Gómez}{Gomez}-Rubio, V., (2008) \href{http://www.asdar-book.org/}{Applied Spatial Data Analysis with R}. Springer, 378 p. 
\item Pebesma, E., (2003) \href{http://www.gstat.org/gstat.pdf}{Gstat user's manual}. Dept. of Physical Geography, Utrecht University, p. 100, www.gstat.org 
}
}
\seealso{ \code{\link{write.data}}, \code{\link{fit.gstatModel}}, \code{gstat::krige} }
\examples{
\dontrun{
library(sp)
library(gstat)

# Meuse data:
demo(meuse, echo=FALSE)
# fit a model:
omm <- fit.gstatModel(observations = meuse, formulaString = om~dist, 
  family = gaussian(log), covariates = meuse.grid)
str(omm@vgmModel)
# write the regression matrix to GeoEAS:
meuse$log_om <- log1p(meuse$om)
write.data(obj=meuse, covariates=meuse.grid["dist"], 
    outfile="meuse.eas", methodid="log_om")
writeGDAL(meuse.grid["dist"], "dist.rst", drivername="RST", mvFlag="-99999")
makeGstatCmd(log_om~dist, vgmModel=omm@vgmModel, 
    outfile="meuse_om_sims.cmd", easfile="meuse.eas", 
    nsim=50, nmin=20, nmax=40, radius=1500)
# compare the processing times:
system.time(system("gstat meuse_om_sims.cmd"))
vgmModel = omm@vgmModel
class(vgmModel) <- c("variogramModel", "data.frame")
system.time(om.rk <- krige(log_om~dist, meuse[!is.na(meuse$log_om),], 
    meuse.grid, nmin=20, nmax=40, model=vgmModel, nsim=50))
}
}
\keyword{methods}
