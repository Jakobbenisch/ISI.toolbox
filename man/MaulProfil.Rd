% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Manning.R
\name{MaulProfil}
\alias{MaulProfil}
\title{A Geometry Defining Function}
\usage{
MaulProfil(h, r)
}
\arguments{
\item{h}{Decimal [m]. Representing the y coordinate to which the width will be returned.}

\item{r}{Decimal [m]. Radius of the Maulprofil.}
}
\value{
Returns width of profile corresponding to input height [m].
}
\description{
This function analytically defines the relation between height and width of a standard Maulprofil.\cr
Mind the units!\cr
If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
}
\examples{
width=MaulProfil(.5,1)

}
\seealso{
\link[ISI.Toolbox]{ManningMaul}
}
\keyword{Manning}
