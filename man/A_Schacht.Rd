% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Manning.R
\name{A_Schacht}
\alias{A_Schacht}
\title{An Area Defining Function (Manhole)}
\usage{
A_Schacht(h, r = 0.6, b = 1.47, hlimit = 0.21)
}
\arguments{
\item{h}{Decimal[m]. Water levels.}

\item{r}{Decimal [m]. Radius of circular bottom part of manhole.}

\item{b}{Decimal [m]. Width of manhole.}

\item{hlimit}{Decimal [m]. Height of transition between circular shaped bottom and straight wall.}
}
\value{
Returns the area of flow cross section [m^2].
}
\description{
This function analytically calculates the area of a flow cross section for a water level.\cr
Mind the units!!\cr
If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)\cr
}
\examples{
A=A_Schacht(h=1,r=0.6,b=1.47,hlimit=0.21)

}
\seealso{
\link[ISI.Toolbox]{ManningSchacht}
}
\keyword{Manning}
