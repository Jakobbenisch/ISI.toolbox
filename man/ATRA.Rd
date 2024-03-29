% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Manning.R
\name{ATRA}
\alias{ATRA}
\title{An Area Defining Function (Trapezoid Profile)}
\usage{
ATRA(w, b, alpha, beta)
}
\arguments{
\item{w}{Decimal or numeric vector [m].}

\item{b}{Decimal [m]. Width of bottom part of the trapezoid profile.}

\item{alpha}{Decimal in rad. Angle of slope -> right (mathematically negative). Alpha=0 results in a vertical line.}

\item{beta}{Decimal in rad. Angle of slope -> left (mathematically positive). Beta=0 results in a vertical line.}
}
\value{
Returns the area of flow cross section [m^2].
}
\description{
This function analytically calculates the area of a flow cross section for a water level.\cr
Mind the units!!\cr
If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)\cr
Also the angle should be between 0 and pi/2. Excluding pi/2!
}
\examples{
A=ATRA(w=1,b=1,alpha=pi/4,beta=pi/5)

}
\seealso{
\link[ISI.Toolbox]{ManningTRA}
}
\keyword{Manning}
