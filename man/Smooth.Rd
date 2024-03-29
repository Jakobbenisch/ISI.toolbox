% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Filters.R
\name{Smooth}
\alias{Smooth}
\title{Function for smoothing an xts vector}
\usage{
Smooth(XTS, window_width = 7)
}
\arguments{
\item{XTS}{xts vector to be smoothened.}

\item{window_width}{The width of the moving window. (Ideally this is a uneven number.)}
}
\value{
Smoothed xts vector.
}
\description{
The function uses a moving window. For each moving window a mean value is calculated and substitutes the unsmoothed one.
}
\examples{
x=DummyTS(Csd=1.3)
y=Smooth(x,window_width=9)
DynPlot(cbind(x,y))
}
\seealso{
\link[ISI.Toolbox]{OneEuroFilter}
}
\keyword{Filter}
\keyword{smooth,}
