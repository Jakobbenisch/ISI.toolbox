% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OtherScripts.R
\name{Rain.Intensity}
\alias{Rain.Intensity}
\title{Rain Analysing Function}
\usage{
Rain.Intensity(Rain_XTS, ignore.gap.length = 60)
}
\arguments{
\item{Rain_XTS}{Xts vector containing the rain events. Ideally no NA values are present and the time step is even. (NA could be replaced with zeros if present. To unify the the time step the Add.Na function could be used.)}

\item{ignore.gap.length}{As a rain event does not actually end once a single rainheight within the event is zero, this value can be specified to set a minimum number of consecutive values that have to be zero so that the event will have ended.
So for a 1 min time step in the rainfall xts object, the value 60 would resemble an hour.}
}
\value{
A list.\cr
[[1]] -> Numerical vector containing the intensities.\cr
[[2]] -> Numerical vector containing the total precipitation.\cr
[[3]] -> Numerical vector containing the start of the event.\cr
[[4]] -> Numerical vector containing the end of the event.\cr
\cr
Ps:\cr
The start and the end must be converted to a datetime object as follows.\cr
start_as_datetime=as.POSIXct(start_numeric,origin = "1970-01-01")
}
\description{
This function calculates the intensity of rain events.
}
\examples{
RainTimeseries=RainGen()
Intensity_Height=Rain.Intensity(RainTimeseries,ignore.gap.length=20)

}
