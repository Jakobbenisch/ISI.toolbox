% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Regression.R
\name{ThinOutData}
\alias{ThinOutData}
\title{A function for thinning out data.}
\usage{
ThinOutData(data, thin.to = 300, style = c("even", "tight top",
  "tight bottom"), median_mode = c("min", "max"))
}
\arguments{
\item{thin.to}{Integer, that specifies the number of ranges that will be filled 1 data point each, if data is available.}

\item{style}{String, that specifies the distribution of the ranges. Available options are "even","tight top" and "tight bottom". Default is "even".
Even - all ranges have the same width. Tight top - lower ranges are wider, larger one are tighter. Thight bottom - the opposite of tight top.}

\item{median_mode}{String, that specifies which datapoint to take when 2 points are available. This is relevant due to the fact, that the median of a even number of samples is not computed with the mean of the 2 inner samples. Available options are "min" and "max". Default is "min".}

\item{TS}{Xts. A time series to be thinned out with dim = n,2!}
}
\value{
Data Frame with dim = n,2.
}
\description{
This function sets up an controled number of evenly or unevenly distributed spaces of a range. For each space all data points within
will be summarized using mean or median. This funtion is mainly useful for regression analysis.
}
\examples{
# generating sample data
H=DummyTS(day=50,Anual_Days_of_P =60)
h=as.numeric(H)/100
Q=40*(h)^1.4
Q1=Q
Q=Q+rnorm(length(Q),seq(0,max(Q)/10,length.out = length(Q)),median(Q)/6) # add noise
print(all(Q==Q1))
hQ=xts(cbind(h,Q),order.by=time(H))
# get rid of negative flow (produced by random sampling)
hQ=hQ[hQ[,2]>0]
h=as.numeric(hQ[,1])
Q=as.numeric(hQ[,2])

# plot original data and fit
Fitall(h,Q)
# thinn out data
hQ_thinn=ThinOutData(hQ,thin.to=50,style="tight top")
# plotting new fit
h_thinn=as.numeric(hQ_thinn[,1])
Q_thinn=as.numeric(hQ_thinn[,2])
Fitall(h_thinn,Q_thinn)
}
