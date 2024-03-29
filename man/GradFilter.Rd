% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Filters.R
\name{GradFilter}
\alias{GradFilter}
\title{Function for deleting outliers in an xts vector}
\usage{
GradFilter(TS, abs_grad = 9999999, acceptableError = 0.01,
  UnplausibleGapsize = NULL, Analyse = FALSE, GradZero = FALSE,
  Automatic = FALSE, OnlyDoubleGradientError = FALSE, Prob = 0.98,
  NormalizedGrad = FALSE, windowlengthNG = 4)
}
\arguments{
\item{TS}{Xts vector to be freed of outliers.}

\item{abs_grad}{The absolute value of the gradient limit.}

\item{acceptableError}{The allowed number of Gradient errors in percent. (Number of Errors / Number of Observations *100)}

\item{UnplausibleGapsize}{This parameter allows to filter noisy signals to some extend. If set to a number >=1 the filter will also delete oberservations between two outlier if the distance is smaller or equal to that value.}

\item{Analyse}{Boolean. If TRUE a summary of the gradients distribution and position will be plotted.}

\item{GradZero}{Boolean. If TRUE values that result in gradients of 0 will be deleted.}

\item{Automatic}{Boolean. If TRUE the abs_grad is set automatically using R's quantile function with a default probability of 0.98.}

\item{OnlyDoubleGradientError}{Boolean. If TRUE the filter will only delete values that result in both a negative and positive gradient error.}

\item{Prob}{The probability of for the Automatic modus.}

\item{NormalizedGrad}{Boolean. If True each gradient will take the its neighbours into account to decide, whether it is a oultier or a valid value.\cr
Each transformed gradient G_i will be computed as follows: \cr
G_{i}=g_{i}/mean(g_(i-l/2)+g_(i-l/2+1)+...g_(i+l/2-1)+g_(i+l/2-1))\cr
g_{i} is the ith gradient\cr
l is the number of neighbour gradients to be considered\cr
G_{i} is the transformed gradient and will be used for filtering}

\item{windowlengthNG}{Integer. Specifies the number of neighbour gradients to be considered for the NormalizedGrad method.This value must be even.}
}
\value{
xts vector without outliers.
}
\description{
The function filters gradients greater than a set absolute value of an xts vector. The filter repeats its operation until the acceptableError is reached.
Use the Analyse option (set value to TRUE) to get some statistics about the distribution of gradients in your time series.
Once you know how to set your gradient limits turn the Analyse function to FALSE and rerun the Gradfilter. Limit exceeding values will then be deleted.
}
\examples{
A = DummyTS(days=50)
B = Outliers(A,5,4)
GradFilter(B,Analyse = T,NormalizedGrad = TRUE)
Bf=GradFilter(B,abs_grad = 3,NormalizedGrad = TRUE)
DynPlot(cbind(B,Bf),Labels=c("TS_with_Outliers","TS_without_Outliers"))

}
\seealso{
\link[ISI.Toolbox]{NoiseFilter}
}
\keyword{Filter}
\keyword{Outlier,}
