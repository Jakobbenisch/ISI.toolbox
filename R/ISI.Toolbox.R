#' R-Tools for the Institute of Urban Hydrology of the Technical University Dresden.
#'
#' This package contains numerous function for working with data. These functions my be categorized as follows.\cr\cr
#' \enumerate{
#'   \item Regression Tools
#'   \item Filter Tools
#'   \item Generating Tools
#'   \item Manning Strickler Tools
#'   \item Import Tools
#'   \item Plot Tools
#'   \item General Purpose Tools
#'   \item SN Transformation Algorithm
#'   }
#'
#'
#'
#' The \pkg{ISI.Toolbox} package contains these functions:\cr\cr
#' \describe{
#'   \item{\strong{Regression Tools}}{
#'   \describe{
#'   \item{\link{POLRegression}}{Fits a ploynomial function: y=a_0+a_1*x+a_2*x^2+...a_n*x^n (order=n)}
#'   \item{\link{LOGRegression}}{Fits a logarithmic function: y=a+b*ln(x).}
#'   \item{\link{EXPRegression}}{Fits a exponential function: y=a*exp(b*x).
#' Note that in order to use linear regression a transformation was used. The returned intercept value needs to be transformed back. use exp() or the function Par_Exp_Pot.}
#'   \item{\link{POTRegression}}{Fits a potency function: y=a*x^b.\cr
#' Note that in order to use linear regression a transformation was used. The returned intercept value needs to be transformed back. use exp() or the function Par_Exp_Pot.\cr
#' ( That's the regression, which is recommended for a h-Q-relation )}
#'   \item{\link{hQ}}{Function calculates a discharge time series using a water level time series and a previously calculated regressionobject (POT=POTRegression(x,y) -> POT is a regressionobject).(mind the units!)}
#'   \item{\link{hQConf}}{Function calculates 2 discharge time series (lower and upper confidence bands) using a waterlevel time series and Regressionobject.(mind the units!)}
#'   \item{\link{FourierRegression}}{Fits Fourier Series to xts vector. As the period is set to one day, this function is mainly to calculate diurnal patterns.
#' Fourier Series: y=a_0 +a_1*cos(p*1*t)+b_1*sin(p*1*t) +a_2*cos(p*2*t)+b_2*sin(p*2*t) +... +a_n*cos(p*n*t)+b_n*sin(p*n*t)\cr
#' (orderoftransformation=n, periode = p --> is fix, set to 1 day)}
#'   \item{\link{FourierPattern}}{Calculates an xts vector containing the Fourier Pattern using xts vector and Fourier regressionobject.}
#'   \item{\link{Par_Poly_Log}}{Extracts the coefficients of a regressionobject from POLRegression or LOGRegression.}
#'   \item{\link{Par_Exp_Pot}}{Extracts the coefficients of a regressionobject from POTRegression or EXPRegression.}
#'   \item{\link{ConL_Poly_Log}}{Extracts the confidence intervals of coefficients of a regressionobject from POLRegression or LOGRegression.}
#'   \item{\link{ConL_Exp_Pot}}{Extracts the confidence intervals of coefficients of a regressionobject from POTRegression or EXPRegression}
#'   \item{\link{Fitall}}{Fits all regression types except of fourier and plots the results. This function aims to give an overview to find the best fit.}
#'   \item{\link{ThinOutData}}{Reduces number of data points to minimize problems with overfitting of regression functions in areas with high data point density.}
#'   }}
#'   \item{\strong{Filter Tools}}{
#'   \describe{
#'   \item{\link{GradFilter}}{The function filters gradients greater than a set absolute value of an xts vector. The filter repeats its operation until the acceptableError is reached.
#' Use the Analyse option (set value to TRUE) to get some statistics about the distribution of gradients in your time series.
#' Once you know how to set your gradient limits turn the Analyse function to FALSE and rerun the Gradfilter. Limit exceeding values will then be deleted.}
#'   \item{\link{NoiseFilter}}{The function filters noise using std within a moving window. The Regression parameter allows to fit polynoms in each moving window and the sd is calculated of the residuals.
#' This allows for reduced influence of trends that might be in the xts vector. Use the Analyse function to get an idea about the sd distribution of the time series and thereby setting the s_max values.}
#'   \item{\link{OneEuroFilter}}{The function is for smoothing an xts or numeric vector. It's based on the publication by Géry Casiez.
#' For detailed description on how to use the filter please read the paper.}
#'   \item{\link{Smooth}}{The function uses a moving window. For each moving window a mean value is calculated and substitutes the unsmoothed one.}
#'   \item{\link{Fourierperday}}{The function uses a linear model to fit a Fourier Series for each day of an xts vector.}
#'   \item{\link{FFT_XTS_Filter}}{The function uses fft to transform the xts vector to the magnitude frequency domain. Limits across frequency and mangitudes can be set to filter certain parts of the signals.}
#'   \item{\link{Decompose4TXS}}{The function uses the standard decompose tool of the stats package and applies it to an xts object.
#' For a more detailed description please see the decompose function.}
#'   \item{\link{Separate.To.HyS.HyW}}{This function seperates a time series into hydrological winter and summer time series.}
#'   \item{\link{FixDrift}}{This function deducts linear drift and returns an xts without drift.}
#'   }}
#'   \item{\strong{Generating Tools}}{
#'   \describe{
#'   \item{\link{DummyTS}}{Generates an xts with an desired Dynamik. The default setting will create a fairly realistic discharge pattern.
#' The function achieves this buy using a Fourier Series of order 5 (the period is set to 2 pi), the Raingenerator function and multiple linear storage models.}
#'   \item{\link{RainGen}}{The function disaggregates a diurnal pattern to a desired timestep. The disaggregation uses the beta parameter to determine the probabilities of rain for each disaggregation step. The first step breaks 24 hours into 3 times 8 hours. The following steps bisect each of the previous time intervals.
#' If the beta parameter contains 4 values the disaggregation will result in hourly values. 6 values will result in 15 min time steps.
#' Any desired sampletime below 15 min will be calculated by simply spreading it out out evenly.
#' The disaggragation matches a probability to each disaggregation step time interval. These will be multiplied to finally determine the rain probility (portion) of the final time steps.}
#'   \item{\link{Noise}}{The function inserts noise into a given xts. (dim(xts) = c(n,1))}
#'   \item{\link{Outliers}}{The function inserts outlier into a given xts. (dim(xts) = c(n,1))}
#'   }}
#'   \item{\strong{Manning Strickler Tools}}{
#'   \describe{
#'   \item{\link{MaulProfil}}{This function analytically defines the relation between height and width of a standard Maulprofil.}
#'   \item{\link{IterMaul}}{This function numerically calculates the discharge for a single waterlevel. Using the MaulProfil function.}
#'   \item{\link{ManningMaul}}{This function numerically calculates the discharge for a vector of water levels. Using the MaulProfil and the IterMaul functions.}
#'   \item{\link{AEGG}}{This function analytically calculates the area of a flow cross section for a water level.(Egg Profile)}
#'   \item{\link{LuEGG}}{This function analytically calculates the wetted perimeter of a flow cross section for a water level.(Egg Profile)}
#'   \item{\link{ManningEGG}}{This function analytically calculates the discharge for a vector of water levels. Using the AEGG and the LuEGG functions.(Egg Profile)}
#'   \item{\link{ACIR}}{This function analytically calculates the area of a flow cross section for a water level.(Circular Profile)}
#'   \item{\link{LuCIR}}{This function analytically calculates the wetted perimeter of a flow cross section for a water level.(Circular Profile)}
#'   \item{\link{ManningCIR}}{This function analytically calculates the discharge for a vector of water levels. Using the ACIR and the LuCIR functions.(Circular Profile)}
#'   \item{\link{AREC}}{This function analytically calculates the area of a flow cross section for a water level.(Rectangular Profile)}
#'   \item{\link{LuREC}}{This function analytically calculates the wetted perimeter of a flow cross section for a water level.(Rectangular Profile)}
#'   \item{\link{ManningREC}}{This function analytically calculates the discharge for a vector of water levels. Using the AREC and the LuREC functions.(Rectangular Profile)}
#'   \item{\link{ATRA}}{This function analytically calculates the area of a flow cross section for a water level.(Trapezoid Profile)}
#'   \item{\link{LuTRA}}{This function analytically calculates the wetted perimeter of a flow cross section for a water level.(Trapezoid Profile)}
#'   \item{\link{ManningTRA}}{This function analytically calculates the discharge for a vector of water levels. Using the ATRA and the LuTRA functions.(Trapezoid Profile)}
#'   \item{\link{A_Schacht}}{This function analytically calculates the area of a flow cross section for a water level.(For a Manhole)}
#'   \item{\link{Lu_Schacht}}{This function analytically calculates the wetted perimeter of a flow cross section for a water level.(For a Manhole)}
#'   \item{\link{ManningSchacht}}{This function analytically calculates the discharge for a vector of water levels. Using the A_Schacht and the Lu_Schacht functions.(For a Manhole)}
#'   \item{\link{FitManning}}{A function for fitting Manning Strickler to a h-Q-realtion.}
#'   }}
#'   \item{\strong{Import Tools}}{
#'   \describe{
#'   \item{\link{ReadTabletoXTS}}{This function allows you to import data using the read.table function and converts to the xts format.}
#'   \item{\link{WriteXTStoTXT}}{This function allows to export an xts object to a file. It's based on the write.zoo function. Please use that one if you want more input options.}
#'   }}
#'   \item{\strong{Plot Tools}}{
#'   \describe{
#'   \item{\link{DynPlot}}{Plots data with less code than in dygraph with limited options.}
#'   \item{\link{PlotKDH2D}}{The function plots the results of the KDH function in a levelplot.}
#'   }}
#'   \item{\strong{General Purpose Tools}}{
#'   \describe{
#'   \item{\link{DeleteNAs}}{Deletes all NA - values of an xts or numeric (is equivalent to na.omit but prints the amount of NA's deleted)}
#'   \item{\link{DeleteNATimes}}{Deletes all NA - values that occur in the timestamp of an xts.}
#'   \item{\link{DeleteDuplicteTime}}{Deletes all duplicate timestamps that occur in an xts.}
#'   \item{\link{DelZeroNegValues}}{Deletes all zero and/or negative values of an xts or numeric.}
#'   \item{\link{AllFilters}}{The function utilizes the DeleteNAs(), DeleteNATimes() and DeleteDuplicteTime() functions in that order.}
#'   \item{\link{ExtendFalse}}{Extends the FALSE sections of a logical vector by a desired amount.}
#'   \item{\link{EraseTrue}}{Erases the TRUE values of a logical vector if the TRUE sections are <= than a given value.}
#'   \item{\link{CutTimeSeries}}{The function cuts an xts to a new time window.}
#'   \item{\link{CutTimeSeriesDel}}{The function deletes a time window of an xts vector.}
#'   \item{\link{FFT_XTS}}{The function uses fft to transform the xts vector to the magnitude frequency domain. Note that the first magnitude in of a n fft corresponds to the average offset of the entire signal. The plot does not contain this value, the returned vectors do.}
#'   \item{\link{FixPeriodicity}}{Fixes periodicity necessary for plotting with dygraphs by deleting NA values in time stamp and add na values for an even time step.}
#'   \item{\link{timetonumber}}{Converts a time to number.}
#'   \item{\link{Histo}}{Calculates the duration and quantity of events exceeding the limit-value.}
#'   \item{\link{KDH}}{Calculates the occurance of exceedance for given concentration limit values and temporal (duration) classes.}
#'   \item{\link{Add.NA}}{The function unifies the time step by adding NA values to gaps.}
#'   \item{\link{Rain.Intensity}}{This function calculates the intensity of rain events.}
#'   \item{\link{FillNaValues}}{This function offers a variety of options to fill gaps in the data.}
#'   \item{\link{Check_Fill_Quality}}{This function illustrates the fill quality of the nearest_neighbour method of the FillNaValues. A day at the time.}
#'   \item{\link{TimeZoneConverter}}{A function changing Timezones.}
#'   \item{\link{EventDetector}}{A function for detection discharge events.}
#'   \item{\link{BaseFlow}}{Calculates base flow Index of a flow hydrograph.}
#'   \item{\link{EventIndex}}{Hysteresis index and first flush index function.}
#'   }}
#'   \item{\strong{SN Transformation}}{
#'   \describe{
#'   \item{\link{XY2SN}}{A function for transforming XY to SN coordinates.}
#'   \item{\link{SN2XY}}{A function for transforming SN to XY coordinates.}
#'   \item{\link{Matrix.to.Raster}}{A function for transforming a matrix to a grid.}
#'
#'   }}
#' }
#'
"_PACKAGE"
#> [1] "_PACKAGE"
