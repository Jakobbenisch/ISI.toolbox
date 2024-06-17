#' Linear Regression for Polynomial Function
#'
#' Fits a ploynomial function: y=a_0+a_1*x+a_2*x^2+...a_n*x^n (order=n)
#' @param x Numeric vector.
#' @param y Numeric vector. Must be same length like x.
#' @param order Integer. The order of the polynomial function to be fitted. Must be smaller than length of x.
#' @param summary Boolean. If TRUE summary will be plotted.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{Fitall}
#' @return Returns lm object (regressionobject).
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' POL=POLRegression(x,y)
#' plot(x,y)
#' lines(x,predict(POL),col=2)
POLRegression=function(x,y,order=2,summary=TRUE){
  ## Regression (lm = linear model)
  #writes string that contains the regression information (if unclear --> print string)
  res=lm(y~polym(x,degree=order,raw=TRUE))

  if(summary){
    print(summary(res))
    print("Function: y=a_0+a_1*x+a_2*x^2+...a_n*x^n (order=n)")
    print("---Fitted Parameters---")
    print(res[1])
    print("---Confidence Levels---")
    con=confint(res, level = 0.95)
    print(con)}
  else{}
  return(res)
}


#' Linear Regression for Logarithmic Function
#'
#' Fits a logarithmic function: y=a+b*ln(x).
#' @param x Numeric vector.
#' @param y Numeric vector. Must be same length like x.
#' @param summary Boolean. If TRUE summary will be plotted.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{Fitall}
#' @return Returns lm object (regressionobject).
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' LOG=LOGRegression(x,y)
#' plot(x,y)
#' lines(x,predict(LOG),col=2)
LOGRegression=function(x,y,summary=TRUE){
  # Regression (lm = linear model)
  res = lm(y ~log(x))

  if(summary){
    print(summary(res))
    print("Function: y=a+b*ln(x)")
    print("---Fitted Parameters---")
    print(res[1])
    print("---Confidence Levels---")
    con=confint(res, level = 0.95)
    print(con)}
  else{}
  return( res)
}

#' Linear Regression for Exponential Function
#'
#' Fits a exponential function: y=a*exp(b*x).\cr
#' Note that in order to use linear regression a transformation was used. The returned intercept value needs to be transformed back. use exp() or the function Par_Exp_Pot.
#' @param x Numeric vector.
#' @param y Numeric vector. Must be same length like x.
#' @param summary Boolean. If TRUE summary will be plotted.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{Fitall}, \link[ISI.Toolbox]{Par_Exp_Pot}
#' @return Returns lm object (regressionobject).
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' EXP=EXPRegression(x,y)
#' plot(x,y)
#' lines(x,exp(predict(EXP)),col=2)
EXPRegression=function(x,y,summary=TRUE){
  # Regression (lm = linear model)
  res = lm(log(y) ~x)

  if(summary){
    print(summary(res))
    print("Function: y=a*exp(b*x)")
    print("---Fitted Parameters---")
    print(res[1])
    print("---Confidence Levels---")
    con=confint(res, level = 0.95)
    print(con) }
  else{}
  return(res)
}

#' Linear Regression for Potency Function.
#'
#' Fits a potency function: y=a*x^b.\cr
#' Note that in order to use linear regression a transformation was used. The returned intercept value needs to be transformed back. use exp() or the function Par_Exp_Pot.\cr
#' ( That's the regression, which is recommended for a h-Q-relation )
#' @param x Numeric vector.
#' @param y Numeric vector. Must be same length like x.
#' @param summary Boolean. If TRUE summary will be plotted.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{Fitall}, \link[ISI.Toolbox]{Par_Exp_Pot}
#' @return Returns lm object (regressionobject).
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' POT=POTRegression(x,y)
#' plot(x,y)
#' lines(x,exp(predict(POT)),col=2)
POTRegression=function(x,y,summary=TRUE){
  # Regression (lm = linear model)
  res = lm(log(y) ~ log(x))

  if(summary){
    print(summary(res))
    print("Function: y=a*x^b")
    print("---Fitted Parameters---")
    print(res[1])
    print("---Confidence Levels---")
    con=confint(res, level = 0.95)
    print(con) }
  else{}
  return( res)
}


#' Waterlevel Discharge Relations via Regression.
#'
#' Function calculates a discharge time series using a water level time series and a previously calculated regressionobject (POT=POTRegression(x,y) -> POT is a regressionobject).(mind the units!)
#' @param TS xts vector containing waterlevels.
#' @param regressionobject A Regressionobject.
#' @param Type String. Type of the Regressionobject. Possible inputs are "POL", "POT", "EXP" or "LOG".
#' @return Returns xts vector containing discharge values.
#' @keywords lm, Regression
#' @seealso \link[ISI.Toolbox]{hQConf}
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' TS=xts(x,order.by=seq(strptime("19.08.2016 10:00","%d.%m.%Y %H:%M",tz=""), strptime("19.08.2016 10:09","%d.%m.%Y %H:%M",tz=""), "mins"))
#' P=POTRegression(x,y,summary=FALSE)
#' Q=hQ(TS,P,Type='POT')
hQ=function(TS,regressionobject,Type="POL"){
  t=time(TS)
  x=as.numeric(TS)
  if(Type=="POL"){
    abc=Par_Poly_Log(regressionobject)
    string='Y=abc[1]'
    for(i in 1:(length(abc)-1)){
      string=paste(string,paste("+ abc[",toString(i+1),"]*x^",toString(i),sep=""),sep="")}
    eval(parse(text=string))
  }
  else if(Type=="LOG"){
    ab=Par_Poly_Log(regressionobject)
    a=ab[1]
    b=ab[2]
    Y=a+log(x)*b}
  else if(Type=="EXP"){
    ab=Par_Exp_Pot(regressionobject)
    a=ab[1]
    b=ab[2]
    Y=a*exp(b*x)}
  else if(Type=="POT"){
    ab=Par_Exp_Pot(regressionobject)
    a=ab[1]
    b=ab[2]
    Y=a*x^b}
  else{warning("Wrong input Type! Valid inputs are: 'POL', 'POT', 'EXP', 'LOG'")}
  TS=xts(Y,order.by=t)
  return(TS)
}


#' Confidence Bands for Waterlevel Discharge Relations via Regression.
#'
#' Function calculates 2 discharge time series (lower and upper confidence bands) using a waterlevel time series and Regressionobject.(mind the units!)
#' @param TS xts vector containing waterlevels.
#' @param regressionobject A Regressionobject.
#' @param Type String. Type of the Regressionobject. Possible inputs are "POL", "POT", "EXP" or "LOG".
#' @return Returns list of length 2. \cr
#' [[1]] -> xts vector containing discharge values for lower confidence band.\cr
#' [[2]] -> xts vector containing discharge values for upper confidence band.
#' @keywords lm, Regression
#' @seealso \link[ISI.Toolbox]{hQ}
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' TS=xts(x,order.by=seq(strptime("19.08.2016 10:00","%d.%m.%Y %H:%M",tz=""), strptime("19.08.2016 10:09","%d.%m.%Y %H:%M",tz=""), "mins"))
#' P=POTRegression(x,y,summary=FALSE)
#' hQC=hQConf(TS,P,Type='POT')
hQConf=function(TS,regressionobject,Type="POL"){
  t=time(TS)
  x=as.numeric(TS)
  if(Type=="POL"|Type=="LOG"){
    Pre=predict(regressionobject,interval="confidence",newdata=data.frame(x=x))
    Y_l=Pre[,2]
    Y_u=Pre[,3]}
  else if(Type=="EXP"|Type=="POT"){
    Pre=exp(predict(regressionobject,interval="confidence",newdata=data.frame(x=x)))
    Y_l=Pre[,2]
    Y_u=Pre[,3]}
  else{warning("Wrong input Type! Valid inputs are: 'POL', 'POT', 'EXP', 'LOG'")}
  TS_l=xts(Y_l,order.by=t)
  TS_u=xts(Y_u,order.by=t)
  TS=cbind(TS_l,TS_u)
  names(TS)=c("lower_Conf_level","upper_Conf_level")
  return(TS)
}


#############################################################################
#####   Regression for periodic time series analysis (Diurnal Patter)   #####
#############################################################################
#Regression for Fourier type function:
#y=a_0 +a_1*sin(p*1*t)+b_1*cos(p*1*t) +a_2*sin(p*2*t)+b_2*cos(p*2*t) +... +a_n*sin(p*n*t)+b_n*cos(p*n*t)
#(orderoftransformation=n, periode = p --> is fix, set to 1 day)
#Input:   Time series, order of transformation (n), print the summary (TRUE or FALSE)
#Output:  R - regression-object  --> contains all information about regression
#Example: source("C:/Work/R/Filters.R")
#         D2=DummyTS(days=14,Rain_off=TRUE)
#         F=FourierRegression(D2)

#' Linear Regression for Fourier Series.
#'
#' Fits Fourier Series to xts vector. As the period is set to one day, this function is mainly to calculate diurnal patterns.
#' Fourier Series: y=a_0 +a_1*cos(p*1*t)+b_1*sin(p*1*t) +a_2*cos(p*2*t)+b_2*sin(p*2*t) +... +a_n*cos(p*n*t)+b_n*sin(p*n*t)
#' (orderoftransformation=n, periode = p --> is fix, set to 1 day)
#' @param timeseries xts vector.
#' @param orderoftransformation Integer. Order of Fourier Series (n).
#' @param summary Boolean. If TRUE summary will be plotted.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{FourierPattern}
#' @return Returns lm object (regressionobject).
#' @keywords lm, Regression, Fourier Series, Fourier, diurnal pattern
#' @export
#' @examples
#' D2=DummyTS(days=14,Rain_off=TRUE)
#' F=FourierRegression(D2)
FourierRegression = function(timeseries, orderoftransformation=5,summary=TRUE){
  per=2*pi
  y=as.numeric(timeseries)
  #converts datetimes to numbers
  x=as.numeric(time(timeseries))/(24*3600)
  ## Regression (lm = linear model)
  #writes string that contains the regression information (if unclear --> print string)
  string="y ~ cos(per*x)+sin(per*x)"
  for( i in 2:orderoftransformation){
    string=paste(string,paste("cos(per*x*",toString(i),")+sin(per*x*",toString(i),")",sep=""), sep="+")
  }
  #executes string
  mod=lm(eval(parse(text = string)))

  if(summary){
    print(summary(mod))
    print("---Confidence Levels---")
    con=confint(mod, level = 0.95)
    print(con)}
  else{}
  return(mod)
}


#' Fourier Pattern from Regressionobject.
#'
#' Calculates an xts vector containing the Fourier Pattern using xts vector and Fourier regressionobject.
#' @param timeseries xts vector.
#' @param Regressionobject Regressionobject from FourierRegression.
#' @param summary Boolean. If TRUE summary will be plotted.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{FourierRegression}
#' @return Returns xts vector with Fourier Pattern.
#' @keywords lm, Regression, Fourier Series, Fourier, diurnal pattern
#' @export
#' @examples
#' D2=DummyTS(days=14,Rain_off=TRUE)
#' F=FourierRegression(D2)
#' DP=FourierPattern(D2,F)
#' DynPlot(cbind(D2,DP))
FourierPattern = function(timeseries, Regressionobject){
  x=as.numeric(time(timeseries))/(24*3600)
  order=(length(Regressionobject[[1]])-1)/2
  per=2*pi
  X=rep(1,length(x))
  for(i in 1:order){
    cos=paste("cos(per*x*",toString(i),")",sep="")
    sin=paste("sin(per*x*",toString(i),")",sep="")
    a=eval(parse(text=cos))
    b=eval(parse(text=sin))
    X=cbind(X,a,b)
  }
  coeff=as.matrix(as.numeric(Regressionobject[[1]]))
  Pattern=X%*%coeff
  Pattern=xts(Pattern,order.by=time(timeseries))
  return(Pattern)
}

#' Extraction Function.
#'
#' Extracts the coefficients of a regressionobject from POLRegression or LOGRegression.
#' @param regressionobject Regressionobject from POLRegression or LOGRegression.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{POLRegression}, \link[ISI.Toolbox]{LOGRegression}
#' @return Returns coefficients from regressionobject without any transformation.
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' POL=POLRegression(x,y)
#' print(Par_Poly_Log(POL))
Par_Poly_Log=function(regressionobject){
  par=c()
  for( i in 1:length( regressionobject[[1]])){
    par=append(par,as.numeric(regressionobject[[1]][i]))}
  return(par)}

#' Extraction Function.
#'
#' Extracts the coefficients of a regressionobject from POTRegression or EXPRegression.
#' @param regressionobject Regressionobject from POTRegression or EXPRegression.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{POTRegression},\link[ISI.Toolbox]{EXPRegression}
#' @return Returns coefficients from regressionobject with transformation.
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' POT=POTRegression(x,y)
#' print(Par_Exp_Pot(POT))
Par_Exp_Pot=function(regressionobject){
  par=as.numeric(c(exp(regressionobject[[1]][1]),regressionobject[[1]][2]))
  return(par)}

#' Extraction Function.
#'
#' Extracts the confidence intervals of coefficients of a regressionobject from POLRegression or LOGRegression.
#' @param regressionobject Regressionobject from POLRegression or LOGRegression.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{POLRegression},\link[ISI.Toolbox]{LOGRegression}
#' @return Returns confidence intervals of coefficients from regressionobject without any transformation.
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' POL=POLRegression(x,y)
#' print(ConL_Poly_Log(POL))
ConL_Poly_Log=function(regressionobject){
  con=confint(regressionobject, level = 0.95)
  low=c()
  up=c()
  count=length(con)/2
  for( i in 1:count){
    low=append(low,as.numeric(con[i]))
    up=append(up,as.numeric(con[(i+count)]))}
  return(list(low,up))}

#' Coefficient Extraction from Regressionobject.
#'
#' Extracts the confidence intervals of coefficients of a regressionobject from POTRegression or EXPRegression
#' @param regressionobject Regressionobject from POLRegression or LOGRegression.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{POTRegression},\link[ISI.Toolbox]{EXPRegression}
#' @return Returns confidence intervals coefficients from regressionobject with transformation.
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' POT=POTRegression(x,y)
#' print(ConL_Exp_Pot(POT))
ConL_Exp_Pot=function(regressionobject){
  con=confint(regressionobject, level = 0.95)
  low=as.numeric(c(exp(con[1]),con[2]))
  up=as.numeric(c(exp(con[3]),con[4]))
  return(list(low,up))}



#Fits all regression types except of fourier and plots the results
#Input:   X-array, y-array      (without NA's)
#Output:  Plot
#Example: x=seq(1,10,1)
#         y=exp(x)+seq(10,100,10)
#         TS=xts(x,order.by=seq(strptime("19.08.2016 10:00","%d.%m.%Y %H:%M",tz=""), strptime("19.08.2016 10:09","%d.%m.%Y %H:%M",tz=""), "mins"))
#         Fitall(x,y)

#' Linear Regression for Log,Pot,Exp,Pol Functions.
#'
#' Fits all regression types except of fourier and plots the results. This function aims to give an overview to find the best fit.
#' @param X Numeric vector.
#' @param Y Numeric vector. Must be same length like x.
#' @param poly_order Integer. Order of polynomial function.
#' @seealso \link[stats]{lm}, \link[ISI.Toolbox]{POLRegression},\link[ISI.Toolbox]{POTRegression},\link[ISI.Toolbox]{LOGRegression},\link[ISI.Toolbox]{EXPRegression}
#' @return Returns a plot.
#' @keywords lm, Regression
#' @export
#' @examples
#' x=seq(1,10,1)
#' y=exp(x)+seq(10,100,10)
#' Fitall(x,y,poly_order=4)
Fitall=function(X,Y,poly_order=2){
  #orders data
  DataF=data.frame(X,Y)
  DataF=DataF[order(DataF[, 1]),]
  x=DataF[,1]
  y=DataF[,2]

  POL=POLRegression(x,y,order=poly_order)
  LOG=LOGRegression(x,y)
  EXP=EXPRegression(x,y)
  POT=POTRegression(x,y)

  plot(x,predict(POL),type='l',col='blue')
  lines(x,predict(LOG),col='red')
  lines(x,exp(predict(EXP)),col="green")
  lines(x,exp(predict(POT)),col="black")

  #generates labelnames with the coefficient of determination
  R2=function(regressionobject){
    R=summary(regressionobject)$r.squared
    name=deparse(substitute(regressionobject))
    Ri=paste(name," R^2=",toString(round(R,3)),sep="")
    return(Ri)}
  legend("topleft",c(R2(POL),R2(LOG),R2(EXP),R2(POT)),lwd=c(2.5,2.5,2.5,2.5),col=c("blue","red","green","black"))
  points(x,y)
}

#' A function for thinning out data.
#'
#' This function sets up an controled number of evenly or unevenly distributed spaces of a range. For each space all data points within
#' will be summarized using mean or median. This funtion is mainly useful for regression analysis.
#' @param TS Xts. A time series to be thinned out with dim = n,2!
#' @param thin.to Integer, that specifies the number of ranges that will be filled 1 data point each, if data is available.
#' @param median_mode String, that specifies which datapoint to take when 2 points are available. This is relevant due to the fact, that the median of a even number of samples is not computed with the mean of the 2 inner samples. Available options are "min" and "max". Default is "min".
#' @param style String, that specifies the distribution of the ranges. Available options are "even","tight top" and "tight bottom". Default is "even".
#' Even - all ranges have the same width. Tight top - lower ranges are wider, larger one are tighter. Thight bottom - the opposite of tight top.
#' @return Data Frame with dim = n,2.
#' @export
#' @examples
#' # generating sample data
#' H=DummyTS(day=50,Anual_Days_of_P =60)
#' h=as.numeric(H)/100
#' Q=40*(h)^1.4
#' Q1=Q
#' Q=Q+rnorm(length(Q),seq(0,max(Q)/10,length.out = length(Q)),median(Q)/6) # add noise
#' print(all(Q==Q1))
#' hQ=xts(cbind(h,Q),order.by=time(H))
#' # get rid of negative flow (produced by random sampling)
#' hQ=hQ[hQ[,2]>0]
#' h=as.numeric(hQ[,1])
#' Q=as.numeric(hQ[,2])
#'
#' # plot original data and fit
#' Fitall(h,Q)
#' # thinn out data
#' hQ_thinn=ThinOutData(hQ,thin.to=50,style="tight top")
#' # plotting new fit
#' h_thinn=as.numeric(hQ_thinn[,1])
#' Q_thinn=as.numeric(hQ_thinn[,2])
#' Fitall(h_thinn,Q_thinn)
ThinOutData=function(data,thin.to=300,style=c("even","tight top","tight bottom"),median_mode=c("min","max")) {
  if (length(style)==3){style="even"}
  if (length(median_mode)==2){median_mode="min"}
  if(median_mode!="max" & median_mode!="min"){stop("Median_mode must be either max or min.")}
  Normalize = function(array, x, y){

    # Normalize to [0, 1]:
    m = min(array)
    range = max(array) - m
    array = (array - m) / range

    # Then scale to [x,y]:
    range2 = y - x
    normalized = (array*range2) + x
    return(normalized)
  }
  if (is.xts(data)){
    XTS_clean = na.omit(data)
    DataF = data.frame(XTS_clean[, 1], XTS_clean[,2])}
  else if (is.data.frame(data)){DataF=data}
  else{stop("The input must be a dataframe or an xts.")}
  DataF = DataF[order(DataF[, 1]), ]
  x = as.numeric(DataF[, 1])
  y = as.numeric(DataF[, 2])
  min = min(x, na.rm = TRUE)
  max = max(x, na.rm = TRUE)
  X = seq(min, max, length.out=thin.to)
  if (style=="even"){}
  else if (style=="tight top"){
    X=log(X)
    X=Normalize(X,min,max)
  }
  else if (style=="tight bottom"){
    X=exp(X)
    X=Normalize(X,min,max)
  }
  else{stop("Incorrect style. Possible inputs are 'even','tight top' and 'tight bottom'.")}

  Median=function(x,choose="max"){
    x_sorted=sort(x)
    len=length(x_sorted)
    if (len%%2==0){
      ind=as.integer(len/2)
      if(choose=="max"){
        ind=ind+1
      }
    }
    else{
      ind=is.integer(len/2)+1
    }
    return(x_sorted[ind])
  }
  newx = c()
  newy = c()
  for (i in 2:length(X)) {
    logi = x >= X[i - 1] & x <= X[i]
    if (!all(!logi)) {
      med=Median(x[logi],choose=median_mode)
      newx = c(newx, med)
      yn=y[x==med]
      if (length(yn>1)){
        if(median_mode=="min"){yn=min(yn)}
        else{yn=max(yn)}
      }
      newy = c(newy,yn )

    }
  }

  return(data.frame(newx = newx, newy = newy))
}






