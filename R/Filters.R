#' A NA Deletion Function
#'
#' Deletes all NA - values of an xts or numeric (is equivalent to na.omit but prints the amount of NA's deleted)
#' @param vector numeric or xts to be freed of NA's.
#' @param Print Boolean. Print the number of NA's deleted?
#' @param return_index Boolean. If TRUE returns a list. First position is the filtered time series. Second position is a vector containing the indexes of the filtered values.
#' @seealso \link[stats]{na.omit}
#' @return Vector or xts without NA - values or depending on the return_index paramter a list.
#' @keywords delete, NA
#' @export
#' @examples
#' TimeSeries=DeleteNAs(TimeSeries)
#'
DeleteNAs = function(vector,Print=TRUE,return_index=FALSE){
  index=is.na(vector)
  vector=vector[!index]
  outputstring=paste(toString(sum(as.integer(index))),"NA values have been deleted!", sep = " ")
  if(Print){print(outputstring)}
  if(return_index){vector=list(vector,which(index));names(vector)=c("filtered xts, indexes of filtered values")}
  return(vector)
}

#' A NA Deletion Function
#'
#' Deletes all NA - values that occur in the timestamp of an xts
#' @param vector xts to be freed of NA - timestamps.
#' @param Print Boolean. Print the number of NA - timestamps deleted?
#' @param return_index Boolean. If TRUE returns a list. First position is the filtered time series. Second position is a vector containing the indexes of the filtered values.
#' @keywords delete, NA
#' @return Vector or xts without NA - values in timestamp or depending on the return_index paramter a list.
#' @export
#' @examples
#' TimeSeries=DeleteNATimes(TimeSeries)
#'
DeleteNATimes = function(vector, Print=TRUE,return_index=FALSE){
  index=is.na(time(vector))
  vector=vector[!index]
  if(Print){
    outputstring=paste(toString(sum(as.integer(index))),"NA datetime/s have been deleted!", sep = " ")
  print(outputstring)
  }
  if(return_index){vector=list(vector,which(index));names(vector)=c("filtered xts, indexes of filtered values")}
  return(vector)
}

#' A Deletion Function
#'
#' Deletes all duplicate timestamps that occur in an xts
#' @param vector xts to be freed of duplicate values in timestamps.
#' @param Print Boolean. Print the number of duplicate timestamps deleted?
#' @param return_index Boolean. If TRUE returns a list. First position is the filtered time series. Second position is a vector containing the indexes of the filtered values.
#' @keywords delete, NA
#' @return Vector or xts without duplicate values in timestamp or depending on the return_index paramter a list.
#' @export
#' @examples
#' TimeSeries=DeleteDuplicteTime(TimeSeries)
#'
DeleteDuplicteTime=function(vector,Print=TRUE,return_index=FALSE){
  index=duplicated(time(vector))
  vector=vector[!index]
  if(Print){
    outputstring=paste(toString(sum(as.integer(index))),"duplicate datetime/s have been deleted!", sep = " ")
    print(outputstring)
  }
  if(return_index){vector=list(vector,which(index));names(vector)=c("filtered xts, indexes of filtered values")}
  return(vector)
}

#' A Function for manupulation of logical vectors
#'
#' Extends the FALSE sections of a logical vector by a desired amount.
#' @param vec logical vector.
#' @param amount amount of extention of FALSE section.
#' @keywords logical
#' @return manipulated logical vector
#' @export
#' @examples
#' a=c(TRUE,TRUE, FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE)
#' b=ExtendFalse(a,1)
#' c=ExtendFalse(a,2)
#'
ExtendFalse = function(vec,amount=1){
  vector2=vec
  if (amount==0){vector2=vec}
  else{
    ind=which(vec)
    lo=diff(ind)>1
    LowEdge=ind[lo]
    TopEdge=diff(ind)[lo]+LowEdge
    for(i in 0:(amount-1)){
      LowEdge=union(LowEdge,LowEdge-i)
      TopEdge=union(TopEdge,TopEdge+i)
    }
    indexes=union(LowEdge,TopEdge)
    indexes[indexes>0&indexes<length(vec)]

    vector2[indexes]=logical(length(indexes))}
  return(vector2)
}

#' A Function for manupulation of logical vectors
#'
#' Erases the TRUE values of a logical vector if the TRUE sections are <= than a given value
#' @param LogiVec logical vector.
#' @param UnplausibleGapsize value up to which the TRUE sections will be replaced with FALSE sections
#' @keywords logical
#' @return manipulated logical vector
#' @export
#' @examples
#' a=c(TRUE,TRUE, FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE)
#' b=EraseTrue(a,1)
#'
EraseTrue=function(LogiVec,UnplausibleGapsize=2){
  index=which(!LogiVec)
  distance=diff(index)-1
  logi=which(distance<=UnplausibleGapsize&distance>0)
  for(i in logi){
    vec=seq(index[i],index[i+1])
    index=union(index,vec)}
  LogiVec[index]=logical(length(index))
  return(LogiVec)}

#' A Deletion Function
#'
#' Deletes all zero and/or negative values of an xts or numeric
#' @param vector xts or numeric vector.
#' @param sign determines the relation that will be used for deleting. Possible inputs are "<0" or "<=0".
#' @param extendGap setting this paramter to a value greater 0 will utilize the ExtendFalse() function, thereby extending the resulting gaps in the data.
#' @keywords delete
#' @return The cleaned xts vector.
#' @export
#' @seealso  \link[ISI.Toolbox]{DummyTS},\link[ISI.Toolbox]{DynPlot},\link[ISI.Toolbox]{ExtendFalse}
#' @examples
#' TS=DummyTS()-15
#' TS2=DelZeroNegValues(TS,"<=0",1)
#' DynPlot(cbind(TS,TS2))
#'
DelZeroNegValues = function(vector,sign="<0",extendGap=0){
  if(sign=="<0"){
    index=!vector<0.0
  }
  else if(sign=="<=0"){
    index=vector>0.0
  }
  else{print("Wrong input string!")}
  index=ExtendFalse(index,1)
  vector=vector[index]
  return(vector)
}

#' A NA Deletion Function
#'
#' The function utilizes the DeleteNAs(), DeleteNATimes() and DeleteDuplicteTime() functions in that order
#' @param XTS xts vector to be freed of NA's.
#' @keywords delete
#' @return The cleaned xts.
#' @export
#' @seealso  \link[ISI.Toolbox]{DeleteNATimes},\link[ISI.Toolbox]{DeleteDuplicteTime}
#' @examples
#' TS=AllFilters(TS)
#'
AllFilters = function(XTS){
  XTS=DeleteNAs(XTS)
  XTS=DeleteNATimes(XTS)
  XTS=DeleteDuplicteTime(XTS)
  return(XTS)
}

#' Function for cutting an xts vector
#'
#' The function cuts an xts to a new time window
#' @param TimeSeries xts vector to be cut.
#' @param start The new start of the xts (format: "\%Y-\%m-\%d \%H:\%M")
#' @param end The new end of the xts (format: "\%Y-\%m-\%d \%H:\%M")
#' @keywords cut
#' @return xts vector cut to the desired time window.
#' @export
#' @seealso  \link[ISI.Toolbox]{CutTimeSeriesDel}
#' @examples
#' XTS = CutTimeSeries(XTS, "2016-06-16 12:00", "2016-08-05 00:00")
#'
CutTimeSeries = function (TimeSeries, start, end)
{
  TZ = tzone(TimeSeries)
  ind1 = time(TimeSeries) >= as.POSIXct(start, tz = TZ)
  ind2 = time(TimeSeries) <= as.POSIXct(end, tz = TZ)
  if(as.POSIXct(end, tz = TZ)<as.POSIXct(start, tz = TZ)){print("care: end of time series is earlier than its start!")}
  TimeSeries = TimeSeries[ind1 & ind2]
  return(TimeSeries)
}
#' Function for cutting an xts vector
#'
#' The function deletes a time window of an xts vector
#' @param TimeSeries xts vector to be cut.
#' @param start The start of the window to be deleted (format: "\%Y-\%m-\%d \%H:\%M")
#' @param end The end of the window to be deleted (format: "\%Y-\%m-\%d \%H:\%M")
#' @keywords cut, delete
#' @return xts vector without the cut time window.
#' @export
#' @seealso  \link[ISI.Toolbox]{CutTimeSeries}
#' @examples
#' XTS = CutTimeSeriesDel(XTS, "2016-06-16 12:00", "2016-08-05 00:00")
#'
CutTimeSeriesDel=function(TimeSeries,start,end){
  TZ = tzone(TimeSeries)
  ind1 = time(TimeSeries) >= as.POSIXct(start, tz = TZ)
  ind2 = time(TimeSeries) <= as.POSIXct(end, tz = TZ)
  if(as.POSIXct(end, tz = TZ)<as.POSIXct(start, tz = TZ)){print("care: end of time series is earlier than its start!")}
  TimeSeries = TimeSeries[!(ind1&ind2)]
  return(TimeSeries)}

#' Function for deleting outliers in an xts vector
#'
#' The function filters gradients greater than a set absolute value of an xts vector. The filter repeats its operation until the acceptableError is reached.
#' Use the Analyse option (set value to TRUE) to get some statistics about the distribution of gradients in your time series.
#' Once you know how to set your gradient limits turn the Analyse function to FALSE and rerun the Gradfilter. Limit exceeding values will then be deleted.
#' @param TS Xts vector to be freed of outliers.
#' @param abs_grad The absolute value of the gradient limit.
#' @param acceptableError The allowed number of Gradient errors in percent. (Number of Errors / Number of Observations *100)
#' @param UnplausibleGapsize This parameter allows to filter noisy signals to some extend. If set to a number >=1 the filter will also delete oberservations between two outlier if the distance is smaller or equal to that value.
#' @param Analyse Boolean. If TRUE a summary of the gradients distribution and position will be plotted.
#' @param GradZero Boolean. If TRUE values that result in gradients of 0 will be deleted.
#' @param Automatic Boolean. If TRUE the abs_grad is set automatically using R's quantile function with a default probability of 0.98.
#' @param OnlyDoubleGradientError Boolean. If TRUE the filter will only delete values that result in both a negative and positive gradient error.
#' @param Prob The probability of for the Automatic modus.
#' @param NormalizedGrad Boolean. If True each gradient will take the its neighbours into account to decide, whether it is a oultier or a valid value.\cr
#' Each transformed gradient G_i will be computed as follows: \cr
#' G_{i}=g_{i}/mean(g_(i-l/2)+g_(i-l/2+1)+...g_(i+l/2-1)+g_(i+l/2-1))\cr
#' g_{i} is the ith gradient\cr
#' l is the number of neighbour gradients to be considered\cr
#' G_{i} is the transformed gradient and will be used for filtering
#' @param windowlengthNG Integer. Specifies the number of neighbour gradients to be considered for the NormalizedGrad method.This value must be even.
#' @keywords Outlier, Filter
#' @return xts vector without outliers.
#' @export
#' @seealso  \link[ISI.Toolbox]{NoiseFilter}
#' @examples
#' A = DummyTS(days=50)
#' B = Outliers(A,5,4)
#' GradFilter(B,Analyse = T,NormalizedGrad = TRUE)
#' Bf=GradFilter(B,abs_grad = 3,NormalizedGrad = TRUE)
#' DynPlot(cbind(B,Bf),Labels=c("TS_with_Outliers","TS_without_Outliers"))
#'
GradFilter=function(TS, abs_grad=9999999,acceptableError=0.01, UnplausibleGapsize=NULL, Analyse=FALSE, GradZero=FALSE, Automatic=FALSE,OnlyDoubleGradientError=FALSE,Prob=0.98,NormalizedGrad=FALSE,windowlengthNG=4){
  if(OnlyDoubleGradientError & NormalizedGrad){stop("The 2 options OnlyDoubleGradientError and NormalizedGrad can not be combined. Please turn one off.")}
  if(NormalizedGrad& windowlengthNG%%2!=0){stop("The parameter windowlengthNG must be even!")}
  #seperating datetime from numeric values
  TStime=as.numeric(time(TS))
  TSvalue=as.numeric(TS)
  #calculating the gradients
  deltaX=diff(TStime)
  deltaY=diff(TSvalue)
  grad=deltaY/deltaX
  if(NormalizedGrad){
    NormGrad=function(i){return(grad[i]/(mean(grad[(i-windowlengthNG/2):(i+windowlengthNG/2)])))}
    index=(1+windowlengthNG/2):(length(grad)-windowlengthNG/2)
    Grad=sapply(index,FUN=NormGrad)
    grad=Grad
  }
  if(Analyse){
    grad_abs=abs(grad)
    #printing and plotting gradient statistics
    cat(paste("mean grad:","\t",sprintf("%.2e",mean(grad_abs)),"\n",
              "min grad:","\t",sprintf("%.2e",min(grad_abs)),"\n",
              "max grad:","\t",sprintf("%.2e",max(grad_abs)),"\n",
              "median grad:","\t",sprintf("%.2e",median(grad_abs)),"\n",sep=""))
    par(mfrow=c(1,2))
    hist(grad, main="Histogram",xlab = "Gradient",breaks = 20)
    plot(index(grad),grad,xlab = "Index",ylab = "Gradient",main = "Location of Outlier")
    par(mfrow=c(1,1))}
  else{
    if(Automatic){abs_grad=quantile(abs(grad), probs=Prob)
    print(abs_grad)}

    #creating logical vector that marks gradient errors wich FALSE
    if(OnlyDoubleGradientError){ind=!((c(FALSE,grad>abs_grad)&c(grad<(-abs_grad),FALSE))|(c(grad>abs_grad,FALSE)&c(FALSE,grad<(-abs_grad))))}
    else{ if(NormalizedGrad){ind=c(TRUE,!logical(windowlengthNG/2),grad<abs_grad&grad>-abs_grad,!logical(windowlengthNG/2))}
      else{ind=c(TRUE,grad<abs_grad&grad>-abs_grad)}
    }

    #filters 0-gradients, if turned on
    if(GradZero){index0=!grad==0;ind=ind&c(TRUE,index0)}
    #if plausible gap size is a number, the set amount equals non-gradient-error values between gradient-errors will also be filtered
    if(is.null(UnplausibleGapsize)==FALSE){ind=EraseTrue(ind,UnplausibleGapsize)}
    #cleaning the timeseries
    TS=TS[ind]

    #running a loop (repeating the filtering) until the acceptable error of gradients within the time series is reached.
    while((length(ind)-table(ind)["TRUE"])/length(ind)>acceptableError/100){
      #seperating datetime from other values
      TStime=as.numeric(time(TS))
      TSvalue=as.numeric(TS)
      #calculating the gradients
      deltaX=diff(TStime)
      deltaY=diff(TSvalue)
      grad=deltaY/deltaX
      if(NormalizedGrad){
        index=(1+windowlengthNG/2):(length(grad)-windowlengthNG/2)
        Grad=sapply(index,FUN=NormGrad)
        grad=Grad
      }
      if(OnlyDoubleGradientError){ind=!((c(FALSE,grad>abs_grad)&c(grad<(-abs_grad),FALSE))|(c(grad>abs_grad,FALSE)&c(FALSE,grad<(-abs_grad))))}
      else{ if(NormalizedGrad){ind=c(TRUE,!logical(windowlengthNG/2),grad<abs_grad&grad>-abs_grad,!logical(windowlengthNG/2))}
        else{ind=c(TRUE,grad<abs_grad&grad>-abs_grad)}
      }
      TS=TS[ind]}
    return(TS)}}
#' Function for removing noise of an xts vector.
#'
#' The function filters noise using std within a moving window. The Regression parameter allows to fit polynoms in each moving window and the sd is calculated of the residuals.
#' This allows for reduced influence of trends that might be in the xts vector. Use the Analyse function to get an idea about the sd distribution of the time series and thereby setting the s_max values.
#' @param TS xts vector to be cut.
#' @param window_width The width of the moving window.
#' @param s_max The maximal Standard Deviation that is allowed.
#' @param lag The lag from one moving window to the next.
#' @param Regression Boolean. If TRUE linear regression will be used for each moving window.
#' @param Polydegree The order of the polynomial function that will be fitted in each moving window.
#' @param Analyse Boolean. If TRUE a summary of the sd's of all moving windows will be plotted.
#' @keywords Filter, Noise
#' @return Xts vector without noise.
#' @export
#' @seealso  \link[ISI.Toolbox]{GradFilter}
#' @examples
#' #Creating xts with noise
#' D=DummyTS(days=7, SampleTime_min = 1, Rain_off = T)
#' D=Noise(D,15,1.5)
#'
#' #Analysing
#' Dclean=NoiseFilter(D,window_width=20,Analyse=TRUE)
#'
#' #Filtering
#' Dclean=NoiseFilter(D,s_max=0.6,window_width=20)
#'
#' #Plotting
#' DynPlot(cbind(D,Dclean),Labels=c("TS_with_Noise","TS_without_Noise"))
#'
NoiseFilter= function(TS,window_width,s_max=99999999,lag=1,Regression=FALSE,Polydegree=1,Analyse=FALSE){
  #seperating datetime from other values
  Times=time(TS)
  Values=as.numeric(TS)
  Index=seq(1,length(Values),1)
  #fit a polynomial function to each moving window and calculates the std of the residuals (thereby eliminating the influence of trends)
  if(Regression){
    order=Polydegree
    PolyN=function(X){
      y=as.numeric(X)
      x=(as.numeric(time(X))-as.numeric(time(X[1])))/60
      if(order==1){string="fit =lm(y ~ x)"}
      else{
        string="fit =lm(y ~ x"
        for(i in 2:order){
          string=paste(string,paste("+ I(x^",toString(i),")",sep=""),sep="")}
        string=paste(string,")",sep="")}
      eval(parse(text=string))
      residuals=resid(fit)
      return(sd(residuals))}
    #applys the moving window function
    SDVEC=rollapply(TS, width = window_width, by=lag, FUN = PolyN, align = "left")
    SDVEC=SDVEC[!is.na(SDVEC)]
  }
  #calculates the std without any regression
  else{
    SDVEC=rollapply(Values, width = window_width, by=lag, FUN = sd, align = "left")}
  #illustrates some std statistics
  if(Analyse){cat(paste("mean std:\t", sprintf("%.2e",mean(SDVEC)),"\n",
                        "min std:\t", sprintf("%.2e",min(SDVEC)),"\n",
                        "max std:\t", sprintf("%.2e",max(SDVEC)),"\n",
                        "median std:\t", sprintf("%.2e",median(SDVEC)),"\n",sep = ""))
    par(mfrow=c(1,2))
    hist(SDVEC, xlab = "Standard Deviation",main="Histogram")
    plot(index(SDVEC),SDVEC,xlab = "Moving Window",ylab = "Standard Deviation",main = "Location of Noise")
    par(mfrow = c(1, 1))}
  else{
    #calculates the indexes of all moving windows
    LowIndexes=rollapply(Index, width = window_width, by=lag, FUN = min, align = "left")
    HighIndexes=LowIndexes+window_width-1
    #identifies the moving window indexes of noise
    if(length(s_max)==1){SDIndex=which(SDVEC>s_max)}
    else{SDIndex=which(SDVEC>s_max[1]&SDVEC<s_max[2])}
    #sets up a logical vector an changes noise sections to false using a loop
    if(length(SDIndex)==0){warning("S_max is too high. No Data has been filtered!")}
    else{lo=!logical(length(Values))
    for(i in 1:length(SDIndex)){
      del=c(LowIndexes[SDIndex[i]]:HighIndexes[SDIndex[i]])
      lo[del]=logical(length(del))}
    Times=Times[lo]
    Values=Values[lo]}
    TS=xts(Values,order.by = Times)
    return(TS)}}



#' Function for smoothing an xts vector
#'
#' The function is for smoothing an xts or numeric vector. It's based on the publication by Géry Casiez.
#' For detailed description on how to use the filter please read the paper.
#' @source \url{http://cristal.univ-lille.fr/~casiez/1euro/}
#' @param x xts vector to be smoothed.
#' @param freq A measure on how many data point are in one cycle. For instance samples per day. (1440 is dt = 1 minute)
#' @param mincutoff The greater mincutoff the less lag and more sensitive (less memory).
#' @param beta The smaller beta the smoother.
#' @param dcutoff The greater mincutoff the less lag and more sensitive (less memory).
#' @keywords cut, delete
#' @return xts vector without the cut time window.
#' @export
#' @seealso  \link[ISI.Toolbox]{Smooth}
#' @examples
#' TS=DummyTS(Csd = 1.2)
#' TS2=OneEuroFilter(TS,freq=1440)
#' DynPlot(cbind(TS,TS2))
OneEuroFilter=function(x,freq=120,mincutoff=1,beta=1,dcutoff=1){

  XTS=FALSE
  if(is.xts(x)){time=time(x);x=as.numeric(x);XTS=TRUE}
  else if(!is.numeric(x)){}

  Alpha=function(freq,cutoff){
    te = 1.0 / freq
    tau = 1.0 / (2*pi*cutoff)
    return(1.0 / (1.0 + tau/te))}

  x_LP_Filter=function(i,xi,alpha){
    if(i==1){x_hatxprev<-xi}
    x_hatx<-alpha*xi+(1-alpha)*x_hatxprev
    x_hatxprev<<-x_hatx
    return(x_hatx)}

  dx_LP_Filter=function(i,xi,alpha){
    if(i==1){dx_hatxprev<-xi}
    dx_hatx<-alpha*xi+(1-alpha)*dx_hatxprev
    dx_hatxprev<<-dx_hatx
    return(dx_hatx)}

  Inner=function(i){
    if(i==1){dx<-0}
    else{dx<-(x[i]-x_hatxprev)*freq} #dx<-(x[i]-x_hatxprev)*freq
    edx<-dx_LP_Filter(i=i,xi=dx,alpha<-Alpha(freq,dcutoff))
    cutoff<-mincutoff+beta*abs(edx)
    Result=x_LP_Filter(i=i,xi=x[i],alpha<-Alpha(freq,cutoff))
    return(Result)}

  vec=sapply(index(x),FUN = Inner)

  if(XTS){vec=xts(vec,order.by = time)}
  return(vec)}

#' Function for smoothing an xts vector
#'
#' The function uses a moving window. For each moving window a mean value is calculated and substitutes the unsmoothed one.
#' @param XTS xts vector to be smoothened.
#' @param window_width The width of the moving window. (Ideally this is a uneven number.)
#' @keywords smooth, Filter
#' @return Smoothed xts vector.
#' @export
#' @seealso  \link[ISI.Toolbox]{OneEuroFilter}
#' @examples
#' x=DummyTS(Csd=1.3)
#' y=Smooth(x,window_width=9)
#' DynPlot(cbind(x,y))

Smooth= function(XTS,window_width=7){
  Time=time(XTS)
  Values=as.numeric(XTS)
  mean=rollapply(Values, width = window_width, by=1, FUN = function(x) mean(x, na.rm=T), align = "center")
  start=as.integer((window_width-1)/2)+1
  XTS=xts(mean,order.by = Time[start:(length(mean)+start-1)])
  return(XTS)}


#' Function for calculating an average Fourier Sieres for each day
#'
#' The function uses a linear model to fit a Fourier Series for each day of an xts vector.
#' @param TS Xts vector to be analysed.
#' @param order The order of the Fourier Series.
#' @param accepltabledev To be set in percent. When this function is used, it initially calculates a Fourier Series for the entire xts vector (resembling the average daily pattern). To ensure that the daily fitted patterns are close to the average pattern this value compares
#' the daily fitted series to the average. If the deviation is greater than the set one the calculated pattern won't be used.
#' @param MinDataPointsDaily To be able to estimate a representative daily pattern enough data points must be provided for each day. If the set value (minimum observations per day) isn't met the day will be left out.
#' @keywords smooth, Filter, Fourier
#' @return The function returns a list. Its first position contains the daily fitted Fourier Patterns while the second entry contains the diurnal residuals. The third position contains the regression object of the average Fourier Series.
#' @export
#' @seealso  \link[ISI.Toolbox]{FourierRegression}
#' @examples
#' Y=DummyTS(days=20,Rain_off=TRUE)
#' Y=Outliers(Y,5,4,keyword="long",maxlenoutlier=30)
#' D=Fourierperday(Y,4)
#' DD=D[[1]];WD=D[[2]]
#' DynPlot(cbind(Y,DD,WD),Labels=c("Data","Fourier","Residuals"))
Fourierperday=function(TS,order=4,accepltabledev=25, MinDataPointsDaily=400){
  TS=na.omit(TS)
  Time=time(TS)
  FT=FourierRegression(TS,order,summary = F)
  DAY=function(x){lubridate::date(time(x[1]))}
  days= lubridate::date(xts::apply.daily(TS,FUN = DAY))
  Resid=TS[1]
  Four=TS[1]
  for(i in 2:length(days)){
    start=as.POSIXct(days[i-1],tz="",format="%Y-%m-%d") -3600
    stop=as.POSIXct(days[i],tz="",format="%Y-%m-%d") -3600
    TS2=TS[time(TS)>start&time(TS)<stop]
    if(length(TS2)<MinDataPointsDaily){}
    else{
      FR=FourierRegression(TS2,order,summary = F)
      Fou=FourierPattern(TS2,FR)
      FouAll=FourierPattern(TS2,FT)
      Error=2*(sum(abs(as.numeric(FouAll)-as.numeric(Fou))))/(sum((as.numeric(FouAll)+as.numeric(Fou))))
      if(Error*100>accepltabledev){}
      else{
        RES=xts(resid(FR),order.by = time(TS2))
        Resid=rbind(Resid,RES)
        Four=rbind(Four,Fou)}}}
  Resid=Resid[-1]
  Four=Four[-1]
  return(list(Four,Resid,FT))}




#' Function for Fast-Fourier-Transformation of an xts vector
#'
#' The function uses fft to transform the xts vector to the magnitude frequency domain. Note that the first magnitude in of a n fft corresponds to the average offset of the entire signal. The plot does not contain this value, the returned vectors do.
#' @param TS Xts vector to be transformed.
#' @param dt The times tep in days. (One Minute time step -> 1/24/60)
#' @param Return_FreqMag_Sin_Cos Boolean. If TRUE the frequency vector and both the cosine and sine magnitudes seperately will be returned.
#' @param Return_Freq_Mag_Phase Boolean. If TRUE the frequency vector, magnitude vector and phase vector will be returned.
#' @param Plot_Freq_Mag Boolean. Returns the Frequency Magnitude plot.
#' @param Ylim Vector of length 2. Set the plot's y axis limits.
#' @param Xlim Vector of length 2. Set the plot's x axis limits.
#' @keywords Fourier, fft
#' @return Depending on the set parameters Return_FreqMag_Sin_Cos, Return_Freq_Mag_Phase and Plot_Freq_Mag a list will be returned and possible a plot will show.\cr\cr
#' If Return_FreqMag_Sin_Cos and Return_Freq_Mag_Phase are TRUE:\cr
#' The list will contain 5 entires.\cr
#' [[1]] -> frequency vector\cr
#' [[2]] -> cosine magnitude vector\cr
#' [[3]] -> sine magnitude vector\cr
#' [[4]] -> total magnitude vector\cr
#' [[5]] -> phase vector\cr\cr
#' If Return_FreqMag_Sin_Cos is TRUE and Return_Freq_Mag_Phase is FALSE:\cr
#' The list will contain 3 entires.\cr
#' [[1]] -> frequency vector\cr
#' [[2]] -> cosine magnitude vector\cr
#' [[3]] -> sine magnitude vector\cr\cr
#' If Return_FreqMag_Sin_Cos is FALSE and Return_Freq_Mag_Phase is TRUE:\cr
#' The list will contain 3 entires.\cr
#' [[1]] -> frequency vector\cr
#' [[2]] -> total magnitude vector\cr
#' [[3]] -> phase vector\cr\cr
#' If Return_FreqMag_Sin_Cos and Return_Freq_Mag_Phase are FALSE:\cr
#' An empty list will be returned.\cr\cr
#' If Plot_Freq_Mag is TRUE the frequency magnitude plot will be shown.\cr
#' @export
#' @seealso  \link[ISI.Toolbox]{OneEuroFilter}
#' @examples
#' x=DummyTS()
#' FFT_XTS(x)
FFT_XTS=function(TS,dt = 1/24/60, Return_FreqMag_Sin_Cos=FALSE,Return_Freq_Mag_Phase=FALSE,Plot_Freq_Mag=TRUE,Ylim=NULL,Xlim=NULL){
  ToReturn=list()
  #CREATE OUR TIME SERIES DATA
  y <- as.numeric(TS)
  time= time(TS)

  #Domain setup
  T <- length(time)*dt

  n <- length(time)

  #CREATE OUR FREQUENCY ARRAY
  f <- 1:length(time)/T

  ende=length(f)/2

  #FOURIER TRANSFORM WORK
  Y <- stats::fft(y)
  mag <- sqrt(Re(Y)^2+Im(Y)^2)*2/n
  phase <- atan(Im(Y)/Re(Y))
  Yr <- Re(Y)*2/n
  Yi <- Im(Y)*2/n*-1
  if(Return_FreqMag_Sin_Cos){
    ToReturn[[1]]=f[1:ende]
    ToReturn[[2]]=Yr[1:ende]
    ToReturn[[3]]=Yi[1:ende]
    names(ToReturn)=c("Frequency","Cosine Values","Sine Values")}
  if(Return_Freq_Mag_Phase){
    if(length(ToReturn)==3){ToReturn[[4]]=mag[1:ende]
    ToReturn[[5]]=phase[1:ende]
    names(ToReturn)=c("Frequency","Cosine Values","Sine Values","Magnitude","Phase")}
    else{ToReturn[[1]]=f[1:ende]
    ToReturn[[2]]=mag[1:ende]
    ToReturn[[3]]=phase[1:ende]
    names(ToReturn)=c("Frequency","Magnitude","Phase")}}
  if(Plot_Freq_Mag){
    mag2=mag
    mag2[1]=0
    plot(f[1:ende],mag2[1:ende],type="l",xlim = Xlim,ylim = Ylim,ylab = "Magnitude [Unit of TS]",xlab = "Frequency [1/d]")}
  return(ToReturn)}


#' Function for Fast-Fourier-Transformation based filtering of an xts vector
#'
#' The function uses fft to transform the xts vector to the magnitude frequency domain. Limits across frequency and mangitudes can be set to filter certain parts of the signals.
#' @param TS Xts vector to be filtered
#' @param dt The time step in days. (One Minute time step -> 1/24/60)
#' @param freq.range A vector of length m=2*n (n is element of positive integers). Specifiy a frequency range, for which max.mag, min.mag or both should be applied. Multiple ranges can be set at once.
#' @param max.mag A number or a vector of length m/2. Magnitudes higher than max.mag will be erased from the signal.
#' @param min.mag A number or a vector of length m/2. Magnitudes lower than min.mag will be erased from the signal.
#' @param Ylim Vector of length 2. Set the plot's y axis limits.
#' @param Xlim Vector of length 2. Set the plot's x axis limits.
#' @keywords Fourier, fft, filter
#' @return Returns the filtered xts object.
#' @export
#' @examples
#' TS1=DummyTS(days = 100)
#' #This command doesnt filter only plot. It will return the unfiltered xts object.
#' A=FFT_XTS_Filter(x,dt = 1/24/60,Xlim = c(0,5))
#' #Different filter examples. The purpose of these examples is only to show possible filter settings, without any deeper context.
#' A=FFT_XTS_Filter(x,dt = 1/24/60,Xlim = c(0,5),max.mag = 9)
#' A=FFT_XTS_Filter(x,dt = 1/24/60,Xlim = c(0,5),min.mag = 9)
#' A=FFT_XTS_Filter(TS1,dt = 1/24/60,Xlim = c(0,1),freq.range = c(0.3,0.4),max.mag = 7)
#' A=FFT_XTS_Filter(TS1,dt = 1/24/60,Xlim = c(0,1),freq.range = c(0.2,0.4),min.mag = 7)
#' A=FFT_XTS_Filter(TS1,dt = 1/24/60,Xlim = c(0,2),freq.range = c(0,0.4,0.9,1.1),min.mag = c(0,2),max.mag = c(4,99))
#'
FFT_XTS_Filter=function(TS,dt = 1/24/60,max.mag=NULL,min.mag=NULL,freq.range=NULL,Xlim=NULL,Ylim=NULL){
  #CREATE OUR TIME SERIES DATA
  y <- as.numeric(TS)
  time= time(TS)
  index=logical(length(y))
  #Domain setup
  T <- length(time)*dt

  n <- length(time)

  #CREATE OUR FREQUENCY ARRAY
  f <- 1:length(time)/T

  #FOURIER TRANSFORM WORK
  Y <- stats::fft(y)
  mag <- sqrt(Re(Y)^2+Im(Y)^2)*2/n
  offset=Y[1]

  if(!is.null(max.mag)&is.null(min.mag)&is.null(freq.range)){
    index=mag>max.mag
  }
  else if(is.null(max.mag)&!is.null(min.mag)&is.null(freq.range)){
    index=mag<min.mag
  }
  else if(is.null(max.mag)&!is.null(min.mag)&!is.null(freq.range)){
    index1=logical(length(Y))
    index2=logical(length(Y))
    index=logical(length(Y))
    for(i in 1:length(min.mag)){
      index1=mag<min.mag[i]
      index2=f>=freq.range[(2*i-1)]&f<=freq.range[2*i]
      index=index|(index1&index2)
    }
  }
  else if(!is.null(max.mag)&is.null(min.mag)&!is.null(freq.range)){
    index1=logical(length(Y))
    index2=logical(length(Y))
    index=logical(length(Y))
    for(i in 1:length(max.mag)){
      index1=mag>max.mag[i]
      index2=f>=freq.range[(2*i-1)]&f<=freq.range[2*i]
      index=index|(index1&index2)
    }
  }
  else if(!is.null(max.mag)&!is.null(min.mag)&!is.null(freq.range)){
    index1=logical(length(Y))
    index2=logical(length(Y))
    index=logical(length(Y))
    for(i in 1:length(max.mag)){
      index1=!(mag>=min.mag[i]&mag<=max.mag[i])
      index2=f>=freq.range[(2*i-1)]&f<=freq.range[2*i]
      index=index|(index1&index2)
    }
  }
  else if(is.null(max.mag)&is.null(min.mag)&!is.null(freq.range)){
    warning("In order to filter the Frequency Range at least one of both magnitude limits must be set!")
    stop()
  }

  mag[which(index)]=0
  mag[1]=0
  plot(f[1:length(f)/2],mag[1:length(f)/2],type="l",xlim = Xlim,ylim = Ylim,ylab = "Magnitude",xlab = "Frequency")
  Y[which(index)]=0
  Y[1]=offset
  Y_=Re(fft(Y,inverse = T))/length(Y)
  FilteredTS=xts(Y_,order.by = time(TS))
  return(FilteredTS)}


#' Function for decomposing an xts vector
#'
#' The function uses the standard decompose tool of the stats package and applies it to an xts object.
#' For a more detailed description please see the decompose function.
#' @param XTS xts vector to be decomposed.
#' @param decomp A string. "daily","weekly" or "monthly" are possible inputs. Default is daily.
#' @param sampletime A string. "1-min","5-min", "15-min" or "1-hour" are possible inputs. Default is 1-min.
#' @param typ A string. "additive" or "multiplicative" or valid inputs. Default is "additive".
#' @keywords smooth, Filter, decompose
#' @return Decompose object.
#' @export
#' @seealso  \link[stats]{decompose}
#' @examples
#' x=DummyTS(Csd=1.3)
#' A=Decompose4TXS(x)
#' plot(A)
Decompose4TXS=function(XTS,decomp=c("daily","weekly","monthly"),sampletime=c("1-min","5-min", "15-min","1-hour"),typ=c("additive","multiplicative")){
  if(length(decomp)>1){decomp=decomp[1]}
  if(length(sampletime)>1){sampletime=sampletime[1]}
  if(length(typ)>1){typ=typ[1]}

  if(decomp=="daily"){fq=24}
  else if(decomp=="weekly"){fq=24*7}
  else if(decomp=="monthly"){fq=24*30.5}
  else{warning("decomp input is wrong! Choose one of the following: 'daily', 'weekly', 'monthly'!")}

  if(sampletime=="1-min"){st=60}
  else if(sampletime=="5-min"){st=12}
  else if(sampletime=="15-min"){st=4}
  else if(sampletime=="1-hour"){st=1}
  else{warning("sampletime input is wrong! Choose one of the following: '1-min', '5-min', '15-min', '1-hour'!")}
  DC=decompose(ts(XTS,frequency = fq*st),type=typ)
  return(DC)}


#' Function for fixing linear drift in xts data
#'
#' This function calculates linear trends according to the given data. These trends are used to deduct the drift.
#' @param XTS_to_fix xts vector to be corrected from drift.
#' @param Fix_Info an xts vector that specifies the drift times and values. It consists of n+1 lines for n being the number of drift sections.
#' Each line contains a datetime value followed by the drift value and the value it will be corrected to. The first line contains the begin of the first drift and contains the two values 0,0 as there is no drift present at this time.
#' the next line marks the end of the first drift and the two values for correcting it. Further lines use the previous datetime as start for its drift and the current one as end of the drift.
#' @return Xts vector without drift.
#' @export
#' @examples
#' x=DummyTS(Csd=1.3)
#' A=Decompose4TXS(x)
#' plot(A)
#' t=seq(as.POSIXct("2016-01-01 00:00:00",tz="UTC"),as.POSIXct("2016-01-02 00:00:00",tz="UTC"),by="min")
#' v=0.5*sin(seq(0,30*pi,length.out = length(t)))+1
#' drift=seq(0,2,length.out = 720)
#' v[361:1080]=v[361:1080]+drift
#' TS=xts(v,order.by = t)
#' #Fixinfo is probably more efficient to be read out of a txt file.
#' Fixinfo=rbind(cbind(xts(0,order.by = as.POSIXct("2016-01-01 06:00:00")),xts(0,order.by = as.POSIXct("2016-01-01 06:00:00"))),
#'               cbind(xts(3.5,order.by = as.POSIXct("2016-01-01 17:59:00")),xts(1.5,order.by = as.POSIXct("2016-01-01 17:59:00"))))
#' print(Fixinfo)
#'
#' DynPlot(cbind(TS,FixDrift(TS,Fixinfo)),Labels = c("TS with Drift","TS without Drift"))
FixDrift=function(XTS_to_fix,Fix_Info){
  TZXTS=tzone(XTS_to_fix)
  #sets up linear function
  Lin_Func=function(t,t1,t2,y1,y2){
    t=as.numeric(t)
    t1=as.numeric(t1)
    t2=as.numeric(t2)
    y1=as.numeric(y1)
    y2=as.numeric(y2)
    m=(y2-y1)/(t2-t1)
    f_x=m*(t-t1)+y1
    return(f_x)}

  #cycles through Fix_Info
  for(i in 2:dim(Fix_Info)[1]){
    FixInfo_i=as.POSIXct(toString(time(Fix_Info[i,])),tz=TZXTS)
    FixInfo_i_1=as.POSIXct(toString(time(Fix_Info[i-1,])),tz=TZXTS)
    Logical=(time(XTS_to_fix)>=FixInfo_i_1)&(time(XTS_to_fix)<=FixInfo_i)
    y2=Fix_Info[i,][,1]-Fix_Info[i,][,2]
    Corr=Lin_Func(time(XTS_to_fix)[Logical],FixInfo_i_1,FixInfo_i,0,y2)
    XTS_to_fix[Logical]=XTS_to_fix[Logical]-Corr
  }
  return(XTS_to_fix)
}



#' A function for sepereating hydrological summer and winter.
#'
#' The dates are divided by 01.05. and 01.11..
#' @param TS Xts. A time series to be seperated.
#' @return A list. First position is hydrological summer xts. Second postion is hydrological winter xts.
#' @export
#' @examples
#' TS=DummyTS(days=800,Rain_off = TRUE)
#' HH=Separate.To.HyS.HyW(TS)
#' HyS=HH[[1]]
#' HyW=HH[[2]]
#'
Separate.To.HyS.HyW=function(TS){
  Time=as.POSIXct(strftime(time(TS),format = "%m-%d %H:%M:%S"),format = "%m-%d %H:%M:%S",tz = tzone(TS))
  start=as.POSIXct("01.05. 00:00:00",tz=tzone(TS),format="%d.%m. %H:%M:%S")
  end=as.POSIXct("01.11. 00:00:00",tz=tzone(TS),format="%d.%m. %H:%M:%S")
  HyS=Time>start&Time<end
  Separated=list(TS[HyS],TS[!HyS])
  names(Separated)=c("Hydrological Summer","Hydrological Winter")
  return(Separated)}

#' A function for correcting offsets.
#'
#' A certain time period can be selected and corrected.
#' @param TS Xts. A time series to be corrected.
#' @return TS Xts. A time series with the corrections.
#' @export
#' @examples
#' TS=DummyTS(days=30,Rain_off = TRUE)
#' TS=rbind(CutTimeSeries(TS,start(TS)+24*3600,start(TS)+48*3600)+30,CutTimeSeriesDel(TS,start(TS)+24*3600,start(TS)+200*3600))
#' TS_corr=CorrectOffset(TS,start(TS)+24*3600,start(TS)+200*3600,-30)
#' DynPlot(cbind(TS,TS_corr))
#'
CorrectOffset=function(xts,start,end,offset){
  xts.new=CutTimeSeries(xts,start,end)
  xts.new=xts.new+offset
  xts.old=CutTimeSeriesDel(xts,start,end)
  xts.final<-rbind(xts.new,xts.old)
  return(xts.final)
}

