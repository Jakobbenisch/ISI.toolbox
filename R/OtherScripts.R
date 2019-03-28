
#' Repair Function
#'
#' Fixes periodicity necessary for plotting with dygraphs by deleting NA values in time stamp and add na values for an even time step.
#' @param vec Xts vector to be ´fixed.
#' @return Xts vector with uniform time step.
#' @keywords fix, periodicity
#' @export
#' @examples
#' TimeSeries=FixPeriodicity(TimeSeries)
#'
FixPeriodicity = function(TS){
  Time=time(TS)
  Time=na.omit(Time)
  time=as.numeric(Time)
  difftime=diff(time)
  per=min(difftime)
  if (per==0){per=min(difftime[difftime>0]);warning("The time series includes repeating datetimes!")}
  tout=seq(Time[1],Time[length(Time)],by=per)
  TS2=xts(rep(NA,length(tout)),order.by = tout)
  new=cbind(TS,TS2)
  TS=subset(new,select=-dim(new)[2])
  return(TS)
}

#' Time Conversion Function
#'
#' Converts a time to number.
#' @param datetime Xts vector or POSIXct vector.
#' @return Time in seconds from the 01.01.1970 01:00 as numeric vector.
#' @export
#' @examples
#' TimeIndexinNumbers=timetonumber(TimeSeries)
#'
timetonumber = function(datetime){
  if(class(datetime)[1]=="xts"){
    num=as.numeric(time(datetime))}
  else if(class(datetime)[1]=="POSIXct"){
    num=as.numeric(datetime)}
  else{print("Wrong input class for timetonumber function!")}
  return(num)
}


#' Event Duration Count Function
#'
#' Calculates the duration and quantity of events exceeding the limit-value.
#' @param TS Xts vector to be analysed.
#' @param limit Decimal. Value to be exceeded.
#' @param maxgapsize Maximal tolerable size of temporal gap in data (in hours). If the gap is larger the event will considered to have ended.
#' @return Returns numerical vector containing duration (in seconds) of exceedance events.
#' @export
#' @examples
#' TS=DummyTS()
#' vec=Histo(TS,limit=30)
#'
Histo = function(TS,limit,maxgapsize=24){
  switch=0
  indstart=0
  indstop=0
  D=c()
  maxgapsize=maxgapsize*3600
  Values=as.numeric(TS)
  Times=time(TS)
  TRFA=Values>limit
  indexes=which(TRFA %in% c(TRUE))
  indexes=append(indexes,(indexes[length(indexes)]+10))
  if(length(indexes)==0){D=0;return(D)}
  for(i in 1:(length(indexes)-1)){
    if(indexes[i+1]-indexes[i]==1 & switch==0){indstart=indexes[i]; switch=1}
    else if((indexes[i+1]-indexes[i])==1 & switch==1){}
    else if(switch==1 & (indexes[i+1]-indexes[i])!=1 ){indstop=indexes[i];
    D=append(D,(timetonumber(Times[indstop])-timetonumber(Times[indstart])));
    switch=0}
    else if(switch==1 & (timetonumber(Times[indexes[i+1]])-timetonumber(Times[indexes[i]]))>maxgapsize ){indstop=indexes[i];
    D=append(D,(timetonumber(Times[indstop])-timetonumber(Times[indstart])));
    switch=0}
    else{}}
  return(D)
}

#' Concentration Duration Occurance Function
#'
#' Calculates the occurance of exceedance for given concentration limit values and temporal (duration) classes.
#' @param TS Xts vector to be analysed.
#' @param limits Numeric vector. Values to be exceeded.
#' @param timeclasses Numeric vector. Values to be set in days.
#' @param maxgapsize Maximal tolerable size of temporal gap in data (in hours). If the gap is larger the event will considered to have ended.
#' @param old Boolean. If TRUE the occurance as an event will be returned. If FALSE the occurance time based on the respective time class will be returned.
#' @return Returns matrix (top line time in seconds, first row intensity, other occurance).
#' @export
#' @examples
#' TS=DummyTS()
#' M=KDH(TS,limits=c(10,20,30,40,50,60),timeclasses=c(0.5/48,1/48.,1/24.,1/12.,1/8.,1/6.,1/3.,1))
#' print(M)
KDH =function(TS,limits,timeclasses,maxgapsize=24,old=FALSE){
  timeclasses=timeclasses*3600*24
  X=matrix(NA, (length(limits)+1), (length(timeclasses)+1))
  for(i in 1:length(limits)){d=Histo(TS,limits[i],maxgapsize)
    X[(i+1),][1]=limits[i]
    for( p in 1:length(timeclasses)){
      X[1,][p+1]=timeclasses[p]
      if(old){
        X[(i+1),][p+1]=sum(as.numeric(d>timeclasses[p]))
      }
      else{
        amount=d[d>timeclasses[p]]
        if(length(amount)==0){
          freq=0
        }
        else{
          freq=c()
          for(q in amount){
            freq=append(freq,q/timeclasses[p])
          }
        }
        X[(i+1),][p+1]=sum(freq)
      }
    }
  }
  return(X)
}



#' Add NA Values Function
#'
#' The function unifies the time step by adding NA values to gaps.
#' @param VEC Xts vector.
#' @param by Time step of xts vector in seconds. Default is one minute (60 s).
#' @return Xts vector with uniform time step.
#' @export
#' @examples
#' TS=DummyTS()
#' TS2=CutTimeSeriesDel(TS,start="2016-01-02 12:00:00",end="2016-01-02 22:00:00")
#' TS2=Add.NA(TS2)
#' DynPlot(cbind(TS,TS2))
#'
Add.NA=function(VEC,by=60){
  # erstellt einen Zeitvektor ohne Lücken
  start=min(time(VEC))
  end=max(time(VEC))
  time_vec=seq(start,end,by=by)
  # Dummy xts wird erstellt
  dummy=xts(rep(0,length(time_vec)),order.by = time_vec)
  # Dummy und Xts werden zusammengeführt. (Dabei der NA's in die Lücken geschrieben)
  VEC=cbind(dummy,VEC)
  # Jetzt wird der Dummyvektor wieder gelöscht
  VEC=VEC[,2:(dim(VEC)[2])]
  return(VEC)}






#' Rain Analysing Function
#'
#' This function calculates the intensity of rain events.
#' @param Rain_XTS Xts vector containing the rain events. Ideally no NA values are present and the time step is even. (NA could be replaced with zeros if present. To unify the the time step the Add.Na function could be used.)
#' @param ignore.gap.length As a rain event does not actually end once a single rainheight within the event is zero, this value can be specified to set a minimum number of consecutive values that have to be zero so that the event will have ended.
#' So for a 1 min time step in the rainfall xts object, the value 60 would resemble an hour.
#' @return A list.\cr
#' [[1]] -> Numerical vector containing the intensities.\cr
#' [[2]] -> Numerical vector containing the total precipitation.\cr
#' [[3]] -> Numerical vector containing the start of the event.\cr
#' [[4]] -> Numerical vector containing the end of the event.\cr
#' \cr
#' Ps:\cr
#' The start and the end must be converted to a datetime object as follows.\cr
#' start_as_datetime=as.POSIXct(start_numeric,origin = "1970-01-01")
#'
#' @export
#' @examples
#' RainTimeseries=RainGen()
#' Intensity_Height=Rain.Intensity(RainTimeseries,ignore.gap.length=20)
#'
Rain.Intensity=function(Rain_XTS,ignore.gap.length=60){
  Dur<<-c()
  Height<<-c()
  Starttime<<-c()
  Endtime<<-c()
  event<<-FALSE
  sumheigth<<-c()
  Value<-as.numeric(Rain_XTS)
  Temp<-time(Rain_XTS)
  start<-0
  Count=function(i){
    Time<-Temp[i]
    x<-Value[i]
    eventprev<-event
    if(x>0){event<<-TRUE; sumheigth<<-c(sumheigth,x)}
    else{if(i<(length(Value)-ignore.gap.length)){if(!any(Value[i:(i+ignore.gap.length)]>0)){event<<-FALSE}}
      else{event<<-FALSE}}
    if(!eventprev&event){sumheigth<<-c(sumheigth[length(sumheigth)]);start<<-i}
    else if(eventprev&!event){Dur<<-c(Dur,(as.numeric(Time)-as.numeric(Temp[start]))/3600); Height<<-c(Height,sum(sumheigth)); Starttime<<-c(Starttime,Temp[start]); Endtime<<-c(Endtime,Time)}}
  sapply(1:length(Rain_XTS),FUN=Count)
  return(list(Height/Dur,Height,Starttime,Endtime))}


#' A function for dealing with NA values in data.
#'
#' This function offers a variety of options to fill gaps in the data.
#' It incoorporates linear interpolation and spline interpolation which can be applied to xts data sets of any dimension.\cr
#' The regression option is designed for xts data sets of dimension (m,2), where the first column represents x in the second y.
#' A polynomial, potency, exponential or logarithmic function will be fitted to x and y and will be used to fill the gaps in y.\cr
#' The option diurnal_pattern requires an xts data set of dimension (m,1) for which a Fourier Series will be fitted und used to fill the gaps.\cr
#' best_fit is another method available. It compares days with gaps to days without gaps and uses the best fit, which is then equalized, to fill the gap. (The algorithms uses apply.daily)
#' Xts data sets may have the dimension (m,n), depending on the cn_exclude_last_row paramter n or (n-1) columns will be used for comparison (by correlation).
#' Either way the last columns gaps will be filled.
#' @param XTS Xts vetor. Data set with NA values to be filled. Depending in the use parameter different dimensions are required (see description).
#' @param use String. Specifies the method for filling the NA values. \cr
#' Possible inputs are: "linear", "spline", "regression", "diurnal_pattern" or "best_fit"
#' @param maxgap Integer. Used by the linear and by the spline method. Specifies the maximum gap size that will be filled.
#' @param Regkey String. Used by the regression method. Specifies the function that will be fitted. \cr
#' Possible inputs are: "POT", "POL", "EXP" or "LOG"
#' @param Poldegree Integer. Specifies the degree of the polynomial function to be fitted. If Regkey=="POL" and use=="regression".
#' @param Fourierdegree Integer. Specifies the degree of the Fourier Series to be fitted. If use=="diurnal_pattern".
#' @param Gaprange Numeric vector of length two. Only used by the best_fit method. The vector specifies the size of gaps that will be filled. The 2 values represent gap size percentages of a day. \cr
#' E.g.: c(0,60) Gaps between the sizes of 0\% and 60\% of a day will be filled.
#' @param Allowed_NAs_in_Fill_Day Number. Only used by the best_fit method. Specifies the allowed percent of NA values in the fill day.
#' @param Min_R Number. Only used by the best_fit method. This values represents the minimal quality of the fit. If it is not reached a day won't be considered for filling a gap. \cr
#' If the dimension of XTS is greater than (m,2) the number will be compared to the average fit quality.
#' @param cn_exclude_last_row Boolean. Only used by the best_fit method. If TRUE the last row won't be used for comparison. (Must be TRUE if dim(XTS)==(m,1))
#' @param Regression100Points Boolean. Only used by the regression method. If TRUE the regression model will be calculated with 100 points only instead of all available data.\cr
#' The advantage is that the fit will be a better representation for areas with little measurement point density.
#' @seealso \link[zoo]{na.approx}, \link[zoo]{na.spline}, \link[ISI.Toolbox]{Check_Fill_Quality}
#' @return If use != "best_fit" the na - filled - xts vector will be returned. \cr
#' Otherwise a list will be returned:
#' [[1]] -> na - filled - xts vector
#' [[2]] -> a numerical vector containing the information of the fill quality
#' [[3]] -> a list containing numericial vectors with the day to be filled and the fill day.
#' @export
#' @examples
#'TS=DummyTS(days=150)
#'for(i in 1:20){
#'  start=as.integer(runif(1,1,length(TS)))
#'  len=as.integer(runif(1,1,1440))
#'  TS[start:(start+len)]=NA
#'}
#'TS_filled=FillNaValues(TS,use = "best_fit",Gaprange = c(0,80),cn_exclude_last_row = FALSE)
#'
#Gaprange of non_NA_values of a day that will be filled by nearest neighbour in %, Allowed_NAs_in_Fill_Day in %
FillNaValues=function(XTS,use="linear",maxgap=4,Regkey="POT",Poldegree=2,Fourierdegree=5,Gaprange=c(0,60),Allowed_NAs_in_Fill_Day=0,Min_R=0.85,cn_exclude_last_row=TRUE,Regression100Points=FALSE){
  #input must be single column
  if(use=="linear"){XTS=na.approx(XTS,maxgap = maxgap)}
  else if(use=="spline"){XTS=na.spline(XTS,maxgap = maxgap)}

  #The input xts will have 2 columns the first will be complete (x) and the second will have gaps (y) (for example a datetime,h,Q)
  else if(use=="regression"){
    if(Regression100Points){
      FixedPoints=function(TS){
        XTS_clean=na.omit(TS)
        DataF=data.frame(XTS_clean[,1],XTS_clean[,2])
        DataF=DataF[order(DataF[1]),]
        x=as.numeric(DataF[,1])
        y=as.numeric(DataF[,2])
        min=min(x,na.rm = TRUE)
        max=max(x,na.rm=TRUE)

        X=seq(min,max,by=(max-min)/100)
        newx=c()
        newy=c()
        for(i in 2:length(X)){
          logi=x>X[i-1]&x<X[i]
          if(!all(!logi)){
            newx=c(newx,mean(x[logi]))
            newy=c(newy,mean(y[logi]))}}
        return(data.frame(newx=newx,newy=newy))
      }
      XTS_clean=FixedPoints(XTS)
    }
    else{
      XTS_clean=na.omit(XTS)
    }
    x=as.numeric(XTS_clean[,1]);y=as.numeric(XTS_clean[,2])
    if(Regkey=="POL"){
      RegObj=POLRegression(x,y,order=Poldegree,summary = FALSE)
      plot(x,y);lines(x,predict(RegObj),col=2)
      gaps_logical=is.na(XTS[,2])
      x_gaps_logical=!is.na(XTS[,1][gaps_logical])
      Fill_values=as.numeric(hQ(XTS[,1][gaps_logical][x_gaps_logical],RegObj,Type=Regkey))
      XTS[,2][which(gaps_logical)][which(x_gaps_logical)]<-Fill_values}
    else if(Regkey=="POT"){
      RegObj=POTRegression(x,y,summary = FALSE)
      plot(as.numeric(XTS[,1]),as.numeric(XTS[,2]));points(x,y,col="blue");lines(x,exp(predict(RegObj)),col=2)
      gaps_logical=is.na(XTS[,2])
      x_gaps_logical=!is.na(XTS[,1][gaps_logical])
      Fill_values=as.numeric(hQ(XTS[,1][gaps_logical][x_gaps_logical],RegObj,Type=Regkey))
      XTS[,2][which(gaps_logical)][which(x_gaps_logical)]<-Fill_values}
    else if(Regkey=="EXP"){
      RegObj=EXPRegression(x,y,summary = FALSE)
      plot(x,y);lines(x,exp(predict(RegObj)),col=2)
      gaps_logical=is.na(XTS[,2])
      x_gaps_logical=!is.na(XTS[,1][gaps_logical])
      Fill_values=as.numeric(hQ(XTS[,1][gaps_logical][x_gaps_logical],RegObj,Type=Regkey))
      XTS[,2][which(gaps_logical)][which(x_gaps_logical)]<-Fill_values}
    else if(Regkey=="LOG"){
      RegObj=LOGRegression(x,y,summary = FALSE)
      plot(x,y);lines(x,predict(RegObj),col=2)
      gaps_logical=is.na(XTS[,2])
      x_gaps_logical=!is.na(XTS[,1][gaps_logical])
      Fill_values=as.numeric(hQ(XTS[,1][gaps_logical][x_gaps_logical],RegObj,Type=Regkey))
      XTS[,2][which(gaps_logical)][which(x_gaps_logical)]<-Fill_values}
    else{warning("Wrong input for parameter Regkey! Possible inputs: c('POL','POT','EXP','LOG')")}}
  else if(use=="diurnal_pattern"){
    XTS_clean=na.omit(XTS)
    RegObj=FourierRegression(XTS_clean,orderoftransformation=Fourierdegree,summary = FALSE)
    Fill_values=FourierPattern(XTS,RegObj)
    gaps_logical=is.na(XTS)
    XTS[which(gaps_logical)]<-Fill_values[which(gaps_logical)]}
  else if(use=="best_fit"){
    Info_R=numeric()
    Info_days=list()
    XTS2=XTS
    Dates=lubridate::date(XTS2)

    BestMatch=function(day_with_NA){

      Bool_Multiline=function(XTS,Boolean){
        for(i in 1:dim(XTS)[2]){
          XTS[,i][Boolean[,i]]=NA}
        return(na.omit(XTS))}

        Comp=function(day_without_NA,day_with_NA){
          fillvec=day_without_NA[,dim(day_without_NA)[2]]
          if(cn_exclude_last_row){day_without_NA=subset(day_without_NA,select=-dim(day_without_NA)[2])}
          else{}
          if(all(dim(day_without_NA)==dim(day_with_NA))){
            lo1=!is.na(day_with_NA)
            lo2=!is.na(day_without_NA)
            NA_index=lo1&lo2
            if(identical(day_with_NA,day_without_NA)){R=-1}
            else if(sum(is.na(fillvec))/length(fillvec)>(Allowed_NAs_in_Fill_Day/100)){R=-1}
            else{R=cor(Bool_Multiline(day_without_NA,!NA_index),Bool_Multiline(day_with_NA,!NA_index))}}
          else{R=-1}

          if(!is.null(dim(R))){R=mean(diag(R))}
          return(R)}

      R=as.numeric(apply.daily(XTS,FUN = Comp,day_with_NA=day_with_NA))
      R_max=max(R,na.rm=T)

      if(R_max>Min_R){
        Info_R<<-append(Info_R,R_max)
        out=na.omit(unique(Dates)[R==R_max])}
      else{out=R_max}

      return(out)}

    Match_Curves=function(XTS,XTS_gapless){
      #x=as.numeric(time(XTS_gapless));x=x-min(x)
      y=as.numeric(XTS_gapless)-as.numeric(XTS)
      Difference=xts(y,order.by=time(XTS_gapless))

      #wenn vorn na und hinten ni und so weiter
      if(is.na(Difference[length(Difference)])){
        Ind=is.na(Difference)
        p=which(!Ind)
        pos=p[length(p)]
        Difference[which(Ind)]=rep(Difference[pos],sum(Ind))}
      else if(is.na(Difference[1])){
        Ind=is.na(Difference)
        p=which(!Ind)
        pos=p[1]
        Difference[which(Ind)]=rep(Difference[pos],sum(Ind))
      }
      else{D=na.approx(Difference)}
      D=na.approx(Difference)
      return(D)}

    FillGap=function(current_day,match_day){
      #simple replacement
      repl_index=Dates==match_day
      faulty_day_index=Dates==current_day
      only_NA_index=is.na(XTS[faulty_day_index][,dim(XTS)[2]])
      XTS_repl=XTS[repl_index][,dim(XTS)[2]]
      Difference=Match_Curves(XTS2[,dim(XTS)[2]][faulty_day_index],XTS_repl)
      Replace=as.numeric(XTS_repl[only_NA_index])-as.numeric(Difference[only_NA_index]);Replace[Replace<0]=0
      XTS2[,dim(XTS)[2]][which(faulty_day_index)][which(only_NA_index)]<<-Replace}

    Gap_Yes_or_No=function(day_of_XTS){
      if(any(is.na(day_of_XTS))){
        per_gap=table(is.na(day_of_XTS))["TRUE"]/length(day_of_XTS)*100
        if(Gaprange[1]<per_gap&per_gap<Gaprange[2]){
          current_Day=unique(lubridate::date(day_of_XTS))
          if(cn_exclude_last_row){match_Day=BestMatch(subset(day_of_XTS,select=-dim(day_of_XTS)[2]))}
          else{match_Day=BestMatch(day_of_XTS)}
          if(is.Date(match_Day)){
            FillGap(current_Day,match_Day)
            daypair=c(current_Day,match_Day)
            Info_days<<-append(Info_days,list(daypair))}
          else{}
        }
      }
    }

    apply.daily(XTS,FUN = Gap_Yes_or_No)
    XTS=list(XTS2,Info_R,Info_days)
  }
  else{warning("Wrong input for parameter use! Possible inputs: c('linear','spline','regression','diurnal_pattern','best_fit')")}
  return(XTS)
}


#' A function for checking the FillNaValues Output.
#'
#' This function illustrates the fill quality of the best_fit method of the FillNaValues. A day at the time.
#' @param Fillcount Integer. This value indexes all filled days.
#' @param Before_Fill Xts Vector. Original data, before using the fill function.
#' @param FillReturnList List (returned from the FillNaValues function).
#' @return A dynamic plot and the fit quality expressed as average coefficient of determination.
#' @seealso \link[ISI.Toolbox]{FillNaValues}
#' @export
#' @examples
#'TS=DummyTS(days=150)
#'for(i in 1:20){
#'  start=as.integer(runif(1,1,length(TS)))
#'  len=as.integer(runif(1,1,1440))
#'  TS[start:(start+len)]=NA
#'}
#'TS_filled=FillNaValues(TS,use = "best_fit",Gaprange = c(0,80),cn_exclude_last_row = FALSE)
#'Check_Fill_Quality(1,TS,TS_filled)
#'
Check_Fill_Quality=function(Fillcount,Before_Fill,FillReturnList){
  originday=format(FillReturnList[[3]][[Fillcount]][1],"%Y-%m-%d %H:%M")
  originday_p1=format(FillReturnList[[3]][[Fillcount]][1]+1,"%Y-%m-%d %H:%M")
  fillday=format(FillReturnList[[3]][[Fillcount]][2],"%Y-%m-%d %H:%M")
  fillday_p1=format(FillReturnList[[3]][[Fillcount]][2]+1,"%Y-%m-%d %H:%M")
  P1=CutTimeSeries(Before_Fill,originday,originday_p1)
  P2=CutTimeSeries(FillReturnList[[1]],originday,originday_p1)
  P3=CutTimeSeries(FillReturnList[[1]],fillday,fillday_p1)
  P3=xts(as.data.frame(P3),order.by = time(P2))
  row=dim(P3)[2]
  Plt=DynPlot(cbind(P2[,row],P1[,row],P3[,row]),Labels = c("Filled","Original_Data","Fill_Day"))
  print(paste("Average R:",toString(FillReturnList[[2]][[Fillcount]])))
  return(Plt)}

#' A function changing Timezones.
#'
#' This function enables you to change timezones. Thereby getting rid of winter and summer daylight savings.
#' @param XTS Xts vector. The timezone does not matter. It will be adjusted to the TZorigin.
#' @param TZorigin TimeZone string. Tzone(XTS) will be adjusted to the TZorigin.
#' @param TZdest TimeZone string. TimeZone that the xts will be converted to.
#' @return Xts with the desired timezone.
#' @export
#' @examples
#' TS=DummyTS(days=365,Rain_off = TRUE,timezone = "Europe/Berlin")
#' print(tzone(TS))
#' TS=TimeZoneConverter(TS)
#' print(tzone(TS))
TimeZoneConverter=function(XTS,TZorigin="Europe/Berlin",TZdest="UTC"){
  Time=as.character(time(XTS))
  dat=data.frame(XTS,row.names=1:dim(XTS)[1])
  XTS=xts(dat,order.by = as.POSIXct(Time,tz=TZorigin))
  return(as.xts(XTS,tz=TZdest))
}

#' A function for detection Discharge Events.
#'
#' This functionwas designed to detect increase in flow. It finds events via 2 criteria, the exceedance of a limit value and the duration of exceedance. Both should be provided by the user.
#' @param q_xts Xts vector. Timeseries containing flow signal.
#' @param quantile_limit Numeric value. Is used for autosetting the limit value for event detection via the quantile of the flow distribution. The quantile_limit must be between 1 and 99.
#' @param fixed_limit Numeric value. If desired a fixed flow limit can be set for event detection. (quantile_limit must be set to NULL in that case.)
#' @param min_duration_in_time_steps Integer. The value speciefies the minimal number of timesteps that the limit value needs to be exceeded so that it will be recognized as event.
#' @param padding_in_time_steps Integer. Padding for detected events, so that the begiing and end will also be captured.
#' @param padding_ratio Numeric value. Must be between 0.05 and 0.95. Specifies the padding backward in time (padding_in_time_steps x padding_ratio) and forward in time (padding_in_time_steps x (1-padding_ratio)).
#' @param inpute_NA Boolean. If True the original time series containing NA values for all non event time steps will be returned. Else the original time series will be cut so that it only contains events.
#' @param na_approx_max_gap Integer. Gaps in flow signal up to size na_approx_max_gap will be interpolated linearly. Greater gaps will be filled with zero.
#' @return Xts with the desired timezone.
#' @export
#' @examples
#' TS=DummyTS(days=50)
#' Events=EventDetector(TS,fixed_limit = 35,quantile_limit = NULL)
#' DynPlot(cbind(TS,Events))
EventDetector=function(q_xts,quantile_limit=90,fixed_limit=NULL,min_duration_in_time_steps=180,padding_in_time_steps=600,padding_ratio=1/3,inpute_NA=TRUE,na_approx_max_gap=5){
  df_numeric=as.numeric(q_xts)
  if (!is.null(quantile_limit) & is.null(fixed_limit)){
    if(quantile_limit<1 | quantile_limit>99){stop("The quantile_limit must be between 1 and 99.")}
    q_lim=as.numeric(quantile(df_numeric,quantile_limit/100))
    print(paste("Determined quantile limit:",toString(q_lim)))
  }else if(is.null(quantile_limit) & !is.null(fixed_limit)){
    q_lim=fixed_limit
  }else if( !is.null(quantile_limit) & !is.null(fixed_limit)){
    stop("Please only provide quantile_limit or fixed_limit.")
  }else if(is.null(quantile_limit) & is.null(fixed_limit)){
    stop("Please provide either quantile_limit or fixed_limit.")
  }

  time_vec=time(q_xts)
  min_ts_in_min=min(diff(time_vec))
  df_uni=Add.NA( q_xts,by=60*min_ts_in_min)
  #gaps up to size na_approx_max_gap will be interpolated linearly
  df_uni=na.approx(df_uni,maxgap = na_approx_max_gap)
  #greater gaps will be filled with zero
  df_uni=na.fill(df_uni,0)

  boo=df_numeric>q_lim
  ev_boo=rollapply(boo,width=min_duration_in_time_steps,FUN=all)
  ev_boo=c(ev_boo,rep(FALSE,min_duration_in_time_steps-1))

  extend_true=function(boolean,up=FALSE){
    len=length(boolean)
    if(up){return(c(boolean[2:len],FALSE)|boolean)}
    else{return(c(FALSE,boolean[1:(len-1)])|boolean)}
  }

  for (i in 1:(min_duration_in_time_steps-1)){
    ev_boo=extend_true(ev_boo)
  }

  for (up in 1:as.integer(min_duration_in_time_steps/2+padding_in_time_steps*padding_ratio)){
    ev_boo=extend_true(ev_boo,up=TRUE)
  }

  for (down in 1:as.integer(min_duration_in_time_steps/2+padding_in_time_steps*(1-padding_ratio))){
    ev_boo=extend_true(ev_boo)
  }

  if (inpute_NA){
    df_uni[!ev_boo]=NA
    return(df_uni)
  }
  else{
    return (df_uni[!ev_boo])
  }
}
