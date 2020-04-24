#' Dissolved Oxygen Saturation Function
#'
#' This function uses the APHA-Method to calculate maximum dissolved oxygen concentration in a open freshwater, from given water-temperature and water-surface-level atmospheric pressure.\cr
#' For more details see: American Public Health Association Eaton, A.D., 2006. Standard methods for the examination of water & wastewater, 21st ed. 2005, centennial ed., [Nachdr.]. ed. American Public Health Association, Washington, DC.\cr
#' @param Temp Water-Temperature in °C
#' @param Pressure Water-Surface-Level Atmospheric Pressure in kPa 
#' @return A numeric vector\cr
#' @keywords APHA, Oxygen Saturation
#' @examples
#'pressure=rnorm(n = 21,mean = 101.32,sd = 2)
#'temperature=sin(x=seq(4,10,by = 0.3))*5+15
#'cs=CS(temperature,pressure)
#' @export
CS=function (Temp,Pressure){
    atm=Pressure/101.325 
    degK <- (Temp + 273.15)
    cs=exp(((-139.34411) + ((157570.1)/degK) - ((66423080)/degK^2) + 
              ((1.2438e+10)/degK^3) - ((862194900000)/degK^4)))
    pwv=exp(11.8571-(3840.7/7)-(21696/degK^2))
    theta=0.000975-(1.426+((10^-5)*Temp))+(6.436*((10^-8)*Temp^2))
    cs*atm*(((1-pwv/atm)*(1-theta*atm))/((1-pwv)*(1-theta)))
}


#Nighttime-Slope-Methode

#'This function calculates Ecosystem Respiration (ER) [mgO2/L*d] and Reaeration Coefficient (kO2) [1/d] from a timeseries of measured DO an calculated DO-Saturation on a daily basis by Nighttime-Slope Method.\cr
#'First step is calculation of sunset and sunrise. If your river is shaded true sunset and sunrise (solar irradiance close to zero) might differ from astronomical times. Use the setdel and risedel parameter to adjust delay (in seconds).\cr
#'It is recommended to determine true sunset by modeling or measurement to improve accuracy of calculation.\cr
#'Next Step is to identify nighttime DO-Minima.\cr
#'Then temporal DO-Gradient (??DO/??t mgO2/L*d) and Saturationdeficit (DOsat-DO) in the time between sunset and DO-Minima is fitted by a linear function. Slope results kO2 and Intersept ER.\cr
#'The timeseries has to be continuous, evenly spaced and without NA's.\cr 
#'Function uses NOAA sunrise sunset calculator.\cr
#'For more details see: Hornberger, G.M., Kelly, M.G., 1975. Atmospheric reaeration in a river using productivity analysis. Journal of the Environmental Engineering Division 101, 729-739.
#'@param data xts-timeseries of DO (mg/L) DOsat (mg/L) and Temperature (°C) (data must be in given order)
#'@param Start start time as "YYYY/MM/DD"
#'@param End end time as "YYYY/MM/DD"
#'@param tstep timestep of the timeseries (minutes)
#'@param tz timezone (sunrise sunset calculation) as "UTC"
#'@param lat latitude in decimal degrees
#'@param long longitude in decimal degrees
#'@param plot boolean to plot resulting linear function (creates a plot for every day)
#'@param plot_path path to save the plots, name is generated from date, as "C:/Users/Test/"
#'@param setdel timedelay of sunset (seconds)
#'@param risedel timedelay of sunrise (seconds)
#'@return A dataframe containing INT (ER), Slp (kO2), Rsquared of the linear function, Date, Mean Temperature for timeframe, and number of datapoints used for calculation\cr
#'@keywords APHA, Oxygen Saturation
#'@examples
#'First create a sample timeseries
#'datetime is our index 
#'datetime=seq(from=as.POSIXct("2019-06-14 12:00"), by=60*5, to=as.POSIXct("2019-06-17 12:00"))
#'Dirunal DO, and Temperature is represented by a Sine wave
#'do=sin(x=seq(1,19,by = (18/(length(datetime)-1))))*3+10
#'temperature=sin(x=seq(1,19,by = (18/(length(datetime)-1))))*5+12
#'Atmospheric Pressure kPa is a normal distributed function with standard pressure (101.32 kPa) as mean and sd=1.
#'pressure=rnorm(length(datetime),101.32,1)
#'Calculate DO-Saturation from Temperature and Atmospheric Pressure 
#'cs=CS(temperature,pressure)
#'Combine the index (datetime), DO, DOSat(cs) and Temperature to a xts Timeseries
#'ts=xts(x = cbind(do,cs,temperature),order.by = datetime)
#'Calculate ER and kO2 with nighttime slope method. 
#'nighttime.slope(data = ts,Start = "2019/06/14",End = "2019/06/16",tstep = 5,tz = "UTC-1",lat=51.050407,long= 13.737262)
#' @export
nighttime.slope=function(data,setdel,risedel,Start,End,tstep,tz,lat,long,plot,plot_path){
    start <- as.Date(Start,format="%Y/%m/%d")
    end   <- as.Date(End,format="%Y/%m/%d")
    Date <- start
    df=data.frame()
    while (Date < end){
      setrise=sunrise.set(lat=lat,long=long,date=Date,timezone = tz,num.days = 2)
      set=setrise[1,2]+setdel
      rise=setrise[2,1]+risedel
      o2nsm=CutTimeSeries(data[,c(1,2)],start=set,end=rise)
      dfo2nsm=data.frame(time=index(o2nsm),coredata(o2nsm))
      tmin=dfo2nsm[(dfo2nsm[,2]==min(dfo2nsm[,2],na.rm=T)),1]
      o2nsm=CutTimeSeries(o2nsm,start=set,end=tmin)
      o2nsm=period.apply(o2nsm,endpoints(o2nsm,on = "mins",k=tstep),mean)
      o2nsm$diff=o2nsm[,2]-o2nsm[,1]
      o2nsm$grad=(diff(o2nsm[,1],lag=1))*(60*24)/tstep
      
      a=as.numeric(coredata(o2nsm[,3]))
      b=as.numeric(coredata(o2nsm[,4]))
      
      if (is.na(all(a==b))) {
        Date <- Date + 1
        next
      }
      
      ablm=(lm(b~a))
      
      temp=mean(CutTimeSeries(data[,3],start=set,end=tmin))
      
      tmp=cbind(as.numeric(coef(ablm)[1]),
                as.numeric(coef(ablm)[2]),
                as.numeric(summary(ablm)$r.squared),
                Date,
                temp,
                length(a))
      df=rbind(df,tmp)
      Date <- Date + 1
      
      if (plot==TRUE) {
        emf(file=paste(c(plot_path,Date,".emf"),collapse = ""),width = 5,height = 5)
        plot(a,b,xlab="Sättigungsdefizit [mg/L]",ylab="??DO/??t [mg O2/L*d]")
        abline(ablm)
        mtext(text = paste("R=",round(summary(ablm)$r.squared,digits=3),"Int(ER)=",round(coef(ablm)[1],digits=3),"mgO2/L*d","Slp(kO2)=",round(coef(ablm)[2],digits=3),"1/d"),adj = 1,side = 3)
        dev.off()
      }
      
       
      
    }
    names(df)=c("Int","Slp","Rsquared","Date","Temp","Length")
    df[,4]=as.POSIXct.Date(as.Date(df[,4]))
    return(df)
  }



datetime=seq(from=as.POSIXct("2019-06-14 12:00"), by=60*5, to=as.POSIXct("2019-06-17 12:00"))
do=sin(x=seq(1,19,by = (18/(length(datetime)-1))))*3+10
temperature=sin(x=seq(1,19,by = (18/(length(datetime)-1))))*5+12
pressure=rnorm(length(datetime),101.32,1)
cs=CS(temperature,pressure)
ts=xts(x = cbind(do,cs,temperature),order.by = datetime)
autoplot.zoo(ts)
t=nighttime.slope(data = ts,Start = "2019/06/14",End = "2019/06/16",tstep = 5,
                tz = "UTC-1",lat=51.050407,long= 13.737262,plot=TRUE,
                plot_path="C:/Users/Tobias/Documents/Universität/Masterarbeit/TEST/")
