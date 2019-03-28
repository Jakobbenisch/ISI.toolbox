
#' A Generating Function
#'
#' Generates an xts with an desired Dynamik. The default setting will create a fairly realistic discharge pattern.
#' The function achieves this buy using a Fourier Series of order 5 (the period is set to 2 pi), the Raingenerator function and multiple linear storage models.
#' For future reference slow catchments refer to the ones that do have a long rainfall concentration time, thus fast catchments do have a short rainfall concentration time.
#' @param days The number of days to be generated
#' @param SampleTime_min The logging interval in minutes. (has been tested for 1,2 and 5 minutes)
#' @param sd_days The standard deviation from day to day (the higher the greater the difference from day to day)
#' @param Cunif The uniform distributed noise to be set on the signal.
#' @param Csd The normal distributed noise to be set on the signal.
#' @param FouP A vector containing 11 values for the Fourier Series starting with a_0. The next coefficients are the amplitudes of cosine and sine terms.(every even index number is a cosine amplitude...)
#' @param AnualPrecipitation_mm The anual rainfall height in mm.
#' @param Anual_Days_of_P The number of days in rains in a year.
#' @param Beta A list containing 4 to 6 values of vectors that describe the probabilities of rain, for the rain's disaggregation. For a more detailed description see the Raingenerator() function.
#' @param Lin_Stor_Catchment_Area The total area used by the linear storage model. It is spread according to the relation between the slow and fast catchments set in the Lin_Stor_Number_of_Catchments_slow_fast vector.
#' @param Lin_Stor_Number_of_Catchments_slow_fast Vector of length 2. The first position represents the number of slow catchments and the second the fast ones.
#' @param Lin_Stor_Average_Lambda_slow_fast Vector of length 2. Specifying the average retardation coefficient of the slow and fast linear storage models. (First slow then fast.)
#' @param Rain_off Logical value. If TRUE the Raingenerator and the linear storage models won't be utilized.
#' @param timezone Timezone of returned xts.
#' @seealso \link[ISI.Toolbox]{Raingenerator}
#' @return Returns a generated xts vector.
#' @keywords generate, create, synthetic
#' @export
#' @examples
#' TS=DummyTS(days=30)
#' DynPlot(TS)
#'
DummyTS=function(days=7,SampleTime_min=1,sd_days=0.25,Cunif=0.01,Csd=0.03,FouP=c(8.95652,-2.42785,-0.58782,-0.97698,-2.67130,0.34104,-0.14011,0.10532,0.57853,-0.15287,-0.15902),AnualPrecipitation_mm=600,Anual_Days_of_P=100,Beta=list(0.3,0.5,0.5,c(0.5,0.5),c(2,2),c(2,2)),Lin_Stor_Catchment_Area=2800,Lin_Stor_Number_of_Catchments_slow_fast=c(3,1),Lin_Stor_Average_Lambda_slow_fast = c(1000,100),Rain_off=FALSE,timezone="UTC"){
  x=seq(1/(SampleTime_min*60*24),1,1/(60*24)*SampleTime_min)
  per=2*pi
  ##dry pattern
  SingleDay=function(i){
    #random daily pattern
    addend=rnorm(n = 11,mean=0,sd=sd_days)
    xx=(FouP[1]+addend[1])+cos(per * x)*(FouP[2]+addend[2])+sin(per * x)*(FouP[3]+addend[3])+cos(per * x * 2)*(FouP[4]+addend[4])+sin(per * x * 2)*(FouP[5]+addend[5])+cos(per * x * 3)*(FouP[6]+addend[6])+sin(per * x * 3)*(FouP[7]+addend[7])+cos(per * x * 4)*(FouP[8]+addend[8])+sin(per * x * 4)*(FouP[9]+addend[9])+cos(per * x * 5)*(FouP[10]+addend[10])+sin(per * x * 5)*(FouP[11]+addend[11])
    #fix offset
    deviation=xx[1]-FouP[1]/2
    xx=xx-deviation
    y<<-c(y,xx)}

  y<-c()
  sapply(1:days,FUN=SingleDay)

  if(!Rain_off){
    ##Rain

    P=RainGen(Num_Days=days,AnualP_mm=AnualPrecipitation_mm,Anual_Days_of_Rain=Anual_Days_of_P,sampletime_min=SampleTime_min,beta=Beta)
    P=as.numeric(P)


    LinStor=function(P_i,Area,lambda){
      dt=1/(60*24)*SampleTime_min
      lambda=lambda*dt
      V_i=(Area*(exp(dt/lambda)-1)*P_i+V_i_1)*exp(-dt/lambda)
      V_i_1<<-V_i
      return(V_i)}

    MultipleLSs=function(i){

      As=abs(rnorm(1,mean=slowA,sd=slowA*0.5))
      Ls=abs(rnorm(1,mean=Ls ,sd=Ls*0.5))
      V_i_1<<-0
      RP=sapply(P,FUN=LinStor,Area=As,lambda=Ls)
      RAINPATTERN<<-RAINPATTERN+RP
    }
    MultipleLSf=function(i){

      As=abs(rnorm(1,mean=fastA,sd=fastA*0.5))
      Ls=abs(rnorm(1,mean=Lf ,sd=Lf*0.5))
      V_i_1<<-0
      RP=sapply(P,FUN=LinStor,Area=As,lambda=Ls)
      RAINPATTERN<<-RAINPATTERN+RP
    }

    Ls=Lin_Stor_Average_Lambda_slow_fast[1]
    Lf=Lin_Stor_Average_Lambda_slow_fast[2]

    slow=Lin_Stor_Number_of_Catchments_slow_fast[1]
    fast=Lin_Stor_Number_of_Catchments_slow_fast[2]
    p=Lin_Stor_Number_of_Catchments_slow_fast[1]/(Lin_Stor_Number_of_Catchments_slow_fast[1]+Lin_Stor_Number_of_Catchments_slow_fast[2])
    slowA=Lin_Stor_Catchment_Area*p
    fastA=Lin_Stor_Catchment_Area*(1-p)


    RAINPATTERN<<-rep(0,length(y))
    if(slow==0 & fast!=0){sapply(1:slow,FUN=MultipleLSs)}
    else if(slow!=0 & fast==0){sapply(1:fast,FUN=MultipleLSf)}
    else{
      sapply(1:slow,FUN=MultipleLSs)
      sapply(1:fast,FUN=MultipleLSf)}

    y=y+RAINPATTERN+runif(length(y),-Cunif,Cunif)+rnorm(length(y),0,Csd)+14}
  else{y=y+runif(length(y),-Cunif,Cunif)+rnorm(length(y),0,Csd)+14}
  start=strptime("01.01.2016 00:00",format="%d.%m.%Y %H:%M",tz=timezone)
  t=seq(start+SampleTime_min*60,start+(24*3600*days),SampleTime_min*60)
  TS=xts(y,order.by = t,tzone = timezone)
  return(TS)
}



#' A Generating Function
#'
#' The function disaggregates a diurnal pattern to a desired timestep. The disaggregation uses the beta parameter to determine the probabilities of rain for each disaggregation step. The first step breaks 24 hours into 3 times 8 hours. The following steps bisect each of the previous time intervals.
#' If the beta parameter contains 4 values the disaggregation will result in hourly values. 6 values will result in 15 min time steps.
#' Any desired sampletime below 15 min will be calculated by simply spreading it out out evenly.
#' The disaggragation matches a probability to each disaggregation step time interval. These will be multiplied to finally determine the rain probility (portion) of the final time steps.
#'
#' @param Num_Days The number of days to be generated.
#' @param AnualP_mm Anual precipitation in mm.
#' @param Anual_Days_of_Rain The number of days in rains in a year.
#' @param sampletime_min The logging interval of the resulting xts vector in minutes.
#' @param beta This list may have 4 to 6 positions. The values at each postion my be a number or a vector of 2.
#' The number represents the probability of each disaggreation slice to contain a rain event. By the same logic the vector contains
#' the 2 shape parameters of a beta distribution which will be sampled to determine the probility.
#'
#' Example:
#' beta=list(0.33,0.5,c(2,2),c(1,5))
#' The disaggregation will be performed down to hourly time steps (24/3 -> 3*8/2 -> 6*4/2 -> 12* 2/2).
#' In the first disaggregation step the probability of each of the 3 8 hour steps to contain a rain event is 0.33.
#' In the second disarggregation step the probability of each of the 6 4 hour steps to contain a rain event is 0.5.
#' The following two disaggregation steps set probabilities by random sampling a beta distribution.
#'
#' @param Rain_per_day_xts If NULL (default) th function will use the Num_Days, AnualP_mm and Anual_Days_of_Rain parameters to set up a vector containing daily rain events. Alternatively an xts with actual daily precipitation values may be set as input, which will be used for disaggregation instead.
#' @return Returns a generated precipitation signal (xts vector).
#' @keywords generate, create, rain
#' @export
#' @examples
#' Rain_xts=RainGen()
#'
RainGen=function(Num_Days=7,AnualP_mm=600,Anual_Days_of_Rain=121,sampletime_min=1,beta=list(0.3,0.5,0.5,c(0.5,0.5),c(2,2),c(2,2)),Rain_per_day_xts=NULL){

  x=seq(1/(sampletime_min*60*24),1,1/(60*24)*sampletime_min)


  if(is.null(Rain_per_day_xts)){
    Rain_or_None=rbinom(Num_Days, 1, Anual_Days_of_Rain/365)
    WRoN=which(Rain_or_None!=0)
    Rain_per_Day=rbeta(n = length(WRoN),shape1 =0.4,shape2 = 2)+abs(rnorm(n = length(WRoN),mean=1.5,sd=0.7))*rbinom(n = length(WRoN), 1, 0.1)++abs(rnorm(n = length(WRoN),mean=3,sd=3))*rbinom(n = length(WRoN), 1, 0.01)
    Rain_per_Day=Rain_per_Day*AnualP_mm/sum(Rain_per_Day)*Num_Days/365
    Rain_or_None[WRoN]=Rain_or_None[WRoN]*Rain_per_Day
    start=strptime("01.01.2016 00:00",format="%d.%m.%Y %H:%M",tz="")}
  else{Num_Days=dim(Rain_per_day_xts)[1]
  Rain_or_None=Rain_per_day_xts
  start=time(Rain_per_day_xts)[1]}

  Cascade=function(beta){
    #sample Alpha
    #8,4,2,1,0.5,0.25 h Kaskadenunterteilung

    if(length(beta)<4){warning("Beta need to be at least 4 long. (first divide is to 8 hours then 4 then 2 then 1)!");stop()}
    else if(length(beta)>6){warning("Beta cant be longer than 4. (first divide is to 8 hours then 4 then 2 then 1 then 0.5 then 0.25)!");stop()}

    Branches=function(i,beta,prev=NULL){
      if(i==1){
        if(length(beta[[i]])==1){
          if(i==1){sampsize=3}
          else{sampsize=2}
          P1=rbinom(sampsize, 1, beta[[i]])
          if(!any(P1==1)){P1[sample(1:sampsize,1)]=1}
          if(sum(P1)==3){P1=c(1/3,1/3,1/3)}
          else if(sum(P1)==2){P1[which(P1==1)]=c(0.5,0.5)}}
        else{Alpha=rbeta(1,beta[[i]][1],beta[[i]][2])
        P1=c(Alpha,1-Alpha)}}
      else{
        Split=function(prev){
          if(prev==0){P1<<-c(P1,0,0)}
          else{
            if(length(beta[[i]])==1){
              if(i==1){sampsize=3}
              else{sampsize=2}
              P2=rbinom(sampsize, 1, beta[[i]])
              if(!any(P2==1)){P2[sample(1:sampsize,1)]=1}
              if(sum(P2)==3){P2=c(0.3,0.3,0.3)}
              else if(sum(P2)==2){P2[which(P2==1)]=c(0.5,0.5)}}
            else{Alpha=rbeta(1,beta[[i]][1],beta[[i]][2])
            P2=c(Alpha,1-Alpha)}
            P1<<-c(P1,P2*prev)}
          return(P1)}

        P1<<-c()
        sapply(prev,FUN=Split)}
      return(P1)}

    #sort
    EventSORT=function(vec){

      sorter=function(vec,num=3){
        index=seq(num,length(vec),by=num)
        vec2=vec[index]
        vec[index]=rep(100,length(index))
        vec=sort(vec)
        vec2=sort(vec2,decreasing = T)
        vec[which(vec==100)]=vec2
        return(rev(vec))}

      ind=which(vec!=0)
      ind2=which(diff(ind)>1)
      start=ind[1]
      if(length(ind2)==0){vec[ind[1]:ind[length(ind)]]<-sorter(vec[ind[1]:ind[length(ind)]])}
      else{
        for(i in 1:(length(ind2)+1)){
          if(i<(length(ind2)+1)){num=ind2[i]
          end=ind[num]}
          else{end=ind[length(ind)]}
          vec[start:end]<-sorter(vec[start:end])
          if(i<(length(ind2)+1)){start=ind[num+1]}}}
      return(vec)}

    AdjustTime=function(vec){
      len=length(beta)
      if(len==4){Minuten=60}
      else if(len==5){Minuten=30}
      else if(len==6){Minuten=15}
      ratio=Minuten/sampletime_min
      adjust=function(i){
        start=1+ratio*(i-1)
        stop=ratio*(i)
        vec2[start:stop]<<-rep(vec[i]/ratio,ratio)}
      vec2<-rep(0,length(x))
      sapply(1:length(vec),FUN=adjust)
      return(vec2)}

    prev=c()
    for(i in 1:length(beta)){
      vec=Branches(i,beta=beta,prev = prev)
      prev=vec}
    vec=EventSORT(vec)
    vec=AdjustTime(vec)
    return(vec)}

  Zeros=function(len){return(rep(0,len))}

  Final_func=function(count){
    if(Rain_or_None[count]==0){
      A=Zeros(length(x))}
    else{A=Cascade(beta = beta)*Rain_or_None[count]}
    add=length(x)*(count-1)
    ARRAY[(1+add):(length(x)+add)]<<-A}

  ARRAY<<-rep(0,length(x)*Num_Days)
  sapply(c(1:length(Rain_or_None)),FUN = Final_func)


  t=seq(start+sampletime_min*60,start+(24*3600*Num_Days),sampletime_min*60)
  TS=xts(ARRAY,order.by = t)

  return(TS)}


#' A Generating Function
#'
#' The function inserts noise into a given xts. (dim(xts) = c(n,1))
#' @param TS Xts vector in which the noise is to be inserted.
#' @param PercentageNoise The percent of the signal that should inherit noise.
#' @param Noisespread Standard Deviation of the noise.
#' @param Overlap Boolean that allows (TRUE) different noise sections to overlap or not(FALSE).
#' For PercentageNoise >50 Overlap should be set to TRUE, as the computation will otherwise take quite long.
#' @param NoiseSections How many noisy section should be created.
#' @return Returns xts vector with noise.
#' @keywords generate, create, noise, outlier
#' @seealso \link[ISI.Toolbox]{DynPlot}
#' @export
#' @examples
#' A = DummyTS(Rain_off = T)
#' B = Noise(A,15,1.5,NoiseSections=3)
#' DynPlot(cbind(B,A),Labels=c("TS_with_Noise","TS_without_Noise"))
#'
Noise=function(TS,PercentageNoise, Noisespread, Overlap=FALSE, NoiseSections=1){
  Times=time(TS)
  Values=as.numeric(TS)
  numnoise=(length(Values)*PercentageNoise/100)
  #function speeds up the partition of the numnoise vector
  Factor=function(num){
    count=0
    while(num>999){
      num=num/10
      count=count+1}
    return(c(num,count))}

  NewNum.Factor=Factor(numnoise)
  lenvec=t(restrictedparts(NewNum.Factor[1], NoiseSections))[as.integer(runif(1,1,R(NoiseSections,NewNum.Factor[1]))),]
  lenvec=lenvec*10^NewNum.Factor[2]
  indexes=sample(1:length(Values), NoiseSections, replace=FALSE)
  if(Overlap){while(any(indexes+lenvec-1>length(Values))){
    indexes=sample(1:length(Values), NoiseSections, replace=FALSE)}}
  else{
    while(any(diff(c(indexes,length(Values)))<=lenvec)){
      indexes=sample(1:length(Values), NoiseSections, replace=FALSE)}}
  for(i in 1:length(lenvec)){
    start=indexes[i]
    stop=indexes[i]+lenvec[i]-1
    Values[start:stop]=Values[start:stop]+rnorm(lenvec[i],0,Noisespread)}
  TS=xts(Values, order.by=Times)
  return(TS)}



#' A Generating Function
#'
#' The function inserts outlier into a given xts. (dim(xts) = c(n,1))
#' @param TS Xts vector in which the outliers are to be inserted.
#' @param PercentageError The percent of the signal that should be outliers.
#' @param Outlierspread Standard Deviation of the Outliers
#' @param longspreadfactor Defines the std within an offset (several Outliers in a row) as fraction of Outlierspread.
#' @param keyword String "short" or "long". If "short" simple Outliers are created. If "long" Outliers are set in a sequence so that a Offset is generated.
#' @param maxlenoutlier Restricts the max length of offsets.
#' @return Returns xts vector with outliers
#' @keywords generate, create, noise, outlier
#' @seealso \link[ISI.Toolbox]{DynPlot}
#' @export
#' @examples
#' A = DummyTS(Rain_off = T)
#' B = Outliers(A,5,4)
#' DynPlot(cbind(B,A),Labels=c("TS_with_Outliers","TS_without_Outliers"))
#'
Outliers =function(TS, PercentageError, Outlierspread, keyword="short", maxlenoutlier=10, longspreadfactor=0.04){
  len=length(TS)
  indexTRFA=logical(len)
  numerror=(len*PercentageError/100)
  indexes=floor(runif(n=floor(numerror), min=2, max=(len-1)))
  indexes=indexes[!duplicated(indexes)]
  if(keyword=="short"){
    indexTRFA[indexes]=!logical(length(indexes))
    noise=rnorm(length(indexes),0,Outlierspread)
    TS[indexTRFA]=TS[indexTRFA]+noise}
  else if(keyword=="long"){
    count=0
    i=1
    while(count<length(indexes)){
      lenout=floor(runif(1, min=1, max=maxlenoutlier))
      alloffset=rnorm(1,0,Outlierspread)
      if(indexes[i]+lenout>len){lenout=len-indexes[i]}
      else{}
      singleoffset=rnorm(lenout,0,(Outlierspread*longspreadfactor))
      indexTRFA=logical(len)
      indexTRFA[seq(indexes[i],(indexes[i]+lenout-1),1)]=!logical(lenout)
      noise=(rep(alloffset,lenout)+singleoffset)
      TS[indexTRFA]=TS[indexTRFA]+noise
      count=count+lenout
      i=i+1}}
  else{print("Wrong input String!")}
  return(TS)
}





