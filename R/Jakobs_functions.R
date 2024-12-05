#' Base Flow Index Function
#'
#' Calculates Base Flow Index of a flow hydrograph (Gustard et al. 1992)
#' @param Flow_XTS Xts vector of flow hydrograph with minimum duration of five days. Ideally no NA values are present and the time step is even. (NA could be replaced with zeros if present. To unify the the time step the Add.Na function could be used.)
#' @param n number of days (numeric) for the calculation of local minimum (default is five days, can be chosen smaller for flushy hydrographs and high resoluted data)
#' @keywords baseflow
#' @return A list.\cr
#' [[1]] -> Base Flow Index (BFI) as numeric\cr
#' [[2]] -> daily mean, five day minima and baseflow caluclated with by difference with the mean daily-Q or min-daily Q as xts\cr
#' @examples
#' TS=DummyTS(days=30)
#' BaseFlowTS=BaseFlowGustard(TS)
#' @export
BaseFlowGustard=function (Flow_XTS,n=5) 
{
  daily = aggregate(Flow_XTS, list(Index = cut(index(Flow_XTS), 
                                               breaks = "day")), mean, na.rm = TRUE)
  daily2 = aggregate(Flow_XTS, list(Index = cut(index(Flow_XTS), 
                                                breaks = "day")), min, na.rm = TRUE)#
  n
  fiveday = aggregate(daily, list(rep(1:(nrow(daily)%/%n), 
                                      each = n, len = nrow(daily))), min, na.rm = T)
  f1 = as.numeric(fiveday)
  f2 = 0.9 * f1
  f3 = NULL
  f4 = NULL
  i = NULL
  n1 = length(f1)
  for (i in seq(2, n1 - 1, 1)) {
    f3[1] = f1[1]
    f3[n1] = f1[n1]
    if (f2[i] < f1[i - 1] & f2[i] < f1[i + 1]) {
      f3[i] = f1[i]
    }
    else {
      f3[i] = NA
    }
  }
  f3 = na.approx(f3, na.rm = FALSE)
  
  c1 = as.data.frame(daily)
  d1 = as.data.frame(daily2)
  
  c2 = as.data.frame(fiveday)
  c3 = as.data.frame(ifelse(c1[, 1] %in% c2[, 1], c1[, 1], 
                            NA))
  c3 = as.data.frame(cbind(c1, c3)) 
  n = nrow(c3)
  c3[1, 2] = c3[1, 1]
  c3[n, 2] = c3[n, 1]
  c3[, 2] = FillNaValues(c3[, 2], maxgap = 10)
  c3[, 3] = ifelse(c3[, 1] %in% f3, c3[, 1], NA)
  c3[1, 3] = c3[1, 1]
  c3[n, 3] = c3[n, 1]
  c3[, 3] = FillNaValues(c3[, 3], maxgap = 100)
  c3[, 3]<-ifelse(c3[,3]>c3[,1],c3[,1],c3[,3])
  c3[,4]<-c3[,3]
  c3[,5]<-d1
  c3[, 4]<-ifelse(c3[,3]>c3[,5],c3[,5],c3[,3])
  c3=c3[,c(1,2,3,4)]
  colnames(c3) = c("daily", "fivedayminima", "baseflow","baseflow_day_min")
  c3=xts(c3,order.by = as.POSIXct(rownames(c3),tz="UTC"))
  c3[, 3] = FillNaValues(c3[, 3], maxgap = 100)
  c4 = sum(c3[, 3])
  c5 = sum(c3[, 1])
  BFI = c4/c5
  result = list(BFI, c3)
  names(result) = c("BFI", "baseflow")
  return(result)
}

#' Base Flow Index Function
#'
#' Calculates Base Flow Index of a flow hydrograph (according to Duncan 2019, https://doi.org/10.1016/j.jhydrol.2019.05.040). Calibration guidelines 
#' 1)The master recession should show a step up during significant rain (except shortly after a previous event, where event runoff is still present), since significant rain should recharge groundwater and
#' thus increase baseflow. If it does not, a steeper recession (lower k) may be more appropriate.
#' 2)The master recession should not lie much below total flow in the absence of rain (except shortly after an event, where event runoff is still present). If it does, the modelled baseflow is receding more 
#' slowly than the observed flows. A steeper recession (lower k) may be more appropriate.
#' 3)The master recession should not cling tightly to total flow in the absence of rain and be continuously held down by the 'not greater than total flow' condition. If it does, the modelled recession is
#' steeper than that of the observed flows. A flatter recession (higher k) may be more appropriate.
#' 4)If zero flows may occur, the constant c must be a small negative flow, typically a few percent of mean flow at the site. If zero flows never occur, c may be initially set to zero. The constant c should be 
#' calibrated using observed flows close to zero.
#' 5)Setting the filter parameter in the second pass equal to the recession constant in the first pass appears to be satisfactory in each case examined.
#' @param Flow_XTS Xts vector of flow hydrograph. Ideally no NA values are present and the time step is even. (NA could be replaced with zeros if present. To unify the the time step the Add.Na function could be used.)
#' @param k recession constant (numeric) ------k=0.925 for daily data from all, btw. 0.91 to 0.98 catchments has been proposed; If the daily time step is now subdivided into 24 hourly sub-steps, the recession parameter in each sub-step must be the 24th root of the original k, to reach the same value at the end of the day. Hence the daily recession parameter of 0.93 at Little Stringybark Creek becomes 0.9970 at hourly sub-steps, and 0.99970 at six-minute sub-steps.
#' @param c is a constant flow added to the exponential decay component. see point 4) in the description above
#' @keywords baseflow
#' @return A xts vector.\cr
#' [1] -> Base Flow after a master recession curve (for checking)\cr
#' [2] -> Base Flow after a master recession curve and a digital filter (final baseflow)\cr
#' @examples
#' TS=DummyTS(days=30)
#' BaseFlowTS=BaseFlowDuncan(TS)
#' @export
BaseFlowDuncan <- function(Flow_XTS,k=.93,c=0) {
  ####establish master recession curve
  M <- rep(0, length(Flow_XTS))
  QB<-rep(0,length(Flow_XTS))
  QD<-rep(0,length(Flow_XTS))
  M[length(Flow_XTS)]<-Flow_XTS[length(Flow_XTS)]
  for (i in length(Flow_XTS):2) {
    M[i-1]<-(M[i]-c)/k+c
    if(M[i-1]>=Flow_XTS[i-1])    {M[i-1]<-Flow_XTS[i-1]  }
  }  
  ####establish digital model
  QD[length(M)]<-M[length(M)]
  for (i in 2:length(Flow_XTS)) {
    QD[i]<-k*QD[i-1]+(M[i]-M[i-1])*(1+k)/2
    if(QD[i]>M[i]){
      QD[i]<-M[i]
    }
    if(QD[i]<0){
      QD[i]<-0
    }
    QB[i]<-M[i]-QD[i]   }
  together=cbind(M,QB)
  colnames(together) = c("master recession curve","master recession curve & digital filter")
  return(as.xts(together,order.by=index(Flow_XTS)))}

#' Gradient filter for Events
#'
#' Calculates Base Flow Index of a flow hydrograph (according to Duncan 2019, https://doi.org/10.1016/j.jhydrol.2019.05.040). Calibration guidelines 
#' 1)The master recession should show a step up during significant rain (except shortly after a previous event, where event runoff is still present), since significant rain should recharge groundwater and
#' thus increase baseflow. If it does not, a steeper recession (lower k) may be more appropriate.
#' 2)The master recession should not lie much below total flow in the absence of rain (except shortly after an event, where event runoff is still present). If it does, the modelled baseflow is receding more 
#' slowly than the observed flows. A steeper recession (lower k) may be more appropriate.
#' 3)The master recession should not cling tightly to total flow in the absence of rain and be continuously held down by the 'not greater than total flow' condition. If it does, the modelled recession is
#' steeper than that of the observed flows. A flatter recession (higher k) may be more appropriate.
#' 4)If zero flows may occur, the constant c must be a small negative flow, typically a few percent of mean flow at the site. If zero flows never occur, c may be initially set to zero. The constant c should be 
#' calibrated using observed flows close to zero.
#' 5)Setting the filter parameter in the second pass equal to the recession constant in the first pass appears to be satisfactory in each case examined.
#' @param Flow_XTS Xts vector of flow hydrograph. Ideally no NA values are present and the time step is even. (NA could be replaced with zeros if present. To unify the the time step the Add.Na function could be used.)
#' @param k recession constant (numeric) ------k=0.925 for daily data from all, btw. 0.91 to 0.98 catchments has been proposed; If the daily time step is now subdivided into 24 hourly sub-steps, the recession parameter in each sub-step must be the 24th root of the original k, to reach the same value at the end of the day. Hence the daily recession parameter of 0.93 at Little Stringybark Creek becomes 0.9970 at hourly sub-steps, and 0.99970 at six-minute sub-steps.
#' @param c is a constant flow added to the exponential decay component. see point 4) in the description above
#' @keywords baseflow
#' @return A xts vector.\cr
#' [1] -> Base Flow after a master recession curve (for checking)\cr
#' [2] -> Base Flow after a master recession curve and a digital filter (final baseflow)\cr
#' @examples
#' TS=DummyTS(days=30)
#' BaseFlowTS=BaseFlowDuncan(TS)
#' @export

Grad_Event=function(ts,ts_waterlevel,show_hy_summer_winter=FALSE,Analyse=FALSE,check_result= F,show_event_selection=F,
grad_event=9999,grade_baseflow=9999,
fixed_limit_summer = .15,fixed_limt_winter = .17,time_buffer=20*60,
min_duration_in_time_steps = 5, padding_in_time_steps = 60){
  
  ### beginn with checking fixed limits (just guess and see with show_hy_summer_winter=T)
  ### ts=oxy/ph/wl etc
  ### accept homogenising of timesteps
  ###make sure you update the EventDetector function
  ###takes 10 min before and after into the filter
  
  sep_ts=Separate.To.HyS.HyW(ts_waterlevel)
  
  hydro_summer <- sep_ts[[1]]
  hydro_winter <- sep_ts[[2]]
  
  if(show_hy_summer_winter){
    print(DynPlot(cbind(hydro_summer,hydro_winter),Labels = c("hy_summer","hy_winter")))}
  else{
    
    checktimes = c(NA, difftime(index(ts_waterlevel[-1]), index(ts_waterlevel[-nrow(ts_waterlevel)]), 
                                units = "secs"))
    checktimestepsdiff = hist(checktimes, plot = F)
    if (length(checktimestepsdiff$counts) > 2) {
      
      print("summer events")
      events_summer=EventDetector(hydro_summer, quantile_limit = NULL, fixed_limit = fixed_limit_summer,
                                  min_duration_in_time_steps = min_duration_in_time_steps, 
                                  padding_in_time_steps = padding_in_time_steps,
                                  single_events_in_list = T)
      print("winter events")
      events_winter=EventDetector(hydro_winter, quantile_limit = NULL, fixed_limit = fixed_limt_winter,
                                  min_duration_in_time_steps = min_duration_in_time_steps, 
                                  padding_in_time_steps = padding_in_time_steps,
                                  single_events_in_list = T)
    }else{
      events_summer=EventDetector(hydro_summer, quantile_limit = NULL, fixed_limit = fixed_limit_summer,
                                  min_duration_in_time_steps = min_duration_in_time_steps, 
                                  padding_in_time_steps = padding_in_time_steps,
                                  single_events_in_list = T)
      events_winter=EventDetector(hydro_winter, quantile_limit = NULL, fixed_limit = fixed_limt_winter,
                                  min_duration_in_time_steps = min_duration_in_time_steps, 
                                  padding_in_time_steps = padding_in_time_steps,
                                  single_events_in_list = T)
    }
    
    events=rbind(events_summer[[2]],events_winter[[2]])
    
    if(show_event_selection)
    {print(DynPlot(cbind(ts_waterlevel,rbind(events_summer[[2]]),rbind(events_winter[[2]])),Labels = c("entire series","summer ts","winter ts")))
      stop("are you happy with the selected events?")}
    
    print(paste0("Number of detected events: ",length(events_summer[[1]])+length(events_winter[[1]])))
    
    starts=c(lapply(events_summer[[1]], function(x) as.character(format(index(x)[1], "%Y-%m-%d %X"))),lapply(events_winter[[1]], function(x) as.character(format(index(x)[1], "%Y-%m-%d %X"))))
    ends=c(lapply(events_summer[[1]], function(x) as.character(format(index(x)[length(x)], "%Y-%m-%d %X"))),lapply(events_winter[[1]], function(x) as.character(format(index(x)[length(x)], "%Y-%m-%d %X"))))
    char_ev_start=unlist(c(starts))
    char_ev_end=unlist(c(ends))
    char_ev=rbind(char_ev_start,char_ev_end)
    char_ev=unique(char_ev)#double events were found...
    
    #print(DynPlot(cbind(ts_waterlevel,events_summer[[2]],events_winter[[2]]),Labels = c("original","summer","winter"),Events = char_ev))
    ts_coutout=c()
    ts_temp_winter=list()
    ts_output=c()
    
    time_before=-time_buffer
    time_after=time_buffer
    ##summer
    for(i in 1:length(events_winter[[1]])){
      ts_cutout_temp=CutTimeSeries(ts,start(events_winter[[1]][[i]])+time_before,end(events_winter[[1]][[i]])+time_after)
      ts_coutout=rbind(ts_coutout,ts_cutout_temp)
      
      ts_temp_winter[[i]]=CutTimeSeriesDel(ts,start(events_winter[[1]][[i]])+time_before,end(events_winter[[1]][[i]])+time_after)
    }
    ts_base_winter=na.omit(do.call("cbind",ts_temp_winter))[,1]
    
    #winter
    ts_temp_summer=list()
    
    for(i in 1:length(events_summer[[1]])){
      ts_cutout_temp=CutTimeSeries(ts,start(events_summer[[1]][[i]])+time_before,end(events_summer[[1]][[i]])+time_after)
      ts_coutout=rbind(ts_coutout,ts_cutout_temp)
      ts_temp_summer[[i]]=CutTimeSeriesDel(ts,start(events_summer[[1]][[i]])+time_before,end(events_summer[[1]][[i]])+time_after)
    }
    
    ts_base=na.omit(cbind(ts_base_winter,do.call("cbind",ts_temp_summer)))[,1]
    
    #baseflow=ts_waterlevel[!(index(ts_waterlevel)%in%index(events))]
    #ts_event=na.omit(cbind(ts,events))[,1]
    #ts_baseflow=na.omit(cbind(ts,baseflow))[,1]
    
    if(Analyse){
      print("event")
      GradFilter(ts_coutout,Analyse = T)
      print("baseflow")
      GradFilter(ts_base,Analyse = T)
    } 
    else{ 
      event=GradFilter(ts_coutout,abs_grad = grad_event)
      baseflow=GradFilter(ts_base,abs_grad = grade_baseflow)
      baseflow_smooth=Smooth(baseflow,window_width = 5)
      ts_output=rbind(event,baseflow_smooth)}
    
    if(check_result){
      print(DynPlot(cbind(ts,ts_output),Events = char_ev))
    }else{}
  }
  return(ts_output)
}






