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
