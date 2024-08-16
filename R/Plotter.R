#' A Simple Dynamic Plotter based on Dygraphs
#'
#' Plots data with less code than in dygraph with limited options.
#' @param X The xts to plot. For multiple use X=cbind(first,second,third)
#' @param Title String. Title of graph.
#' @param Ylab String. Title of left y-axis.
#' @param Y2lab String. Title of right y-axis.
#' @param Labels Vector of Strings. Length must be equal to dim(X)[2].
#' @param Dotsize Decimal.
#' @param YRange Vector of length 2. Setting the limits for the left y-axis.
#' @param Y2Range Vector of length 2. Setting the limits for the right y-axis.
#' @param colset A color palette. Possible inputs are: 1, 2 or 3.
#' @param Axis A vector. Length must be equal to dim(X)[2]. The vector specifies which y axis the xts should be plotted on. 1 equals the left axis, 2 equals the right axis.
#' @param Events A vector of datetime strings (format="\%Y-\%m-\%d \%H-\%M-\%S"). For each entry a vertical line will be added to the plot.
#' @param Eventnames A vector of strings. Must have same length as Events vector. Lables each event.
#' @param legendwidth Integer. Width of the legend shown in the dynamic plot. Default is 200.
#' @seealso \link[dygraphs]{dygraph}
#' @return Dynamic Plot.
#' @keywords dygraph, dynamic plot, plot
#' @export
#' @examples
#' a=DummyTS(days=14)
#' b=DummyTS(days=14)
#' c=DummyTS(days=14)
#' DynPlot(cbind(a,b,c),Axis=c(1,2,2),Labels=c("First","Second","Third"),Events = c("2016-01-02 01:00:00"),Eventnames = c("Eventname"))
DynPlot = function(X, Title = NULL, Ylab = NULL, Y2lab = NULL, 
                   Labels = seq(1, dim(X)[2], 1), Dotsize = 1, 
                   YRange = NULL, Y2Range = NULL, colset = 2, 
                   Axis = rep(1, dim(X)[2]), Events = NULL, 
                   Eventnames = NULL, legendwidth = 200) {
  
  # Assign labels to the columns if not provided
  if (is.null(names(X))) {
    names(X) <- Labels
  }
  
  # Define the axis types
  axis <- c("y", "y2")
  
  # Predefined colors
  COLS <- c("rgb(240,163,255)", "rgb(0,117,220)", "rgb(153,63,0)", 
            "rgb(76,0,92)", "rgb(25,25,25)", "rgb(0,92,49)", "rgb(43,206,72)", 
            "rgb(255,204,153)", "rgb(128,128,128)", "rgb(148,255,181)", 
            "rgb(143,124,0)", "rgb(157,204,0)", "rgb(194,0,136)", 
            "rgb(0,51,128)", "rgb(255,164,5)", "rgb(255,168,187)", 
            "rgb(66,102,0)", "rgb(255,0,16)", "rgb(94,241,242)", 
            "rgb(0,153,143)", "rgb(224,255,102)", "rgb(116,10,255)", 
            "rgb(153,0,0)", "rgb(255,255,128)", "rgb(255,255,0)", 
            "rgb(255,80,5)")
  
  # Start building the dygraph
  d <- dygraph(X, main = Title, ylab = Ylab, group = "1") %>%
    dyOptions(useDataTimezone = TRUE, drawPoints = TRUE, 
              connectSeparatedPoints = TRUE, pointSize = Dotsize, 
              colors = randomColor(dim(X)[2]))
  
  # Loop through columns and add dySeries to the dygraph
  for (i in 1:dim(X)[2]) {
    d <- d %>% dySeries(names(X)[i], axis = axis[Axis[i]], color = COLS[i])
  }
  
  # Add events if provided
  if (!is.null(Events)) {
    for (k in 1:length(Events)) {
      d <- d %>% dyEvent(Events[k], Eventnames[k], labelLoc = "top", color = "black", strokePattern = "dashed")
    }
  }
  
  # Add axes settings and other options
  d <- d %>%
    dyAxis("y", valueRange = YRange) %>%
    dyAxis("y2", valueRange = Y2Range, label = Y2lab) %>%
    dyRangeSelector() %>%
    dyLegend(width = legendwidth, show = "auto", showZeroValues = TRUE, labelsDiv = NULL, hideOnMouseOut = TRUE)
  
  return(d)
}

#basic 2D plot of the time frequency intensity analysis
#Input: Matrix generated with KDH function, directory and filename
#Output: Plot, pdf
#Example: PlotKDH2D(Matrix,"C:/Work/R/KDH plots/QMS6.pdf")
#' Plot Function
#'
#' The function plots the results of the KDH function in a levelplot.
#' @param M Matrix from KDH function.
#' @param savetofolder String. To specify folder and name where plot should be saved.
#' @param ColBarRange Vector of length 2. Range for the levelplot colorbar.
#' @param LogScale Boolean. If TRUE the colorbar (and occurances) will be log scaled.
#' @return Returns Levelplot.
#' @seealso \link[ISI.Toolbox]{KDH}
#' @export
#' @examples
#' TS=DummyTS()
#' M=KDH(TS,limits=c(10,20,30,40,50,60),timeclasses=c(0.5/48,1/48.,1/24.,1/12.,1/8.,1/6.,1/3.,1))
#' PlotKDH2D(M)
PlotKDH2D = function(M,savetofolder=NULL,ColBarRange=NULL,LogScale=FALSE){
  #format Matrix
  Xn=M[2:(length(M[,1])),1]
  Yn=M[1,][2:length(M[1,])]
  Z=M[2:(length(M[,2])),2:(length(M[2,]))]
  if(LogScale){Z=log10(Z)}
  else{}
  if(is.null(ColBarRange)){SEQ=seq(0,max(Z))}
  else{SEQ=seq(ColBarRange[1],ColBarRange[2])}
  YnL= paste(fractions(Yn/3600))
  XnL=paste(round(Xn,1))
  if(!is.null(savetofolder)){
  pdf(savetofolder)
  ABC=levelplot(Z, at = SEQ,
                col.regions = colorRampPalette(c("yellow", "red")),
                ylab="t [h]",ylim=YnL,xlim=XnL ,xlab="Q [l/s]")
  print(ABC)
  dev.off()}
  else{
    ABC=levelplot(Z, at = seq(0,max(Z)),
                  col.regions = colorRampPalette(c("yellow", "red")),
                  ylab="t [h]",ylim=YnL,xlim=XnL ,xlab="Q [l/s]")
  }
  return(ABC)
}



