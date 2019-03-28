#' An Import Function
#'
#' This function allows you to import data using the read.table function and converts to the xts format.
#' @param filename File path and name.
#' @param sep Seperator of rows.
#' @param format Specifictation of the date time format for the conversion to xts.
#' @param header Does the document you intent to read from have a (1-line) header?
#' @param dec Decimal seperator of numbers.
#' @param skip How many lines should be skipped when reading the document?
#' @param linestoread Do you love cats?
#' @param SpacesbeforeDateTime How many spaces are in front the date-time expression?
#' @param TimeZone The timezone variable of the xts conversion.
#' @param fill Boolean. Should gaps in be filled with NA? (For reading documents with a differing number of rows.)
#' @param comment.char Character. Specifying what strings are to be ignored. (All characters in one line after the specified comment.char will be ignored). Default is turned off comment.char="").
#' @param na.strings Character vector. Specifyfing all characters strings that will be turned to NA's.
#'
#' @seealso \link[utils]{read.table},\link[base]{as.POSIXct},\link[xts]{xts},\link[base]{strptime}
#' @export
#' @examples
#' TimeSeries=ReadTabletoXTS("C:/file.txt",header=TRUE)
#'
ReadTabletoXTS = function(filename, sep=",", format="%Y-%m-%d %X", header=FALSE, dec = ".", skip=0,linestoread=-1,SpacesbeforeDateTime=0,TimeZone="UTC",fill=FALSE,comment.char="",na.strings=c("#-1")){
  data.4 <-read.table(filename, header=header, sep=sep,dec=dec,skip=skip, nrow = linestoread,fill=fill,comment.char = comment.char,na.strings=na.strings)
  x=strsplit(format,split="%")[[1]][4]
  DateTimeSep=substring(x,2,nchar(x))
  if(sep!=DateTimeSep){
    time.4 <- as.POSIXct(strptime(data.4[,1+SpacesbeforeDateTime],format=format, tz = TimeZone))
    test.1 <- subset(data.4,select=-(1+SpacesbeforeDateTime))}
  else{
    time.4 <- as.POSIXct(strptime(paste(data.4[,1+SpacesbeforeDateTime],data.4[,2+SpacesbeforeDateTime],sep=DateTimeSep),format=format, tz = TimeZone))
    test.1 <- subset(data.4,select=c(-(1+SpacesbeforeDateTime),-(2+SpacesbeforeDateTime)))}
  for(i in 1:dim(test.1)[2]){if(is.factor(test.1[,i])){test.1[,i]=as.numeric(as.character(test.1[,i]))}
    else{test.1[,i]=as.numeric(test.1[,i])}}
  PCM <- xts(test.1,order.by=time.4,tzone = TimeZone)
  return(PCM)
}


#' An Import Function
#'
#' This function allows to export an xts object to a file. It's based on the write.zoo function. Please use that one if you want more input options.
#' @param TS Xts to be written to a txt.
#' @param DirName Destination and name of the file to be safed.
#' @seealso \link[zoo]{write.zoo}
#' @export
#' @examples
#' WriteXTStoTXT(TS,"C:/Work/R/output.txt")
#'
WriteXTStoTXT = function(TS,DirName="C:/Work/R/output.txt"){
  write.zoo(TS,file=DirName,col.names = T,sep=",", quote=FALSE, na="")}


