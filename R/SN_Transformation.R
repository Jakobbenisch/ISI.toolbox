SN.GetSNInfo=function(CenterLineCoo,DepthCoo){
  ReturnList=list()
  Index=SN.GetIndex(CenterLineCoo,DepthCoo)

  if(Index==1){Index=2}
  else if(is.na(CenterLineCoo[,1][(Index+1)])){Index=Index-1}
  A=CenterLineCoo[(Index-1),][,1:2]
  B=CenterLineCoo[Index,][,1:2]
  C=CenterLineCoo[(Index+1),][,1:2]
  P=DepthCoo

  #All 5 distances for the 2 triangles
  L_AP=SN.EukDist(A,P)
  L_AB=SN.EukDist(A,B)
  L_BP=SN.EukDist(B,P)
  L_BC=SN.EukDist(B,C)
  L_PC=SN.EukDist(P,C)

  #The 4 relevant angles
  #left triangle
  Angle_PAB=acos((L_AP^2+L_AB^2-L_BP^2)/(2*L_AP*L_AB))
  Angle_ABP=acos((L_BP^2+L_AB^2-L_AP^2)/(2*L_BP*L_AB))
  #right triangle
  Angle_PBC=acos((L_BP^2+L_BC^2-L_PC^2)/(2*L_BC*L_BP))
  Angle_BCP=acos((L_PC^2+L_BC^2-L_BP^2)/(2*L_PC*L_BC))

  #print(c(Angle_PAB,Angle_ABP,Angle_PBC,Angle_BCP)*180/pi)

  #case e)
  if((Angle_PAB == Angle_BCP) | (Angle_PAB <= pi/2 & Angle_BCP <= pi/2 & Angle_PBC > pi/2 & Angle_ABP > pi/2) | (Angle_PAB > pi/2 & Angle_BCP > pi/2 & Angle_PBC <= pi/2 & Angle_ABP <= pi/2)){
    M=B
    P_AM = SN.EukDist(P,M)*SN.PointAboveOrBelow(P,A,M)
    ReturnList[[1]]=SN.EukDist(A,M)
    ReturnList[[2]]=P_AM
    ReturnList[[3]]="M"}
  #case d)
  else if(Angle_PAB <= pi/2 & Angle_ABP <= pi/2 & Angle_PBC <= pi/2 & Angle_BCP <= pi/2){
    P_AB = L_AP*sin(Angle_PAB)*SN.PointAboveOrBelow(P,A,B)
    P_BC = L_BP*sin(Angle_PBC)*SN.PointAboveOrBelow(P,B,C)

    if(abs(P_AB)<abs(P_BC)){
      S_PAB = L_AP*cos(Angle_PAB)
      ReturnList[[1]]=S_PAB
      ReturnList[[2]]=P_AB
      ReturnList[[3]]="AB"}
    else{
      S_PBC = L_BP*cos(Angle_PBC)
      ReturnList[[1]]=S_PBC
      ReturnList[[2]]=P_BC
      ReturnList[[3]]="BC"}}
  #cases b)
  else if(Angle_PAB <= pi/2 & Angle_ABP <= pi/2){
    P_AB = L_AP*sin(Angle_PAB)*SN.PointAboveOrBelow(P,A,B)
    S_PAB = L_AP*cos(Angle_PAB)
    ReturnList[[1]]=S_PAB
    ReturnList[[2]]=P_AB
    ReturnList[[3]]="AB"}
  else if(Angle_PBC <= pi/2 & Angle_BCP <= pi/2){
    P_BC = L_BP*sin(Angle_PBC)*SN.PointAboveOrBelow(P,B,C)
    S_PBC = L_BP*cos(Angle_PBC)
    ReturnList[[1]]=S_PBC
    ReturnList[[2]]=P_BC
    ReturnList[[3]]="BC"}
  else{
    warning("Something unanticipated happened!")

    stop()}

  ReturnList[[4]]=Index-1
  names(ReturnList)=c("S_part","N","Closest_Part","Index")
  return(ReturnList)}

SN.EukDist=function(XY1,XY2){
  return(sqrt((XY1[1]-XY2[1])^2+(XY1[2]-XY2[2])^2))}

SN.EasyArcLength=function(Matrix){
  return(sum(sqrt(diff(Matrix[,1])^2+diff(Matrix[,2])^2)))}

SN.PointAboveOrBelow=function(P,A,B){
  v=sign((B[2] - A[2]) * (P[1] - A[1]) - (B[1] - A[1]) * (P[2] - A[2]))
  return(v)
}

SN.GetIndex=function(CenterLineCoo,DepthCoo){
  Index=which.min(abs((as.numeric(CenterLineCoo[,1])-as.numeric(DepthCoo[1]))^2+(as.numeric(CenterLineCoo[,2])-as.numeric(DepthCoo[2]))^2))
  return(Index)
}

SN.Transformation2SN=function(LoopIndex,CenterlineCoordinates,DepthCoordinates){

  ORIGIN=1
  SNI_List=SN.GetSNInfo(CenterlineCoordinates,DepthCoordinates[LoopIndex,])
  if(SNI_List[[3]]=="M" | SNI_List[[3]]=="AB"){END=SNI_List[[4]]}
  else{END=SNI_List[[4]]+1}

  S=SN.EasyArcLength(as.matrix(CenterlineCoordinates[ORIGIN:END,]))

  newDepthCoordinates=c(S+SNI_List[[1]],SNI_List[[2]],DepthCoordinates[LoopIndex,][3])
  return(newDepthCoordinates)
}


SN.TransPoints=function(X,Alpha,N){
  #X-> H,R (y,x)
  dx=sin(Alpha)*abs(N)
  dy=cos(Alpha)*abs(N)
  if(N>0){
    if(X[2]>0 &X[1]>0){
      P=X+c(-dy,dx)
    }else if(X[2]>0 &X[1]<0){
      P=X+c(-dy,-dx)
    }else if(X[2]<0 &X[1]>0){
      P=X+c(dy,dx)
    }else if(X[2]<0 &X[1]>0){
      P=X+c(dy,-dx)
    }else{warning("Undefined Point");stop()}

  }else if(N<0){
    if(X[2]>0 &X[1]>0){
      P=X+c(dy,-dx)
    }else if(X[2]>0 &X[1]<0){
      P=X+c(dy,dx)
    }else if(X[2]<0 &X[1]>0){
      P=X+c(-dy,-dx)
    }else if(X[2]<0 &X[1]>0){
      P=X+c(-dy,dx)
    }else{warning("Undefined Point");stop()}

  }else{P=X}

  return(P)
}

SN.Transformation2XY=function(LoopIndex,CenterlineCoordinates,S,N,Z,ThalwegCumSum){

  #sorts DF after S coordinates
  LenDiff=ThalwegCumSum-S[LoopIndex]
  MinDiff=min(LenDiff[LenDiff>=0])
  Ind=which(LenDiff==MinDiff)
  if(MinDiff==0){Ind=Ind+1}

  P1=as.numeric(CenterlineCoordinates[Ind,])
  P2=as.numeric(CenterlineCoordinates[(Ind+1),])

  if(MinDiff==0){
    X=P1
    P1=as.numeric(CenterlineCoordinates[(Ind-1),])
    Alph=abs(atan((P2[1]-P1[1])/(P2[2]-P1[2])))

  }else{
    X=c(P2[1]-P1[1],P2[2]-P1[2])
    X=X/sqrt(X[1]^2+X[2]^2)*MinDiff+c(P1[1],P1[2])
    Alph=abs(atan(X[2]/X[1]))
  }

  XY<<-list()
  for(i in 1:length(N)){
    XY[[i]]<<-SN.TransPoints(X,Alph,N[i])}

  # XY=foreach(i = 1:length(N),.combine = rbind)%do%
  #   SN.TransPoints(X,Alph,N[i])

  XY<<-matrix(unlist(XY), ncol = 2,byrow = TRUE)

  newDepthCoordinates=list(XY[,1],XY[,2],Z[LoopIndex,])
  return(newDepthCoordinates)
}


#' A function for transforming XY to SN coordinates.
#'
#' The transformation needs 2 dataframes as input. One resembling the centerline, the other containing measured depth information and the riversides.\cr
#' The measured points must be within the riversides!
#' @param ThalwegCoordinates_HR Dataframe with dim = n,2. First postion is the Hochwert (y-coordinate) and second ist the Rechtswert (x-coordinate). This dataset can be estimated by the centerline oy your river.
#' @param DepthCoordiantes_HRZ Dataframe with dim = n,2.  First postion is the Hochwert (y-coordinate), second ist the Rechtswert (x-coordinate) and third is the Depth (z-coordinate). This dataset contains the measured depth as well as the riverside coordinates with a corresponding depth of 0.
#' @return A matrix containing the S-coordinates [,1], N-coordinates [,2] and the depth [,3].
#' @seealso \link[ISI.Toolbox]{SN2XY}
#' @export
#' @examples
#' ### This example incoorporates an entire SN transformation ###
#' ##In order to use the example you need to install and load the akima package as we need it for our grid interpolation##
#' # Create Testdata for this example
#' R=seq(362000,364000)
#' X=R-363000
#' H=X^2*10/80000+20*sin(X/50/pi)+5666100
#' H2=X^2*10/80000+20*sin(X/50/pi)-20+5666100
#' CL=0.5*(H2+H)
#'
#' DPR=runif(n = 1500,min = 362010,max=363990)
#' DPH=c()
#' DPZ=c()
#' for(i in 1:length(DPR)){
#'   index=which.min(abs(R-DPR[i]))
#'   value=rnorm(1,mean=0,sd=5)+CL[index]
#'   d=abs(value-CL[index])
#'   depth=5+d*-5/10
#'   if(value>H[index]){value=H[index]-0.5}
#'   else if(value<H2[index]){value=H2[index]+0.5}
#'   d=abs(value-CL[index])
#'   depth=5+d*-5/10
#'   DPH=c(DPH,value)
#'   DPZ=c(DPZ,depth)
#' }
#'
#' UB=data.frame(H=H,R=R,Z=rep(0,length(R)))
#' LB=data.frame(H=H2,R=R,Z=rep(0,length(R)))
#' Data=data.frame(H=DPH,R=DPR,Z=DPZ)
#' DP=rbind(UB[10:(dim(UB)[1]-10),],LB[10:(dim(LB)[1]-10),],Data)
#' CLP=data.frame(H=CL,R=R)
#'
#' # plotting data before transformation and interpolation
#' XY=DP
#' R01=XY[,3]/max(XY[,3])
#' plot(XY[,2],XY[,1],col=rgb(1,1-R01,0))
#'
#' # transform 2 SN
#' # CLP [H,R]; DP [H,R,Z]
#' SN_DP=XY2SN(CLP,DP)
#'
#' # interpolate using interp of the akima package
#' mesh=interp(SN_DP[,1], SN_DP[,2], SN_DP[,3], duplicate = "strip", nx = 500, ny = 50)
#' # mesh list(S,N,Z)
#'
#' # transform to Gauss Kruger
#' XY=SN2XY(CLP,mesh)
#' # XY matrix [,1]->H; [,2]->R; [,3]->Z
#'
#' # getting rid of NA values that were produced by the interp function
#' XY=na.omit(XY)
#'
#' # plotting result
#' R01=XY[,3]/max(XY[,3])
#' plot(XY[,2],XY[,1],col=rgb(1,1-R01,0))
#'
#' # convert Matrix to Raster to have further plotting options
#' XYZ=Matrix.to.Raster(XY,h.res=50,r.res=500)
#' levelplot(t(XYZ[[3]]),col.regions = colorRampPalette(c("blue","yellow", "red")))
XY2SN=function(ThalwegCoordinates_HR,DepthCoordiantes_HRZ){

  SNZ<<-list()
  for(ind in 1:dim(DepthCoordiantes_HRZ)[1]){
    SNZ[[ind]]<<-SN.Transformation2SN(ind,ThalwegCoordinates_HR,DepthCoordiantes_HRZ)}

  # SNZ=foreach(ind=1:dim(DepthCoordiantes_HRZ)[1],.combine = rbind)%do%
  #   SN.Transformation2SN(ind,ThalwegCoordinates_HR,DepthCoordiantes_HRZ)

  return(matrix(unlist(SNZ), ncol = 3,byrow = TRUE))
}


#' A function for transforming SN to XY coordinates.
#'
#' The transformation needs a dataframe (ThalwegCoordinates_HR) and list containing the grid information (DepthGrid_SNZ) as input.
#' @param ThalwegCoordinates_HR Dataframe with dim = n,2. First postion is the Hochwert (y-coordinate) and second ist the Rechtswert (x-coordinate). This dataset can be estimated by the centerline oy your river.
#' @param DepthGrid_SNZ List containing the grid information. First position is S, second position is N and third position is Z (depth).
#' @return A matrix containing the H-coordinates (y) [,1], R-coordinates (x) [,2] and the depth (z) [,3].
#' @seealso \link[ISI.Toolbox]{XY2SN}
#' @export
#' @examples
#' ### This example incoorporates an entire SN transformation ###
#' ##In order to use the example you need to install and load the akima package as we need it for our grid interpolation##
#' # Create Testdata for this example
#' R=seq(362000,364000)
#' X=R-363000
#' H=X^2*10/80000+20*sin(X/50/pi)+5666100
#' H2=X^2*10/80000+20*sin(X/50/pi)-20+5666100
#' CL=0.5*(H2+H)
#'
#' DPR=runif(n = 1500,min = 362010,max=363990)
#' DPH=c()
#' DPZ=c()
#' for(i in 1:length(DPR)){
#'   index=which.min(abs(R-DPR[i]))
#'   value=rnorm(1,mean=0,sd=5)+CL[index]
#'   d=abs(value-CL[index])
#'   depth=5+d*-5/10
#'   if(value>H[index]){value=H[index]-0.5}
#'   else if(value<H2[index]){value=H2[index]+0.5}
#'   d=abs(value-CL[index])
#'   depth=5+d*-5/10
#'   DPH=c(DPH,value)
#'   DPZ=c(DPZ,depth)
#' }
#'
#' UB=data.frame(H=H,R=R,Z=rep(0,length(R)))
#' LB=data.frame(H=H2,R=R,Z=rep(0,length(R)))
#' Data=data.frame(H=DPH,R=DPR,Z=DPZ)
#' DP=rbind(UB[10:(dim(UB)[1]-10),],LB[10:(dim(LB)[1]-10),],Data)
#' CLP=data.frame(H=CL,R=R)
#'
#' # plotting data before transformation and interpolation
#' XY=DP
#' R01=XY[,3]/max(XY[,3])
#' plot(XY[,2],XY[,1],col=rgb(1,1-R01,0))
#'
#' # transform 2 SN
#' # CLP [H,R]; DP [H,R,Z]
#' SN_DP=XY2SN(CLP,DP)
#'
#' # interpolate using interp of the akima package
#' mesh=interp(SN_DP[,1], SN_DP[,2], SN_DP[,3], duplicate = "strip", nx = 500, ny = 50)
#' # mesh list(S,N,Z)
#'
#' # transform to Gauss Kruger
#' XY=SN2XY(CLP,mesh)
#' # XY matrix [,1]->H; [,2]->R; [,3]->Z
#'
#' # getting rid of NA values that were produced by the interp function
#' XY=na.omit(XY)
#'
#' # plotting result
#' R01=XY[,3]/max(XY[,3])
#' plot(XY[,2],XY[,1],col=rgb(1,1-R01,0))
#'
#' # convert Matrix to Raster to have further plotting options
#' XYZ=Matrix.to.Raster(XY,h.res=50,r.res=500)
#' levelplot(t(XYZ[[3]]),col.regions = colorRampPalette(c("blue","yellow", "red")))
SN2XY=function(ThalwegCoordinates_HR,DepthGrid_SNZ){

  ThalwegCumSum=cumsum(sqrt(diff(ThalwegCoordinates_HR[,1])^2+diff(ThalwegCoordinates_HR[,2])^2))

  S=DepthGrid_SNZ[[1]]
  N=DepthGrid_SNZ[[2]]
  Z=DepthGrid_SNZ[[3]]
  Num_SN=dim(Z)

  HH<<-list()
  RR<<-list()
  ZZ<<-list()
  for(ind in 1:Num_SN[1]){
    XYZ=SN.Transformation2XY(ind,ThalwegCoordinates_HR,S,N,Z,ThalwegCumSum)
    HH[[ind]]<<-XYZ[[1]]
    RR[[ind]]<<-XYZ[[2]]
    ZZ[[ind]]<<-XYZ[[3]]}

  # XYZ=foreach(ind=1:Num_SN[1],.combine = rbind)%do%
  #   SN.Transformation2XY(ind,ThalwegCoordinates_HR,S,N,Z,ThalwegCumSum)

  return(matrix(c(unlist(HH),unlist(RR),unlist(ZZ)), ncol = 3))
}


#' A function for transforming a matrix to a grid.
#'
#' Some plotting tools require a grid. This tool converts a Matrix to a grid (list).
#' @param MA A matrix containing the H-coordinates (y) [,1], R-coordinates (x) [,2] and the depth (z) [,3].
#' @param h.res Integer secifying the resolution of the H coordinate.
#' @param r.res Integer secifying the resolution of the R coordinate.
#' @return A grid (list).
#' @export
#' @examples
#' ### This example incoorporates an entire SN transformation ###
#' ##In order to use the example you need to install and load the akima package as we need it for our grid interpolation##
#' # Create Testdata for this example
#' R=seq(362000,364000)
#' X=R-363000
#' H=X^2*10/80000+20*sin(X/50/pi)+5666100
#' H2=X^2*10/80000+20*sin(X/50/pi)-20+5666100
#' CL=0.5*(H2+H)
#'
#' DPR=runif(n = 1500,min = 362010,max=363990)
#' DPH=c()
#' DPZ=c()
#' for(i in 1:length(DPR)){
#'   index=which.min(abs(R-DPR[i]))
#'   value=rnorm(1,mean=0,sd=5)+CL[index]
#'   d=abs(value-CL[index])
#'   depth=5+d*-5/10
#'   if(value>H[index]){value=H[index]-0.5}
#'   else if(value<H2[index]){value=H2[index]+0.5}
#'   d=abs(value-CL[index])
#'   depth=5+d*-5/10
#'   DPH=c(DPH,value)
#'   DPZ=c(DPZ,depth)
#' }
#'
#' UB=data.frame(H=H,R=R,Z=rep(0,length(R)))
#' LB=data.frame(H=H2,R=R,Z=rep(0,length(R)))
#' Data=data.frame(H=DPH,R=DPR,Z=DPZ)
#' DP=rbind(UB[10:(dim(UB)[1]-10),],LB[10:(dim(LB)[1]-10),],Data)
#' CLP=data.frame(H=CL,R=R)
#'
#' # plotting data before transformation and interpolation
#' XY=DP
#' R01=XY[,3]/max(XY[,3])
#' plot(XY[,2],XY[,1],col=rgb(1,1-R01,0))
#'
#' # transform 2 SN
#' # CLP [H,R]; DP [H,R,Z]
#' SN_DP=XY2SN(CLP,DP)
#'
#' # interpolate using interp of the akima package
#' mesh=interp(SN_DP[,1], SN_DP[,2], SN_DP[,3], duplicate = "strip", nx = 500, ny = 50)
#' # mesh list(S,N,Z)
#'
#' # transform to Gauss Kruger
#' XY=SN2XY(CLP,mesh)
#' # XY matrix [,1]->H; [,2]->R; [,3]->Z
#'
#' # getting rid of NA values that were produced by the interp function
#' XY=na.omit(XY)
#'
#' # plotting result
#' R01=XY[,3]/max(XY[,3])
#' plot(XY[,2],XY[,1],col=rgb(1,1-R01,0))
#'
#' # convert Matrix to Raster to have further plotting options
#' XYZ=Matrix.to.Raster(XY,h.res=50,r.res=500)
#' levelplot(t(XYZ[[3]]),col.regions = colorRampPalette(c("blue","yellow", "red")))
Matrix.to.Raster=function(MA,h.res=250,r.res=250){
  MA=na.omit(MA)
  MA_H=seq(min(MA[,1]),max(MA[,1]),length.out = h.res)
  MA_R=seq(min(MA[,2]),max(MA[,2]),length.out = r.res)
  Z<<-matrix(NA, h.res, r.res)
  i<<-0
  Loop =function(h){
    i<<-i+1
    if(i%%(h.res/100)==0){cat(paste(toString(i/h.res*100),"%\n"))}
    R_h<<-MA_H[(h:(h+1))]
    Loop2=function(r){
      R_r=MA_R[(r:(r+1))]
      bool_h=MA[,1]>=R_h[1]&MA[,1]<=R_h[2]
      bool_r=MA[,2]>=R_r[1]&MA[,2]<=R_r[2]
      boo=bool_h&bool_r
      if(any(boo)){
        Z[h,r]<<-mean(MA[,3][boo])
      }
    }
    sapply(1:(r.res-1),FUN=Loop2)
  }
  sapply(1:(h.res-1),FUN=Loop)
  ma=list(MA_H,MA_R,Z)
  names(ma)=c("y","x","z")
  return(ma)
}

