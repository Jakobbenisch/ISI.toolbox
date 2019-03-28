
#' A Geometry Defining Function
#'
#' This function analytically defines the relation between height and width of a standard Maulprofil.\cr
#' Mind the units!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h Decimal [m]. Representing the y coordinate to which the width will be returned.
#' @param r Decimal [m]. Radius of the Maulprofil.
#' @seealso \link[ISI.Toolbox]{ManningMaul}
#' @return Returns width of profile corresponding to input height [m].
#' @keywords Manning
#' @export
#' @examples
#' width=MaulProfil(.5,1)
#'
MaulProfil=function(h,r){
  if(h<=1/4*r){
    x=sqrt(4*r^2-(h-2*r)^2)}
  else if(h<=1.5*r){
    x=sqrt(r^2-(h-0.5*r)^2)}
  else{warning("The given height exceeds the pipe's heigt")}
  return(x)}


#' An Flow Calculation Function
#'
#' This function numerically calculates the discharge for a single waterlevel. Using the MaulProfil function.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h Decimal [m]. Water level.
#' @param r Decimal [m]. Radius of the Maulprofil.
#' @param I Decimal [-]. Slope of flow cross section.
#' @param kst Decimal [m^{1/3}/s]. Friction coefficient according to Manning Strickler. (Rauheitsbeiwert nach Strickler).
#' @param increment Decimal [m]. This value defines the height of the trapezoids used for approximating the Maulprofil. The result's precision rises the smaller the increment value. 0.001 m is the default value.
#' @seealso \link[ISI.Toolbox]{ManningMaul}
#' @return Returns discharge corresponding to input height [m^3/s].
#' @keywords Manning
#' @export
#' @examples
#' Q=IterMaul(h=0.5,r=1,I=0.001,kst=100)
#'
IterMaul=function(h,I,kst,r,increment=0.001){
  h_vec=seq(0,h,increment)
  b_vec=sapply(h_vec,FUN=function(x) MaulProfil(x,r))
  A=2*sum(rollapply(b_vec,width=2,by=1,FUN=function(x){increment*(x[2]+x[1])/2}))
  L_u=2*sum(rollapply(b_vec,width=2,by=1,FUN=function(x){sqrt((x[2]-x[1])^2+increment^2)}))
  Q=A*sqrt(I)*(A/L_u)^(2/3)*kst
  return(Q)}

#' A Flow Calculation Function (Maulprofil)
#'
#' This function numerically calculates the discharge for a vector of water levels. Using the MaulProfil and the IterMaul functions.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h numeric vector containing water levels [m].
#' @param r Decimal [m]. Radius of the Maulprofil.
#' @param I Decimal [-]. Slope of flow cross section.
#' @param kst Decimal [m^{1/3}/s]. Friction coefficient according to Manning Strickler. (Rauheitsbeiwert nach Strickler).
#' @param increment Decimal [m]. This value defines the height of the trapezoids used for approximating the Maulprofil. The result's precision rises the smaller the increment value. 0.001 m is the default value.
#' @seealso \link[ISI.Toolbox]{ManningEGG},\link[ISI.Toolbox]{ManningCIR},\link[ISI.Toolbox]{ManningREC},\link[ISI.Toolbox]{ManningTRA},\link[ISI.Toolbox]{ManningSchacht}
#' @return Returns numerical vector containing discharge values [m^3/s].
#' @keywords Manning
#' @export
#' @examples
#' H=runif(100,0,1.5)
#' Q=ManningMaul(h=H,r=1,I=0.001,kst=100)
#'
ManningMaul=function(h,I,kst,r,increment=0.001){
  Q=sapply(as.numeric(h),FUN = function(x) IterMaul(h=x,I=I,kst=kst,r=r,increment = increment))
  return(Q)
}


#' An Area Defining Function (Egg Profile)
#'
#' This function analytically calculates the area of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param w Decimal [m]. Water level of the egg profile.
#' @param r Decimal [m]. Radius of the egg profile.
#' @seealso \link[ISI.Toolbox]{ManningEGG}
#' @return Returns the area of flow cross section [m^2].
#' @keywords Manning
#' @export
#' @examples
#' A=AEGG(w=1,r=1)
#'
AEGG=function(w,r){
  w=round(w,3)
  if(w<=0.2*r){
    Area=(((w-2.*r)+3/2.*r)*((r/2.)^2.-((w-2.*r)+3/2.*r)^2. )^0.5+((r/2.)^2.)*asin(((w-2.*r)+3/2.*r)/(r/2.))+(pi/2.)*(r/2.)^2.)}
  else if(w<=2*r){
    Area=((w-2.*r)*((3.*r)^2.-(w-2.*r)^2.)^0.5+((3.*r)^2.)*asin((w-2.*r)/(3.*r))-4.*r*(w-2.*r)+9/5.*r*((3.*r)^2.-(9/5.*r)^2. )^0.5-((3.*r)^2.)*asin(-3/5.)-36/5.*(r^2.) + ((0.2*r-2.*r)+3/2.*r)*((r/2.)^2.-((0.2*r-2.*r)+3/2.*r)^2. )^0.5+((r/2.)^2.)*asin(((0.2*r-2.*r)+3/2.*r)/(r/2.))+(pi/2.)*(r/2.)^2.)}
  else{
    Area=((w-2.*r)*(r^2.-(w-2.*r)^2.)^0.5+(r^2.)*asin((w-2.*r)/r)+(2.*r-2.*r)*((3.*r)^2.-(2.*r-2.*r)^2.)^0.5+((3.*r)^2.)*asin((2.*r-2.*r)/(3.*r))-4.*r*(2.*r-2.*r)+9/5.*r*((3.*r)^2.-(9/5.*r)^2. )^0.5-((3.*r)^2.)*asin(-3/5.)-36/5.*(r^2.)+((0.2*r-2.*r)+3/2.*r)*((r/2.)^2.-((0.2*r-2.*r)+3/2.*r)^2. )^0.5+((r/2.)^2.)*asin(((0.2*r-2.*r)+3/2.*r)/(r/2.))+(pi/2.)*(r/2.)^2.)}
  return(Area)}

#' A Wetted Perimeter Defining Function (Egg Profile)
#'
#' This function analytically calculates the wetted perimeter of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param w Decimal [m]. Water level of the egg profile.
#' @param r Decimal [m]. Radius of the egg profile.
#' @seealso \link[ISI.Toolbox]{ManningEGG}
#' @return Returns the wetted perimeter of flow cross section [m].
#' @keywords Manning
#' @export
#' @examples
#' Lu=LuEGG(w=1,r=1)
#'
LuEGG=function(w,r){
  w=round(w,3)
  if(w<=0.2*r){
    L=(2.*(asin(2.*(w-2.*r)/r+3.)+pi/2.)*r/2.)}
  else if(w<=2*r){
    L=(2.*(asin((w-2.*r)/(3.*r))-asin(-3/5.))*3.*r)+(2.*(asin(2.*(0.2*r-2.*r)/r+3.)+pi/2.)*r/2.)}
  else{
    L=(2.*asin((w-2.*r)/r)*r)+(2.*(asin((2.*r-2.*r)/(3.*r))-asin(-3/5.))*3.*r)+(2.*(asin(2.*(0.2*r-2.*r)/r+3.)+pi/2.)*r/2.)}
  return(L)}


#' An Area Defining Function (Circular Profile)
#'
#' This function analytically calculates the area of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param w Decimal or numeric vector [m]. Water level of the circular profile..
#' @param r Decimal [m]. Radius of the circular profile.
#' @seealso \link[ISI.Toolbox]{ManningCIR}
#' @return Returns the area of flow cross section [m^2].
#' @keywords Manning
#' @export
#' @examples
#' A=ACIR(w=1,r=1)
#'
ACIR=function(w,r){
  h=as.numeric(w)
  alpha=acos((r-h)/r)
  A= r^2*alpha - (r-h)*sin(alpha)*r
  return(A)}

#' A Wetted Perimeter Defining Function (Circular Profile)
#'
#' This function analytically calculates the wetted perimeter of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param w Decimal [m]. Water level of the circular profile.
#' @param r Decimal [m]. Radius of the circular profile.
#' @seealso \link[ISI.Toolbox]{ManningCIR}
#' @return Returns the wetted perimeter of flow cross section [m].
#' @keywords Manning
#' @export
#' @examples
#' Lu=LuEGG(w=1,r=1)
#'
LuCIR=function(w,r){
  w=as.numeric(w)
  alpha=acos((r-w)/r)
  L=alpha*r*2
  return(L)}


#' An Area Defining Function (Rectangular Profile)
#'
#' This function analytically calculates the area of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param w Decimal or numeric vector [m]. Water level of the rectangular profile.
#' @param b Decimal [m]. Width of the rectangular profile.
#' @seealso \link[ISI.Toolbox]{ManningREC}
#' @return Returns the area of flow cross section [m^2].
#' @keywords Manning
#' @export
#' @examples
#' A=AREC(w=1,b=1)
#'
AREC=function(w,b){
  w=as.numeric(w)
  A=w*b
  return(A)}

#' A Wetted Perimeter Defining Function (Rectangular Profile)
#'
#' This function analytically calculates the wetted perimeter of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param w Decimal [m]. Water level of the rectangular profile.
#' @param b Decimal [m]. Width of the rectangular profile.
#' @seealso \link[ISI.Toolbox]{ManningREC}
#' @return Returns the wetted perimeter of flow cross section [m].
#' @keywords Manning
#' @export
#' @examples
#'Lu=LuREC(w=1,b=1)
#'
LuREC=function(w,b){
  w=as.numeric(w)
  L=2*w+b
  return(L)}


#' An Area Defining Function (Trapezoid Profile)
#'
#' This function analytically calculates the area of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)\cr
#' Also the angle should be between 0 and pi/2. Excluding pi/2!
#' @param w Decimal or numeric vector [m].
#' @param b Decimal [m]. Width of bottom part of the trapezoid profile.
#' @param alpha Decimal in rad. Angle of slope -> right (mathematically negative). Alpha=0 results in a vertical line.
#' @param beta Decimal in rad. Angle of slope -> left (mathematically positive). Beta=0 results in a vertical line.
#' @seealso \link[ISI.Toolbox]{ManningTRA}
#' @return Returns the area of flow cross section [m^2].
#' @keywords Manning
#' @export
#' @examples
#' A=ATRA(w=1,b=1,alpha=pi/4,beta=pi/5)
#'
ATRA=function(w,b,alpha,beta){
  w=as.numeric(w)
  A=b*w+0.5*tan(alpha)*w^2+0.5*tan(beta)*w^2
  return(A)}

#' A Wetted Perimeter Defining Function (Trapezoid Profile)
#'
#' This function analytically calculates the wetted perimeter of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param w Decimal or numeric vector [m].
#' @param b Decimal [m]. Width of bottom part of the trapezoid profile.
#' @param alpha Decimal in rad. Angle of slope -> right (mathematically negative). Alpha=0 results in a vertical line.
#' @param beta Decimal in rad. Angle of slope -> left (mathematically positive). Beta=0 results in a vertical line.
#' @seealso \link[ISI.Toolbox]{ManningTRA}
#' @return Returns the wetted perimeter of flow cross section [m].
#' @keywords Manning
#' @export
#' @examples
#' Lu=LuTRA(w=1,b=1,alpha=pi/4,beta=pi/5)
#'
LuTRA=function(w,b,alpha,beta){
  w=as.numeric(w)
  L=w/cos(alpha)+b+w/cos(beta)
  return(L)}



#' A Flow Calculation Function (Egg Profile)
#'
#' This function analytically calculates the discharge for a vector of water levels. Using the AEGG and the LuEGG functions.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h numeric vector containing water levels [m].
#' @param r Decimal [m]. Radius of the egg profile.
#' @param I Decimal [-]. Slope of flow cross section.
#' @param kst Decimal [m^{1/3}/s]. Friction coefficient according to Manning Strickler. (Rauheitsbeiwert nach Strickler).
#' @seealso \link[ISI.Toolbox]{ManningMaul},\link[ISI.Toolbox]{ManningCIR},\link[ISI.Toolbox]{ManningREC},\link[ISI.Toolbox]{ManningTRA},\link[ISI.Toolbox]{ManningSchacht}
#' @return Returns numerical vector containing discharge values [m^3/s].
#' @keywords Manning
#' @export
#' @examples
#' H=runif(100,0,2.1)
#' Q=ManningEGG(h=H,r=1,I=0.001,kst=100)
#'
ManningEGG=function(h,r,kst,I){
  A=sapply(h,FUN=function(x) AEGG(x,r))
  L_u=sapply(h,FUN=function(x) LuEGG(x,r))
  Rhy=A/L_u
  Q=A*kst*(Rhy)^(2/3.)*(I^(1/2.))
  return(Q)}

#' A Flow Calculation Function (Circular Profile)
#'
#' This function analytically calculates the discharge for a vector of water levels. Using the ACIR and the LuCIR functions.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h numeric vector containing water levels [m].
#' @param r Decimal [m]. Radius of the Maulprofil
#' @param I Decimal [-]. Slope of flow cross section.
#' @param kst Decimal [m^{1/3}/s]. Friction coefficient according to Manning Strickler. (Rauheitsbeiwert nach Strickler).
#' @seealso \link[ISI.Toolbox]{ManningMaul},\link[ISI.Toolbox]{ManningEGG},\link[ISI.Toolbox]{ManningREC},\link[ISI.Toolbox]{ManningTRA},\link[ISI.Toolbox]{ManningSchacht}
#' @return Returns numerical vector containing discharge values [m^3/s].
#' @keywords Manning
#' @export
#' @examples
#' H=runif(100,0,2)
#' Q=ManningCIR(h=H,r=1,I=0.001,kst=100)
#'
ManningCIR=function(h,r,kst,I){
  Rhy=ACIR(h,r)/LuCIR(h,r)
  Q=(ACIR(h,r)*kst*(Rhy)^(2/3.)*(I^(1/2.)))
  return(Q)}

#' A Flow Calculation Function (Rectangular  Profile)
#'
#' This function analytically calculates the discharge for a vector of water levels. Using the AREC and the LuREC functions.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h numeric vector containing water levels [m].
#' @param b Decimal [m]. Width of rectangular profile.
#' @param I Decimal [-]. Slope of flow cross section.
#' @param kst Decimal [m^{1/3}/s]. Friction coefficient according to Manning Strickler. (Rauheitsbeiwert nach Strickler).
#' @seealso \link[ISI.Toolbox]{ManningMaul},\link[ISI.Toolbox]{ManningEGG},\link[ISI.Toolbox]{ManningCIR},\link[ISI.Toolbox]{ManningTRA},\link[ISI.Toolbox]{ManningSchacht}
#' @return Returns numerical vector containing discharge values [m^3/s].
#' @keywords Manning
#' @export
#' @examples
#' H=runif(100,0,1)
#' Q=ManningREC(h=H,b=2,I=0.001,kst=100)
#'
ManningREC=function(h,b,kst,I){
  Rhy=AREC(h,b)/LuREC(h,b)
  Q=(AREC(h,b)*kst*(Rhy)^(2/3.)*(I^(1/2.)))
  return(Q)}


#' A Flow Calculation Function (Trapezoid  Profile)
#'
#' This function analytically calculates the discharge for a vector of water levels. Using the ATRA and the LuTRA functions.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)\cr
#' Also the angle should be between 0 and pi/2. Excluding pi/2!
#' @param h numeric vector containing water levels [m].
#' @param b Decimal [m]. Width of bottom part of the trapezoid profile.
#' @param alpha Decimal in rad. Angle of slope -> right (mathematically negative). Alpha=0 results in a vertical line.
#' @param beta Decimal in rad. Angle of slope -> left (mathematically positive). Beta=0 results in a vertical line.
#' @param I Decimal [-]. Slope of flow cross section.
#' @param kst Decimal [m^{1/3}/s]. Friction coefficient according to Manning Strickler. (Rauheitsbeiwert nach Strickler).
#' @seealso \link[ISI.Toolbox]{ManningMaul},\link[ISI.Toolbox]{ManningEGG},\link[ISI.Toolbox]{ManningCIR},\link[ISI.Toolbox]{ManningREC},\link[ISI.Toolbox]{ManningSchacht}
#' @return Returns numerical vector containing discharge values [m^3/s].
#' @keywords Manning
#' @export
#' @examples
#' H=runif(100,0,1)
#' Q=ManningTRA(h=H,b=2,I=0.001,kst=100,alpha=pi/3,beta=pi/3)
#'
ManningTRA=function(h,b,alpha,beta,kst,I){
  Rhy=ATRA(h,b,alpha,beta)/LuTRA(h,b,alpha,beta)
  Q=(ATRA(h,b,alpha,beta)*kst*(Rhy)^(2/3.)*(I^(1/2.)))
  return(Q)}



#' An Area Defining Function (Manhole)
#'
#' This function analytically calculates the area of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)\cr
#' @param h Decimal[m]. Water levels.
#' @param b Decimal [m]. Width of manhole.
#' @param r Decimal [m]. Radius of circular bottom part of manhole.
#' @param hlimit Decimal [m]. Height of transition between circular shaped bottom and straight wall.
#' @seealso \link[ISI.Toolbox]{ManningSchacht}
#' @return Returns the area of flow cross section [m^2].
#' @keywords Manning
#' @export
#' @examples
#' A=A_Schacht(h=1,r=0.6,b=1.47,hlimit=0.21)
#'
A_Schacht=function(h,r=0.6,b=1.47,hlimit=0.21) {
  if (h<=hlimit){
    A=ACIR(w=h,r=r)}
  else{
    A=ACIR(w=hlimit,r=r)+(h-hlimit)*b}
  return(A)
}

#' A Wetted Perimeter Defining Function (Manhole)
#'
#' This function analytically calculates the wetted perimeter of a flow cross section for a water level.\cr
#' Mind the units!!\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h Decimal[m]. Water levels.
#' @param b Decimal [m]. Width of manhole.
#' @param r Decimal [m]. Radius of circular bottom part of manhole.
#' @param hlimit Decimal [m]. Height of transition between circular shaped bottom and straight wall.
#' @seealso \link[ISI.Toolbox]{ManningSchacht}
#' @return Returns the wetted perimeter of flow cross section [m].
#' @keywords Manning
#' @export
#' @examples
#' Lu=Lu_Schacht(h=1,r=0.6,b=1.47,hlimit=0.21)
#'
Lu_Schacht=function(h,r=0.6,b=1.47,hlimit=0.21) {
  if (h<=hlimit){
    Lu=LuCIR(w=h,r=r)}
  else{
    winkel=acos((r-hlimit)/r)
    Lu=LuCIR(w=hlimit,r=r)+b-2*(sin(winkel)*r)+2*(h-hlimit)}
  return(Lu)
}

#' A Flow Calculation Function (For a Manhole)
#'
#' This function analytically calculates the discharge for a vector of water levels. Using the A_Schacht and the Lu_Schacht functions.\cr
#' Mind the units!!\cr
#' As backwater effects are likely for this cross section (depending on the following) the results will be less accurate with rising water level.\cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param h numeric vector containing water levels [m].
#' @param b Decimal [m]. Width of manhole.
#' @param r Decimal [m]. Radius of circular bottom part of manhole.
#' @param hlimit Decimal [m]. Height of transition between circular shaped bottom and straight wall.
#' @param I Decimal [-]. Slope of flow cross section.
#' @param kst Decimal [m^{1/3}/s]. Friction coefficient according to Manning Strickler. (Rauheitsbeiwert nach Strickler).
#' @seealso \link[ISI.Toolbox]{ManningMaul},\link[ISI.Toolbox]{ManningEGG},\link[ISI.Toolbox]{ManningCIR},\link[ISI.Toolbox]{ManningREC},\link[ISI.Toolbox]{ManningTRA}
#' @return Returns numerical vector containing discharge values [m^3/s].
#' @keywords Manning
#' @export
#' @examples
#' H=runif(100,0,1.6)
#' Q=ManningSchacht(h=H,b=2,r=0.8,hlimit=0.21,I=0.001,kst=100)
#'
ManningSchacht=function(h,kst,I,r=0.6,b=1.47,hlimit=0.21){
  DO=function(h){
    Rhy=A_Schacht(h=h,r=r,b=b,hlimit=hlimit)/Lu_Schacht(h=h,r=r,b=b,hlimit=hlimit)
    Q=(A_Schacht(h=h,r=r,b=b,hlimit=hlimit)*kst*(Rhy)^(2/3.)*(I^(1/2.)))
    return(Q)}

  Q=sapply(h,FUN=DO)

  return(Q)}


#' A function for fitting Manning Strickler to a h-Q-realtion.
#'
#' This tool uses the FME modFit function to fit Manning Strickler to h,Q data. There are different fitting methods available (\link[FME]{modFit}).\cr
#' Mind the units, when using this function!!\cr
#' In order to get a realistic fit for the kst and the I parameters of Manning please give plausible lower and upper boundries. \cr
#' If you get errors using this function please make sure no NA values are within the water level vector and the waterlevels are plausible. (No negative water levels or water levels higher than the profile.)
#' @param Geometry String that specifies the channels geometry. Valid options are 'EGG','CIR','REC' or 'TRA'.Default is "EGG".
#' @param kst_I_start Numerical vector of length 2. The fitting algorithm requires all parameters to be in a vector. First position kst value [m^{1/3}/s] and second position I [-] (Slope).
#' @param kst_I_lwr Numerical vector of length 2. Lower boudries for kst and I.
#' @param kst_I_upr Numerical vector of length 2. Upper boudries for kst and I.
#' @param modelfitmethod String specifying the fitting algorithm to use. Options can be found under \link[FME]{modFit} for the parameter method.
#' @param h_vec Numerical vector containing the water level [m].
#' @param Q_vec Numerical vector containing the flow [m^{3}/s].
#' @param width Number. Width of geometry. Only Relevant for "REC" and "TRA".
#' @param radius Number. Radius of geometry. Only Relevant for "EGG" and "CIR".
#' @param alpha Decimal in rad. Angle of slope -> right (mathematically negative). Alpha=0 results in a vertical line. Only relevant for "TRA".
#' @param beta Decimal in rad. Angle of slope -> left (mathematically positive). Beta=0 results in a vertical line. Only relevant for "TRA".
#' @seealso \link[FME]{modFit},\link[ISI.Toolbox]{ManningMaul},\link[ISI.Toolbox]{ManningEGG},\link[ISI.Toolbox]{ManningCIR},\link[ISI.Toolbox]{ManningREC},\link[ISI.Toolbox]{ManningTRA}
#' @return Returns the modFit object from the FME package. Use $par to extract the fitted parameters. You can also use indexing, depending on the modelfitmethod the first or the second postions contains the parameters.
#' @keywords Manning, Regression, h-Q
#' @export
#' @examples
#' # creating H,Q with set values for I=0.001 and kst=100
#' H=runif(100,0,1.6)
#' Q=ManningEGG(h=H,r=0.8,I=0.001,kst=100)
#'
#' # fitting Manning to data
#' Mod=FitManning(Geometry = "EGG",kst_I_start = c(80,0.0008), kst_I_lwr = c(70,0.0005), kst_I_upr = c(110,0.003),h_vec = H,Q_vec = Q,r=0.8)
#'
#' # extraction the kst and the I value:
#' kstI=Mod$par
#'
#' # sorting the waterlevel values in order to get a nice line plot
#' Hfit=sort(H)

#' # calculating Q with the fitted parameters
#' Qfit=ManningEGG(h=Hfit,r=0.8,I=kstI[2],kst=kstI[1])
#'
#' # plotting original data (black dots) and the fitted data (red line)
#' plot(H,Q)
#' lines(Hfit,Qfit,col=2)
#'
#' # As shown in this example, there are multiple solutions for the kst and I paramters.
#' # Although the values fitted aren't the same as the ones given for this example, the function lies within the data points perfectly.
FitManning=function(Geometry=c("EGG","CIR","REC","TRA"),kst_I_start=c(80,0.001),kst_I_lwr=c(50,0.0001),kst_I_upr=c(150,0.01),modelfitmethod="BFGS",h_vec,Q_vec,width,radius,alpha,beta){
  #set default
  if(length(Geometry)==4){Geometry="EGG"}

  if(Geometry=="EGG"){
    ManningEGGSolveRES=function(kstI,h,r,Q){
      q=ManningEGG(kst=kstI[1],I=kstI[2],h=h,r=r)
      return((Q-q))}

    Model=modFit( ManningEGGSolveRES,kst_I_start, h=h_vec,r=radius, Q=Q_vec, method = modelfitmethod,lower = kst_I_lwr,upper=kst_I_upr)

  }else if(Geometry=="CIR"){
    ManningCIRSolveRES=function(kstI,h,r,Q){
      q=ManningCIR(kst=kstI[1],I=kstI[2],h=h,r=r)
      return((Q-q))}

    Model=modFit( ManningCIRSolveRES,kst_I_start, h=h_vec,r=radius,Q=Q_vec,method = modelfitmethod,lower = kst_I_lwr,upper=kst_I_upr)

  }else if(Geometry=="REC"){
    ManningRECSolveRES=function(kstI,h,b,Q){
      q=ManningREC(kst=kstI[1],I=kstI[2],h=h,b=b)
      return((Q-q))}

    Model=modFit( ManningRECSolveRES,kst_I_start, h=h_vec,b=width,Q=Q_vec,method = modelfitmethod,lower = kst_I_lwr,upper=kst_I_upr)

  }else if(Geometry=="TRA"){
    ManningTRASolveRES=function(kstI,h,b,alph,bet,Q){
      q=ManningTRA(kst=kstI[1],I=kstI[2],h=h,b=b,alpha = alph,beta = bet)
      return((Q-q))}

    Model=modFit( ManningTRASolveRES,kst_I_start, h=h_vec,b=width,alph=alpha,bet=beta,Q=Q_vec,method = modelfitmethod,lower = kst_I_lwr,upper=kst_I_upr)

  }else{
    stop("The Geometry you have chosen is not implemented. Valid options are 'EGG','CIR','REC' or 'TRA'.")
    Model=NULL
  }

  return(Model)
}
