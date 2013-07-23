#---------------------------------------------------------------------------#
#                                                                           #
# Cellular Dune dynamics model                                              #
#                                                                           #
# Copyright (C) 2013 Fernando Barrio Parra & Inmaculada Rodriguez Santalla  #
#                                                                           #
# This program is free software; you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with GNU gv; see the file COPYING.  If not, see                     #
# <http://www.gnu.org/licenses/>.                                           #
#                                                                           #
#   Author: Fernando Barrio Parra      Biology and Geology department       #
#    email: fernando.barrio@urjc.es    Universidad Rey Juan Carlos          #
#                                      Tulipán s/n 28933, Móstoles          #           
#                                      Madrid (Spain)                       #
#                                                                           #
#                                                                           #
#---------------------------------------------------------------------------#





##Librarys to load before work
library(lattice)
library(caTools)



###############################################-User data entry-########################################################

#Enter the path were dune DEMs and wind data are and were you want to save your results

  Data_Path<-"C:\\Users\\fernando\\Documents\\ModelData" 
  Results_Path<-"C:\\Users\\fernando\\Documents\\ModelResults"


#Enter the names of the files in wich are the initial and final dune DEMs and the decimal indicator

  Initial_Dune_file<-"20120415ascii.txt"
  Final_Dune_file<-"20120418ascii.txt"
  Decimal_Dune_Files<-"," #Choose between "," or "."

#Enter the X and Y cell dimensions (in meters), the basement height and the Repose Angle (RA) in degrees

  Xcellsize<-0.7 
  Ycellsize<-0.7
  Basement_Height<-0.6
  RA<-35

#Enter the name of the Wind data file, decimal indicator and columns were wind intensity (m·s^-1) and direction (º) are

  Wind_data_file<-"WindData.txt"
  Decimal_Wind_File<-"." #Choose between "," or "."
  Wind_intensity_column<-2
  Wind_direction_column<-3

#Enter the wind correlation function between the meteorological station and the dune field site, and the wind aeolian transport correlation
  
  Wind_Corr_Function<-function(V){1.435*V-1.264} #Sánchez-García (2008) 
  Transport_Function<-function(V){0.111*exp(abs(V)*0.434)}#Sánchez-García (2008) (Kg·m-1·h-1)


#Enter the number of records to average the wind data by. In our case, each wind measurement is representative of 10 min

  n<-10

#Enter the transport wind velocity threshold Vt(m·s^(-1))

  Vt<-6.45 

#Enter sand parameters

  Rs<-2650              #Sand grains density (kg/m3)
  Ps<-0.4               #Sand porosity (m3 air/m3 soil)

#Enter here the list of phenomenological variables to enter in the simulations and the number of the starting simulation. 
#Simulations will be stored in files wich names will be startSim+(n-1), where n is the position of the phenomenological variable
#list position

alist<-c(0.5)
blist<-c(1)
clist<-c(0)
dlist<-c(2)
elist<-c(2)
flist<-c(0)

startSim<-1

#Especify the kind of Trapping

Trapping<-0 #Choose between 0=No trapping, 1=Traping in the wind ward side, 2=Trapping in the lee side, 3= Trapping in all the dune
            #4= Adding slabs in the north and west edges
###############################################-Data preprocessing-###############################################


##We set the directory were the data is

  setwd(Data_Path)

#####-Dune DEMs preprocesing-#####

#We Read the DEMs files

  MAT<-as.matrix(read.table(Initial_Dune_file,dec=Decimal_Dune_Files,skip=6));for (i in 1:length(MAT)){ if(MAT[i]<0) MAT[i]<-Basement_Height}
  MATfin<-as.matrix(read.table(Final_Dune_file,dec=Decimal_Dune_Files,skip=6));for (i in 1:length(MATfin)){ if(MATfin[i]<0) MATfin[i]<-Basement_Height}
  AsciiHeader<-as.matrix(readLines(Initial_Dune_file,n=6))

#We obtain the array X and Y dimensions

  numcol<-dim(MAT)[1] #Number of columns
  numrow<-dim(MAT)[2] #Number of rows

#We correct the DEMs height by the basement height
  
  MAT<-MAT-Basement_Height
  MATfin<-MATfin-Basement_Height

#We estimate the Initial and Final sediment volume (m^3), and the net sediment budget
  
  Vini<-sum(MAT)*Xcellsize*Ycellsize     #Initial dune volume
  Vfin<-sum(MATfin)*Xcellsize*Ycellsize  #Final dune volume
  Bal<-Vfin-Vini                         #Sediment net budget

#We recalculate the repose angle (RA) in Radians                               

  RA<-RA*pi/180                         

#####-Wind data preprocesing-####

#We Read the Wind data file, extract the intensity and direction, and recalculate it to make it match withe reference 
#system. We apply the correlation function 

  Wind_Data<-read.table(Wind_data_file,dec=Decimal_Wind_File,h=T)
  Wind_Intensity<-Wind_Corr_Function(Wind_Data[,Wind_intensity_column])
  Wind_Direction<-Wind_Data[,Wind_direction_column]
  Wind_Directiona<-ifelse((Wind_Direction-180)>0,Wind_Direction-180,Wind_Direction-180+360) 
  
##We apply the wind averaging algorithm 

#The auxiliar variables are  
  s<-seq(1,length(Wind_Intensity),n) 
  WIa<-numeric(length(s))
  WDa<-numeric(length(s))
  WD<-numeric(length(s))

#The averaging loop algorithm is


    for (i in 1:length(s))
    {
    WIa[i]<-mean(Wind_Intensity[s[i]:(s[i]+9)])
    WDa[i]<-mean(Wind_Directiona[s[i]:(s[i]+9)])
    WD[i]<-mean(Wind_Direction[s[i]:(s[i]+9)])
    }
  WIa[length(s)]<-mean(Wind_Intensity[max(s):length(s)])
  WDa[length(s)]<-mean(Wind_Directiona[max(s):length(s)])
  WD[length(s)]<-mean(Wind_Direction[max(s):length(s)])

#We join the wind intensity and direction and filter the records by the transport wind velocity threshold

  wind<-cbind(WIa,WDa,WD)
  wind<-wind[(wind[,1])>=Vt,]

#we stablish the model time steps as the number of wind records bigger than Vt

  steps<-length(wind[,1])

#We estimate the aeolian transport in terms of cell height this way.

  q<-Transport_Function(wind[,1])   
  q<-q/Rs                                                      
  q<-q/(1-Ps)                                                  
  q<-q*Xcellsize                        
  q<-n*q/6                              
  q<-q/(Xcellsize*Ycellsize)            

#Here we estimate the saltation lenght this way
l0<-wind[,1]/Vt
l0x<-l0*sin(wind[,2]*pi/180) #X axis saltation length, postivive to the east 
l0y<--l0*cos(wind[,2]*pi/180)#Y axis saltation length, postivive to the south

#Parameters needed for hillshade calculations

altitude<-0
Zenith_rad<-(90-altitude)*pi/180 #Radian angle indicating the position of the illumination measured from the zenith (perpendicular to the surface)
Azimuth<-ifelse(wind[,3]+180>360,wind[,3]-180,wind[,3]+180)  #Geographical angle from "wind illumination" comes from. 
Azimuth_math<-360-Azimuth+90
Azimuth_rad<-ifelse(Azimuth_math>=360,Azimuth_math-360,Azimuth_math)
Azimuth_rad<-Azimuth_math*pi/180

###############################################-Auxiliary variables definition-###############################################

#Saltation algorithm auxiliar variables

T<-seq(1:steps)                                              #Time steps vector
DEM<-array(0,dim=c(numcol,numrow,length(T)))                 #Array to store the DEMs from each time step
DEM[,,1]<-MAT                                                #Here we introduce the DEM initial condition
DEMz1<-array(0,dim=c(numcol,numrow,length(T)))               #In this array we subtract the eroded volume to each cell in the DEM
DZ2<-DEMz1                                                   #In this array we displace the eroded cells downwind
qs<-DEMz1                                                    #In this array we estimate the volume to be ereoded in each cell
qd<-qs                                                       #In this array we estimate the trapped sand
lsx<-DEMz1                                                   #Here we store the saltation lenght in the X axis
lsy<-DEMz1                                                   #Here we store the saltation lenght in the Y axis
DEMFin<-array(0,dim=c(numcol,numrow,1))                      #Array to store the surveyed final DEM
DEMFin[,,1]<-MATfin                                          #Here we introduce the DEM final condition

#Slope, aspect and hillshade algorithms auxiliar variables

Z<-array(min(MAT),dim=c(numcol+2,numrow+2,length(T)))        #Array to store DEMs with suitable dimensions for this algorithms
DZDX<-array(0,dim=c(numcol,numrow,length(T)))                #Array to store the local slope in the X axis
DZDY<-DZDX                                                   #Array to store the local slope in the Y axis
aspect<-array(0,dim=c(numcol,numrow,length(T)))              #Array to store the aspect
slope<-array(0,dim=c(numcol,numrow,length(T)))               #Array to store the slope
hillshade<-array(0,dim=c(numcol,numrow,length(T)))           #Array to store the hillshade
Cond_Hillshade<-array(0,dim=c(numcol,numrow,length(T)))      #Array to store the decisions about cell exposition

###Avalanche auxiliar variables

jj<-c(-1,0,1,-1,0,1,-1,0,1)
ii<-c(-1,-1,-1,0,0,0,1,1,1)

###Color ramp for results representation
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                                 "yellow", "#FF7F00", "red", "#7F0000"))

##############################-Slope, Aspect and hillshade application to the initial condition-###################################################
  #Generation of a progress bar



  pb <- winProgressBar(title = "Progress bar", min = 0, max = (numrow), width = 300)

  
  
  
  
for (i in 1:numrow){
  for (j in 1:numcol){
    
    #Progress bar actualization
    Sys.sleep(0.00001)
    setWinProgressBar(pb, i, title=paste("Initial Aspect and slope progress (", round(i/(numrow)*100, 0), "% done) please wait"))
    
    
    #Previous estimations, comon to aspect and slope calculations
    
    Z[2:(numcol+1),2:(numrow+1),1]<-DEM[,,1]
    
    DZDX[j,i,1]<-((Z[j,i+2,1]+2*Z[j+1,i+2,1]+Z[j+2,i+2,1])-(Z[j,i,1]+2*Z[j+1,i,1]+Z[j+2,i,1]))/(8*Xcellsize)
    DZDY[j,i,1]<-((Z[j+2,i,1]+2*Z[j+2,i+1,1]+Z[j+2,i+2,1])-(Z[j,i,1]+2*Z[j,i+1,1]+Z[j,i+2,1]))/(8*Ycellsize)
    
    #Aspect calculation
    aspect[j,i,1]<-57.29578*atan2(DZDY[j,i,1],-DZDX[j,i,1])
    if (aspect[j,i,1]<=90) aspect[j,i,1]<-90-aspect[j,i,1] else{
      if (aspect[j,i,1]>90) aspect[j,i,1]<-360-aspect[j,i,1]+90 }
    
    #Slope calculation
    slope[j,i,1]<-atan(sqrt(DZDY[j,i,1]^2+DZDX[j,i,1]^2)) #Radians
    
  }} 
  close(pb)
    #hillshade calculation
    hillshade[,,1]<-255*((cos(Zenith_rad)*cos(slope[,,1]))+
   (sin(Zenith_rad)*sin(slope[,,1])*cos(Azimuth_rad[1]-aspect[,,1]*pi/180)))

setwd(Results_Path)

######################################################-MODEL-##########################################################

#######################################-Erosion and saltation function definition-#####################################

for(S in 1:length(alist)){

a<-alist[S]
b<-blist[S]
c<-clist[S]

qH<-function(H){a+b*H+c*H^2}                        #This is the erosion function.

d<-dlist[S]
e<-elist[S]
f<-flist[S]

salt<-function(H){d+e*H+f*H^2}                       #This is the saltation function

qs_function<-function(T){                            #This is the function to apply the erosion function to a DEM
  
    
  q[T]*qH(DEM[,,T])                                  
  
}


L_function<-function(T,L0){                          #This is the function to apply the erosion function to a DEM
  
    
  round(L0*(1+salt(DEM[,,T])))
}
#Generation of a progress bar

simnumber<-startSim+S-1

pb <- winProgressBar(title = "Progress bar", min = 0, max = (max(T)-1), width = 300)
##########################################-Saltation algorithm-######################################################

for(t in 1:(max(T)-1)){
    
      
      
      #Progress bar actualization
      Sys.sleep(0.00001)
      setWinProgressBar(pb, t, title=paste("Simulation number ",simnumber," (", round(t/(max(T)-1)*100, 0), "% done) please wait"))
                                            
  
      #Eroded volume estimation
      
      qs[,,t]<-qs_function(t) 
      
      #Saltation length estimation
      
      lsx[,,t]<-L_function(T=t,L0=l0x[t])
      lsy[,,t]<-L_function(T=t,L0=l0y[t])
      
      #Cell exposition, deposition and erosion
      
      Cond_Hillshade[,,t]<-ifelse(hillshade[,,t]>=0,1,0) #1 if the cell is to be eroded
      
      if (Trapping==0){DEM[,,t]<-DEM[,,t]}
      if (Trapping==1){DEM[,,t]<-DEM[,,t]+Cond_Hillshade[,,t]*ifelse(DEM[,,t]>0,1,0)*q[t]}
      if (Trapping==2){DEM[,,t]<-DEM[,,t]+(1-Cond_Hillshade[,,t])*ifelse(DEM[,,t]>0,1,0)*q[t]}
      if (Trapping==3){DEM[,,t]<-DEM[,,t]+ifelse(DEM[,,t]>0,1,0)*q[t]}
      if (Trapping==4){DEM[1,,t]<-DEM[1,,t]+q[t];DEM[,2:numrow,t]<-DEM[,2:numrow,t]+q[t]}
      
      
      
      
      DEMz1[,,t]<-ifelse(DEM[,,t]>qs[,,t],DEM[,,t]-qs[,,t],0)*Cond_Hillshade[,,t]
      
      #Eroded cell deposition algorithm
      
      for (i in 1:numrow){
        for (j in 1:numcol){
      
      if 
      ((j<=(numcol-lsy[j,i,t]))&(i<=(numrow-lsx[j,i,t])))
        
      {
        DZ2[j+lsy[j,i,t],i+lsx[j,i,t],t]<-DEM[j,i,t]*Cond_Hillshade[j,i,t]-DEMz1[j,i,t]
      } 
      
        }}
      
      #DEM state after saltation
      
      DEM[,,t+1]<-DEMz1[,,t]+DZ2[,,t]+DEM[,,t]*(1-Cond_Hillshade[,,t])
      
    
  
##########################################-Avalanche algorithm-###################################################
  
  
  for (i in 1:numrow){
    for (j in 1:numcol){
      
      
            if (((j-1)&(i-1)>0)&((j+1)<numcol)&((i+1)<numrow)){  
            CentralCell<-DEM[j,i,t+1]
            
            Neighbors<-c(                                               #Here we look for the neighbors
              
              DEM[j-1,i-1,t+1],  DEM[j,i-1,t+1],   DEM[j+1,i-1,t+1],
              DEM[j-1,i,t+1],    CentralCell,     DEM[j+1,i,t+1],
              DEM[j-1,i+1,t+1],  DEM[j,i+1,t+1],   DEM[j+1,i+1,t+1]
            )
            
            localslope<-(CentralCell-Neighbors)/Xcellsize               #We estimate the local slope
            
            #iteration<-0
            
            while ((max(abs(localslope))>tan(RA))){
              
                            
              for (X in c(1,2,3,4,6,7,8,9)){
              
              if ((abs(localslope[X]))>tan(RA)){
              
              TarjetHigh<-Neighbors[X]
              
              if (CentralCell>TarjetHigh){
              
              qa<-0.5*(CentralCell-TarjetHigh-Xcellsize*tan(RA))        #we stimate here the volume to be desplaced
              
              DEM[j,i,t+1]<-CentralCell-qa
              
              DEM[j+jj[X],i+ii[X],t+1]<-DEM[j+jj[X],i+ii[X],t+1]+qa
              
              } else{
                
                qa<-0.5*(TarjetHigh-CentralCell-Xcellsize*tan(RA)) 
                
                DEM[j,i,t+1]<-CentralCell+qa
                
                DEM[j+jj[X],i+ii[X],t+1]<-DEM[j+jj[X],i+ii[X],t+1]-qa
              }
              Neighbors<-c(
                
                DEM[j-1,i-1,t+1],  DEM[j,i-1,t+1],   DEM[j+1,i-1,t+1],
                DEM[j-1,i,t+1],    DEM[j,i,t+1],     DEM[j+1,i,t+1],
                DEM[j-1,i+1,t+1],  DEM[j,i+1,t+1],   DEM[j+1,i+1,t+1]
              )
              
              
              localslope<-(DEM[j,i,t]-Neighbors)/Xcellsize
              
              
              #iteration<-iteration+1
              #Sys.sleep(0.1)
              #plot(j,i,main=paste("Avalanching position ",j,i),xlim=c(0,numcol),ylim=c(0,numrow),xlab="Y",ylab="X")
              #text(j+5*Xcellsize,i+5*Xcellsize,paste("Iteration ",iteration))
              
              }} 
              if (max(abs(localslope)>=tan(RA))) {break}
              }}}}
              
              
         
      
############################-Slope, Aspect and hillshade application to each time step-#################################
  
  #See coments in line 186
  for (i in 1:numrow){
        for (j in 1:numcol){
      
      Z[2:(numcol+1),2:(numrow+1),t+1]<-DEM[,,t+1]
      
      DZDX[j,i,t+1]<-((Z[j,i+2,t+1]+2*Z[j+1,i+2,t+1]+Z[j+2,i+2,t+1])-(Z[j,i,t+1]+2*Z[j+1,i,t+1]+Z[j+2,i,t+1]))/(8*Xcellsize)
      DZDY[j,i,t+1]<-((Z[j+2,i,t+1]+2*Z[j+2,i+1,t+1]+Z[j+2,i+2,t+1])-(Z[j,i,t+1]+2*Z[j,i+1,t+1]+Z[j,i+2,t+1]))/(8*Ycellsize)
      
      slope[j,i,t+1]<-atan(sqrt(DZDY[j,i,t+1]^2+DZDX[j,i,t+1]^2)) 
      aspect[j,i,t+1]<-57.29578*atan2(DZDY[j,i,t+1],-DZDX[j,i,t+1])
        if (aspect[j,i,t+1]<=90) aspect[j,i,t+1]<-90-aspect[j,i,t+1] else{
          if (aspect[j,i,t+1]>90) aspect[j,i,t+1]<-360-aspect[j,i,t+1]+90 }
        
        }}
      
        hillshade[,,t+1]<-255*((cos(Zenith_rad)*cos(slope[,,t+1]))+
          (sin(Zenith_rad)*sin(slope[,,t+1])*cos(Azimuth_rad[t+1]-aspect[,,t+1]*pi/180)))

###############################################-Error estimation-############################################################
}   
MASK<-ifelse(DEMFin[,,1]>0,1,0)
ERRORS<-DEMFin[,,1]-DEM[,,steps]              #An array showing the diference between model and surveyed dune final state
MSE<-mean((DEMFin[,,1]-DEM[,,steps]*MASK)^2) #Mean Square Error

###############################################-Simulation save-###########################################################

      

save.image(paste(as.character(simnumber),".RData",sep="") )
write.gif(DEM, paste(as.character(simnumber),".gif",sep=""), col=jet.colors, delay=1,scale="always")
close(pb)
}    

    




#############################################-Results post processing-##############################################

#Results visualization
levelplot(DEM[,,steps],at=(seq(from=0.001,to=max(DEM), by=(max(DEM)/60))),col.regions=jet.colors,xlab="Y",ylab="X")
levelplot(ERRORS, by=(max(DEM)/60),col.regions=jet.colors,xlab="Y",ylab="X")

##Mass balance
Vout<-Vini-sum(DEM[,,steps]*MASK*Xcellsize*Ycellsize) #Volume that leaves the system
Vin<-Vfin-sum(DEM[,,steps]*MASK*Xcellsize*Ycellsize)  #Volume that remains in the system

