###Date 1 September 2022 Venice###
#
##Begin of license
#Copyright (c) 2022 Sebastiano Trevisani (strevisani@iuav.it)
#MIT type license
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#1)The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#2)When using this software, in particular for publications, please cite the related papers (if then accepted refer to that one)
#Trevisani, S. and Teza, G. and Guth, P., A Simplified Geostatistical Approach for Characterizing Key Aspects of Short-Range Roughness. Available at SSRN: https://ssrn.com/abstract=4223135 or http://dx.doi.org/10.2139/ssrn.4223135
#Trevisani S., Rocca M.,MAD: robust image texture analysis for applications
#in high resolution geomorphometry. Computer & Geosciences, 2015 https://doi.org/10.1016/j.cageo.2015.04.003
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.
##End of License
#
##Main functions implementing geostatistical-based
#surface/image texture indexes,using Terra package. 
#With these functions you can
#compute classical geostatistical indexes
#(e.g., variogram and madogram) as well as the robust
#version MAD based on the median of absolute directional differences.
#
#You will find also some functions for computing roughness based
#on dispersion of normal vectors to surface (see paper for details).
#
#In this public version of the code, I did not included the functions for
#creating the kernels for computing the
#directional differences on the fly for any direction and lag.
#I provide a basic set of pre-computed kernels
#for calculating directional differences
#for basic lags (1,2,4,6,and 8 pixels) in four directions. I also provide the kernels for
#computing directional differences of order 2, permitting to calculate the new MADk2 based metrics, 
#which do not require detrending and have some improvement in presence of sharp  morphological transitions.
#However, consider that in the development of ad-hoc kernels,
#not necessarily limited to bivariate indexes, there is a lot
#of potential for detecting interesting aspects of surface/image texture.
#Other potential relies in the various approaches for the decomposition of 
#trend and residuals, and in multiscale smoothing approaches.
#Finally, when you need to study long-range distances you can resample
#the original DEM/image and the resampling should be function of the lag distance and/or the wavelengths you need to highlight.
#
#These functions are developed as a starting point, for promoting
#creativity in surface/image texture analysis, and eventually implementing
#these in other software environments.
#

#load the required library Terra
library(terra)

###Utilities###
#conversion from geographical directions to mathematical ones
matDeg<-function(alpha){
  eta=450-alpha
  if(eta>=360) eta=eta-360 else eta=eta
}

#conversion to radians
rad<-function (degree) {
  (degree * pi)/180
}

#Load the kernels.
#I use square kernels, however in N-S and W-E directions
#you could use a vector (i.e., 1D kernel)
#The kernels are furnished with the source code.
load("basicKernels.RData")
###End Utilities###
   
######Surface texture functions######

KernelCircular=function(radius){
  #This function builds a simple circular kernel
  #with a given radius expressed in pixels
  size=radius*2+1
  center=c((size+1)/2,(size+1)/2)
  kernel=matrix(nrow=size,ncol=size,NA)
  for (i in 1:size) {
    for (j in 1:size){
      if (sqrt((j-center[1])**2+(i-center[1])**2)<=radius){kernel[i,j]=1}
    }
  }
  return(kernel)
}

KernelRectangular=function(lenx,leny){
  #This function builds a simple rectangular kernel
  #with given sizes expressed in pixels
  kernel=matrix(nrow=leny,ncol=lenx,1)
  return(kernel)
}

CalcMedians=function(deltas,w){
  #compute the medians of directional
  #absolute differences from a list
  #and returns a list of rasters.
  #Deltas-> list with rasters containing directional differences (of any order...)
  #w-> the search window (kernel, e.g., w=KernelCircular(3)).
  #For madogram you need to use the mean instead of the median in the
  # focal function, then divide by two (if you need it). If you need the variogram you should consider the
  #square of directional differences, calculate the mean, and divide by 2. Below you will find the functions for conventional geostatistical estimators
  
  medians=list()
  nlayers=nlyr(deltas)
  #Instead of a for cycle you could
  #use sapp...but generally we have
  #very few iterations
  for (i in 1:nlayers){
    medians[[i]]=focal(abs(deltas[[i]]),w,median,na.rm=FALSE,expand=F)
  }
  return(rast(medians))
}

anisoDir=function(N,NE,E,SE){
  #AnisoDir calculates the direction of maximum continuity
  #using a circular statistics approach and using four directions
  #along N, NE, E and SE (in this precise order!).
  #Returns a raster.
  #So the input is MAD (or other spatial variability indexes) 
  #calculated in the four directions
  180-(57.2957*0.5*atan2((NE-SE),(E-N)))
  #If you need more directions you should define a new function
  #taking as argument direction and modulus, quite easy.
}

anisoDirL=function(x){
  #AnisoDirL calculates the direction of maximum continuity
  #using a circular statistics approach and using a list of four
  #rasters of directional differences stored in the following order of
  #directions:N, NE, E and SE.
  ##Returns a raster.
  #See function anisoDir() for details.
  #anisoDir(x[[1]],x[[2]],x[[3]],x[[4]])
  lapp(x, anisoDir)
}

anisoR=function(N,NE,E,SE){
  #Standardized resultant length. This is used as anisotropy index
  #Use four rasters with directional differences: N, NE,E, SE
  #Returns a raster.
  sqrt((NE-SE)**2+(E-N)**2)/(N+NE+E+SE)
}


anisoRL=function(x){
  #Standardized resultant length. This is used as anisotropy index
  #Use a list of four rasters with directional differences: N, NE,E, SE
  #anisoR(x[[1]],x[[2]],x[[3]],x[[4]])
  #Returns a raster.
  lapp(x, anisoR)
} 

Madscan<-function(inRaster,kernels,w){
  #Calculate MAD basic indexes based on 4 directions in this
  #order N,NE,SE,S.
  #Returns 3 rasters: 1)isotropic roughness; 2) direction of anisotropy;
  #3)index of anisotropy.
  #If you need more directions you need to generalize
  #the functions for anisotropy.
  #With very large files and few space on disk I use a modified version
  #that works a little differently. But I do not insert it in this public library.
  #inRaster->input raster, depending from the kernels
  #it may be a detrended version (i.e., high pass filtered) or directly the DTM/image.
  #kernels->a list of kernels (e.g.,myKernels=list(N2c,NE2c,E2c,SE2c)).
  #w->search window (e.g., w=KernelCircular(3)).
  deltas=list()
  #instead of a for loop you could use lapply
  for (i in 1:length(kernels)){
    deltas[[i]]=focal(inRaster, w=data.matrix(kernels[[i]]), na.rm=FALSE,expand=F)
  }
  #directional MADs
  deltas=rast(deltas)
  dirMad=CalcMedians(deltas,w)
  rm(deltas)
  #isotropic mad
  madIso=app(dirMad,fun=mean)
  #direction of maximum continuity
  anisoDirection=anisoDirL(dirMad)
  #anisotropy computed with circular statistics
  #standardized resultant length
  anisoR=anisoRL(dirMad)
  result=c(madIso,anisoDirection,anisoR)
  result
}

###Less robust geostatistical indexes###

CalcMeans=function(deltas,w,exponent){
  #Compute the means of directional
  #absolute differences elevated at an exponent.
  #Returns a raster.
  #With this you can compute variogram and madogram (but remember that for
  #classical geostatistical indexes you need to divide by 2!)
  #Deltas-> list with rasters containing directional differences (of any order...)
  #w-> the search window (e.g., w=KernelCircular(3))
  #Exponent->the exponent to consider (e.g., 2 for variogram and 1 for madogram;
  #but other exponents are possible if you need).
  
  means=list()
  nlayers=nlyr(deltas)
  for (i in 1:nlayers){
    means[[i]]=focal((abs(deltas[[i]]))^exponent,w,mean,na.rm=FALSE,expand=F)
  }
  return(rast(means))
  #
}

Meanscan<-function(inRaster,kernels,w,exponent){
  #this function permits to calculate Variogram (put exponent=2) and madogram (put exponent=1),
  #and divide by two the isotropic roughness to get the usual geostatistical formulation.
  #calculate basic indexes based on 4 directions in this
  #order N,NE,SE,S
  #if you need more directions you need to generalize
  #the functions for anisotropy.
  #Returns 3 rasters: 1)isotropic roughness; 2) direction of anisotropy;
  #3)index of anisotropy.
  #inRaster->input raster, depending from the kernels
  #it may be a detrended version (i.e., high pass filtered) or directly the DTM.
  #kernels->a list of kernels (e.g.,myKernels=list(N2c,NE2c,E2c,SE2c))
  #w->search window (e.g., w=KernelCircular(3))
  deltas=list()
  #instead of a for loop you could use lapply
  for (i in 1:length(kernels)){
    deltas[[i]]=focal(inRaster, w=data.matrix(kernels[[i]]), na.rm=FALSE,expand=F)
  }
  #directional MADs
  deltas=rast(deltas)
  dirMad=CalcMeans(deltas,w,exponent)
  rm(deltas)
  #isotropic variability
  madIso=app(dirMad,fun=mean)
  #direction of maximum continuity
  anisoDirection=anisoDirL(dirMad)
  #anisotropy computed with circular statistics
  #standardized resultant length
  anisoR=anisoRL(dirMad)
  result=c(madIso,anisoDirection,anisoR)
  result
  #
}
###End Less robust geostatistical indexes###


###Other roughness indexes### 

#Here some roughness indexes related to Vector dispersion of normals to surface

circularDispersionGV=function(inraster,window){
  #Circular dispersion of gradient vectors (analogous to circular variance of aspect)
  #using the mean resultant length approach.
  #working with angles in the geographical convention
  #window-> the search window/kernel (e.g., window=KernelCircular(3))
  slope=terrain(inraster,v="slope")
  aspect=terrain(inraster,v="aspect")
  #convert in radians!
  rad<-function (degree) {
    (degree * pi)/180
  }
  x=cos(rad(slope))*cos(rad(aspect))
  y=cos(rad(slope))*sin(rad(aspect))
  z=sin(rad(slope))
  w=window
  X=focal(x, w,fun=sum,expand=F,na.rm=F)
  Y=focal(y, w,fun=sum,expand=F,na.rm=F)
  Z=focal(z, w,fun=sum,expand=F,na.rm=F)
  R=sqrt(X^2+Y^2+Z^2)
  return(1-(R/sum(w,na.rm=T)))
}

###
circularDispersionNV=function(inraster,window){#to test
  #Circular dispersion of normal vectors using the mean resultant length approach.
  #Analogous to VRM in Saga
  #working with angles in the geographical convention
  #window-> the search window/kernel (e.g., window=KernelCircular(3))
  slope=terrain(inraster,v="slope")
  aspect=terrain(inraster,v="aspect")
  #convert in radians!
  #rad<-function (degree) {
  #  (degree * pi)/180
  #}
  #respect to the formulas of Davis book
  # we consider that sin(90-slope)=cos(slope)
  # with 90-slope the angle of the normal vector respect to the horizontal plane
  x=sin(rad(slope))*cos(rad(aspect))
  y=sin(rad(slope))*sin(rad(aspect))
  z=cos(rad(slope))
  w=window
  X=focal(x, w,fun=sum,expand=F,na.rm=F)
  Y=focal(y, w,fun=sum,expand=F,na.rm=F)
  Z=focal(z, w,fun=sum,expand=F,na.rm=F)
  R=sqrt(X^2+Y^2+Z^2)
  return(1-(R/sum(w,na.rm=T)))
}

#vector dispersion using eigenvalues
roory<-function(x){
  #function for using the eigenvalues approach
  x=matrix(x,nrow=3,ncol=3)
  if (anyNA(x)){
    smoothness=NA
  }
  else{
    ei=eigen(x)
    l1=ei$values[1]
    l2=ei$values[2]
    smoothness=log(l1/l2)
  }    
  return(smoothness)
}

circularEigenNV=function(inraster,window){
  #normal vector dispersion using the eigenvalues approach.
  #inraster->input raster (DTM, image, etc.)
  #window-> search window/kernel
  #for example w=KernelCircular(3).
  #it is not very efficient...and with topographic data, if using 2.5D representation,
  #gives analogous results to one based on circular dispersion via resultant length (after considering the log).
  #It could be coded more elegantly and efficiently, I coded it only for testing purposes.
  #With small variations of the code you can also calculate anisotropy, considering other eigenvalues
  #s2 and s3.
  
  slope=terrain(inraster,v="slope")
  aspect=terrain(inraster,v="aspect")
  #convert in radians!
  rad<-function (degree) {
    (degree * pi)/180
  }
 
  #respect to the formulas of Davis book
  # we consider that sin(90-slope)=cos(slope)
  # with 90-slope the angle of the normal vector respect to the horizontal plane
  x=sin(rad(slope))*cos(rad(aspect))
  y=sin(rad(slope))*sin(rad(aspect))
  z=cos(rad(slope))
  w=window
  #Build the matrix of cross products
  x2=x^2
  y2=y^2
  z2=z^2
  xy=x*y
  xz=x*z
  yz=y*z
  X2=focal(x2, w,fun=sum,expand=F,na.rm=F)
  Y2=focal(y2, w,fun=sum,expand=F,na.rm=F)
  Z2=focal(z2, w,fun=sum,expand=F,na.rm=F)
  XY=focal(xy, w,fun=sum,expand=F,na.rm=F)
  XZ=focal(xz, w,fun=sum,expand=F,na.rm=F)
  YZ=focal(yz, w,fun=sum,expand=F,na.rm=F)
  cp<-c(X2,XY,XZ,XY,Y2,YZ,XZ,YZ,Z2)
  cp=values(cp)
  #colnames(cp)<-(c("X2","XY","XZ","XY","Y2","YZ","XZ","YZ","Z2"))
  eigen=apply(cp,1,roory)
  #print ("done")
  rooryz=X2
  return(init(rooryz,eigen))
}


###End other roughness indexes### 


######End Surface texture functions######