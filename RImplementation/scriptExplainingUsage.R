#Here an example of script that shows how to derive basic short-range surface roughness indexes
#Uncomment the following command for cleaning things
#rm(list=ls())
#ls()
# set the directory containing source code with functions (you need to set yours)
setwd("F:/backmad/CodeProfiling/CodeAndDataSharing")
dir()
#library(terra)
#load the functions for geostatistical analysis
source("GeostTextureLibraryShared.R")
#This will upload automatically the "Terra" package and the
#precomputed basic kernels for directional differences of order 1 and 2
#i.e.
#load("basicKernels.RData")
#
#The naming convention for kernels
#of order 1 is "kNc" in which "N" is the lag distance in pixels (i.e.,  K1c->1 pixel, K2c->2 pixels, etc.)
#For differences of order 2 you have the kernels K05ck2->half pixel, K1ck2->1 pixel and k2ck2->2 pixels.
#I suggest to use K1ck2 if you don't need measures of relative roughness.
#Each kernel object is a list of four directional kernels for directions
#NS, NE-SW, E-W, SE-NW


#Let's see an example with the usual geostatistical approach
#using the residual DEM, in this case a residual DTM.
#Approach described in paper:
#Trevisani, S. & Rocca, M., 2015. MAD: Robust image texture analysis for applications in high resolution geomorphometry. 
#Computers and Geosciences, vol. 81, pp. 78-92.

#The residual DTM has been derived by means of a simple two pass moving average approach,
#but you can derive it via other approaches, e.g. gaussian kernels.

#In this example we need to use kernels of order 1 (i.e., directional differences)
#load the residual DTM

#Set directory with data (define your path)
setwd("F:/backmad/CodeProfiling/CodeAndDataSharing/DATA")
dir()
#Load residual DTM
res=rast("residualDem.tif")
res
plot(res)

#if you want to compute short-range roughness indexes for a lag of 2 pixels
#you need to use the kernel k2C
#and set the search window,
#for example a circular one with a radius of 3 pixels
w=KernelCircular(3)
#w
#for a rectangular use KernelRectangular(x,y)

#Calcualate basic short-range roughness indexes
mad2c=Madscan(res,k2c,w)
#the object mad2c contains 3 rasters: isotropic roughness, anisotropy direction and anisotropy strength
#If you need to save on disk for importing in GIS
#writeRaster(mad2c[[1]],"rIso2c.tif",NAflag =-3.40282346639e+038,overwrite=T)
#writeRaster(mad2c[[2]],"rAnisoDir2c.tif",NAflag =-3.40282346639e+038,overwrite=T)
#writeRaster(mad2c[[3]],"rAnisoR2c.tif",NAflag =-3.40282346639e+038,overwrite=T)
#if you need other lag distances just change the kernel, e.g. mad4c=Madscan(res,K4c,w).

#Now lest see the calculation of roughness
#using the approach based on differences of order 2, i.e. differences of differences.
#This approach is applied directly to the DEM/image, it bypass detrending
#The approach is fully described in the paper(currently a preprint):
#Trevisani, S., Teza, G. and Guth, P., A Simplified Geostatistical Approach for Characterizing Key Aspects of Short-Range Roughness. 
#Available at SSRN: https://ssrn.com/abstract=4223135 or http://dx.doi.org/10.2139/ssrn.4223135
#
#Load the DTM from Trentino Province
trento=rast("trentoDEM.tif")
trento
#plot(trento)
#calculate short-range roughness with kernels of order 2 and a lag of 1 pixel
mad1ck2=Madscan(trento,k1ck2,w)
#same as above, 3 rasters: isotropic roughness, anisotropy direction and anisotropy strength
#writeRaster(mad1ck2[[1]],"rIso1ck2.tif",NAflag =-3.40282346639e+038,overwrite=T)
##writeRaster(mad1ck2[[2]],"rAnisoDir1ck2.tif",NAflag =-3.40282346639e+038,overwrite=T)
#writeRaster(mad1ck2[[3]],"rAnisoR1ck2.tif",NAflag =-3.40282346639e+038,overwrite=T)
###

#For understanding how things work or for testing purposes you can
#calculate directional differences with the various kernels, e.g.:
#
####deltasN2c=focal(res, w=data.matrix(k2c[[1]]), na.rm=FALSE,expand=F)
###deltasNE2c=focal(res, w=data.matrix(k2c[[2]]), na.rm=FALSE,expand=F)
##deltasE2c=focal(res, w=data.matrix(k2c[[3]]), na.rm=FALSE,expand=F)
#d#eltasSE2c=focal(res, w=data.matrix(k2c[[4]]), na.rm=FALSE,expand=F)


