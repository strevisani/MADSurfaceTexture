#Upgraded on December 2021 with functions for computing Variogram and Madogram
#
#Upgraded on August 2021 from the old version
#also changed the costant for float conversion to 10000
#in medfloat() function.
#####
#Venice, June 15, 2014
#After the (possible) acceptance of the paper, we think to deliver the code,
#including the kernels used for the derivation
#of directional differences, according
#to the open MIT license. In the meanwhile please keep 
#it reserved and limit the use for the review process.

#MIT license to be used (plus citation and disclaimer)
###Date 1 September 2022 Venice###
#
##Begin of license 1 September 2022
#Copyright (c) 2015 and 2022 Sebastiano Trevisani (strevisani@iuav.it)
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
#2)When using this software, in particular for publications, please cite the related papers:
#Trevisani, S. and Teza, G. and Guth, P., A Simplified Geostatistical Approach for Characterizing Key Aspects of Short-Range Roughness. Available at SSRN: https://ssrn.com/abstract=4223135 or http://dx.doi.org/10.2139/ssrn.4223135
#Trevisani S., Rocca M.,MAD: robust image texture analysis for applications
#in high resolution geomorphometry. Computer & Geosciences, 2015 https://doi.org/10.1016/j.cageo.2015.04.003
	
#
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
#

#In this code we define the core functions of MAD based indexes presented
#in the paper "Trevisani S., Rocca M.,MAD: robust image texture analysis for applications
#in high resolution geomorphometry. Submitted to Computer & Geosciences, 2014".
#These functions are written in order to be easy to understand and
#to modify or re implement in other softwares. Where relevant, we suggest possible modifications,
#generalizations and special use of functions. We also inserted some "classical" geostatistical
#functions such as variogram and madogram for comparison purposes (but we use our modified sampling
#approach, i.e. the search window increases with the lag size).
#We avoided to be too much pythonic in programming styles so as to make
#clear the different steps of algorithms. Moreover we opted to code the program
#with "static" kernels files in order to facilitate their use in other softwares,
#i.e. just copy the weights. However it is clear that these kernels can be created
#on the fly in relation to specific needs and also more directions can be calculated.

#Important note: in ArcGIS is not possible to calculate
#the median with focalstastics on a float raster (!!!).
#So, to do that we need to convert temporary the float raster to an integer
# using a multiplying factor (1000)see the function "medFloat()"


#These functions make use of spatial analyst extension
#Anyway, in all those software codes implementing custom kernels 
# these function can be easily implemented
from arcpy.sa import *

###Directional differences
#The basic step for MAD as well as other bivariate spatial continuity indexes (variogram, madogram, etc.)
#is to calculate quickly and accurately differences between points pairs. This is accomplished via ad-hoc
#defined kernels (in the folder kernels) that permit via focal statistic function to perform bilinear interpolation
# and the difference between interpolated points in one step. We mainly defined kernels for different lag distances
# considering four directions (N-S,SW-NE,W-E, NW-SE). Clearly, other kernels can be defined,
#in order to perform the calculation in more directions (if necessary).
#The calculation of directional differences, at least
#for not large kernels is not computationally demanding.

##Core directional differences functions

#This is the core function for performing multiscale directional analysis of surface texture
#Use these directional differences for calculating anisotropy parameters and basic surface texture indexes
def calcDelta(inRaster,kernel):
	"""
	Calculate directional differences on a raster "inRaster" using a kernel file
	"kernel". calcDelta returns a raster with directional differences with the same
	resolution of input raster. The "kernel" should report the full path to kernel file.
	"""
	myNbrWeight=NbrWeight(kernel)
	delta=FocalStatistics(inRaster,myNbrWeight,"SUM","NODATA")
	return delta
#End calcDelta.

#A special function for playing with increments of order-K
def calcDeltak2(inRaster,kernel):
	"""
	This function calculate difference of differences, i.e. increments of order 2, and can be
	easily expanded to higher order k. Applying this function directly on a DTM, without detrending,
	permits to highlight fine scale morphologies.
	"""
        deltak2=calcDelta(calcDelta(inRaster,kernel),kernel)
	return deltak2
#End calcDelta.

#This functions is the most used for deriving directional differences for surface texture
#analysis
#For a given list of kernels for directional differences and a given input raster (likely a residual DTM)
# return a list of rasters with the corresponding directional differences
def calcAllDelta(inRaster,kernels):
	"""
	This function calculates for an input raster "inRaster" all the directional differences
	reported in a list of kernels "kernels", reporting the full path to kernels. 
	"""
	return [calcDelta(inRaster,X) for X in kernels]
#End calcAllDelta.         

##End Core directional differences functions


##Special (Absolute) directional differences functions for lag of 1 pixel and 1.44 pixels.
#these functions are useful for relative roughness calculation and for special applications.
#These are not intended for the calculation of anisotropy parameters.
##Roughness computed with these kernels can be also useful if you need
#to filter out some fine scale artifacts from the DTM still mantaining
#high detail.

#Lag 1 pixel
#This function calculates the absolute value of directional differences for lag of 1 cell.
#Given that for a specific direction we have two opposite directional differences (e.g. in NS direction
# we have the difference between central pixel and the pixel one step at N and the difference between
# the central pixel and the pixel at one step at S) we need two average their absolute value to have a
#symmetric value to associate to the central pixel and to be used for MAD calculation. So the results of this
# kernels can be used in analogy to output of the basic kernels, but consider some smoothing of roughness.
#This kind of kernel can be also modified to compute non-symmetrical surface texture indexes i.e. where Index(h) is
#different from index(-h).
def calcAbsDelta1c(inRaster,kernelDir):
	"""
	Function for calculating the absolute value of directional differences for lag o 1 pixel.
	It take 2 inputs: inRaster (the input raster) and kernelDir (the directory where kernels are stored).
	The output, even if we have already derived the absolute value, can be used with MAD functions
	as the Basic kernels. We don't suggest the use of these kernels for the calculation of anisotropy indexes.
	These are conceived mainly for the calculation of relative roughness or other specific needs.
	"""
	asym1c=["N1cAsym","NE1cAsym","E1cAsym","SE1cAsym","S1cAsym","SW1cAsym","W1cAsym","NW1cAsym"]
	dirAsym=fullPath(kernelDir,asym1c,".txt")
	deltas1c=calcAllDelta(inRaster,dirAsym)
	#number of differences
	nDeltas=len(deltas1c)
	symAbsDelta1c=[(Abs(deltas1c[i])+Abs(deltas1c[i+nDeltas/2]))/2 for i in range(nDeltas/2)]
	return symAbsDelta1c
#End calcAbsDelta1c.
##End lag 1 pixel


#Lag 1.4142 Pixel
#These set of functions are used for calculating directional
#differences using the shortest lag along diagonals. This function is a variant
#for calculating directional differences and can be useful in some circumstances.
#The first step is to use the kernels:
#myKernels=((WeightDir+"NE1.txt"),(myWeightDir+"SE1.txt"))
#these are 2x2 kernels, that calculate differences along diagonals and the output value is stored on the
#NW corner, but geometrically these differences are representative of the center on the kernel 2x2.
#These differences are then manipulated to derive cell centered absolute directional differences.
#Using a kernel 2x2 we can use simple trigonometry to derive directional differences in any direction
#starting from the 2 calculated differences along diagonals. In this case we consider simply NS and EW directions
#to be compatible with the output of basic kernels (the value 0.7071 come from the sin or cos of 45):
#It is easy to generalize this function for calculating directional differences in any direction.
#This function can be useful for special needs and eventually for calculation of relative roughness.
#Also in this case some smoothing of roughness is expected.
def deltaDiag(inRaster,kernels):
	"""
	Utility function.
	Calculates directional differences in four directions (N-S, NE-SW, E-W, SE-NW) using
	2x2 kernels calculating differences along diagonals. The value is stored on the NW pixel.
	The kernels to be used are "NE1Diag.txt" and "SE1Diag.txt"
	"""
        myNbrWeight = NbrWeight(kernels[0])
        deltaR2NE=FocalStatistics(inRaster,myNbrWeight,"SUM","NODATA")
        myNbrWeight = NbrWeight(kernels[1])
        deltaR2SE=FocalStatistics(inRaster,myNbrWeight,"SUM","NODATA")
        deltaR2N=(-deltaR2SE+deltaR2NE)*0.7071
        deltaR2E=(deltaR2SE+deltaR2NE)*0.7071
        return deltaR2N,deltaR2NE,deltaR2E,deltaR2SE
#end deltaDiag 


def centerDiag(inDelta,shifted):
	"""
	Utility function:
	We use the mean of the Abs of four deltas to re center on the cell
	using kernel "shifted", see kernel "shifted3x3.txt".
	"""
        myNbrWeight = NbrWeight(shifted)
	delta=FocalStatistics(Abs(inDelta),myNbrWeight,"SUM","NODATA")/4
	return delta
#End absDiag.


def absDiagCentered(inRaster,kernelDir):
	"""
	Function for calculating the absolute value of directional differences for lag o 1.4142 pixel.
	It take 2 inputs: inRaster (the input raster) and kernelDir (the directory where kernels are stored).
	The output, even if we have already derived the absolute value, can be used with MAD functions
	as the Basic kernels. We don't suggest the use of these kernels for the calculation of anisotropy indexes.
	These are conceived for the calculation of surface roughness, relative roughness and in those situations
	where we need to calculate directional differences quickly (and approximately) in any possible direction
	(in this case this function and the function deltaDiag should be modified to take as input also the desidered
	angle direction(s)). 
	"""
	kernels=((kernelDir+"NE1Diag.txt"),(kernelDir+"SE1Diag.txt"))
        shifted=kernelDir+"shifted3X3.txt"
	return [centerDiag(X,shifted) for X in deltaDiag(inRaster,kernels)]
#End absDiagCentered.
#End Lag 1.4142 Pixel

##End Special (Absolute) directional differences functions for lag of 1 pixel and 1.4142 pixels.

###End Directional differences


###Calculation of indexes

#Utility function for calculating the median of a float raster with focalstatistics in ArcGIS 10.2 (
#which doesn't permit to calculate with focal statistics the median of a float!!!).
#We are thinking to DTM so using as multiplying factor of 1000 is ok.
#You should change this maybe with other kind of data.
#August 2021
#with the new versions of arcgis you do not need it
#we change 1000 to 10000 for better accuracy
def medFloat(floatRaster,window):
	"""
	Function for performing median calculation on float raster.This is necessary in ArcGis10.2.
	The input parameters "window" is a search window defined according
	to spatial analyst terminology, e.g.:
	radius=3
	window=NbrCircle(radius,"CELL")
	"""
	med=Float(FocalStatistics(Int(floatRaster*10000),window,"MEDIAN","NODATA"))/10000
	return med
#End medFloat.

#This function calculates MAD from directional differences. Using the median as estimator
#we can also start the calculation with squared differences and then convert to absolute
#values after computing the median (the ordering of values doesn't change!)
def medAbs(inDelta,window):
	"""
	Calculate MAD i.e. the median absolute directional differences 
	from differences "inDelta" and using a search window "window".
	Given the way in which directional differences are calculated and stored
	the effective search area is linked to lag size. In order to use other sampling strategies ad-hoc
	kernels can be defined.
	"""
	medabs=medFloat(Abs(inDelta),window)
	return medabs
#end Medabs.


def medAbsAll(myDeltas,window):
	"""
	For all furnished deltas, for example coming from calcAlldelta(),
        calc all median absolute differences with a given search window "window"
	"""
        return [medAbs(X,window) for X in myDeltas]
#end medAbsList.

## Anisotropy parameters
#The multiplication to 57.2957 is related to radiant conversion
#and for 0.5 because of we double the angles for not considering
#the orientation. See Davis, J.C., 2002. Statistics and Data Analysis in Geology: 
#John Wiley & Sons Inc., New York(here we consider also the modulo of vectors).
#This function could be generalized to consider more directions.
def anisoDir(N,NE,E,SE):
	"""
	AnisoDir calculates the direction of maximum continuity
	using a circular statistics approach and using four directions
	along N, NE, E and SE (in this precise order!).
	"""
	direction=180-(57.2957*0.5*ATan2((NE-SE),(E-N)))
	return direction
#end anisoDir.

#The same of anisoDir but working with a list, so this will be the most used.
def anisoDirL (x):
	"""
	AnisoDirL calculates the direction of maximum continuity
	using a circular statistics approach and using a list of four
	rasters of directional differences stored in the following order of
	directions:N, NE, E and SE.
	See function anisoDir() for details.
	"""
        direction=anisoDir(x[0],x[1],x[2],x[3])
        return direction
#end anisoDirL


#Standardized resultant length. This is used as anisotropy index
def anisoR(N,NE,E,SE):
	"""
	Standardized resultant length. This is used as anisotropy index
	Use four rasters with directional differences: N, NE,E, SE
	"""
        meanR=SquareRoot(Square(NE-SE)+Square(E-N))/(N+NE+E+SE)
        return meanR
#end anisoR


#Standardized resultant length
def anisoRL(x):
	"""
	Standardized resultant length. This is used as anisotropy index
	Use a list of four rasters with directional differences: N, NE,E, SE
	"""
        meanR=anisoR(x[0],x[1],x[2],x[3])
        return meanR
#end anisoRL 

#####Compounds functions##############
#These functions makes with one step all the steps required for basic surface texture indexes calculation

#This function is implemented using circular windows
def madscan(inRaster,kernels,windowRadius):
	"""
	this function calculates directly isoMad, anisotropy direction and anisotropy R
	having 3 inputs: 1) input raster, 2)list of kernels (in the 4 directions), 3) Search window radius
	"""
	#calculate deltas in 4 directions
	deltas=calcAllDelta(inRaster,kernels)
	#define search window
	radius=windowRadius
 	window=NbrCircle(radius,"CELL")
 	#calculate MAD
	myAbs=medAbsAll(deltas,window)
	#calculate desired indexes (these are the outputs of the function)
	isoMad=CellStatistics(myAbs,"MEAN","NODATA")
	anisoMad=anisoDirL(myAbs)
 	anisoRMad=anisoRL(myAbs)
 	return isoMad,anisoMad,anisoRMad
###End madscan


#This function calculate isotropic MAD for a lag of 1 pixel
def madscan1c(inRaster,kernelDir,windowRadius):
	deltasAbs1c=calcAbsDelta1c(inRaster,kernelDir)
	radius=windowRadius
	window=NbrCircle(radius,"CELL")
	MAD1c=medAbsAll(deltasAbs1c,window)
	MAD1cIso=CellStatistics(MAD1c,"MEAN","NODATA")
	return MAD1cIso
#End mad1c


#if you need to use spatial analyst aggregate tools, to get
#a "low resolution" (i.e. non overlapping moving windows) result use these functions
def medAbsAllAg(myDeltas,cellFactor):
	"""
	For all furnished deltas, for example coming from calcAlldelta(),
        calc all median absolute differences and aggregate according to a "cellfactor"
	"""
        return [Aggregate(Abs(X),cellFactor, "MEDIAN", "TRUNCATE", "NODATA") for X in myDeltas]
#end medAbsAllAg.

#madscan with use aggregate spatial analyst tool. For a "low resolution" result (i.e. non overlapping
#moving windows)
def madscanAg(inRaster,kernels,cellFactor):
	"""
	this functions calculates directly isoMad, anisotropy direction and anisotropy R
	having 3 inputs: 1) input raster, 2)list of kernels (in the 4 directions), 3) Search window radius
	"""
	#calculate deltas in 4 directions
	deltas=calcAllDelta(inRaster,kernels)
 	#calculate MAD
	myAbs=medAbsAllAg(deltas,cellFactor)
	#calculate desired indexes (these are the outputs of the function)
	isoMad=CellStatistics(myAbs,"MEAN","NODATA")
	anisoMad=anisoDirL(myAbs)
 	anisoRMad=anisoRL(myAbs)
 	return isoMad,anisoMad,anisoRMad
 ###End madscanAg.





### Utilities

## Classical geostatistical estimators
#This function is just used for comparing MAD based indexes with variogram based,
#using the improved sampling approach. To use a classical sampling approach (i.e. fixed
#search window for all lags define custom kernels). 
def gamma(inDelta,window):
	"""
	Just an example of function for computing variogram for specific lags e.g. gamma
	Take care that this function uses the improved sampling approach, i.e.
	search window increasing with lag size.
	"""
	gamma=(FocalStatistics(Square(inDelta),window,"MEAN","NODATA"))/2
	return gamma
#end gamma


def gammaAll(myDeltas,window):
	"""
	For all furnished deltas, for example coming from calcAlldelta(), calculate variogram
	"""
        return [gamma(X,window) for X in myDeltas]
#end gammaAll

#This function is implemented using circular windows
def gammascan(inRaster,kernels,windowRadius):
	"""
	December 2021
	this function calculates directly isotropic variogram, anisotropy direction and anisotropy R
	having 3 inputs: 1) input raster, 2)list of kernels (in the 4 directions), 3) Search window radius
	"""
	#calculate deltas in 4 directions
	deltas=calcAllDelta(inRaster,kernels)
	#define search window
	radius=windowRadius
 	window=NbrCircle(radius,"CELL")
 	#calculate MAD
	myAbs=gammaAll(deltas,window)
	#calculate desired indexes (these are the outputs of the function)
	isoGamma=CellStatistics(myAbs,"MEAN","NODATA")
	anisoGamma=anisoDirL(myAbs)
 	anisoRGamma=anisoRL(myAbs)
 	return isoGamma,anisoGamma,anisoRGamma
###End gammascan

def madogram(inDelta,window):
	"""
	Just an example of function for computing madogram for specific lags.
	Take care that this function uses the improved sampling approach, i.e.
	search window increasing with lag size.
	"""
	mado=(FocalStatistics(Abs(inDelta),window,"MEAN","NODATA"))/2
	return mado
#end madogram

def madogramAll(myDeltas,window):
	"""
	For all furnished deltas, for example coming from calcAlldelta(), calculate madogram
	"""
        return [madogram(X,window) for X in myDeltas]
#end madogramAll

#This function is implemented using circular windows
def madogramscan(inRaster,kernels,windowRadius):
	"""
	December 2021
	this function calculates directly isotropic madpgram, anisotropy direction and anisotropy R
	having 3 inputs: 1) input raster, 2)list of kernels (in the 4 directions), 3) Search window radius
	"""
	#calculate deltas in 4 directions
	deltas=calcAllDelta(inRaster,kernels)
	#define search window
	radius=windowRadius
 	window=NbrCircle(radius,"CELL")
 	#calculate MAD
	myAbs=madogramAll(deltas,window)
	#calculate desired indexes (these are the outputs of the function)
	isoMadogram=CellStatistics(myAbs,"MEAN","NODATA")
	anisoMadogram=anisoDirL(myAbs)
 	anisoRMadogram=anisoRL(myAbs)
 	return isoMadogram,anisoMadogram,anisoRMadogram
###End madogramscan

##End classical geostatistical estimators

#this standardize a raster
def standardize(x):
	return (x-x.minimum)/(x.maximum-x.minimum)
#end standardize

#build the full path for kernels
##an helper function
def fullPath(myDir,myNames,mySuffix):
	fullDir=[]
	for i in range(len(myNames)):fullDir.append(myDir+myNames[i]+mySuffix)
	return fullDir
#End fullPath
































	
