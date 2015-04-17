##Venice, June 15, 2014
#This example script illustrates basic operations
#that can be performed with MAD functions.
#This is intended to be used interactively from
#the phyton console of arcgis
#You can, setting from the script the workspace, also use
#this as a stand alone script.
#Put kernel directory and mad function definitions
#in the same folder of arcgis project.

###Import required modules
import os
import arcpy
from arcpy import env
from arcpy.sa import *
import os
import re
#Import MAD functions
from MADfunctionsV1 import *
#Remember to paste the file with functions on the folder before STARTING ARCGIS.
#I suggest to save scracth rasters in a simple folder not a geodatabase
#it seems faster and safer.
##############################################################
# set overwrite option
#env.overwriteOutput = True
# Check out ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
########################################################
#Take care of international settings use . as decimal separator!
#define workspace
base=os.getcwd()
base=base.replace("\\","/")
env.workspace = base
#to avoid the problem of slashes
#you could use os.path.join function

#if you need to run the script from console you need to define
#the correct path to the workspace and kernels


#define kernels directory (here we suppose in the
#same folder of arcgis project
kernelDir=base+"/kernels/"

###End prepare things.

#I define my output directory, here you should
#define yor correct path and name
outputDir=base+"/finalizing/out/"

#the input residual DTM (or whatever you need to analyze)
inRaster=Raster("sample")
#Now we start with a simple calculation of basic MAD based indexes for a lag of 2 pixel
#We define the list of kernels used for the calculation of directional differences
#the list of available basic kernels
#k2c=["N2c","NE2c","E2c","SE2c"]
#k4c=["N4c","NE4c","E4c","SE4c"]
#k6c=["N6c","NE6c","E6c","SE6c"]
#k8c=["N8c","NE8c","E8c","SE8c"]
k2c=["N2c","NE2c","E2c","SE2c"]
dirK2c=fullPath(kernelDir,k2c,".txt")
#Calculate all deltas. It generate a list of 4 rasters containing
#directional differences along N,NE,E and SE
deltas2c=calcAllDelta(inRaster,dirK2c)
#if you need for validation purposes or for other reason, you can save
#these rasters and load them
deltas2c[0].save(outputDir+"N2c")
deltas2c[1].save(outputDir+"NE2c")
deltas2c[2].save(outputDir+"E2c")
deltas2c[3].save(outputDir+"SE2c")
#
#now let's see the calculation of Mad based indexes.
#You could use directly "madscan()" function, but here we want to show
#all the steps.
#Calculate directional mad for lag 2 for the 4 directions.
#Define a circular window for calculating the median absolute differences
radius=3
window=NbrCircle(radius,"CELL")
#calculate MAD(h) for all deltas
MAD2c=medAbsAll(deltas2c,window)
#MAD2c is  list of raster containing MAD calculated in the 4 directions
#with a lag=2 pixels
#if we need to save and load the results
MAD2c[0].save(outputDir+"MAD2cN")
MAD2c[1].save(outputDir+"MAD2cNE")
MAD2c[2].save(outputDir+"MAD2cE")
MAD2c[3].save(outputDir+"MAD2cSE")
#
#This procedure was followed just to show intermediate results.
#If you need to jump directly to basic surface texture indexes
#you can work directly on the list of rasters.
#Consider that MAD have been calculated with a circular moving window with
#a radius of 3 pixels assuring 29 samples on the neighborhood.
#Given the way in which directional differences are calculated the real area
#sampled increases anisotropically in accordance to lag considered for 
#directional differences calculation. The maximum radius for a directional MAD
#is equal to window radius + lag/2. Take care that other sampling approaches are
#possible.

#Basic surface texture indexes on lag 2 pixels
#Omnidirectional MAD
MAD2cIso=CellStatistics(MAD2c,"MEAN","NODATA")
#Direction of maximum continuity
MAD2cAnD=anisoDirL(MAD2c)
#A measure of anisotropy
MAD2cAnR=anisoRL(MAD2c)
#These are basic usages. Manipulating directional differences
#as well as directional MAD values other indexes can be defined.
#From very simple e.g. MAD2cMin=CellStatistics(MAD2c,"MINIMUM","NODATA")
#to more complex like MAD computed using differences orientated in
#the direction of flow.
#The same procedure can be repeated for other lags.


#Compound functions: directly to indexes
#The 3 last indexes can also computed directly from the raster
#with the function madscan()

inRaster=Raster("sample")
#We already calculated MAD indexes for lag 2 pixels
#If you deleted them uncomment these lines
radius=3
#k2c=["N2c","NE2c","E2c","SE2c"]
#dirK2c=fullPath(kernelDir,k2c,".txt")
#MAD2cIso, MAD2cAnD, MAD2cAnR=madscan(inRaster,dirK2c,radius)
# try with lag= 4 picxels
k4c=["N4c","NE4c","E4c","SE4c"]
dirK4c=fullPath(kernelDir,k4c,".txt")
MAD4cIso, MAD4cAnD, MAD4cAnR=madscan(inRaster,dirK4c,radius)
# try with lag= 6 picxels
k6c=["N6c","NE6c","E6c","SE6c"]
dirK6c=fullPath(kernelDir,k6c,".txt")
MAD6cIso, MAD6cAnD, MAD6cAnR=madscan(inRaster,dirK6c,radius)
# try with lag= 8 picxels
k8c=["N8c","NE8c","E8c","SE8c"]
dirK8c=fullPath(kernelDir,k8c,".txt")
MAD8cIso, MAD8cAnD, MAD8cAnR=madscan(inRaster,dirK8c,radius)
#The distance up to which calculate MAD is dependent
#on the data, the way in which you calculated residual DTM, and the
#target of the study....

#In order to calculate relative roughness as well as for a fine-scale
#calculation of MAD roughness is useful to consider also  a lag of 1 pixels
#
#Roughness computed with madscan1c is also useful if you need
#to filter out some fine scale artifacts from the DTM
radius=3
MAD1cIso=madscan1c(inRaster,kernelDir,radius)
#we can calculate relative roughness
#we define a small number for avoiding division by 0
small=0.001
relRoug=MAD1cIso/(MAD2cIso+small)
#
#We can also compare the results with other surface texture indexes
#for example roughness computed as STD of residual DTM
#we compute with two radius ....just to show that
#defects are visible in both outputs
radius=2
window=NbrCircle(radius,"CELL")
stdR2=FocalStatistics(inRaster,window,"STD","NODATA")
radius=3
window=NbrCircle(radius,"CELL")
stdR3=FocalStatistics(inRaster,window,"STD","NODATA")
#We can also compute variogram and madogram
#Here we consider omnidirectional indexes for lag = 2 cells.
#Calculate differences for lag 2 pixels
k2c=["N2c","NE2c","E2c","SE2c"]
dirK2c=fullPath(kernelDir,k2c,".txt")
deltas2c=calcAllDelta(inRaster,dirK2c)
#we calculate variogram
radius=3
window=NbrCircle(radius,"CELL")
#directional variograms
gamma2c=gammaAll(deltas2c,window)
#isotropic variogram
gam2cIso=CellStatistics(gamma2c,"MEAN","NODATA")
#for comparison purpose we can standardize the square root
#of variogram or at least consider its square root, given
#the strong positive skeweness of variogram
#stGam2c=standardize(SquareRoot(gam2cIso))
#or
#SqRtGam2c=SquareRoot(gam2cIso)
#we can also play with madogram
madogram2c=madogramAll(deltas2c,window)
isoMadogram2c=CellStatistics(madogram2c,"MEAN","NODATA")
###########################################################
