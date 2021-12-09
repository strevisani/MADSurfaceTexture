from arcpy.sa import *
# 
class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Madtool]


class Madtool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mad2c"
        self.description = "Calculates short-range roughness using a lag of 2 pixels and a search radius of 3 pixels"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            # raster for input
            displayName = "Input Raster",
            name = "in_raster",
            datatype = "GPRasterLayer",
            parameterType = "Required",
            direction = "Input")


        params = [param0]
        return params

    def isLicensed(self):
        """This tool can be execute with Spatial Analist."""
        try:
            if arcpy.CheckExtension("Spatial") == "Available":
                arcpy.CheckOutExtension("Spatial")
            else:
                raise Exception
        except:
            return False
            
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        
        """The source code of the tool"""
        ###############################################
        """This code is derived (and is a subset) from the first published version
        in June 2014, made available in GITHUB in relation to the paper:
        Trevisani, S. & Rocca, M. 2015, "MAD: Robust image texture analysis for applications in high resolution geomorphometry", Computers and Geosciences, vol. 81, pp. 78-92
	
	https://github.com/strevisani/MADSurfaceTexture
	
	#MIT license
	Begin of license, limited to the following code.
	
	Copyright (c) 2014 and 2021 Sebastiano Trevisani
	Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
	
	1)The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
	2)When using this software, in particular for publications, please cite the related paper: " Trevisani, S. & Rocca, M. 2015, "MAD: Robust image texture analysis for applications in high resolution geomorphometry", Computers and Geosciences, vol. 81, pp. 78-92."
.
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
	
	End of license
	
	
	#Important note: in Arcmap old versions was not possible to calculate
	#the median with focalstastics on a float raster (!!!).
	#So, to do that we had to convert temporary the float raster to an integer
	# using a multiplying factor (10000)see the function "medFloat()".
	#With the current Arcmap version (e.g.,10.8 you would not need this).
	
	
	#These functions make use of spatial analyst extension
	"""
	
	#Surface/image texture related functions an utilities
	
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
	
	
	###Calculation of indexes
	
	#Utility function for calculating the median of a float raster with focalstatistics in ArcGIS 10.2 (
	#which doesn't permit to calculate with focal statistics the median of a float!!!).
	#We are thinking to DTM so using as multiplying factor of 1000 is ok.
	#You should change this maybe with other kind of data.
	#it seems that with newer version it is not necessary (Agust 2021)
	#IN every case now I increase the number of digits for conversion 1000->10000

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
		this functions calculates directly isoMad, anisotropy direction and anisotropy R
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
	#################################################################################
	###End Surface/image texture related functions an utilities###
	
	###The algorithm of the toolbox###
        import os
        from arcpy import env
        #write the kernels
        #in the working directory
	base=os.getcwd()
	base=base.replace("\\","/")
	env.workspace = base
	kernelDir=base+"/"
	N2c=['3 3\n', '0 1 0\n', '0 0 0\n', '0 -1 0']
	NE2c=['3 3\n', '0 0.207106781186547 0.5\n', '-0.207106781186548 0 0.207106781186547\n', '-0.5 -0.207106781186547 0\n']
	E2c=['3 3\n', '0 0 0\n', '-1 0 1\n', '0 0 0']
	SE2c=['3 3\n', '-0.5 -0.207106781186548 0\n', '-0.207106781186547 0 0.207106781186547\n', '0 0.207106781186548 0.5\n']
	#the name of the kernels for directional differences lag 2 pixels
        k2c=["N2c","NE2c","E2c","SE2c"]
	outfile=open("N2c.txt",'w')
	outfile.writelines(N2c)
	outfile.close()
	outfile=open("NE2c.txt",'w')
	outfile.writelines(NE2c)
	outfile.close()
	outfile=open("E2c.txt",'w')
	outfile.writelines(E2c)
	outfile.close()
	outfile=open("SE2c.txt",'w')
	outfile.writelines(SE2c)
	outfile.close()
	dirK2=fullPath(kernelDir,k2c,".txt")
        toDetrend = parameters[0].valueAsText
        #define moving window
        radius=3
	window=NbrCircle(radius,"CELL")
	#calculate trend
	trend=FocalStatistics(toDetrend,window,"MEAN","NODATA")
	trend.save("trend.tif")
	#calculate residual
	inraster=toDetrend-trend
	inraster.save("res.tif")
	#calculate basic roughness indices
        iso, anisoDir, anisoR = madscan(inraster,dirK2,radius)
        iso.save("RoughIsoc2.tif")
        anisoDir.save("RoughAnisoDirc2.tif")
        anisoR.save("RoughAnisoc2.tif")
        for i in dirK2: os.remove(i)
        return
        ###End of algorithm