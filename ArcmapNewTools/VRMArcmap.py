###############################Roughness based on vector dispersion of normal vectors#
#need spatial analyst:
from arcpy.sa import *

def circularDispNV(inRaster,window,nodes):
	#tested
	# calculates vector roughness (for some "ruggedness") using spherical dispersion
	# of Normal Vectors.
	# The parameter "nodes" is the number of considered pixels (i.e., vectors) in the kernel
	# e.g. for a circular moving window of radius 3 are 29.
	# for a square kernel is sizeX x sizeY
	#now, I have no time to code the automatic derivation of the number of pixels,
	#in R is simpler!
	slope=Slope(inRaster)
	aspect=Aspect(inRaster)
	#deg->rad
	def rad(degree):
	    angle=(degree * math.pi)/180
	    return angle
	#tested
	x=Sin(rad(slope))*Cos(rad(aspect))
	y=Sin(rad(slope))*Sin(rad(aspect))
	z=Cos(rad(slope))
	X=FocalStatistics(x,window,"SUM","NODATA")
	Y=FocalStatistics(y,window,"SUM","NODATA")
	Z=FocalStatistics(z,window,"SUM","NODATA")
	R=SquareRoot(Square(X)+Square(Y)+Square(Z))
	return (1-(R/nodes))
#
