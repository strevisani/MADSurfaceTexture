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
        self.label = "MadScan"
        self.description = "Calculates MAD basic indexes"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            # residual raster for input
            displayName = "Input Raster",
            name = "in_raster",
            datatype = "GPRasterLayer",
            parameterType = "Required",
            direction = "Input")
        
        param1 = arcpy.Parameter(
            # output iso
            displayName = "ISO",
            name = "out_iso",
            datatype = "DERasterDataset",
            parameterType = "Required",
            direction = "Output")

        param2 = arcpy.Parameter(
            # output anisodir
            displayName = "ANISODIR",
            name = "out_anisodir",
            datatype = "DERasterDataset",
            parameterType = "Required",
            direction = "Output")

        param3 = arcpy.Parameter(
            # output anisor
            displayName = "ANISOR",
            name = "out_anisor",
            datatype = "DERasterDataset",
            parameterType = "Required",
            direction = "Output")

        param4 = arcpy.Parameter(
            displayName = "Kernel Dir",
            name = "kerneldir",
            datatype = "DEDiskConnection",
            parameterType = "Required",
            direction = "Input")
        
        param5 = arcpy.Parameter(
            # lag value list
            displayName = "Select Lag value",
            name = "in_lagvalue",
            datatype = "GPLong", #maybe not long
            parameterType = "Required",
            direction = "Input")
        
        param5.filter.type = "ValueList"
        param5.filter.list = [2, 4, 6, 8]

        param6 = arcpy.Parameter(
            # radius value 
            displayName = "Insert Radius value",
            name = "in_radius",
            datatype = "GPDouble",
            parameterType = "Required",
            direction = "Input")

        


        params = [param0, param1, param2, param3, param4, param5, param6]
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
        """The source code of the tool."""
        
        from MADfunctionsV1 import *
        def makeKernelList(lagval,kernbase):
	    """Generate kernel list to be used"""
	    kernlist =[]
	    if lagval == 2:
	        kernlist = [kernbase+"\\N2c.txt",kernbase+"\\NE2c.txt",kernbase+"\\E2c.txt",kernbase+"\\SE2c.txt"]
	    elif lagval == 4:
	        kernlist = [kernbase+"\\N4c.txt",kernbase+"\\NE4c.txt",kernbase+"\\E4c.txt",kernbase+"\\SE4c.txt"]
	    elif lagval == 6:
	        kernlist = [kernbase+"\\N6c.txt",kernbase+"\\NE6c.txt",kernbase+"\\E6c.txt",kernbase+"\\SE6c.txt"]
	    elif lagval == 8:
	        kernlist = [kernbase+"\\N8c.txt",kernbase+"\\NE8c.txt",kernbase+"\\E8c.txt",kernbase+"\\SE8c.txt"]
	    return kernlist
        
        inraster = parameters[0].valueAsText
        out_iso = parameters[1].valueAsText
        out_anisodir = parameters[2].valueAsText
        out_anisor = parameters[3].valueAsText
        kernbase = parameters[4].valueAsText
        lagval = parameters[5].value
        radiusval = parameters[6].value

        dirk = makeKernelList(lagval, kernbase)
        iso, anisoDir, anisoR = madscan(inraster,dirk,radiusval)
        iso.save(out_iso)
        anisoDir.save(out_anisodir)
        anisoR.save(out_anisor)

        
        return
