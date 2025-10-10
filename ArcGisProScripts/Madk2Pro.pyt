from arcpy.sa import *
import arcpy
import os

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "MADk2 lag 1 and 2 pixels"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Madtool]

class Madtool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "MADk2"
        self.description = "Calculates short-range roughness indexes (isotropic roughness, anisotropy direction and strength) using differences of order 2 for user-selectable lag distance and circular window. For details see: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4223135"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            displayName="Input Raster",
            name="in_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
            
        param1 = arcpy.Parameter(
            displayName="Output Workspace",
            name="out_workspace", 
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        param2 = arcpy.Parameter(
            displayName="Search Window Radius (cells)",
            name="window_radius",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        
        # Set default radius value
        param2.value = 3
        
        # Set filter for radius parameter
        param2.filter.type = "Range"
        param2.filter.list = [1, 100]  # Minimum 1 cell, maximum 100 cells

        param3 = arcpy.Parameter(
            displayName="Lag Distance (cells)",
            name="lag_distance",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        # Set lag distance options
        param3.filter.type = "ValueList"
        param3.filter.list = ["1 pixel", "2 pixels"]
        param3.value = "1 pixel"

        params = [param0, param1, param2, param3]
        return params

    def isLicensed(self):
        """This tool can be execute with Spatial Analyst."""
        try:
            if arcpy.CheckExtension("Spatial") == "Available":
                arcpy.CheckOutExtension("Spatial")
                return True
            else:
                return False
        except:
            return False

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        # Validate that radius is a positive integer
        if parameters[2].value is not None:
            if parameters[2].value < 1:
                parameters[2].setErrorMessage("Search window radius must be at least 1 cell")
            elif parameters[2].value > 100:
                parameters[2].setWarningMessage("Large radius values may significantly increase processing time")
        return

    # Core functions moved inside the execute method as local functions
    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        # Core functions as local functions within execute
        def calcDelta(inRaster, kernel):
            """Calculate directional differences on a raster"""
            myNbrWeight = NbrWeight(kernel)
            delta = FocalStatistics(inRaster, myNbrWeight, "SUM", "NODATA")
            return delta

        def calcAllDelta(inRaster, kernels):
            """Calculate all directional differences for given kernels"""
            return [calcDelta(inRaster, X) for X in kernels]

        def medFloat(floatRaster, window):
            """Function for performing median calculation on float raster"""
            # Use higher precision for modern ArcGIS Pro
            med = Float(FocalStatistics(Int(floatRaster * 10000), window, "MEDIAN", "NODATA")) / 10000
            return med

        def medAbs(inDelta, window):
            """Calculate MAD - median absolute directional differences"""
            medabs = medFloat(Abs(inDelta), window)
            return medabs

        def medAbsAll(myDeltas, window):
            """Calculate median absolute differences for all deltas"""
            return [medAbs(X, window) for X in myDeltas]

        def anisoDir(N, NE, E, SE):
            """Calculate direction of maximum continuity"""
            direction = 180 - (57.2957 * 0.5 * ATan2((NE - SE), (E - N)))
            return direction

        def anisoDirL(x):
            """Calculate direction of maximum continuity from list"""
            direction = anisoDir(x[0], x[1], x[2], x[3])
            return direction

        def anisoR(N, NE, E, SE):
            """Standardized resultant length as anisotropy index"""
            meanR = SquareRoot(Square(NE - SE) + Square(E - N)) / (N + NE + E + SE)
            return meanR

        def anisoRL(x):
            """Standardized resultant length from list"""
            meanR = anisoR(x[0], x[1], x[2], x[3])
            return meanR

        def fullPath(myDir, myNames, mySuffix):
            """Build full path for kernels"""
            fullDir = []
            for i in range(len(myNames)):
                fullDir.append(myDir + myNames[i] + mySuffix)
            return fullDir

        def madscan(inRaster, kernels, windowRadius):
            """Calculate isoMad, anisotropy direction and anisotropy R"""
            # Calculate deltas in 4 directions
            deltas = calcAllDelta(inRaster, kernels)
            # Define search window
            radius = windowRadius
            window = NbrCircle(radius, "CELL")
            # Calculate MAD
            myAbs = medAbsAll(deltas, window)
            # Calculate desired indexes
            isoMad = CellStatistics(myAbs, "MEAN", "NODATA")
            anisoMad = anisoDirL(myAbs)
            anisoRMad = anisoRL(myAbs)
            return isoMad, anisoMad, anisoRMad

        # Main execution
        try:
            # Get parameters
            inraster = parameters[0].valueAsText
            out_workspace = parameters[1].valueAsText
            window_radius = parameters[2].value
            lag_distance = parameters[3].valueAsText
            
            # Set workspace
            arcpy.env.workspace = out_workspace
            arcpy.env.overwriteOutput = True
            
            # Create kernel directory in scratch workspace
            kernelDir = arcpy.env.scratchFolder
            if not kernelDir:
                kernelDir = arcpy.env.scratchGDB
            kernelDir = kernelDir + os.sep
            
            messages.addMessage("Starting MAD roughness calculation...")
            messages.addMessage(f"Input raster: {inraster}")
            messages.addMessage(f"Output workspace: {out_workspace}")
            messages.addMessage(f"Search window radius: {window_radius} cells")
            messages.addMessage(f"Lag distance: {lag_distance}")
            
            messages.addMessage("Creating kernel files...")
            
            # Define kernel contents based on lag distance selection
            if lag_distance == "1 pixel":
                # 3x3 kernels for 1-pixel lag
                N1ck2 = ['3 3\n', '0 -1 0\n', '0 2 0\n', '0 -1 0\n']
                NE1ck2 = ['3 3\n', '0 -0.207106781186547 -0.5\n', '-0.207106781186548 1.82842712474619 -0.207106781186547\n', '-0.5 -0.207106781186547 0\n']
                E1ck2 = ['3 3\n', '0 0 0\n', '-1 2 -1\n', '0 0 0\n']
                SE1ck2 = ['3 3\n', '-0.5 -0.207106781186548 0\n', '-0.207106781186547 1.82842712474619 -0.207106781186547\n', '0 -0.207106781186548 -0.5\n']
                kernels = [N1ck2, NE1ck2, E1ck2, SE1ck2]
                kernel_names = ["N1ck2", "NE1ck2", "E1ck2", "SE1ck2"]
                lag_suffix = "c1k2"
            else:
                # 5x5 kernels for 2-pixel lag
                N2ck2 = ['5 5\n', '0 0 -1 0 0\n', '0 0 0 0 0\n', '0 0 2 0 0\n', '0 0 0 0 0\n', '0 0 -1 0 0\n']
                NE2ck2 = ['5 5\n', '0.0000000 0.0000000 0 -0.2426407 -0.1715729\n', '0.0000000 0.0000000 0 -0.3431458 -0.2426407\n', '0.0000000 0.0000000 2 0.0000000 0.0000000\n', '-0.2426407 -0.3431458 0 0.0000000 0.0000000\n', '-0.1715729 -0.2426407 0 0.0000000 0.0000000\n']
                E2ck2 = ['5 5\n', '0 0 0 0 0\n', '0 0 0 0 0\n', '-1 0 2 0 -1\n', '0 0 0 0 0\n', '0 0 0 0 0\n']
                SE2ck2 = ['5 5\n', '-0.1715729 -0.2426407 0 0.0000000 0.0000000\n', '-0.2426407 -0.3431458 0 0.0000000 0.0000000\n', '0.0000000 0.0000000 2 0.0000000 0.0000000\n', '0.0000000 0.0000000 0 -0.3431458 -0.2426407\n', '0.0000000 0.0000000 0 -0.2426407 -0.1715729\n']
                kernels = [N2ck2, NE2ck2, E2ck2, SE2ck2]
                kernel_names = ["N2ck2", "NE2ck2", "E2ck2", "SE2ck2"]
                lag_suffix = "c2k2"
            
            kernel_files = []
            
            # Write kernel files
            for name, content in zip(kernel_names, kernels):
                kernel_path = kernelDir + name + ".txt"
                with open(kernel_path, 'w') as outfile:
                    outfile.writelines(content)
                kernel_files.append(kernel_path)
            
            messages.addMessage("Calculating roughness indexes...")
            
            # Calculate roughness indexes with user-defined radius
            iso, anisoDir, anisoR = madscan(inraster, kernel_files, window_radius)
            
            # Save outputs with lag distance and radius in filename
            iso_output = os.path.join(out_workspace, f"RoughIso{lag_suffix}_r{window_radius}.tif")
            anisoDir_output = os.path.join(out_workspace, f"RoughAnisoDir{lag_suffix}_r{window_radius}.tif") 
            anisoR_output = os.path.join(out_workspace, f"RoughAniso{lag_suffix}_r{window_radius}.tif")
            
            iso.save(iso_output)
            anisoDir.save(anisoDir_output)
            anisoR.save(anisoR_output)
            
            messages.addMessage("Calculation completed successfully!")
            messages.addMessage(f"Results saved to: {out_workspace}")
            messages.addMessage(f" - Isotropic roughness: RoughIso{lag_suffix}_r{window_radius}.tif")
            messages.addMessage(f" - Anisotropy direction: RoughAnisoDir{lag_suffix}_r{window_radius}.tif")
            messages.addMessage(f" - Anisotropy strength: RoughAniso{lag_suffix}_r{window_radius}.tif")
            
            # Clean up kernel files
            for kernel_file in kernel_files:
                try:
                    if os.path.exists(kernel_file):
                        os.remove(kernel_file)
                except:
                    messages.addWarningMessage("Could not delete temporary kernel file: " + kernel_file)
                    
        except Exception as e:
            messages.addErrorMessage("Error executing tool: " + str(e))
            raise
            
        finally:
            # Check in extension
            if arcpy.CheckExtension("Spatial") == "CheckedOut":
                arcpy.CheckInExtension("Spatial")
                
        return
        