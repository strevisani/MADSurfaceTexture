import arcpy
import numpy as np

class RRIFunction(object):
    def __init__(self):
        self.label = "RRI Function"
        self.description = "Calculate Radial Roughness Index, RRI, (5x5 window, second order differences)"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameters
        param0 = arcpy.Parameter(
            name="in_dem",
            displayName="Input DEM Raster",
            direction="Input",
            datatype="GPRasterLayer",
            parameterType="Required"
        )
        
        param1 = arcpy.Parameter(
            name="out_raster",
            displayName="Output Raster",
            direction="Output",
            datatype="DERasterDataset",
            parameterType="Required"
        )
        
        return [param0, param1]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        try:
            # Get parameters
            in_dem = parameters[0].valueAsText
            out_raster = parameters[1].valueAsText
            
            arcpy.AddMessage("Calculating RRI...")
            
            # Convert to numpy array
            dem_array = arcpy.RasterToNumPyArray(in_dem)
            
            # Get raster properties
            desc = arcpy.Describe(in_dem)
            lower_left = desc.extent.lowerLeft
            cell_size = desc.meanCellWidth
            spatial_reference = desc.spatialReference
            
            # Calculate RRI
            rri_array = self.calculate_rri_array(dem_array)
            
            # Convert back to raster
            rri_raster = arcpy.NumPyArrayToRaster(
                rri_array,
                lower_left,
                cell_size,
                cell_size,
                0
            )
            
            # Set spatial reference
            arcpy.management.DefineProjection(rri_raster, spatial_reference)
            
            # Save output
            rri_raster.save(out_raster)
            
            arcpy.AddMessage("RRI calculation completed successfully!")
            
        except Exception as e:
            arcpy.AddError(f"Error: {str(e)}")
            raise

    def calculate_rri_array(self, dem_array):
        """Calculate RRI for entire array"""
        rows, cols = dem_array.shape
        rri_array = np.zeros((rows, cols), dtype=np.float32)
        
        # Add NoData handling
        dem_array = np.where(dem_array == -9999, np.nan, dem_array)
        
        # Process each cell with 5x5 window (skip 2-pixel border)
        for row in range(2, rows-2):
            for col in range(2, cols-2):
                # Extract 5x5 window
                window = dem_array[row-2:row+3, col-2:col+3]
                
                # Skip if any NaN values in window
                if np.any(np.isnan(window)):
                    rri_array[row, col] = np.nan
                    continue
                    
                # Flatten and calculate RRI
                flat_window = window.flatten()
                rri_array[row, col] = self.calculate_rri(flat_window)
        
        return rri_array

    def calculate_rri(self, x):
        """Calculate RRI for a single 5x5 window (25 elements)"""
        # Note: Indexing adjusted for row-major order (ArcGIS Pro)
        # x[12] is the center cell in 5x5 window
        rri = (
            # External differences
            abs(-0.5*x[0]-0.5*x[12]-0.207106781186547*x[1]-0.207106781186547*x[5]-0.207106781186547*x[7]-0.207106781186547*x[11]+1.82842712474619*x[6]) +
            abs(-x[10]+2*x[11]-x[12]) +
            abs(-0.5*x[20]-0.5*x[12]-0.207106781186547*x[15]-0.207106781186547*x[21]-0.207106781186547*x[11]-0.207106781186547*x[17]+1.82842712474619*x[16]) +
            abs(-x[22]+2*x[17]-x[12]) +
            abs(-0.5*x[24]-0.5*x[12]-0.207106781186547*x[19]-0.207106781186547*x[23]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[18]) +
            abs(-x[14]+2*x[13]-x[12]) +
            abs(-0.5*x[4]-0.5*x[12]-0.207106781186547*x[3]-0.207106781186547*x[9]-0.207106781186547*x[7]-0.207106781186547*x[13]+1.82842712474619*x[8]) +
            abs(-x[2]+2*x[7]-x[12]) +
            # Internal directional differences
            abs(-0.5*x[6]-0.5*x[18]-0.207106781186547*x[7]-0.207106781186547*x[11]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[12]) +
            abs(-x[11]+2*x[12]-x[13]) +
            abs(-0.5*x[16]-0.5*x[8]-0.207106781186547*x[7]-0.207106781186547*x[11]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[12]) +
            abs(-x[7]+2*x[12]-x[17])
        ) / 12
        
        return rri

class Toolbox(object):
    def __init__(self):
        self.label = "RRI Toolbox"
        self.alias = "RRI"
        self.description = "Calculate Radial Roughness Index (RRI)"
        
        # List of tool classes associated with this toolbox
        self.tools = [RRIFunction]
        