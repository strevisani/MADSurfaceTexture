import arcpy
import numpy as np

# Define constants for the RRI calculation
# Use float64 for high precision
K_COEFF = 0.20710678118654746 
K_FACTOR = 1.8284271247461901

class RRIFunction(object):
    def __init__(self):
        self.label = "RRI optimized"
        self.description = "Calculate Radial Roughness Index, RRI, (5x5 window, second order differences),optimized to be fast"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # ... (Parameters remain the same) ...
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
        return arcpy.CheckExtension("Spatial") == "Available"

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        arcpy.env.overwriteOutput = True
        
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise arcpy.ExecuteError("Spatial Analyst license is required for this tool.")

        try:
            in_dem = parameters[0].valueAsText
            out_raster = parameters[1].valueAsText
            
            arcpy.AddMessage("Calculating RRI (Vectorized and Stable)...")
            
            # --- Get Raster Properties ---
            desc = arcpy.Describe(in_dem)
            lower_left = desc.extent.lowerLeft
            cell_size = desc.meanCellWidth
            spatial_reference = desc.spatialReference
            
            # Identify the input NoData value
            input_nodata = desc.noDataValue
            if input_nodata is None:
                # If NoData is not explicitly set in metadata, use the common -9999
                input_nodata = -9999

            # --- Convert to NumPy Array with Float64 Precision ---
            # arcpy.RasterToNumPyArray handles NoData values differently depending on the input.
            # Convert to a standard value (e.g., 0) temporarily to convert to float
            dem_array_temp = arcpy.RasterToNumPyArray(in_dem, nodata_to_value=0)
            
            # Convert the array to high-precision float64 for calculation
            dem_array = dem_array_temp.astype(np.float64) 
            
            # Now, explicitly set the NoData value to NaN
            dem_array = np.where(dem_array_temp == 0, np.nan, dem_array)
            # Handle common implicit NoData values
            dem_array = np.where(dem_array == input_nodata, np.nan, dem_array)
            dem_array = np.where(dem_array == -9999, np.nan, dem_array)


            # --- Calculate RRI (Vectorized) ---
            rri_array = self.calculate_rri_array_vectorized(dem_array)
            
            # --- Convert back to raster and save ---
            # Force output data type to float32, which is common for derived raster products
            output_nodata_value = -9999 
            rri_array = np.where(np.isnan(rri_array), output_nodata_value, rri_array).astype(np.float32)

            rri_raster = arcpy.NumPyArrayToRaster(
                rri_array,
                lower_left,
                cell_size,
                cell_size,
                output_nodata_value
            )
            
            arcpy.management.DefineProjection(rri_raster, spatial_reference)
            rri_raster.save(out_raster)
            
            arcpy.AddMessage("RRI calculation completed successfully!")
            
        except Exception as e:
            arcpy.AddError(f"Error: {str(e)}")
            raise
        finally:
            arcpy.CheckInExtension("Spatial")


    def calculate_rri_array_vectorized(self, dem_array):
        """
        Calculates RRI for the entire array using NumPy vectorization with float64.
        """
        rows, cols = dem_array.shape
        # Initialize the output array with NaNs (float32 is sufficient for output)
        rri_array = np.full((rows, cols), np.nan, dtype=np.float32)

        # Define the inner region where the 5x5 calculation is possible
        rri_window = rri_array[2:rows-2, 2:cols-2]
        
        # Define the center cells (P[2, 2])
        center = dem_array[2:rows-2, 2:cols-2]
        
        abs_diff_terms = []

        # --- Vectorized Term Definitions (5x5 Window P[i,j]) ---

        # The indices map P[i, j] to dem_array[r+i-2, c+j-2]
        # Example: P[0,0] is at (r-2, c-2), P[4,4] is at (r+2, c+2)

        # T1 (x[0], x[12], x[1], x[5], x[7], x[11], x[6])
        T1 = abs(
            -0.5*dem_array[0:rows-4, 0:cols-4] - 0.5*center - 
            K_COEFF*dem_array[0:rows-4, 1:cols-3] - K_COEFF*dem_array[1:rows-3, 0:cols-4] - 
            K_COEFF*dem_array[1:rows-3, 2:cols-2] - K_COEFF*dem_array[2:rows-2, 1:cols-3] + 
            K_FACTOR*dem_array[1:rows-3, 1:cols-3]
        )
        abs_diff_terms.append(T1)

        # T2 (-x[10], +2*x[11], -x[12])
        T2 = abs(-dem_array[2:rows-2, 0:cols-4] + 2*dem_array[2:rows-2, 1:cols-3] - center)
        abs_diff_terms.append(T2)
        
        # T3 (x[20], x[12], x[15], x[21], x[11], x[17], x[16])
        T3 = abs(
            -0.5*dem_array[4:rows, 0:cols-4] - 0.5*center - 
            K_COEFF*dem_array[3:rows-1, 0:cols-4] - K_COEFF*dem_array[4:rows, 1:cols-3] - 
            K_COEFF*dem_array[2:rows-2, 1:cols-3] - K_COEFF*dem_array[3:rows-1, 2:cols-2] + 
            K_FACTOR*dem_array[3:rows-1, 1:cols-3]
        )
        abs_diff_terms.append(T3)
        
        # T4 (-x[22], +2*x[17], -x[12])
        T4 = abs(-dem_array[4:rows, 2:cols-2] + 2*dem_array[3:rows-1, 2:cols-2] - center)
        abs_diff_terms.append(T4)
        
        # T5 (x[24], x[12], x[19], x[23], x[13], x[17], x[18])
        T5 = abs(
            -0.5*dem_array[4:rows, 4:cols] - 0.5*center - 
            K_COEFF*dem_array[3:rows-1, 4:cols] - K_COEFF*dem_array[4:rows, 3:cols-1] - 
            K_COEFF*dem_array[2:rows-2, 3:cols-1] - K_COEFF*dem_array[3:rows-1, 2:cols-2] + 
            K_FACTOR*dem_array[3:rows-1, 3:cols-1]
        )
        abs_diff_terms.append(T5)
        
        # T6 (-x[14], +2*x[13], -x[12])
        T6 = abs(-dem_array[2:rows-2, 4:cols] + 2*dem_array[2:rows-2, 3:cols-1] - center)
        abs_diff_terms.append(T6)
        
        # T7 (x[4], x[12], x[3], x[9], x[7], x[13], x[8])
        T7 = abs(
            -0.5*dem_array[0:rows-4, 4:cols] - 0.5*center - 
            K_COEFF*dem_array[0:rows-4, 3:cols-1] - K_COEFF*dem_array[1:rows-3, 4:cols] - 
            K_COEFF*dem_array[1:rows-3, 2:cols-2] - K_COEFF*dem_array[2:rows-2, 3:cols-1] + 
            K_FACTOR*dem_array[1:rows-3, 3:cols-1]
        )
        abs_diff_terms.append(T7)
        
        # T8 (-x[2], +2*x[7], -x[12])
        T8 = abs(-dem_array[0:rows-4, 2:cols-2] + 2*dem_array[1:rows-3, 2:cols-2] - center)
        abs_diff_terms.append(T8)
        
        # T9 (x[6], x[18], x[7], x[11], x[13], x[17], x[12])
        T9 = abs(
            -0.5*dem_array[1:rows-3, 1:cols-3] - 0.5*dem_array[3:rows-1, 3:cols-1] - 
            K_COEFF*dem_array[1:rows-3, 2:cols-2] - K_COEFF*dem_array[2:rows-2, 1:cols-3] - 
            K_COEFF*dem_array[2:rows-2, 3:cols-1] - K_COEFF*dem_array[3:rows-1, 2:cols-2] + 
            K_FACTOR*center
        )
        abs_diff_terms.append(T9)
        
        # T10 (-x[11], +2*x[12], -x[13])
        T10 = abs(-dem_array[2:rows-2, 1:cols-3] + 2*center - dem_array[2:rows-2, 3:cols-1])
        abs_diff_terms.append(T10)
        
        # T11 (x[16], x[8], x[7], x[11], x[13], x[17], x[12])
        T11 = abs(
            -0.5*dem_array[3:rows-1, 1:cols-3] - 0.5*dem_array[1:rows-3, 3:cols-1] - 
            K_COEFF*dem_array[1:rows-3, 2:cols-2] - K_COEFF*dem_array[2:rows-2, 1:cols-3] - 
            K_COEFF*dem_array[2:rows-2, 3:cols-1] - K_COEFF*dem_array[3:rows-1, 2:cols-2] + 
            K_FACTOR*center
        )
        abs_diff_terms.append(T11)
        
        # T12 (-x[7], +2*x[12], -x[17])
        T12 = abs(-dem_array[1:rows-3, 2:cols-2] + 2*center - dem_array[3:rows-1, 2:cols-2])
        abs_diff_terms.append(T12)

        # Sum all the vectorized terms and divide by 12.0
        rri_window[:] = np.sum(abs_diff_terms, axis=0) / 12.0
        
        return rri_array

class Toolbox(object):
    def __init__(self):
        self.label = "RRI Toolbox (optimized)"
        self.alias = "RRIopti"
        self.description = "Calculate Radial Roughness Index (Optimized)"
        self.tools = [RRIFunction]
        