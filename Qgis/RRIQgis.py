from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingParameterRasterLayer
from qgis.core import QgsProcessingParameterRasterDestination
from qgis.core import QgsProcessingException
import processing
import os

class RRITool(QgsProcessingAlgorithm):
    """
    RRI (Radial Roughness Index) Tool
    Calculates surface roughness using Radial Roughness Index
    """
    
    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT,
                'Input raster layer'
            )
        )
        
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT,
                'Output RRI raster'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        input_layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        output_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)
        
        # Ensure output has proper file extension
        if not output_path.lower().endswith('.tif'):
            output_path += '.tif'
        
        feedback.pushInfo(f"Input layer: {input_layer.name()}")
        feedback.pushInfo(f"Output path: {output_path}")
        
        # The complex RRI expression
        expression = (
            "(abs(-0.5*A[-2,-2]-0.5*A[0,0]-0.207106781186547*A[-2,-1]-0.207106781186547*A[-1,-2]"
            "-0.207106781186547*A[-1,0]-0.207106781186547*A[0,-1]+1.82842712474619*A[-1,-1])"
            "+abs(-A[0,-2]+2*A[0,-1]-A[0,0])"
            "+abs(-0.5*A[2,-2]-0.5*A[0,0]-0.207106781186547*A[2,-1]-0.207106781186547*A[1,-2]"
            "-0.207106781186547*A[0,-1]-0.207106781186547*A[1,0]+1.82842712474619*A[1,-1])"
            "+abs(-A[2,0]+2*A[1,0]-A[0,0])"
            "+abs(-0.5*A[2,2]-0.5*A[0,0]-0.207106781186547*A[2,1]-0.207106781186547*A[1,2]"
            "-0.207106781186547*A[0,1]-0.207106781186547*A[1,0]+1.82842712474619*A[1,1])"
            "+abs(-A[0,2]+2*A[0,1]-A[0,0])"
            "+abs(-0.5*A[-2,2]-0.5*A[0,0]-0.207106781186547*A[-2,1]-0.207106781186547*A[-1,2]"
            "-0.207106781186547*A[-1,0]-0.207106781186547*A[0,1]+1.82842712474619*A[-1,1])"
            "+abs(-A[-2,0]+2*A[-1,0]-A[0,0])"
            "+abs(-0.5*A[-1,-1]-0.5*A[1,1]-0.207106781186547*A[-1,0]-0.207106781186547*A[0,-1]"
            "-0.207106781186547*A[0,1]-0.207106781186547*A[1,0]+1.82842712474619*A[0,0])"
            "+abs(-A[0,-1]+2*A[0,0]-A[0,1])"
            "+abs(-0.5*A[1,-1]-0.5*A[-1,1]-0.207106781186547*A[-1,0]-0.207106781186547*A[0,-1]"
            "-0.207106781186547*A[0,1]-0.207106781186547*A[1,0]+1.82842712474619*A[0,0])"
            "+abs(-A[-1,0]+2*A[0,0]-A[1,0]))/12"
        )
        
        # Set up GRASS r.mapcalc.simple parameters
        params = {
            'a': input_layer,
            'expression': expression,
            'output': output_path,
            'GRASS_REGION_PARAMETER': input_layer.extent(),
            'GRASS_REGION_CELLSIZE_PARAMETER': 0,
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_RASTER_FORMAT_META': ''
        }
        
        try:
            # Execute the algorithm
            feedback.pushInfo("Calculating Radial Roughness Index (RRI)...")
            feedback.pushInfo("This may take a while for large rasters...")
            
            result = processing.run("grass7:r.mapcalc.simple", params, context=context, feedback=feedback)
            
            # Return the output
            if 'output' in result:
                feedback.pushInfo("RRI calculation completed successfully!")
                return {self.OUTPUT: result['output']}
            else:
                raise QgsProcessingException("Output file was not created successfully")
                    
        except Exception as e:
            feedback.pushError(f"Error during RRI calculation: {str(e)}")
            raise QgsProcessingException(f"RRI calculation failed: {str(e)}")

    def name(self):
        return 'rritool'

    def displayName(self):
        return 'RRI (Radial Roughness Index)'

    def group(self):
        return 'Raster Analysis'

    def groupId(self):
        return 'rasteranalysis'

    def shortHelpString(self):
        return """
        RRI (Radial Roughness Index) Tool
        
        This tool calculates the Radial Roughness Index using a specialized 
        5x5 convolution kernel that analyzes surface texture in radial directions.
        It is an improvement of TRI (Topographic Ruggedness Index)
        https://doi.org/10.1016/j.geomorph.2023.108838
        
        <b>Input:</b> Any raster layer (DEM or image)
        <b>Output:</b> RRI, same units of input raster, where higher values indicate rougher surfaces
        
        The tool automatically applies the RRI calculation without requiring
        manual expression input.
        """

    def createInstance(self):
        return RRITool()
        