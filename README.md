#Quick update December 2021
This is a quick update to comunicate that in the following weeks I'll upload new code (as soon as I go on with coding and testing I'll add it) 
related to geostatistical based indices for surface/image texture analysis.

I will do some minor updates for the version in phyton for arcmap/arcgis pro (new kernels). In this environments I'm adding a new simplified toolbox and a function for computing
vector ruggedness (i.e., dispersion of vector normals to surface). You will find new arcmap arcgis pro staff in the folder "ArcmapNewTools". 

However, the main change is related to the new coe for R environment, using the facilities of Terra package. Still to upload things.

Regarding the kernels, these will be updated using more significative digits.



########################################


# MADSurfaceTexture
MAD functions for surface/image texture characterization

### Begin of license ########################################################

Copyright (c) 2014 Sebastiano Trevisani

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

- 1)The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
- 2)When using this software, in particular for publications, please cite the related paper:
"S. Trevisani and M. Rocca. MAD: robust image texture analysis for applications in high resolution geomorphometry. Computer & Geosciences (2015), 10.1016/j.cageo.2015.04.003."
(substitute the DOI with the correct volume and pages numbers when available).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

### End of license ##########################################################

## Instructions

In this folder there is the code necessary to test the algorithm and understand its structure (and rebuild the procedure in other software environments).

The software has been tested with ArcGis 10.2, 10.2.1. and 10.3. Take care of how international settings are interpreted by python interpreter in Arcgis 10.41: be sure that python uses correctly the decimal separator (i.e., ".")! Otherwise the kernels with decimal values will not work!



List of items and quick description:
- 1. MADfunctionsV1.py --> Here we find the core functions definition related to MAD indexes. These are used both by "exampleScript.py" as well as by the developed toolbox.
- 2. exampleScript.py --> This script describes the usage of MAD function from the interactive Python shell of ArcGIS
- 3. MADwin	-->	This folder contains the Arctoolbox "madscan" that can be called directly from ArcCatalog. This is a window-based tool (i.e., no scripts). From ArcCatalog, just go in the directory of MADWin and open the interactive tool.
- 4. Kernels --> The set of kernels used by MAD functions. We preferred give these as static files so as to facilitate their re-use in other software.
- 5. Sample --> (compressed folder) Sample data to be used to test the code. The 1 m  Residual HR-DTM furnished in the code has been kindly offered by the Dr. Marco Cavalli, Research Institute CNR-IRPI of Padova, Italy (marco.cavalli@irpi.cnr.it).

The use of function from python scripting is quite easy. A quick way to be sure all works fine is to save function definition in the same folder of your ArcGIS project (*.mxd file) together with kernels directory. Then you
should define an output folder as well as a scratch folder where to save intermediate results (use environment settings of arctoolbox). From our experience we suggest to don't use a geodatabase for storing intermediate outputs but a simple folder. This because of with large files (DTMs on the order of 600 Mb or more) geodatabase can give problem of access on some PCs. Remember that these functions make use of Spatial Analyst ArcGis extension.

## Notes on directional kernels derivation
In this section we add some details regarding the calculation of directional differences that we cannot include in the paper for the restrictions in paper length. The basic steps for MAD as well as other bivariate spatial continuity indexes (variogram, madogram, etc.) is to calculate quickly and accurately differences between points pairs. This is accomplished with ad hoc defined kernels that permit via focal statistic functions to perform in one step bilinear interpolation and to calculate the differences between the interpolated points. As reported in the text the basic kernels for a multiscale directional analysis are based on the scheme represented in the figure 10 of the paper. Under this scheme the directional difference is automatically associated to the central pixel. In the code we defined kernels for different lag distances (2, 4, 6, and 8 pixels) and for four directions (N-S, SW-NE, W-E, NW-SE). Following the same approach, i.e. deriving the weight of bilinear interpolation for a given distribution of the point pairs,  other kernels can be defined in order to perform the calculation in more directions and for larger distances (however, keep in mind that large kernels increase computational cost).
Two special sets of kernels are mentioned in the work (see figure 11). A first kernel permits (figure 11 left) the calculation of directional differences for a lag distance of one cell.  Using the scheme presented in the figure for a specific direction we derive two opposite directional differences. For example, in NS direction we derive the difference between central pixel and the pixel one step at N and the difference between the central pixel and the pixel at one step at S. Consequently we need to average their absolute value to have a symmetric value to associate to the central pixel and to be used for MAD calculation. The results obtained using one cell  kernel can be used in analogy to the output of the basic kernels (but having directly the absolute value); also consider that some smoothing of variability should be expected. It is interesting to note that this kernel-approach can be also modified to compute non-symmetrical surface texture indexes, i.e. where Index(h) is different from index(-h) and also for larger lags. Note the strong similarities to the local binary approach.  
We also defined another kernel-approach (fig. 11, right) that calculates directional differences using the shortest lag along diagonals with a lag of 1.414*cell size. This function can be useful in some circumstances, especially if it should be adapted for the calculation of directional differences in many locally variable directions. The first step is to calculate (figure above) directional differences along the diagonals; we can derive the directional differences, centered on the corners of the central pixel, in any direction using simple trigonometry. In the kernel here used we implemented the NS and EW directions to maintain the same scheme of basic kernels, but a function taking as input the angle direction can be easily coded. After this step, four directional differences are calculated in correspondence of the four corners of the central pixel. Again in order to obtain a cell centered absolute difference we need to average the absolute values of the four directional differences. This function can be useful for special needs and eventually for calculation of relative roughness. Again, in this case some smoothing of roughness is expected.

