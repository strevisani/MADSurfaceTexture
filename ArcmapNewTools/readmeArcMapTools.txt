30th September 2022

Here you will find new developed arcmap tools related to roughness calculation, based on the following papers (one published the other a pre-print).
1) Trevisani, S. and Teza, G. and Guth, P., A Simplified Geostatistical Approach for Characterizing Key Aspects of Short-Range Roughness. 
Available at SSRN: https://ssrn.com/abstract=4223135 or http://dx.doi.org/10.2139/ssrn.4223135
2)Trevisani S., Rocca M.,MAD: robust image texture analysis for applications in high resolution geomorphometry. 
Computer & Geosciences, 2015 https://doi.org/10.1016/j.cageo.2015.04.003

Currently you will find two toolboxes "Mad2cv1" and "Mad1ck2": 
the first calculates short-range roughness with the usual method with detrending (it does a simple detrending based on moving averages);
the latter adopt differences of order 2 and bypass detrending (see http://dx.doi.org/10.2139/ssrn.4223135).


Then there are two folder with kernels of order 1 and of order 2 (respect to the version of 2014 in these the weights are not rounded).
These are useful if you need to use custom scripts and not the toolboxes.

Finally I provide also a function for computing roughness based on the dispersion of normal vectors (VRM).


A really important thing!!!!!!!!!!!!!!!!!!
Arcmap has some trickeries (at least on window systems) in managing the international settings (i.e., comma or point as decimal separator). The
code uses custom weigth kernels that use the point as decimal separator. I suggest to set in international setting as English (US) and 
the point as decimal separator. If you get an error message "invalid weight neighborhood mask" is related to this issue. 

bye

Sebastiano

