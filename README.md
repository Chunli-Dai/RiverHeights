# RiverHeights
River shoreline detection and elevation extraction.

Here are the codes and data files for generating figures in the paper by Dai et 
al. (Estimating river surface elevation from ArcticDEM, GRL, 2018).

Step 1: Shoreline detection using entropy, brightness, and NDWI method: 
maskentropy.m: data processing for Fig. 1.
Input files (orthoimage, metafile, matchtag) are not shared here due to copyright limitation.
tancl.mat: river center line file for winter cases.
plotentropy.gmt, plotentropyb.gmt, plotentropyc.gmt, plotentropyd.gmt: plotting 
script of Fig. 1.
2013MayOr.dat, 2013MayJ.dat, 2013MayMask.dat, 2013MayCoast.dat: data files of Fig. 1, and some of the files are too big for uploading.
multispec.m: NDWI method for shoreline detection. 
Input files: Original multi-spectral imageries are not shared.
Output:shoreline file, coastndwi.dat, not shared here.
plotentropyd2.gmt, plotentropyd2b.gmt: script for plotting Fig. S1.

Step 2: Elevation extraction.
riverprof.m: data processing for extracting river heights. (Some of the subrouti
nes are not shared here due to copyright limitation).
FAirport2.gmt, park.gmt, park2.gmt, road.gmt: control surface files for coregist
ration.

Step 3: Filtering and Fitting.
ProcessTananaFairbanks.m: data processing and plotting of Fig.2.
Preparation:Creat a folder name Elevations in the working directory for storing data files.
tan_cl_Close7.shp: manually drawn river centerlines.
rivp20110804b.mat, rivp20121011a.mat: example data for Fig. 2b. 

Other codes and subroutines see https://github.com/mikedurand/SmoothRiverElevations

Step 4: Plot river height time series and discharge time series.
gageheights.m: plotting the time series of river heights at the USGS gage in Fairbanks, Alaska (Fig 3a).
Preparation: 
gagefo.txt, river height time series at the gage for all seasons.
gageft.txt, same as above but for winter season only.
usgsgage.txt, USGS gage height time series.

stagedischarge.m: Plot Fig.3b for the discharge time series.
Preparation:
uv10to16.txt: data files for stage-discharge rating curve.
legs.m, legsd.m: computation of legendre function for fitting rating curve.
scale4legs.m: preparation of the variables for legendre polynomials, which map variables into range -1 to 1.
Subroutines legs.m, legsd.m, and scale4legs.m are provided by Prof. Michael G. Bevis.

The comparison of river height time series by two different methods.
comparemethods.m


