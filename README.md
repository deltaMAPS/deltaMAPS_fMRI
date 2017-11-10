# deltaMAPS repository (fMRI)

This repository contains the java code for the application of deltaMAPS in FMRI data.
deltaMAPS is described in detail in https://arxiv.org/abs/1602.07249.

The two separate subfolders contain the java source code for (a) identifying domains in fMRI data and (b)
infering the network between the domains. 

The domain identification part of the code accepts as an input a surface file (either for the right or left cortical hemisphere) and a file with the timeseries for each surface voxel (the time series should be preprocessed as described in detail in the deltaMAPS paper). The surface file consists of triangles of connected voxels. Each voxel has a unique id. Similarly, the time series file consists of the id of the voxel followed by the time series of the surface voxel (all values are comma separated). 

Example of two surface files are  also provided. 
