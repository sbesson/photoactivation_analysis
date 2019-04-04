## README

This bundle contains the series of scripts used to produce the analysis
of [Young et al](https://doi.org/10.1242/bio.201410363). All scripts are
written under MATLAB and have been tested using MATLAB R2013b. The following
MATLAB toolboxes are required to run all the components:

- Curve Fitting Toolbox
- Image Processing Toolbox
- Mapping Toolbox
- Optimization Toolbox
- Statistics Toolbox


The bundle contains 3 directories:

-   the [spindleAnalysis](spindleAnalysis) directory contains the scripts used
    to align the spindle width/length measurements. The startup script is
    called [runSpindleAnalysis.m](spindleAnalysis/runSpindleAnalysis.m) and
    looks for XLS files containing time-series of the spindle lengths.

-   the [photoactivation](photoactivation) directory contains the scripts used
    to process the spindle photoactivation experiments into flux/turnover measurements. It assumes the data has been imported onto an OMERO server and makes extensive use of the OMERO.matlab toolbox 

-   the [extern](extern) directory contains third-party MATLAB functions under
    GPL-license used for running the code.
