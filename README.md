# A statistical Analysis of Gaia's data

## Introduction

For this project we use a data file which contains astrometric and photometric data for 401448 stars observed by the Gaia mission, 
which have been identified as belonging to 1229 open star clusters in our galaxy. Besides measuring extremely accurate positions on the sky, Gaia provides 
measurement of the star’s proper motion and parallax. The parallax can be used to directly estimate the distance to the star. 

The data is taken from an analysis of Gaia’s Data Release 2 (DR2), a database of 1.3 billion sources which includes the astrometric data as well as 
photometry in 3 optical bands. The open cluster study used a sophisticated clustering analysis of the astrometric information to identify which objects 
belong with high probability to a known open cluster that is along the line-of-sight.

[Data](https://drive.google.com/file/d/1rcn2CV_0qCfbMOs8UHsMpSdTlnrylA71/view?usp=sharing)

The data columns are described in the additional file `data_description.txt`.


The astrometric data are especially complex, with the method for making the astrometric solution leading to correlations between the different astrometric 
quantities. However, here we are only be focusing on the parallax measurements `plx` and `e_plx`, alongside the photometric apparent magnitude in the 
G band `Gmag` and the photometric colour  `BP-RP`. Note that the magnitude is a logarithmic measure of the flux with a scale runs backwards, 
i.e. brighter  sources have smaller magnitudes. The cluster member data is limited to stars brighter than G magnitude 18. We use the Cluster name which
the star is associated with and the membership probability `PMemb` which is less than 1 for a number of outlier cases that are not certain to be associated
with that cluster.  
