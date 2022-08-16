# An Analysis of Gaia's data, estimating the distance of cluster NGC2506 using Bayesian statistics.

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

## Code Overview
1. We read the entire data file on cluster members into a Pandas dataframe and perform data cleaning. From the cleaned dataframe, we make a new dataframe for the data for cluster NGC2506, removing from it any stars with PMemb<1. 

2. We make a scatter plot of the apparent G magnitude vs BP-RP colour. 

3. We check that there are no flux-dependent biases in the parallax which might affect our results. 

4. We use the NGC2506 parallax data with Bayes’ theorem, to calculate the posterior pdf for the distance 𝑑 (in kpc) to NGC2506, using the formula 𝑑 = 1/𝑝 where 𝑝 is the parallax in milliarcsec (mas). Gaia has a known‘zero-point’ offset - a systematic error – in the parallax, so before we do our calculation we should first add a correction of 0.029 mas to the parallax measurements. We assume that the corrected parallax measurements are normally distributed about the true parallax, with standard deviation given by the errors on the parallax measurements. We plot thr posterior pdf and determine the 1-𝜎 confidence interval on the distance and plot the interval on our 
pdf.

5. Finally, we choose another open cluster in the data set (in this case NGC2168), remove stars with PMemb<1 and obtain the posterior distribution. Then we plot this cluster and NGC2506 on the same colour-magnitude diagram, but using absolute G magnitudes (corrected to a common distance of 10 pc), so that we can compare the diagrams for each cluster. For the purposes of estimating a distance, we assume the best distance for each cluster corresponds to the maximum of the posterior pdf (known as the ‘maximum likelihood estimate’). 
