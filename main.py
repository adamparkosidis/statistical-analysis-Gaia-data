import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sps
import scipy.integrate as spint
import scipy.interpolate as spinterp
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import scipy.optimize as spopt
import os


########################
### Define Functions ###
########################

def corrs(x,y):
    '''
    calculate the Spearman coefficient, rho, and
    the corresponding p-value
    '''
    
    (rhocor, rhopval) = sps.spearmanr(x,y)
    print("For quantities {:s} and {:s}, Spearman's rho and p-value are {:.4f} and {:.4e} respectively \n".format(x.name,y.name,rhocor, rhopval))
    
def constrain_cl(df_name):
    '''
     For a given cluster removes any star that has a probability smaller than 1 of being part of the cluster
     In other words ensures (based on our data) that the remaining stars are trully member of the cluster
    '''
    df = cluster_df[cluster_df['Cluster'] == df_name]
    df = df[df['PMemb'] >= 1]
    return df

def prob(df, d_array):
    '''
    calculate the posterior pdf for the distance d (in kpc) to a given cluster assuming a uniform prior. 
    The prior probability should integrate to 1, although since the  normalisation of the prior 
    divides out of Bayes' formula it could have been arbitrary. Furthermore, I
    assume that the corrected parallax measurements are normally distributed about the true 
    parallax, with standard deviation given by the errors on the parallax measurements.
    '''
    par = df['plx']
    p_errors = df['e_plx']
    prior = 2                                         
    # We need to sum the individual log-likelihoods and also the log(prior):
    loglikel_prior = np.sum(np.log(sps.norm.pdf(par.values.reshape(len(par),1), loc=p_vals, scale=p_errors.values.reshape(len(p_errors), 1))), axis=0)+np.log(prior)
    likel_prior = np.exp(loglikel_prior-np.amax(loglikel_prior))
    likel_prior_int = spint.simpson(likel_prior,d_array,axis=0)
    
    # Now we normalise and we have our posterior pdf for d
    posterior_pdf = likel_prior/likel_prior_int
    
    prob_d = d_array[np.argmax(posterior_pdf)]                  # most likely value from our measurment
    standard_error = prob_d/np.sqrt(len(posterior_pdf))         # the standard error on the value
    
    posterior_cdf = spint.cumulative_trapezoid(posterior_pdf,d_array,initial=0) 
    posterior_ppf = spinterp.interp1d(posterior_cdf,d_array)
    int68 = posterior_ppf([0.16,0.84])                            # 0.5 +- p/2 
    print(r'68% confidence interval on lambda =', int68)
    return posterior_pdf, prob_d, standard_error, int68

#################
### Main Code ###
#################

# Read and explore the data

path =o s.getcwd()

cluster_df = pd.read_csv(path+'/cluster_members.txt', delim_whitespace=True, na_values='---')       # read the data and save it in a DataFrame
cluster_df.info()                            # print info about my data such as type, size ,nan-values etc. 
cluster_df.head()                            # visualize the head of my data to have a better insight

# clean the dat afrom nan-values
cluster_df = cluster_df.dropna()

# Make a new dataframe using data from cluster NGC2506 and remove any stars with PMemb < 1
NGC_2506= constrain_cl('NGC_2506')


# we create a color-magnitude diagram with the apparent G magnitude in the Y-axis and the BP-RP colour in the X-axis

sns.set_style('darkgrid')
plt.figure()
ax = sns.scatterplot(data=NGC_2506, x='BP-RP', y='Gmag')
ax.set_xlabel('BP-RP [mag]', fontsize = 10)
ax.set_ylabel('Gmag [mag]', fontsize = 10)
ax.invert_yaxis()
plt.tight_layout
plt.show()

'''
Analysis: 
Open clusters have a specific age and from the resulting colour-magnitude diagram it is evident that the stars from the top left are missing, since the
bluest, most luminous stars (also shortest-lived) have already reached the end of their lives and are no longer seen.
'''

# We will now create a scatter-plot of the parallax (in milliarcsec) vs. the G magnitude for NGC2506 and try briefly to explain why it looks like that
# and why we would not expect to see this pattern for field stars (i.e. those along the line of sight but not associated with a single cluster).

# we create a scatterplot with the apparent G magnitude in the Y-axis and the parallax in the X-axis for each star of the NGC2506 cluster

plt.figure()
ax = sns.scatterplot(data=NGC_2506, x='plx', y='Gmag', hue='e_plx', palette='coolwarm')
ax.set_xlabel('plx [mas]', fontsize = 10)
ax.set_ylabel('Gmag [mag]', fontsize = 10)
ax.invert_yaxis()
plt.tight_layout
plt.show()

'''
Analysis: 

Let's try to interpret the scatterplot.The x-axis shows the parallax of each star, which is inversly proporsional to the distance of the star, thus as 
parallax becomes smaller we see stars that are further away from us, while as parallax becomes bigger we see stars that are closer to us. The y-axis shows
the apparent magnitude in the G band of each star. In our case, where we know that all the stars belong to the same cluster, the bigger the apparent 
magnitude the dimmer the star is. 

Regarding of the shape of the distribution, it is evident that the majority of the cluster's stars are concentrated towards
bigger $y$ values (dimmer stars). This is normal, because as we previously saw the stars of the cluster have approximatelly the same age and most of the 
bluest, most luminous stars, are already dead. 

Field stars can have significantly different age, significantly different distance and in general can be 
uncorrelated, thus in that case we would also see many luminous stars in our scatterplot, which wouldn't be members of the cluster.

In order to get further insight in the data we can also plot the error in parallax measurements. We can clearly see that the error in the measurements 
becomes bigger for stars with bigger apparent magnitudes. This pattern is logical, because it is more difficult to accurately measure dimmer stars and also
explains why we see a bigger dispersion for the stars with the bigger y-axis values. We can also plot the apparent G magnitude with error parallax to further
pin-point that pattern
'''

plt.figure()
ax = sns.scatterplot(data=NGC_2506, x='e_plx', y='Gmag')
ax.set_xlabel('e_plx [mas]', fontsize = 10)
ax.set_ylabel('Gmag [mag]', fontsize = 10)
ax.invert_yaxis()
plt.tight_layout
plt.show()

# We use an appropriate test on the NGC2506 data to determine whether the parallax depends on the G magnitude.
# From the previous scatterplot we can see that data does not imply a linear relationship and there are also some outliers. Hence, Spearman's method is 
# more appropriate to examine the correlation between the parallax and the G magnitude, even though it seems that there is not a monotonic relationship 
# between them. We assume that the data is **i.i.d.** and we apply Spearman's method:

corrs(NGC_2506['plx'], NGC_2506['Gmag'])

'''
Analysis:
Spearman's ρ is pretty close to zero, thus we can conclude that our initial suspicion is true and the parallax is not monotonic realated to the 
apparent G magnitude.
'''

systematic_error = 0.029
NGC_2506['plx'] = NGC_2506['plx'] + systematic_error   # corrected parallax

d_array = np.linspace(3.05,3.2,1000)  # we define a linespace of d values (where d is the distance)
p_vals = np.array(list(map(lambda x: 1/x,d_array)))  # from distance to parallax
posterior_pdf, dist_NGC_2506, standard_error, int68 = prob(NGC_2506, d_array) 


# And plot it...
plt.figure()
plt.plot(d_array,posterior_pdf)
plt.vlines(int68,ymin=0,ymax=25,color='gray',linestyle='dotted')
plt.xlabel(r'$d \; (kpc)$',fontsize=10)
plt.ylabel(r'posterior density',fontsize=10)
plt.show()
print(r'The area below the integral is {:.2}, the true most probable value of $d$ corresponds to {:.4f} $kpc$ with a standard error of {:.4f}. The dotted lines correspond to 1-σ confidence interval, where 68% of the of the distribution is constraint'.format(np.sum(posterior_pdf)*(d_array[1]-d_array[0]),dist_NGC_2506,standard_error))

# Finally, we chose another open cluster in the data set, remove stars with PMemb<1 and 
# obtain once again the posterior distribution. Then plot this cluster and NGC2506 on the 
# same colour-magnitude diagram, but using absolute G magnitudes (corrected to a common distance 
# of 10 pc), so that we can compare the diagrams for each cluster.
# We select the cluster NGC_2168 and remove any star that has a probability smaller than 1 of being part of cluster NGC_2168
# In other words we ensure (based on our data) that the remaining stars are trully member of the NGC_2168 cluster

NGC_2168= constrain_cl('NGC_2168')

# As mentioned, Gaia has a known ‘zero-point’ offset - a systematic error – in the parallax, thus 
# we first add a correction of 0.029 mas to the parallax measurements.

NGC_2168['plx'] = NGC_2168['plx'] + systematic_error   # corrected parallax


d_array = np.linspace(0.85,0.86,1000)   # we define a linespace of d values (where d is the distance)
p_vals = np.array(list(map(lambda x: 1/x,d_array)))  # from distance to parallax
posterior_pdf, dist_NGC_2168, standard_error, int68 = prob(NGC_2168, d_array)

# And plot it...
plt.figure()
plt.plot(d_array,posterior_pdf)
plt.vlines(int68,ymin=0,ymax=360,color='gray',linestyle='dotted')
plt.xlabel(r'$d \; (kpc)$',fontsize=10)
plt.ylabel(r'posterior density',fontsize=10)
plt.show()
print(r'The area below the integral is {:.2}, the true most probable value of $d$ corresponds to {:.4f} $kpc$ with a standard error of {:.4f}. The dotted lines correspond to 1-σ confidence interval, where 68% of the of the distribution is constraint'.format(np.sum(posterior_pdf)*(d_array[1]-d_array[0]),dist_NGC_2168,standard_error)) 


# we create a new column in the dataframes with the absolute magnitude assuming that the most probable distance for each cluster corresponds

NGC_2506= NGC_2506.assign(Absolute_G=lambda x: x['Gmag'] - 5 * np.log10(dist_NGC_2506) - 10)
NGC_2168= NGC_2168.assign(Absolute_G=lambda x: x['Gmag'] - 5 * np.log10(dist_NGC_2168) - 10)

# we create a color-magnitude diagram with the absolute G magnitude in the Y-axis and the BP minus RP colour in the X-axis 

plt.figure()
ax = sns.scatterplot(data=NGC_2168, x='BP-RP', y='Absolute_G', label='NGC_2168 Cluster')
sns.scatterplot(data=NGC_2506, x='BP-RP', y='Absolute_G', label='NGC_2506 Cluster')
ax.set_xlabel('BP-RP [mag]', fontsize = 10)
ax.set_ylabel('Absolute_G [mag]', fontsize = 10)
ax.invert_yaxis()
plt.tight_layout
plt.show()

'''
Analysis:
The above scatter-plot actually corresponds to an HR diagram and we can see the majority of stars of the two clusters being on the main-sequence. 
Knowing that more massive stars burn brighter and live shorter, we conclude that a cluster has lived only shortly if it still contains bright stars 
along its main sequence and has few stars that have evolved off the main sequence towards later evolutionary stages (such as the red giant or horizontal
giant branch). For NGC2506 cluster, we can clearly see the cut-off (around $Absolute_{G} =2$) on the main-sequence and a number of stars that have already
been evolved off towards later evolutionary stages. On the other hand, cluster NGC2168 still contains hotter, brighter and more massive stars (top left 
corner of the diagram), thus we can conclude that NGC2506 is older than NGC2168.
'''
