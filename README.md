# Temporal Clustering of Fast Coronal Mass Ejections

Minimum Requirements:

    Python 3.5.2

    --> numpy==1.15.0

    --> numpy-indexed==0.3.5

    --> scipy==1.1.0

    --> matplotlib==2.2.2


**Acknowledgements**

The analysis methods in this repository make use of CME data made available made available by [SOHO LASCO](https://cdaw.gsfc.nasa.gov/CME_list/) and sunspot data by [SILSO](http://www.sidc.be/silso/datafiles).

The max spectrum method was developed by [Stoev, S. A., Michailidis, G., & Taqqu, M. S. (2006, 2011), Estimating Heavy–Tail Exponents Through Max Self–Similarity](https://arxiv.org/abs/math/0609163v1). The max spectrum method was applied for the first time by [Ruzmaikin, A., Feynman, J., & Stoev, S. A. (2011), Distribution and Clustering of Fast Coronal Mass Ejections, J. Geophys. Res., 116, A04220](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JA016247) to the analysis of CME velocity time-series for SC 23, for which they adopted the 1000 km/s speed threshold. Some of the additional time-series analyses presented here also follow the latter authors. The results presented here corroborate their analysis.

The '.fits' files in this directory contain information about the solar surface and are [courtesy of NASA/SDO and the AIA, EVE, and HMI science teams](http://soi.stanford.edu/magnetic/index6.html).  APLpy, an open-source plotting package for Python [(Robitaille and Bressert, 2012)](http://adsabs.harvard.edu/abs/2012ascl.soft08017R), was used to convert these .fits files into viewable images that show ARs evolve over the course of solar rotations.

**Synopsis**

The goal of this thesis was to investigate the clustering in time of fast coronal mass ejections (CMEs) during solar cycles SC 23 and SC 24. CMEs are violent eruptions of plasma in the outer atmosphere of the Sun that eject large amounts of particles and magnetic flux into space. Many CMEs originate from active regions (ARs), which are areas in the photosphere of especially strong magnetic fields. CMEs can range in speeds from tens of km/s up to 3500 km/s, though extreme events are the most important in the context of space-weather. Extreme events are characterized as only those CMEs that are as fast or faster than a given speed threshold.

Additional documentation will be developed and elaborated upon. Until then, the main points are covered below. There are additional methods in the PyCodes directory, and additional figures/files in the Results directory. These routines can be used to analyze more solar cycles as more data becomes available. More detailed information can be found [here](http://scholarworks.csun.edu/handle/10211.3/207815); although the methodology is the same, the algorithms have since been updated.

First, we compared SC 23 and SC 24 by looking at the frequencies of sunspots, CMEs, and extreme events - all taken daily, monthly, or yearly. As shown in the figure below, SC 24 had a slightly higher frequency of events than SC 23, but SC 23 contained more extreme events. One can also see that the number of extreme events appears to be directly proportional with the number of observed sunspots.

![Sunspot-CME Correlation](https://i.imgur.com/FHwVsqN.png)    

We then looked at the distribution of CME speeds. For both SC 23 and SC 24, we computed a histogram of all events faster than 20 km/s. These speeds follow a lognormal distribution with heavy tails, allowing us to optimize the fit of these speed distributions. The optimization methods used in this analysis were maximum likelihood estimation and minimum g-test estimation; the code allows for minimum chi square estimation as an alternative to minimum g-test estimation.

![Probability Density of Speeds via MLE](https://i.imgur.com/cH9QPf4.png)    

Then, we compared the distribution of inter-exceedance times (time in-between extreme events) with the distribution of the corresponding exponential distribution; this exponential distribution was obtained by applying the inverse-transform sample method. Both distributions share the same standard deviation. The figure below shows that extreme events tend to cluster at much shorter time-scales and that the distribution of inter-exceedance times decays faster than the corresponding exponential distribution.

![Distribution of Inter-Exceedance Times](https://i.imgur.com/rTgGNI5.png)    


We then examined the events in the tail of the distribution and found that the weighted fit of the averages of the maximum speeds of these extreme events via the max spectrum method obey a power-law over progressively increasing time-scales; in other words, a change in the temporal scale corresponds to a proportional relative change in the average maximum CME speed. 

![Power-Law Errors](https://i.imgur.com/cK2Vhsn.png)    

![Power-Law](https://i.imgur.com/18rBP4M.png)    

This power-law exponent is comparable for both solar cycles despite the fact that SC 23 was much more active.

![Histogram of Power-Law Exponent](https://i.imgur.com/vlOWYS1.png)    


The similarity of the power-law exponent across both SC 23 and SC 24 is reflected in the tail of the speed distributions; note that they appear to decay at comparable rates.

![Decay of Tail of Speed Distribution](https://i.imgur.com/QhObbHz.png)    

We then estimated the extremal index, which corresponds to the intercept parameter of the power-law. This allows us to group temporally dependent extreme events into independent clusters of events. The extremal index ranges from zero (high levels of dependence) to one (all events are independent). This extremal index can serve as a global estimate over all time-scales, or as an estimate at a particular time-scale; both are shown below.

![Extremal Index (global)](https://i.imgur.com/ZfSw4xn.png)    


![Extremal Index (by time-scale)](https://i.imgur.com/muI3xMU.png)    


We then calculated the moment estimator of the extremal index as a function of speed threshold. These moment estimators can be used to calculate the critical time threshold that is used in conjunction with a condition applied to inter-exceedance times in order to group extreme events into clusters. 

![Moment Estimators & Time Threshold](https://i.imgur.com/ywdHHwd.png)    


Once these clusters have been obtained, we can note various cluster statistics. These include relative statistics, such as the average cluster size and duration, as well as histograms of the intra-cluster times (time in-between events within a cluster), intra-cluster durations (duration of cluster), and inter-cluster durations (time in-between clusters). 


![Relative Clustering Statistics](https://i.imgur.com/7QpEnaD.png)    

![Intra-Event Times](https://i.imgur.com/B6rs5Nf.png)    

![Intra-Cluster Durations](https://i.imgur.com/VF7Madz.png)    


[Previous results in the literature](https://arxiv.org/abs/1110.1787) suggest that a single AR cannot produce more than one fast CME within a time-interval of 15 hours and that the average waiting time for same-AR CMEs is 8 hours. Our results are consistent with this observation. However, we also find multiple extreme events that occur within time-scales less than 15 hours and some less than 8 hours; some are only one hour apart! To resolve this, we suggest that these extreme events separated by time-scales less than 15 hours are produced by networks of multiple ARs; this was the case for many of the observed extreme events and has also [been observed independently](https://arxiv.org/abs/1505.01384).





##
