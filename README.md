# Temporal Clustering of Fast Coronal Mass Ejections

Minimum Requirements:

    Python 3.5.2
    
    --> numpy==1.15.0
    
    --> numpy-indexed==0.3.5
    
    --> scipy==1.1.0

    --> matplotlib==2.2.2


**Synopsis**

The goal of my Master's Thesis was to investigate the clustering in time of fast coronal mass ejections (CMEs) during solar cycles SC 23 and SC 24. These extreme events are distinguished by a speed threshold, which we found from the analysis to be 800 km/s. Using data made available by [SOHO LASCO](https://cdaw.gsfc.nasa.gov/CME_list/), we found that the weighted fit of the averages of the maximum speeds of these extreme events via the max spectrum method obey a power-law over progressively increasing time-scales; in other words, a change in the temporal scale corresponds to a proportional relative change in the average maximum CME speed. This power-law exponent is comparable for both solar cycles despite the fact that SC 23 was much more active. 

We also show that these events are not independent, but rather cluster in time. We estimated this level of dependency and the corresponding time threshold that is used to organize these events into independent clusters of events. From these clusters, we obtain inter-exceedance times and cluster durations. It was found that SC 23 exhibited a higher level of dependency (and consequently a lower time threshold) than SC 24. 

[Previous results in the literature](https://arxiv.org/abs/1110.1787) suggest that a single AR cannot produce more than one fast CME within a time-interval of 15 hours and that the average waiting time for same-AR CMEs is 8 hours. Our results are consistent with this observation. However, we also find multiple extreme events that occur within time-scales less than 15 hours and some less than 8 hours; some are only one hour apart! To resolve this, we suggest that these extreme events separated by time-scales less than 15 hours are produced by networks of multiple ARs; this was the case for many of the observed extreme events and has also [been observed independently](https://arxiv.org/abs/1505.01384).
    
**Acknowledgements**

The max spectrum method was developed by [Stoev, S. A., Michailidis, G., & Taqqu, M. S. (2006, 2011), Estimating Heavy–Tail Exponents Through Max Self–Similarity](https://arxiv.org/abs/math/0609163v1). The max spectrum method was applied for the first time by [Ruzmaikin, A., Feynman, J., & Stoev, S. A. (2011), Distribution and Clustering of Fast Coronal Mass Ejections, J. Geophys. Res., 116, A04220](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JA016247) to the analysis of CME velocity time-series for SC 23, for which they adopted the 1000 km/s speed threshold. Some of the additional time-series analyses presented here also follow the latter authors. The results presented here corroborate their analysis.
    
**Supplementary**

ARs can be seen evolving over the course of solar rotations. Images of the Sun (via '.fits' files) taken chronologically (one image per solar rotation) were converted into an animation; see the ['Image Manipulation' repository](https://github.com/allthemikeysaretaken/Image-Manipulation).
    
