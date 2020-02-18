# Temporal Clustering of Fast Coronal Mass Ejections

Minimum Requirements:

    Python 3.5.2

    --> numpy==1.15.0

    --> numpy-indexed==0.3.5

    --> scipy==1.1.0

    --> matplotlib==2.2.2


**Acknowledgements**

The analysis methods in this repository (see `PyCodes`) make use of CME data made available made available by [SOHO LASCO](https://cdaw.gsfc.nasa.gov/CME_list/) and sunspot data by [SILSO](http://www.sidc.be/silso/datafiles).

The max spectrum method was developed by [Stoev, S. A., Michailidis, G., & Taqqu, M. S. (2006, 2011), Estimating Heavy–Tail Exponents Through Max Self–Similarity](https://arxiv.org/abs/math/0609163v1). The max spectrum method was applied for the first time by [Ruzmaikin, A., Feynman, J., & Stoev, S. A. (2011), Distribution and Clustering of Fast Coronal Mass Ejections, J. Geophys. Res., 116, A04220](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JA016247) to the analysis of CME velocity time-series for SC 23, for which they adopted the 1000 km/s speed threshold. Some of the additional time-series analyses presented here also follow the latter authors. The results presented here corroborate their analysis.

The '.fits' files in this repository (see `Data`) contain information about the solar surface and are [courtesy of NASA/SDO and the AIA, EVE, and HMI science teams](http://soi.stanford.edu/magnetic/index6.html).  [APLpy,](https://aplpy.github.io/) an open-source plotting package for Python [(Robitaille and Bressert, 2012)](http://adsabs.harvard.edu/abs/2012ascl.soft08017R), was used to convert these .fits files into viewable images that show ARs evolve over the course of solar rotations.






##
