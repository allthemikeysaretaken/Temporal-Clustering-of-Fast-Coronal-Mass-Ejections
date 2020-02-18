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

**Example - Timelapse of Active Regions**

The following code (via `main_active_region_timelapse.py` in `PyCodes`) can be used to create a timelapse of solar active regions as a function of the Carrington Rotation Number. 

The directory-tree is identical to the directory-tree of this repository. In the first section, make sure to specify the home directory `dir_home` and import the necessary modules from `PyCodes`.

    from active_region_timelapse import *

    ## SPECIFY DIRECTORIES EXPLICITLY
    dir_home = '/Users/.../'
    dir_figs = '{}Figures/'.format(dir_home)
    dir_data = '{}Data/'.format(dir_home)
    dir_fits = '{}MichelsonDopplerImager/'.format(dir_data)
    dir_active_regions = '{}ActiveRegions/'.format(dir_figs)
    ##

In the next section, the .fits files - which correspond to Carrington Rotations 1990 - 2030 in this example - are converted to an image-type file (.png in this case). Each conversion allows for the optional use of a false-color mapping; if not specified, the grayscale image will be saved. The image shown below depicts Carrington Rotation 1990. Once the conversion is done and the images have been saved into `dir_active_regions`, these images can be animated as a timelapse.

    ## GENERATE TIME-LAPSE OF ACTIVE REGIONS
    fps = 5
    img_extension = '.png'
    dpi = 800
    mov_extension = '.mkv'
    for cmap in (None, 'hot'):
        MDI = MichelsonDopplerImager(rdirectory=dir_fits, sdirectory=dir_active_regions)
        MDI.convert_from_fits(cmap, img_extension, dpi)
        MDI.save_timelapse(fps, cmap, mov_extension, codec='mpeg4', search_extension=img_extension)
    ##

![Carrington Rotation 1990](https://i.imgur.com/KhgIuhk.png)  


**Example - Comparing Event Frequencies**

The following code (via `main_frequency_series.py` in `PyCodes`) can be used to compare frequencies of solar events as a function of time. 

The directory-tree is identical to the directory-tree of this repository. In the first section, make sure to specify the home directory `dir_home` and import the necessary modules from `PyCodes`. 

    from frequency_series_methods import *
    
    ## SPECIFY DIRECTORIES EXPLICITLY
    dir_home = '/Users/.../'
    dir_data = '{}Data/'.format(dir_home)
    dir_figs = '{}Figures/'.format(dir_home)
    dir_freq_series = '{}FrequencySeries/'.format(dir_figs)
    ##

In the next section, solar data (files in `Data`) is read and stored. This solar data consists of sunspot frequency, CME frequency, and solar flare frequency - all as a function of a timestep. 

    ## LOAD DATABASE
    DB = DataBase(directory=dir_data)
    DB.load_data(read_soho_lasco=True, read_silso=True, read_hessi=True)
    ##

    ## INITIALIZE FREQUENCY-SERIES
    timestep = 'month' # 'day', 'year'
    FS = FrequencySeries(DB, timestep, directory=dir_freq_series)

Then, one must specify the event type with optional search criteria (for example, CMEs above a certain speed). 

    search_kwargs = dict()
    search_kwargs['search_parameters'] = ('solar cycle', 'solar cycle')
    search_kwargs['search_conditions'] = ('greater than or equal', 'less than or equal')
    search_kwargs['search_values'] = (23, 24)
    FS.load_events('sunspot', **search_kwargs)
    #
    search_kwargs = dict()
    search_kwargs['search_parameters'] = ('solar cycle', 'solar cycle', 'linear speed')
    search_kwargs['search_conditions'] = ('greater than or equal', 'less than or equal', 'greater than')
    search_kwargs['search_values'] = (23, 24, 0)
    FS.load_events('cme', dt_prefix='', **search_kwargs)
    #
    search_kwargs = dict()
    search_kwargs['search_parameters'] = ('solar cycle', 'solar cycle', 'linear speed')
    search_kwargs['search_conditions'] = ('greater than or equal', 'less than or equal', 'greater than')
    search_kwargs['search_values'] = (23, 24, 800)
    FS.load_events('cme', dt_prefix='', **search_kwargs)
    #
    search_kwargs = dict()
    search_kwargs['search_parameters'] = ('start solar cycle', 'start solar cycle')
    search_kwargs['search_conditions'] = ('greater than or equal', 'less than or equal')
    search_kwargs['search_values'] = (23, 24)
    FS.load_events('flare', dt_prefix='start ', **search_kwargs)

When viewing the frequency comparisons between solar events, one has the option of selecting the layout. The image shown below uses `layout='vertical'`.

    for layout in ('overlay', 'vertical'):
        FS.view_frequency_distributions(layout=layout, save=True, figsize=(7,7), sharex=True)
        FS.view_frequency_distributions(layout=layout, save=True, show_solar_cycle_legends=True, figsize=(7,7), sharex=True)
        FS.view_frequency_distributions(layout=layout, save=True, show_solar_cycle_separations=True, figsize=(7,7), sharex=True)
        FS.view_frequency_distributions(layout=layout, save=True, show_solar_cycle_legends=True, show_solar_cycle_separations=True, figsize=(7,7), sharex=True)
    ##

![Frequency Comparison of Solar Events](https://i.imgur.com/2LRs8v8.png)  


**Example - CME Speed Distribution**

The following code (via `main_speed_series.py` in `PyCodes`) can be used to view the distribution of cme speeds. These can be viewed as probability distributions, distributions of observed events, histograms, optimized fits, additional methods to accentuate the distribution tail, and the error-space associated with the desired error metric used to optimize the fits. 

The directory-tree is identical to the directory-tree of this repository. In the first section, make sure to specify the home directory `dir_home` and import the necessary modules from `PyCodes`. 

    from speed_series_methods import *
    
    ## SPECIFY DIRECTORIES EXPLICITLY
    dir_home = '/Users/.../'
    dir_data = '{}Data/'.format(dir_home)
    dir_figs = '{}Figures/'.format(dir_home)
    dir_speed_series = '{}SpeedSeries/'.format(dir_figs)
    ##

In the next section, CME data (files in `Data`) is read and stored. 

    ## LOAD DATABASE
    DB = DataBase(directory=dir_data)
    DB.load_data(read_soho_lasco=True)
    ##

One can select a subset of the stored CME data that satisfies some specified conditions. 

    ## INITIALIZE SPEED-SERIES
    solar_cycles = (23, 24)
    speed_types = ('second order initial speed', 'linear speed')
    nan_policy = 'discard' # 'replace'
    nan_repl = 0
    nan_kwargs = dict(policy=nan_policy, repl=nan_repl)
    SS = SpeedSeries(DB, timestep='day', directory=dir_speed_series)
    
    ## LOAD EVENTS
    for solar_cycle in solar_cycles:
        for speed_type in speed_types:
            nan_kwargs['parameter'] = speed_type
            # SS.load_events(speed_type=speed_type, search_parameters=('solar cycle', speed_type, 'acceleration'), search_conditions=('equal', 'greater than', 'greater than'), search_values=(solar_cycle, 0, 0), nan_kwargs=nan_kwargs)
            SS.load_events(speed_type=speed_type, search_parameters=('solar cycle', speed_type), search_conditions=('equal', 'greater than'), search_values=(solar_cycle, 0), nan_kwargs=nan_kwargs)

Each of the events that are loaded will contain the desired speed data (specified by `distribution_model='lognormal distribution'`) and its natural logarithm (specified by `distribution_model='normal distribution'`). We can specify the binning criteria by which we obtain a histogram of the data. The keyword argument `threshold` corresponds to a bin-threshold; if a bin contains less than `threshold` counts, then the bin is merged with the next bin (next taken to be in the direction of the central peak). These histograms can be used to obtain find the optimized distribution parameters that extremize error functions, which include `'maximum likelihood estimation'`, `'chi square'`, and `'g-test'`. Typically, one selects either Pearson's chi-squared test or the lesser-known [G-Test statistic](https://en.wikipedia.org/wiki/G-test); however, both are shown for the sake of showing functionality. 

    distribution_models = ('lognormal distribution', 'normal distribution')
    SS.load_histograms(distribution_model='lognormal distribution', value=70, criterion='bin width', threshold=5, remove_empty_boundaries=True)
    SS.load_histograms(distribution_model='normal distribution', value=50, criterion='number of bins', threshold=5) #, remove_empty_boundaries=True)
    #
    fit_kwargs = dict()
    fit_kwargs['method'] = 'Nelder-Mead'
    error_metrics = ('maximum likelihood estimation', 'chi square', 'g-test')
    for distribution_model in distribution_models:
        SS.load_optimizations(error_metrics, distribution_model, parameter_guess=None, scale='local', reduce_statistic=True, minimizer_kwargs=None, **fit_kwargs)
        for error_metric in error_metrics:
            if error_metric in ('maximum likelihood estimation',):
                kws = dict(xfrac=0.5, xn=11, xspace_by='number', yfrac=0.5, yn=11, yspace_by='number')
            else:
                x = np.arange(2, 9.5, 0.1)
                y = np.arange(0.3, 1.2, 0.01)
                kws = dict(x=x, y=y)
            SS.load_error_dimensionalizations(error_metric, distribution_model, **kws)

Make sure to specify the layout when viewing a figure. 

    layout = 'vertical'

When viewing the distribution, one has the option of viewing any combination of the following: 
a)  histogram (as bars and/or steps)
b)  the fit via optimized parameters
c)  the confidence interval (obtained from the optimized parameters)
d)  vertical lines denoting distribution statistics (such as mean)
e)  highlighting the tail (must define value and condition for extreme event)
f)  arrow pointing to the tail (must define value and condition for extreme event)

The image shown below shows the fit via optimized parameters (b), the confidence interval (c), and vertical lines (d).

    for distribution_model in distribution_models:
        SS.view_distribution(distribution_model=distribution_model, layout=layout, show_histogram=True, as_bars=True, as_steps=True, save=True, figsize=(7,7), sharex=True, sharey=True)
        SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, save=True, figsize=(7,7), sharex=True, sharey=True)
        SS.view_distribution(distribution_model=distribution_model, layout=layout, confidence_metric='maximum likelihood estimation', confidence_color='darkgreen', save=True, figsize=(7,7), sharex=True, sharey=True)
        SS.view_distribution(distribution_model=distribution_model, layout=layout, confidence_metric='maximum likelihood estimation', confidence_color='darkgreen', stats_metric='maximum likelihood estimation', save=True, figsize=(7,7), sharex=True, sharey=True)
        SS.view_distribution(distribution_model=distribution_model, layout=layout, show_histogram=True, as_steps=True, error_metrics=error_metrics, save=True, figsize=(7,7), sharex=True, sharey=True)
        SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, confidence_metric='maximum likelihood estimation', confidence_color='darkgreen', save=True, figsize=(7,7), sharex=True, sharey=True)
        SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, point_to_tail=True, extreme_value=800, extreme_condition='greater than or equal', save=True, figsize=(7,7), sharex=True, sharey=True)
        SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, tail_metric='maximum likelihood estimation', extreme_value=800, extreme_condition='greater than or equal', save=True, figsize=(7,7), sharex=True, sharey=True)

![Lognormal Speed Distribution with Confidence-Interval and Statistics](https://i.imgur.com/rDuz6I2.png)  

One can also view a log-log plot of the distribution tail.

    SS.view_histogram_tail(layout=layout, extreme_value=800, extreme_condition='greater than or equal', as_steps=True, as_bars=True, save=True, figsize=(7,7), sharex=True, sharey=True)

One can also view the error-space surrounding the extremum of the error function (as evaluated at the optimized parameters). In 2-dimensions, this will appear as a contour mapping. In 3-dimensions, one can view error-surfaces and/or error contours. The image below depicts the dim-2 contour map of the error-space that corresponds to the G-Test statistic.

    SS.view_dim2_errorspace(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, include_colorbar=True, save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_dim3_errorspace(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, show_contours=True, show_surface=True, azim=30, elev=30, include_colorbar=True, save=True, figsize=(7,7))

![2-D Contour Map of G-Test Error-Space](https://i.imgur.com/cBSZZx3.png)  

**Example - Time-Series Analysis**





##
