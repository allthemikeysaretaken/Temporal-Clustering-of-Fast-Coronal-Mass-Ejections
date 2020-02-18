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

    ## RUN & VIEW FREQUENCY-SERIES
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


##
