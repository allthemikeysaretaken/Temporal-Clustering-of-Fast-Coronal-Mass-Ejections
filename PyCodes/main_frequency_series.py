from frequency_series_methods import *

## SPECIFY DIRECTORIES EXPLICITLY
dir_home = '/Users/.../'
dir_data = '{}Data/'.format(dir_home)
dir_figs = '{}Figures/'.format(dir_home)
dir_freq_series = '{}FrequencySeries/'.format(dir_figs)
##

## LOAD DATABASE
DB = DataBase(directory=dir_data)
DB.load_data(read_soho_lasco=True, read_silso=True, read_hessi=True)
##

## RUN & VIEW FREQUENCY-SERIES
timestep = 'month' # 'day', 'year'
FS = FrequencySeries(DB, timestep, directory=dir_freq_series)
#
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
#
for layout in ('overlay', 'vertical'):
    FS.view_frequency_distributions(layout=layout, save=True, figsize=(7,7), sharex=True)
    FS.view_frequency_distributions(layout=layout, save=True, show_solar_cycle_legends=True, figsize=(7,7), sharex=True)
    FS.view_frequency_distributions(layout=layout, save=True, show_solar_cycle_separations=True, figsize=(7,7), sharex=True)
    FS.view_frequency_distributions(layout=layout, save=True, show_solar_cycle_legends=True, show_solar_cycle_separations=True, figsize=(7,7), sharex=True)
##
