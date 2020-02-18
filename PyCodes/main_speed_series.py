from speed_series_methods import *

## SPECIFY DIRECTORIES EXPLICITLY
dir_home = '/Users/.../'
dir_data = '{}Data/'.format(dir_home)
dir_figs = '{}Figures/'.format(dir_home)
dir_speed_series = '{}SpeedSeries/'.format(dir_figs)
##

## LOAD DATABASE
DB = DataBase(directory=dir_data)
DB.load_data(read_soho_lasco=True)
##

## RUN & VIEW SPEED-SERIES
solar_cycles = (23, 24)
speed_types = ('second order initial speed', 'linear speed')
nan_policy = 'discard' # 'replace'
nan_repl = 0
nan_kwargs = dict(policy=nan_policy, repl=nan_repl)
SS = SpeedSeries(DB, timestep='day', directory=dir_speed_series)
for solar_cycle in solar_cycles:
    for speed_type in speed_types:
        nan_kwargs['parameter'] = speed_type
        # SS.load_events(speed_type=speed_type, search_parameters=('solar cycle', speed_type, 'acceleration'), search_conditions=('equal', 'greater than', 'greater than'), search_values=(solar_cycle, 0, 0), nan_kwargs=nan_kwargs)
        SS.load_events(speed_type=speed_type, search_parameters=('solar cycle', speed_type), search_conditions=('equal', 'greater than'), search_values=(solar_cycle, 0), nan_kwargs=nan_kwargs)
#
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
#
layout = 'vertical'
for distribution_model in distribution_models:
    SS.view_distribution(distribution_model=distribution_model, layout=layout, show_histogram=True, as_bars=True, as_steps=True, save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_distribution(distribution_model=distribution_model, layout=layout, confidence_metric='maximum likelihood estimation', confidence_color='darkgreen', save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_distribution(distribution_model=distribution_model, layout=layout, confidence_metric='maximum likelihood estimation', confidence_color='darkgreen', stats_metric='maximum likelihood estimation', save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_distribution(distribution_model=distribution_model, layout=layout, show_histogram=True, as_steps=True, error_metrics=error_metrics, save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, confidence_metric='maximum likelihood estimation', confidence_color='darkgreen', save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, point_to_tail=True, extreme_value=800, extreme_condition='greater than or equal', save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_distribution(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, tail_metric='maximum likelihood estimation', extreme_value=800, extreme_condition='greater than or equal', save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_dim2_errorspace(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, include_colorbar=True, save=True, figsize=(7,7), sharex=True, sharey=True)
    SS.view_dim3_errorspace(distribution_model=distribution_model, layout=layout, error_metrics=error_metrics, show_contours=True, show_surface=True, azim=30, elev=30, include_colorbar=True, save=True, figsize=(7,7))
SS.view_histogram_tail(layout=layout, extreme_value=800, extreme_condition='greater than or equal', as_steps=True, as_bars=True, save=True, figsize=(7,7), sharex=True, sharey=True)
##
