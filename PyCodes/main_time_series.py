from time_series_methods import *

## SPECIFY DIRECTORIES EXPLICITLY
dir_home = '/Users/.../'
dir_data = '{}Data/'.format(dir_home)
dir_figs = '{}Figures/'.format(dir_home)
dir_time_series = '{}TimeSeries/'.format(dir_figs)
##

## LOAD DATABASE
DB = DataBase(directory=dir_data)
DB.load_data(read_soho_lasco=True)
##

## SET RANDOM STATE SEED
np.random.seed(327)
##

## RUN & VIEW TIME-SERIES
timestep = 'hour'
extreme_values = (800, 1000)
all_extreme_values = np.arange(500, 2025, 25).astype(int)
solar_cycles = (23, 24)
speed_types = ('second order initial speed', 'linear speed')
nan_policy = 'replace' # 'discard'
nan_repl = 1
nan_kwargs = dict(policy=nan_policy, repl=nan_repl)
TS = TimeSeries(DB, timestep=timestep, event_type='cme', ref_parameter='speed', directory=dir_time_series)
for solar_cycle in solar_cycles:
    for speed_type in speed_types:
        nan_kwargs['parameter'] = speed_type
        TS.load_events(extreme_parameter=speed_type, extrema='maximum', search_parameters=('solar cycle', 'cycle type'), search_conditions=('equal', 'equal'), search_values=(solar_cycle, 'high-activity'), nan_kwargs=nan_kwargs)
#
TS.load_inter_exceedances(all_extreme_values)
TS.load_inter_exceedance_distributions(extreme_values, bin_width=30) # bin_width=15
TS.view_inter_exceedance_histograms(extreme_values, layout='double-vertical', yspace=50, figsize=(7,7), sharex=True, sharey=True, save=True) # yspace=25
#
jndices = np.arange(4, 14).astype(int)
fit_kwargs = dict()
fit_kwargs['method'] = 'Nelder-Mead'
TS.load_unbiased_estimators(jndices, alpha_model='normal distribution', intercept_model='normal distribution', exclude_original_sample=False, nresamples=1000, nshuffles=3, nan_repl=nan_repl, **fit_kwargs)
statistics = 'mean' # ('skew', 'kurtosis') # None
for attribute in ('alpha', 'intercept', 'theta'):
    TS.view_unbiased_estimators_histograms(attribute, layout='overlay', counts_type='normalized', statistics=statistics, figsize=(7,7), save=True)
    TS.view_unbiased_estimators_histograms(attribute, layout='double-vertical', counts_type='observed', statistics=statistics, figsize=(7,7), sharex=True, sharey=True, save=True)
#
TS.view_max_spectrum(layout='vertical', show_errors=True, facecolors=('gray', 'darkorange', 'steelblue'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='double-vertical', show_errors=True, facecolors=('gray', 'darkorange', 'steelblue'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='vertical', show_fit=True, facecolors=('darkorange', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='double-vertical', show_fit=True, facecolors=('darkorange', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='vertical', show_errors=True, show_fit=True, facecolors=('gray', 'darkorange', 'mediumpurple', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='double-vertical', show_errors=True, show_fit=True, facecolors=('gray', 'darkorange', 'mediumpurple', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
#
TS.view_max_spectrum(layout='vertical', show_errors=True, show_logscale=True, facecolors=('gray', 'darkorange', 'steelblue'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='double-vertical', show_errors=True, show_logscale=True, facecolors=('gray', 'darkorange', 'steelblue'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='vertical', show_fit=True, show_logscale=True, facecolors=('darkorange', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='double-vertical', show_fit=True, show_logscale=True, facecolors=('darkorange', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='vertical', show_errors=True, show_fit=True, show_logscale=True, facecolors=('gray', 'darkorange', 'mediumpurple', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_max_spectrum(layout='double-vertical', show_errors=True, show_fit=True, show_logscale=True, facecolors=('gray', 'darkorange', 'mediumpurple', 'r'), save=True, figsize=(7,7), sharex=True, sharey=True)
#
TS.view_point_estimators(layout='vertical', save=True, figsize=(7,7), sharex=True, sharey=True)
TS.view_point_estimators(layout='double-vertical', save=True, figsize=(7,7), sharex=True, sharey=True)
#
TS.load_moment_estimators_and_time_thresholds(all_extreme_values, baseline=0.5)
TS.view_moment_estimators_with_time_threshold(save=True, figsize=(7,7), sharex=True) # show_max_extreme_value=True,
#
cluster_search_kwargs = dict()
cluster_search_kwargs['search_parameters'] = 'cluster size'
cluster_search_kwargs['search_conditions'] = 'greater than or equal'
cluster_search_kwargs['search_values'] = 2
TS.load_clusters(extreme_values, bias='threshold') # bias='first-order'
#
TS.view_chronological_clusters(parameter='cluster size', extreme_values=extreme_values, save=True, figsize=(7,7), rotation=15, sharex=True, **cluster_search_kwargs)
TS.view_chronological_clusters(parameter='cluster size', statistic='max', extreme_values=extreme_values, save=True, figsize=(7,7), rotation=15, sharex=True, **cluster_search_kwargs)
TS.view_chronological_clusters(parameter='speed', statistic='max', extreme_values=extreme_values, save=True, figsize=(7,7), rotation=15, sharex=True, **cluster_search_kwargs)
TS.view_chronological_clusters(parameter='speed', statistic='mean', extreme_values=extreme_values, save=True, figsize=(7,7), rotation=15, sharex=True, **cluster_search_kwargs)
#
TS.view_relative_cluster_statistics(extreme_values=extreme_values, save=True, figsize=(7,7), **cluster_search_kwargs)
#
TS.view_cluster_duration_histograms(extreme_values=extreme_values, layout='overlay', show_intra_times=True, show_intra_durations=True, show_inter_durations=True, save=True, figsize=(7,7), **cluster_search_kwargs)
TS.view_cluster_duration_histograms(extreme_values=extreme_values, layout='vertical', show_intra_times=True, show_intra_durations=True, show_inter_durations=True, save=True, figsize=(7,7), sharex=True, sharey=True, **cluster_search_kwargs)
TS.view_cluster_duration_histograms(extreme_values=extreme_values, layout='double-vertical', show_intra_times=True, show_intra_durations=True, show_inter_durations=True, save=True, figsize=(7,7), sharex=True, sharey=True, **cluster_search_kwargs)
TS.view_cluster_duration_histograms(extreme_values=extreme_values, layout='overlay', zoom_in=True, show_intra_times=True, show_intra_durations=True, show_inter_durations=True, save=True, figsize=(7,7), **cluster_search_kwargs)
TS.view_cluster_duration_histograms(extreme_values=extreme_values, layout='vertical', zoom_in=True, show_intra_times=True, show_intra_durations=True, show_inter_durations=True, save=True, figsize=(7,7), sharex=True, sharey=True, **cluster_search_kwargs)
TS.view_cluster_duration_histograms(extreme_values=extreme_values, layout='double-vertical', zoom_in=True, show_intra_times=True, show_intra_durations=True, show_inter_durations=True, save=True, figsize=(7,7), sharex=True, sharey=True, **cluster_search_kwargs)
#
for extreme_value in extreme_values:
    TS.view_extreme_event_analysis_table(extreme_value=extreme_value, save=True)
    TS.view_cluster_analysis_table(extreme_value=extreme_value, save=True, **cluster_search_kwargs)
    TS.view_clusters_by_size_table(extreme_value=extreme_value, save=True, **cluster_search_kwargs)
#
TS.view_frechet_distribution(save=True, figsize=(7,7), sharex=True)
TS.view_cluster_example(save=True, figsize=(7,7))
##
