from analysis_viewer import *

## INITIALIZE FILE PATHS & DIRECTORIES
path_cme = '.../Data/CMEs/univ_all.txt'
path_sunspots = '.../Data/Sunspots/SN_d_tot_V2.0.txt'
dir_fits = '.../Data/ARs/'
saveloc = '.../Results/' # None

## INITIALIZE RANDOM STATE (reproduce randomized results)
np.random.seed(327)

## INITIALIZE DATA PROCESSOR
DP = DataProcessing()
DP.store_data(path_soho_lasco=path_cme, path_silso=path_sunspots)

def main(compare_sunspots=False, show_speed_distribution=False, run_analysis=False, create_timelapse=False, view_supplements=False):
    """ """
    if ((compare_sunspots == True) or (show_speed_distribution == True)):
        ## INITIALIZE CME DATA
        sc23_cmes = DP.dispatch_cme_timeseries(speed_type='second order initial speed', solar_cycle=23, subcycle=None, pad_value=None)
        sc24_cmes = DP.dispatch_cme_timeseries(speed_type='second order initial speed', solar_cycle=24, subcycle=None, pad_value=None)
        # sc23_cmes = DP.dispatch_cme_timeseries(speed_type='linear speed', solar_cycle=23, subcycle=None, pad_value=None)
        # sc24_cmes = DP.dispatch_cme_timeseries(speed_type='linear speed', solar_cycle=24, subcycle=None, pad_value=None)
        if compare_sunspots == True:
            ## INITIALIZE SUNSPOT DATA
            sc23_sunspots = DP.dispatch_sunspot_timeseries(solar_cycle=23, subcycle=None)
            sc24_sunspots = DP.dispatch_sunspot_timeseries(solar_cycle=24, subcycle=None)
            ## VIEW CME-SUNSPOT FREQUENCY COMPARISON
            CSSV = CorrelatedSunSpotViewer(saveloc=saveloc)
            CSSV.view_cme_sunspot_frequency_comparison((sc23_sunspots, sc24_sunspots), (sc23_cmes, sc24_cmes), denote_cycles=(23, 24), add_legend=True, desired_timescale='daily', figsize=(7, 6))
            CSSV.view_cme_sunspot_frequency_comparison((sc23_sunspots, sc24_sunspots), (sc23_cmes, sc24_cmes), denote_cycles=(23, 24), add_legend=True, desired_timescale='monthly', figsize=(7, 6))
            CSSV.view_cme_sunspot_frequency_comparison((sc23_sunspots, sc24_sunspots), (sc23_cmes, sc24_cmes), denote_cycles=(23, 24), add_legend=True, desired_timescale='yearly', figsize=(7, 6))
        if show_speed_distribution == True:
            ## INITIALIZE OPTIMIZATION ESTIMATORS
            bin_estimator = 'reduced gtest' # 'reduced chi square', 'gtest', 'chi square'
            alt_estimator = 'maximum likelihood'
            # sc_data = [sc23_cmes]
            # sc_data = [sc24_cmes]
            sc_data = (sc23_cmes, sc24_cmes)
            SDV = SpeedDistributionViewer(sc_data, saveloc=saveloc)
            SDV.initialize_optimizations(w=70, bin_estimator=bin_estimator, alt_estimator=alt_estimator, search_parameters='speed', search_conditions='greater than', search_values=20)
            ## VIEW HISTOGRAM OF CME SPEED DISTRIBUTION TAIL
            SDV.view_tail_histogram(textsize=7)
            ## VIEW COMPARISON OF ESTIMATED FITS
            SDV.view_fit_comparisons(textsize=7) # highlight_tail=True,  point_to_tail=True
            for estimator in (bin_estimator, alt_estimator):
                SDV.view_empirical_fits(estimator=estimator, highlight_tail=True, textsize=7) # point_to_tail=True, show_confidence_interval=True, show_statistics=True
                SDV.view_probability_density_fits(estimator=estimator, point_to_tail=True, show_confidence_interval=True, show_statistics=True, textsize=7) # highlight_tail=True
                SDV.view_errorspace(estimator, show_surface=True, textsize=7, elev='auto', azim='auto', figsize=(7, 6), antialiased=True, shade=True)
                SDV.view_errorspace(estimator, show_contours=True, textsize=7, ticksize=6)
                SDV.view_errorspace(estimator, show_surface=True, show_contours=True, textsize=7, elev='auto', azim='auto', figsize=(7, 6), antialiased=True, shade=True)
    if run_analysis == True:
        ## INITIALIZE CME DATA
        sc23_v20r = DP.dispatch_cme_timeseries(speed_type='second order initial speed', solar_cycle=23, subcycle='high-activity', pad_value=1)
        sc23_vlin = DP.dispatch_cme_timeseries(speed_type='linear speed', solar_cycle=23, subcycle='high-activity', pad_value=1)
        sc24_v20r = DP.dispatch_cme_timeseries(speed_type='second order initial speed', solar_cycle=24, subcycle=None, pad_value=1)
        sc24_vlin = DP.dispatch_cme_timeseries(speed_type='linear speed', solar_cycle=24, subcycle=None, pad_value=1)
        ## INITIALIZE ANALYSIS SPECS
        inter_exceedance_speed_thresholds = np.arange(500, 2001, 10) # np.arange(500, 1501, 10)
        preset_speed_thresholds = (800, 1000)
        nresamples = 1000
        # density_kwargs = None
        # agglomerative_kwargs = None
        density_kwargs = {'method' : 'dbscan', 'distance_metric' : 'manhattan', 'minimum_npoints' : 2}
        agglomerative_kwargs = {'method' : 'single', 'distance_metric' : 'euclidean'}
        ## COLLECT ANALYSIS RESULTS
        analysis_results = []
        # for sc_data in (sc23_v20r, ):
        # for sc_data in (sc24_v20r, ):
        # for sc_data in (sc23_v20r, sc23_vlin):
        # for sc_data in (sc23_v20r, sc24_v20r):
        for sc_data in (sc23_v20r, sc23_vlin, sc24_v20r, sc24_vlin):
            DA = DataAnalysis(sc_data['data'], identifier=sc_data['identifier'])
            DA.run_analysis(preset_speed_thresholds, inter_exceedance_speed_thresholds, nresamples=nresamples, speed_condition='greater than or equal', cluster_bias='threshold', density_kwargs=density_kwargs, agglomerative_kwargs=agglomerative_kwargs)
            analysis_results.append(DA)
        ## VIEW CME ANALYSIS RESULTS
        AV = AnalysisViewer(analysis_results, saveloc)
        AV.save_clusters(speed_threshold=800, search_parameters='cluster size', search_conditions=('greater than or equal', 'less than'), search_values=(2, 4))
        AV.view_tabular_cluster_size_statistics(speed_threshold=800, search_parameters='cluster size', search_conditions=('greater than or equal', 'less than'), search_values=(2, 7), labelsize=7, textsize=8, titlesize=9)
        AV.view_tabular_parameter_results(speed_threshold=800, search_parameters='cluster size', search_conditions='greater than or equal', search_values=2, labelsize=7, textsize=8, titlesize=9)
        AV.view_inter_exceedance_distribution_histogram(speed_threshold=800)
        AV.view_power_law_errors()
        AV.view_power_law()
        AV.view_alpha_hat_histogram(counts_type='observed frequency') # counts_type='probability density'
        AV.view_theta_hat_histogram(counts_type='observed frequency') # counts_type='probability density'
        AV.view_point_estimators()
        AV.view_cluster_parameterizations()
        AV.view_relative_cluster_statistics(ticksize=6, labelsize=7, textsize=7, titlesize=10, speed_threshold=800, search_parameters='cluster size', search_conditions='greater than or equal', search_values=2)
        AV.view_intra_times_histogram(speed_threshold=800, search_parameters='cluster size', search_conditions='greater than or equal', search_values=2)
        AV.view_intra_durations_histogram(speed_threshold=800, search_parameters='cluster size', search_conditions='greater than or equal', search_values=2)
        AV.view_inter_durations_histogram(speed_threshold=800, search_parameters='cluster size', search_conditions='greater than or equal', search_values=2)
        AV.view_chronological_cluster_sizes(ticksize=6, speed_threshold=800, search_parameters='cluster size', search_conditions='greater than or equal', search_values=2)
        AV.save_optimized_clusters(speed_threshold=800)
        AV.view_dendrogram(speed_threshold=800)
        AV.view_density_optimized_clusters(speed_threshold=800, show_clusters=True, show_noise=True, zoom_clusters=True)
        AV.view_density_optimized_clusters(speed_threshold=800, show_clusters=True, show_noise=True, zoom_clusters=False)
    if create_timelapse == True:
        ## SAVE '.fits' AS TIME-LAPSE OF IMAGES
        FI = FitsInterpreter(dir_fits, saveloc)
        FI.save_images(extension='.png', apply_false_color=True, apply_grayscale=True)
        FI.save_timelapse_animation(fps=5, savename='CarringtonRotations_1990-2030__grayscale', include_grayscale=True)
        FI.save_timelapse_animation(fps=5, savename='CarringtonRotations_1990-2030__falsecolor', include_false_color=True)
    if view_supplements == True:
        ## VIEW SUPPLEMENTS
        SV = SupplementaryViewer(saveloc)
        SV.view_cluster_example()
        SV.view_frechet_distribution()

## run all analyses & save figures
main(compare_sunspots=True, show_speed_distribution=True, run_analysis=True, create_timelapse=True, view_supplements=True)

## convert all figures from '.png' to '.eps'
# FO = FigureOptions(saveloc) ## search saveloc recursively for '.png' files
# FO.convert_to_eps('.png')

### ** THINGS TO ADD:

#   centroid-based clustering
#   identifiers for cluster optimization should include datetime and speed
#   save identifiers as textfile; do not label in plot (overcrowded)

#   search events by difference condition (and cumulative condition)
#   --  np.unique(index_pair[0], index_pair[1] for index_pair in indices)

### correlate specific clusters of ARs to specific clusters of CMEs
##
