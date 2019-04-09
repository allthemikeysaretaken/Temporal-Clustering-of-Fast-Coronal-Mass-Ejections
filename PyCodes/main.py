from data_retrieval import *
from timeseries_methods import *
from viewer import *

################################################################################
### *-  SPECIFY FILEPATHS
###     LINE 26     --      PATH TO READ SOHO LASCO CME DATA
###     LINE 55     --      PATH TO READ SILSO SUNSPOT DATA
###     LINE 58     --      DIRECTORY TO SAVE FIGURES
###     LINE 59     --      DIRECTORY TO SAVE TEXT FILES
### *-  SPECIFY DATA SUBTYPES AND ANALYSIS METHODS
###     LINE 31     --      TIME-INTERVAL, SPEED TYPES, AND NAN DATA
###     LINE 36     --      SELECT SUBSETS OF DATA
###     LINE 41     --      SPECIFY SPEED THRESHOLD
###     LINE 45     --      SPECIFY TIMESERIES ANALYSIS METHODS
###     LINE 50     --      SPECIFY KWARGS FOR EACH TIMESERIES ANALYSIS METHOD
### *-  SPECIFY DATA CONDITIONS FOR CLUSTERS (TEXT FILE)
###     LINE 62     --      SPECIFY CONDITIONS TO VIEW CLUSTER DATA
###     LINE 63     --      SPECIFY CONDITIONS TO VIEW CLUSTER STATISTICS
################################################################################

## REPRODUCE RANDOMIZED RESULTS
np.random.seed(327)

## READ DATA FROM FILE & AUTOCORRECT NON-NUMERICAL ENTRIES
path = '.../Data/univ_all.txt'
D = Data(path)
D.initialize()

## ORGANIZE SUBSETS OF DATA
kwargs = {'cycle_type' : 'high-activity', 'nan_policy' : 'replace', 'nan_repl' : 1}
sc23_vinit = D.dispatch('second order initial speed', solar_cycle=23, **kwargs)
sc23_vlinr = D.dispatch('linear speed', solar_cycle=23, **kwargs)
sc24_vinit = D.dispatch('second order initial speed', solar_cycle=24, **kwargs)
sc24_vlinr = D.dispatch('linear speed', solar_cycle=24, **kwargs)
sc_data = (sc23_vinit, sc24_vinit)
# sc_data = (sc23_vinit,)
# sc_data = (sc23_vinit, sc23_vlinr, sc24_vinit, sc24_vlinr)

## SPECIFY SPEED THRESHOLD USED TO DISTINGUISH EXTREME EVENTS
vthreshold = 800
# vthreshold = 1000

## SELECT TIMESERIES ANALYSIS METHODS
keys = ['inter-exceedance decay', 'unbiased estimators', 'extremal index', 'temporal clustering', 'relative cluster statistics']
keys += ['cluster thresholds', 'sunspot correlation', 'speed distribution fit']
keys = tuple(keys)

## SPECIFY PARTICULARS FOR EACH TIMESERIES ANALYSIS METHOD
kwargs = {}
kwargs['unbiased_estimator_kw'] = {'nresamples' : 1000}
kwargs['extremal_index_kw'] = {'n' : 31, 'edges_by' : 'number'}
kwargs['cluster_statistic_kw'] = {'skip_lone_ejecta' : True}
kwargs['speed_distribution_kw'] = {'n' : 20, 'edges_by' : 'width', 'vmin' : 20}
kwargs['sunspot_kw'] = {'path' : '.../Data/SN_d_tot_V2.0.txt'}

## SPECIFY DIRECTORIES TO SAVE FIGURES AND RESULTS
saveloc = '.../Figures/' ## SAVE PLOTS
textloc = '.../Results/' ## SAVE TEXT FILE

## EXPLICIT FUNCTION INPUTS TO SAVE TIMESERIES ANALYSIS DATA & RESULTS
cluster_kwargs = dict(search_parameters='cluster size', search_conditions='exact match', search_values=3, fkeys=None, load_keys=('ith cluster', 'cluster size'), apply_to='all', view_keys=None, use_abs=False, use_halo=None)
stats_kwargs = dict(search_parameters='cluster size', search_conditions='exact match', search_values=3, fkeys=None, load_keys=('ith cluster', 'cluster size'), apply_to='all', view_keys=None, use_abs=False, use_halo=None, skip_lone_ejecta=True)

def main(sc_data, keys, vthreshold, saveloc=None, textloc=None, cluster_kwargs=None, stats_kwargs=None, save_sample_figures=False, **kwargs):
    """
    This function first runs the timeseries analysis methods, which are
    specified by keys and for which vthreshold specifies the speed threshold
    that distinguishes extreme events.

    Then, a text file of the data and/or results from the timeseries
    analysis is saved into the saveloc directory. For each TimeSeries
    instance, the data identifiers and results from the analysis will be
    saved into a text file named 'analysis_results__...'. If cluster_kwargs
    are specified, then the desired clusters are saved into a text file
    named 'cluster_data__...'. If stats_kwargs are specified, then the
    relative statistics and mean duration that correspond to the desired
    clusters are saved into a text file named 'cluster_stats__...'.

    Lastly, the function will save figures that show various data, results,
    and parameters from the analysis.

    *-
    sc_data             :   type <dict>
    *-

        Contains keys:
            'padded data'                    :   type <dict>
            'unique data'                    :   type <dict>
            'original data'                  :   type <dict>

        Each key unlocks a subdictionary, each of which unlocks
        an array. These subkeys and array dtypes are listed below:
            'speed'                         :   type <int / float>
            'datetime'                      :   type <str>
            'dt object'                     :   type <datetime object>
            'year'                          :   type <int>
            'month'                         :   type <int>
            'day'                           :   type <int>
            'hour'                          :   type <int>
            'minute'                        :   type <int>
            'second'                        :   type <int>
            'mass'                          :   type <int / float>
            'kinetic energy'                :   type <int / float>
            'acceleration'                  :   type <int / float>
            'mean position angle'           :   type <int / float>
            'central position angle'        :   type <int / float>
            'halo'                          :   type <bool>

    *-
    keys                :   type <tuple / list / array>
    *-

        Each key corresponds to a timeseries method. After running the
        timeseries methods (by calling TimeSeries.run(...), the corresponding
        results will be available to the TimeSeries instance as an attribute.
        These keys and instance attributes are listed below:

            'inter-exceedance decay'
                unlocks TimeSeries.InterExceedanceDecay

            'unbiased estimators'
                unlocks TimeSeries.UnbiasedEstimators

            'extremal index'
                unlocks TimeSeries.ExtremalIndex

            'temporal clustering'
                unlocks TimeSeries.TemporalClustering

            'cluster thresholds'
                unlocks TimeSeries.cluster_thresholds_by_speed

            'relative cluster statistics'
                unlocks TimeSeries.relative_cluster_statistics

            'sunspot correlation'
                unlocks TimeSeries.SunSpotCorrelator

            'speed distribution fit'
                unlocks TimeSeries.SpeedDistribution

            By default, TimeSeries.InterExceedance is unlocked.

        These methods are detailed in TimeSeriesAnalysis (type <cls>),
        which is located in 'timeseries_methods.py'.

    *-
    vthreshold          :   type <int / float>
    *-

        The value of vthreshold is used in the condition v > vthreshold
        to distinguish extreme events in the calculation of the
        inter-exceedance times.

    *-
    saveloc             :   type <str> or None
    *-

        If specified, saveloc will specify the directory in which figures
        are saved. Otherwise, the figures will be shown (instead of saved).

    *-
    textloc             :   type <str> or None
    *-

        If specified, textloc will specify the directory in which text
        files are saved; these text files contain sample data and/or
        analysis results.

    *-
    cluster_kwargs      :   type <dict> or None
    *-

        If specified, the following key-value pairings should be included:
            'search_parameters'             :   type <tuple / list / array / str>
            'search_conditions'             :   type <tuple / list / array / str>
            'search_values'                 :   type <tuple / list / array / str>

        Optional key-value pairings include:
            fkeys = None
            load_keys = ('ith cluster', 'cluster size')
            apply_to = 'all'
            view_keys = None
            use_abs = False
            use_halo = None

    *-
    stats_kwargs        :   type <dict> or None
    *-

        If specified, the following key-value pairings should be included:
            'search_parameters'             :   type <tuple / list / array / str>
            'search_conditions'             :   type <tuple / list / array / str>
            'search_values'                 :   type <tuple / list / array / str>

        Optional key-value pairings include:
            fkeys = None
            load_keys = ('ith cluster', 'cluster size')
            apply_to = 'all'
            view_keys = None
            use_abs = False
            use_halo = None
            include_lone_ejecta = True

    *-
    save_sample_figures :   type <bool>
    *-

        if True:
            -   saves an example of intra- and inter- cluster durations
            -   saves the pdf and cdf of the Frechet distribution

        if False:
            skip

    *-
    kwargs              :   type <dict> containing type <dict>
    *-

        >>  unbiased_estimator_kw   :   type <dict>

                Contains key-value pairing:
                    'nresamples'                    :   type <int>

        >>  extremal_index_kw       :   type <dict>

                Contains key-value pairing:
                    'n'                             :   type <int / float> or type <array>
                    'edges_by'                      :   'width' or 'number', or 'custom'
                    'bias'                          :   'left' or 'right'

        >>  cluster_statistic_kw    :   type <dict>

                Contains key-value pairing:
                    'skip_lone_ejecta'              :   type <bool>

        >>  speed_distribution_kw   :   type <dict>

                Contains key-value pairing:
                    'n'                             :   type <int / float>, or type <array>
                    'edges_by'                      :   'width' or 'number', or 'custom'

        >>  sunspot_kw              :   type <dict>

                Contains key-value pairing:
                    'path'                          :    type <str>
    """
    ## RUN ANALYSIS
    MultipleTimeSeries = []
    for data in sc_data:
        TimeSeries = TimeSeriesAnalysis(data, keys)
        TimeSeries.run(vthreshold, **kwargs)
        MultipleTimeSeries.append(TimeSeries)
    ## SAVE RESULTS
    if textloc is not None:
        Viewer = FigureVisualizer(MultipleTimeSeries, saveloc=textloc)
        Viewer.save_analysis_results()
        if cluster_kwargs is not None:
            if isinstance(cluster_kwargs, dict):
                Viewer.save_clusters(**cluster_kwargs)
        if stats_kwargs is not None:
            if isinstance(stats_kwargs, dict):
                Viewer.save_cluster_statistics(**stats_kwargs)
    ## VIEW FIGURES
    if saveloc is not None:
        if 'inter-exceedance decay' in keys:
            Viewer = FigureVisualizer(MultipleTimeSeries, 'inter-exceedance distribution decay', saveloc)
            Viewer.dispatch()
        if 'unbiased estimators' in keys:
            Viewer = FigureVisualizer(MultipleTimeSeries, 'alpha-hat histogram', saveloc)
            Viewer.dispatch()
            Viewer = FigureVisualizer(MultipleTimeSeries, 'max spectrum power-law', saveloc)
            Viewer.dispatch()
            Viewer = FigureVisualizer(MultipleTimeSeries, 'unbiased estimator errors', saveloc)
            Viewer.dispatch()
        if 'extremal index' in keys:
            if 'unbiased estimators' not in keys:
                raise ValueError("'unbiased estimators' should be included in keys to run 'extremal index' method")
            Viewer = FigureVisualizer(MultipleTimeSeries, 'point estimators', saveloc)
            Viewer.dispatch()
            Viewer = FigureVisualizer(MultipleTimeSeries, 'extremal index histogram', saveloc)
            Viewer.dispatch()
        if 'cluster thresholds' in keys:
            Viewer = FigureVisualizer(MultipleTimeSeries, 'moment estimates', saveloc)
            Viewer.dispatch()
        if 'sunspot correlation' in keys:
            Viewer = FigureVisualizer(MultipleTimeSeries, 'sunspot correlations', saveloc, **dict(figsize=(15, 6)))
            Viewer.dispatch()
        if 'temporal clustering' in keys:
            Viewer = FigureVisualizer(MultipleTimeSeries, 'chronological cluster size', saveloc, **dict(figsize=(12, 5)))
            Viewer.dispatch()
            Viewer = FigureVisualizer(MultipleTimeSeries, 'durational histograms', saveloc)
            Viewer.dispatch() # **dict(cluster_size=2, search_conditions='greater than')
        if 'relative cluster statistics' in keys:
            if 'temporal clustering' not in keys:
                raise ValueError("'temporal clustering' should be included in keys to run 'relative cluster statistics' method")
            Viewer = FigureVisualizer(MultipleTimeSeries, 'relative cluster statistics', saveloc)
            Viewer.dispatch()
        if 'speed distribution fit' in keys:
            Viewer = FigureVisualizer(MultipleTimeSeries, 'speed tail', saveloc)
            Viewer.dispatch()
            Viewer = FigureVisualizer(MultipleTimeSeries, 'speed distribution', saveloc)
            Viewer.dispatch()
            Viewer = FigureVisualizer(MultipleTimeSeries, 'dim3 speed distribution', saveloc, **dict(figsize=(13, 5)))
            Viewer.dispatch(**dict(err='minimum gtest estimation')) # **dict(err='maximum likelihood estimation')
        ## VIEW SUPPLEMENTARY FIGURES
        if save_sample_figures is True:
            Viewer = FigureVisualizer(MultipleTimeSeries, 'frechet distribution', saveloc)
            Viewer.dispatch()
            Viewer = FigureVisualizer(MultipleTimeSeries, 'illustrative cluster example', saveloc)
            Viewer.dispatch()

if __name__ == '__main__':
    main(sc_data, keys, vthreshold, saveloc, textloc, cluster_kwargs, stats_kwargs, save_sample_figures=True, **kwargs)
