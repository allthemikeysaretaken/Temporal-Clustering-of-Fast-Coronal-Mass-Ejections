from data_processing import *
from analysis_methods import *

class InterExceedance():

    def __init__(self, extreme_events):
        """ """
        self.extreme_events = extreme_events
        self._times = None
        self._distribution = None

    def initialize_times(self):
        self._times = np.diff(self.extreme_events['elapsed hour'])

    @property
    def times(self):
        return self._times

    def initialize_distribution(self):
        IED = InterExceedanceDistribution(self.times)
        IED.initialize_inverse_sample_histogram()
        IED.initialize_time_delta_histogram()
        self._distribution = IED

    @property
    def distribution(self):
        return self._distribution

class DataAnalysis():

    def __init__(self, cme_data, identifier=None):
        self.cme_data = cme_data
        self.identifier = identifier
        self._unbiased_estimators = None
        self._speed_distribution = None
        self.inter_exceedances = {}
        self.cluster_parameters = {}
        self.temporal_clustering = {}
        self._inter_exceedance_speed_thresholds = []

    @property
    def inter_exceedance_speed_thresholds(self):
        return np.unique(self._inter_exceedance_speed_thresholds)

    def initialize_speed_distribution(self, edges=None, n=20, w=None, bin_threshold=5, bias='left'):
        """ """
        SD = SpeedDistribution(self.cme_data['original']['speed'])
        SD.initialize_histograms(edges=edges, n=n, w=w, bin_threshold=bin_threshold, bias=bias)
        # SD.initialize_optimizations()
        self._speed_distribution = SD

    @property
    def speed_distribution(self):
        return self._speed_distribution

    def initialize_unbiased_estimators(self, nresamples=100):
        """ """
        UB = UnbiasedEstimators(self.cme_data['padded']['speed'])
        UB.store_max_spectra(js_subinterval=np.arange(4, 13), nresamples=nresamples)
        alpha_edges = np.arange(3, 4.01, 0.02) # np.arange(2.8, 4.11, 0.02)
        UB.initialize_alpha_hat_histogram(edges=alpha_edges)
        UB.initialize_thetas()
        UB.initialize_point_estimators()
        theta_edges = np.arange(0, 1.1, 0.05) # np.arange(0, 1.1, 0.1)
        UB.initialize_theta_hat_histogram(edges=theta_edges)
        UB.initialize_speed_thresholds()
        self._unbiased_estimators = UB

    @property
    def unbiased_estimators(self):
        return self._unbiased_estimators

    def initialize_inter_exceedances(self, speed_thresholds, search_conditions='greater than or equal'):
        """ """
        if isinstance(speed_thresholds, (int, float)):
            speed_thresholds = (speed_thresholds,)
        Searcher = SearchEvents(self.cme_data['padded'], string_keys=('',))
        search_kwargs = {'search_parameters' : 'speed', 'search_conditions' : search_conditions}
        for speed in speed_thresholds:
            extreme_events = Searcher.search(search_values=speed, **search_kwargs)
            IE = InterExceedance(extreme_events)
            IE.initialize_times()
            self.inter_exceedances[speed] = IE
            self._inter_exceedance_speed_thresholds.append(speed)

    def initialize_inter_exceedance_distribution(self, speed_threshold):
        """ """
        available_speed_thresholds = list(self.inter_exceedances.keys())
        if speed_threshold not in available_speed_thresholds:
            raise ValueError("invalid speed threshold: {}; inter_exceedances calculated for speeds in {}".format(speed_threshold, available_speed_thresholds))
        IE = self.inter_exceedances[speed_threshold]
        IE.initialize_distribution()
        self.inter_exceedances[speed_threshold] = IE

    def initialize_cluster_parameters(self, speed_threshold):
        """ """
        available_speed_thresholds = list(self.inter_exceedances.keys())
        if speed_threshold not in available_speed_thresholds:
            raise ValueError("invalid speed threshold: {}; inter_exceedances calculated for speeds in {}".format(speed_threshold, available_speed_thresholds))
        IE = self.inter_exceedances[speed_threshold]
        CP = ClusterParametrization(IE.times)
        CP.initialize_moment_estimators(baseline=0.5)
        CP.initialize_time_thresholds(baseline=None)
        self.cluster_parameters[speed_threshold] = CP

    def initialize_temporal_clustering(self, speed_threshold, bias='threshold', density_kwargs=None, agglomerative_kwargs=None): #, **search_kwargs):
        """ """
        available_speed_thresholds = list(self.cluster_parameters.keys())
        if speed_threshold not in available_speed_thresholds:
            raise ValueError("invalid speed threshold: {}; cluster_parameters calculated for speeds in {}".format(speed_threshold, available_speed_thresholds))
        IE = self.inter_exceedances[speed_threshold]
        CP = self.cluster_parameters[speed_threshold]
        CP.initialize_bias(bias)
        if density_kwargs is not None:
            CP.initialize_density_optimization(CP.bias, **density_kwargs) # minimum_npoints = 2
        if agglomerative_kwargs is not None:
            CP.initialize_agglomerative_optimization(**agglomerative_kwargs)
        self.cluster_parameters[speed_threshold] = CP
        clusters = CP.get_clusters(IE.extreme_events, bias)
        Searcher = SearchClusters(clusters)
        Searcher.add_scale_parameters()
        self.temporal_clustering[speed_threshold] = Searcher

    def run_analysis(self, preset_speed_thresholds, inter_exceedance_speed_thresholds=None, nresamples=None, speed_condition='greater than or equal', cluster_bias='threshold', density_kwargs=None, agglomerative_kwargs=None):
        """ """
        if inter_exceedance_speed_thresholds is None:
            inter_exceedance_speed_thresholds = np.array(preset_speed_thresholds).copy()
        if not isinstance(preset_speed_thresholds, (tuple, list, np.ndarray)):
            raise ValueError("preset_speed_thresholds = type <tuple / list / array>; not <{}>".format(type(preset_speed_thresholds)))
        if not isinstance(inter_exceedance_speed_thresholds, (tuple, list, np.ndarray)):
            raise ValueError("inter_exceedance_speed_thresholds = type <tuple / list / array>; not <{}>".format(type(inter_exceedance_speed_thresholds)))
        if nresamples is not None:
            self.initialize_unbiased_estimators(nresamples=nresamples)
        self.initialize_inter_exceedances(speed_thresholds=inter_exceedance_speed_thresholds, search_conditions=speed_condition)
        for speed_threshold in preset_speed_thresholds:
            self.initialize_inter_exceedance_distribution(speed_threshold=speed_threshold)
        for speed_threshold in inter_exceedance_speed_thresholds:
            self.initialize_cluster_parameters(speed_threshold=speed_threshold)
        for speed_threshold in preset_speed_thresholds:
            self.initialize_temporal_clustering(speed_threshold=speed_threshold, bias=cluster_bias, density_kwargs=density_kwargs, agglomerative_kwargs=agglomerative_kwargs)









##
