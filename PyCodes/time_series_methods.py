from data_processing import *
from optimization_methods import *
from visual_configuration import *

class TimeSeriesViewer(VisualConfiguration):

    """
    This class is inherited by all classes pertaining to
    timeseries methods. This class inherits methods from
    `VisualConfiguration` to be used as convenience
    functions for plotting.
    """

    def __init__(self, timestep, event_type, ref_parameter, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize):
        """
        timestep:
            type <str>

        event_type:
            type <str>

        ref_parameter:
            type <str>

        directory:
            type <str>

        ticksize:
            type <int / float>

        labelsize:
            type <int / float>

        textsize:
            type <int / float>

        titlesize:
            type <int / float>

        headersize:
            type <int / float>

        cellsize:
            type <int / float>
        """
        super().__init__(directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize)
        self.timestep = timestep
        self.event_type = event_type
        self.ref_parameter = ref_parameter
        self.event_type_label = self.get_event_type_label(event_type)
        self.unit_label = self.available['label'][event_type]['unit'][ref_parameter]

    @staticmethod
    def autocorrect_extreme_values(extreme_values):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None
        """
        if extreme_values is None:
            extreme_values = [None]
        else:
            if isinstance(extreme_values, (int, float)):
                extreme_values = [extreme_values]
            elif not isinstance(extreme_values, (tuple, list, np.ndarray)):
                raise ValueError("invalid type(extreme_values): {}".format(type(extreme_values)))
        return extreme_values

    @staticmethod
    def get_event_type_label(event_type):
        """
        event_type:
            type <dict>
        """
        if event_type == 'cme':
            result = '{}'.format(event_type.upper())
        else:
            result = '{}'.format(event_type.title())
        return result

    @staticmethod
    def get_solar_cycle_label(events):
        """
        events:
            type <dict>
        """
        scycles = []
        for key in list(events['parameter'].keys()):
            if 'solar cycle' in key:
                scycles += events['parameter'][key].tolist()
        scycles = np.unique(scycles)
        if len(scycles) == 1:
            result = 'SC ${}$'.format(scycles[0])
        else:
            if np.all(np.diff(scycles) <= 1):
                result = 'SC ${}$ - ${}$'.format(scycles[0], scycles[-1])
            else:
                args = np.array(['SC ${}$'.format(scycle) for scycle in scycles])
                result = ", ".join(args[:-1].tolist()) + '{}'.format(args[-1])
        return result

    def get_extreme_parameter_label(self, events):
        """
        events:
            type <dict>
        """
        identifiers = events['identifier']
        extreme_parameter = identifiers['extreme parameter']
        ref_key = '{} type'.format(self.ref_parameter)
        label = self.available['label'][self.event_type][ref_key][extreme_parameter]
        return label

    def get_search_label(self, search_kwargs):
        """
        search_kwargs:
            type <dict> or None
        """
        if search_kwargs is None:
            return None
        else:
            _parameters, _conditions, _values = search_kwargs['parameters'], search_kwargs['conditions'], search_kwargs['values']
            result = ''
            f = lambda args : list(args) if isinstance(args, (tuple, list, np.ndarray)) else [args]
            parameters, conditions, values = f(_parameters), f(_conditions), f(_values)
            nps, ncs, nvs = len(parameters), len(conditions), len(values)
            is_same_size = (nps == ncs == nvs)
            if is_same_size == True:
                for idx, (parameter, condition, value) in enumerate(zip(parameters, conditions, values)):
                    sym_condition = self.available['label']['comparison'][condition]
                    if isinstance(value, int):
                        result += '{} ${}$ ${:,}$'.format(parameter.title(), sym_condition, value)
                    else:
                        result += '{} ${}$ ${}$'.format(parameter.title(), sym_condition, value)
                    if parameter == self.ref_parameter:
                        result += '{}'.format(self.unit_label)
                    else:
                        try:
                            unit_label = self.available['label'][self.event_type]['unit'][parameter]
                            result += ' {}'.format(unit_label)
                        except:
                            pass
                    if idx != nps - 1:
                        result += '\n'
            else:
                raise ValueError("{} parameters for {} conditions and {} values".format(nps, ncs, nvs))
            return result


class TemporalClusters():

    def __init__(self, clusters, bias, parameters=None, conditions=None, values=None, apply_to='all', modifiers=None):
        """
        clusters:
            type <dict>

        bias:
            type <str>

        ...
        """
        search_kwargs = dict(parameters=parameters, conditions=conditions, values=values, apply_to=apply_to, modifiers=modifiers)
        S = ClusterSearcher(dict(clusters))
        self.clusters = S.search_clusters(**search_kwargs)[0]
        self.bias = bias
        self.search_kwargs = search_kwargs
        self._relative_statistics = {}
        self._intra_times = {}
        self._intra_durations = {}
        self._inter_durations = {}

    @property
    def relative_statistics(self):
        return self._relative_statistics

    @property
    def intra_times(self):
        return self._intra_times


    @property
    def intra_durations(self):
        return self._intra_durations

    @property
    def inter_durations(self):
        return self._inter_durations

    def store_relative_statistics(self):
        initial_cluster_sizes = np.array([cluster.size for cluster in self.clusters['elapsed']])
        cluster_sizes, size_counts = np.unique(initial_cluster_sizes, return_counts=True)
        nevents = cluster_sizes * size_counts
        rel_prob = nevents / np.sum(nevents)
        result = dict()
        result['cluster size'] = cluster_sizes
        result['number clusters'] = size_counts
        result['number events'] = nevents
        result['relative probability'] = rel_prob
        self._relative_statistics.update(result)

    def store_cluster_times_and_durations(self):
        arr = self.clusters['elapsed']
        intra_times, intra_durations, inter_durations = [], [], []
        first_intra_times = np.diff(arr[0])
        first_intra_duration = arr[0][-1] - arr[0][0]
        intra_times.append(first_intra_times)
        intra_durations.append(first_intra_duration)
        for prev_cluster, curr_cluster in zip(arr[:-1], arr[1:]):
            intra_time = np.diff(curr_cluster)
            intra_dur = curr_cluster[-1] - curr_cluster[0]
            inter_dur = curr_cluster[0] - prev_cluster[-1]
            intra_times.append(intra_time)
            intra_durations.append(intra_dur)
            inter_durations.append(inter_dur)
        self._intra_times['value'] = np.array(intra_times)
        self._intra_durations['value'] = np.array(intra_durations)
        self._inter_durations['value'] = np.array(inter_durations)

    def store_intra_times_histogram(self, bias='left', distribution_model=None):
        """
        bias:
            type <str>

        distribution_model:
            type <str> or None
        """
        _intra_times = np.concatenate(self.intra_times['value'], axis=0)
        edges = np.arange(-0.5, max(_intra_times) + 1.5, 1)
        H = Histogram(_intra_times, distribution_model=distribution_model)
        H.initialize_histogram(edges, criterion='edges', bias=bias)
        self._intra_times['histogram'] = H

    def store_intra_durations_histogram(self, bias='left', distribution_model=None):
        """
        bias:
            type <str>

        distribution_model:
            type <str> or None
        """
        edges = np.arange(-0.5, max(self.intra_durations['value']) + 1.5, 1)
        H = Histogram(self.intra_durations['value'], distribution_model=distribution_model)
        H.initialize_histogram(edges, criterion='edges', bias=bias)
        self._intra_durations['histogram'] = H

    def store_inter_durations_histogram(self, bias='left', distribution_model=None):
        """
        bias:
            type <str>

        distribution_model:
            type <str> or None
        """
        _percentiles = np.array([25, 75])
        iqr_bounds = np.percentile(self.inter_durations['value'], _percentiles)
        iqr_delta = max(iqr_bounds) - min(iqr_bounds)
        wbin = (iqr_delta * 2) // self.inter_durations['value'].size ** (1/3)
        lbin = np.floor(np.min(self.inter_durations['value']))
        rbin = np.ceil(np.max(self.inter_durations['value']))
        nbin = 2 * int((rbin - lbin) // wbin)
        edges = np.linspace(lbin, rbin, nbin, dtype=int)
        H = Histogram(self.inter_durations['value'], distribution_model=distribution_model)
        H.initialize_histogram(edges, criterion='edges', bias=bias)
        self._inter_durations['histogram'] = H

class UnbiasedEstimators():

    """
    This class contains methods to obtain and view the max spectrum
    (and associated parameters) from an extreme-valued distribution.
    """

    def __init__(self, arr, jndices, distribution_model, nresamples=100, nshuffles=3, nan_repl=None, exclude_original_sample=False):
        """
        arr:
            type <array>

        jndices:
            type <array / slice>

        distribution_model:
            type <str>

        nresamples:
            type <int>

        nshuffles:
            type <int>

        nan_repl:
            type <int / float> or None

        exclude_original_sample:
            type <bool>
        """
        super().__init__()
        self.arr = arr
        self.jndices = jndices
        self.model = DistributionModel(distribution_model)
        self.nresamples = nresamples
        self.nshuffles = nshuffles
        self.nan_repl = nan_repl
        self.exclude_original_sample = exclude_original_sample
        self.jmax = int(np.floor(np.log2(arr.size)))
        self.js = np.linspace(1, self.jmax, self.jmax, dtype=int)
        self.xj = self.js[jndices]
        self.function_map = self.load_function_map()
        self._max_spectra = []
        self._alpha = {}
        self._intercept = {}
        self._theta = {}
        self._point_estimators = {}
        self._extreme_value_estimates = {}

    @property
    def max_spectra(self):
        return self._max_spectra

    @property
    def alpha(self):
        return self._alpha

    @property
    def intercept(self):
        return self._intercept

    @property
    def theta(self):
        return self._theta

    @property
    def point_estimators(self):
        return self._point_estimators

    @property
    def extreme_value_estimates(self):
        return self._extreme_value_estimates

    @staticmethod
    def load_function_map():
        fkeys = ('mean', 'median', 'skew', 'kurtosis', 'standard deviation', 'standard error', 'maximum', 'minimum')
        fs = (np.mean, np.median, skew, kurtosis, np.std, sem, np.max, np.min)
        return dict(zip(fkeys, fs))

    @staticmethod
    def get_shapeshifted_data(arr, jprime):
        """
        arr:
            type <array>

        jprime:
            type <tuple>
        """
        desired_size_factor = np.prod([n for n in jprime if n != -1])
        if -1 in jprime:
            desired_size = arr.size // desired_size_factor * desired_size_factor
        else:
            desired_size = desired_size_factor
        return arr.flat[:desired_size].reshape(jprime)

    def get_data_maximums(self, arr):
        """
        arr:
            type <array>
        """
        res = np.max(arr.T, axis=0)
        if self.nan_repl is not None:
            res = res[res != self.nan_repl]
        return res[res > 0]

    def resample_indices(self, indices):
        """
        indices:
            type <array>
        """
        if self.nshuffles > 0:
            for n in range(self.nshuffles):
                np.random.shuffle(indices)
        np.random.shuffle(indices)
        return indices

    def get_full_max_spectrum(self, arr, ddof=0):
        """
        arr:
            type <array>

        ddof:
            type <int>
        """
        keys = ('dlog', 'npad', 'standard deviation', 'standard error', 'yj initial')
        res = {key : [] for key in keys}
        for j in self.js:
            jprime = (-1, 2**j)
            shapeshifted_arr = self.get_shapeshifted_data(arr, jprime)
            partial_spectrum = self.get_data_maximums(shapeshifted_arr)
            dlog = np.log2(partial_spectrum)
            st_dev = np.std(dlog, ddof=ddof)
            st_err = sem(dlog, ddof=ddof)
            yj_init = np.sum(dlog) / dlog.size # np.mean(dlog)
            args = (dlog, dlog.size, st_dev, st_err, yj_init)
            for key, arg in zip(keys, args):
                res[key].append(arg)
        return {key : np.array(value) for key, value in res.items()}

    def get_optimized_max_spectrum(self, max_spectrum, **kwargs):
        """
        max_spectrum:
            type <dict>
        """
        yj = max_spectrum['yj initial'][self.jndices]
        weights = max_spectrum['npad'][self.jndices]
        GLS = GeneralizedLeastSquares([self.xj, yj], self.model)
        optimization_result = GLS.fit(weights=weights, **kwargs)
        slope, intercept = optimization_result['parameters']
        alpha = 1 / slope
        yj_fit = self.model.f((slope, intercept), self.xj)
        return alpha, intercept, yj_fit

    def store_max_spectra(self, ddof=0, **kwargs):
        """
        ddof:
            type <int>
        """
        original_spectrum = self.get_full_max_spectrum(self.arr, ddof)
        alpha, intercept, yj_fit = self.get_optimized_max_spectrum(original_spectrum, **kwargs)
        original_spectrum['yj optimized'] = yj_fit
        self._max_spectra.append(original_spectrum)
        alphas, intercepts = [], []
        alphas.append(alpha)
        intercepts.append(intercept)
        indices = np.arange(self.arr.size, dtype=int)
        for ith_resample in range(self.nresamples):
            indices = self.resample_indices(indices)
            max_spectrum = self.get_full_max_spectrum(self.arr[indices], ddof)
            alpha, intercept, yj_fit = self.get_optimized_max_spectrum(max_spectrum, **kwargs)
            max_spectrum['yj optimized'] = yj_fit
            self._max_spectra.append(max_spectrum)
            alphas.append(alpha)
            intercepts.append(intercept)
        self.alpha['value'] = np.array(alphas)
        self.intercept['value'] = np.array(intercepts)

    def store_extremal_index(self):
        exponents = []
        yj_original_fit = self.max_spectra[0]['yj optimized']
        for max_spectrum, alpha in zip(self.max_spectra[1:], self.alpha['value'][1:]):
            dy = yj_original_fit - max_spectrum['yj optimized']
            exponents.append(dy.T * alpha)
        self._theta['value'] = 2 ** np.array(exponents)

    def store_point_estimators(self):
        for key in ('mean', 'median', 'maximum', 'minimum'):
            f = self.function_map[key]
            self._point_estimators[key] = f(self.theta['value'], axis=0)

    def store_histograms(self, walpha, wintercept, wtheta, alpha_model=None, intercept_model=None, theta_model=None, bias='left'):
        """
        walpha:
            type <int / float>

        wintercept:
            type <int / float>

        wtheta:
            type <float>

        alpha_model:
            type <str> or None

        intercept_model:
            type <str> or None

        theta_model:
            type <str> or None

        bias:
            type <str>
        """
        if self.exclude_original_sample == True:
            alpha_arr, intercept_arr = self.alpha['value'], self.intercept['value']
        else:
            alpha_arr, intercept_arr = np.delete(self.alpha['value'], 0), np.delete(self.intercept['value'], 0)
        theta_arr = self.theta['value'].copy().reshape(-1)
        alpha_hist = Histogram(alpha_arr, alpha_model)
        alpha_hist.initialize_histogram(value=walpha, criterion='bin width')
        intercept_hist = Histogram(intercept_arr, intercept_model)
        intercept_hist.initialize_histogram(value=wintercept, criterion='bin width')
        theta_hist = Histogram(theta_arr, theta_model)
        theta_hist.initialize_histogram(value=np.arange(0, 1 + wtheta, wtheta), criterion='edges')
        self._alpha.update(dict(histogram=alpha_hist))
        self._intercept.update(dict(histogram=intercept_hist))
        self._theta.update(dict(histogram=theta_hist))

    def store_statistics(self, ddof=0):
        """
        ddof:
            type <int>
        """
        if self.exclude_original_sample == True:
            alpha_arr, intercept_arr = self.alpha['value'], self.intercept['value']
        else:
            alpha_arr, intercept_arr = np.delete(self.alpha['value'], 0), np.delete(self.intercept['value'], 0)
        theta_arr = self.theta['value'].copy().reshape(-1)
        alpha_stats, intercept_stats, theta_stats = dict(), dict(), dict()
        for key, f in self.function_map.items():
            if key == 'standard deviation':
                alpha_stats[key] = f(alpha_arr, ddof=ddof)
                intercept_stats[key] = f(intercept_arr, ddof=ddof)
                theta_stats[key] = f(theta_arr, ddof=ddof)
            else:
                alpha_stats[key] = f(alpha_arr)
                intercept_stats[key] = f(intercept_arr)
                theta_stats[key] = f(theta_arr)
        self._alpha.update(dict(statistics=alpha_stats))
        self._intercept.update(dict(statistics=intercept_stats))
        self._theta.update(dict(statistics=theta_stats))

    def store_extreme_value_estimates(self):
        result = {}
        yj = np.array([max_spectra['yj optimized'] for max_spectra in self.max_spectra])
        for key in ('mean', 'median'):
            f = self.function_map[key]
            result[key] = int(round(np.min(2 ** f(yj, axis=0))))
        self._extreme_value_estimates.update(result)

class InterExceedanceDistribution():

    """
    This class contains methods to compare the distribution of
    inter-exceedance times with the corresponding distribution
    obtained from the inverse-transform sampling method.
    """

    def __init__(self, events):
        """
        events:
            type <dict>
        """
        super().__init__()
        self.events = events
        self._histograms = {}
        self._moment_estimators = {}
        self._time_thresholds = {}
        self._clusters = {}
        self._cluster_bias = None

    @property
    def histograms(self):
        return self._histograms

    @property
    def moment_estimators(self):
        return self._moment_estimators

    @property
    def time_thresholds(self):
        return self._time_thresholds

    @property
    def clusters(self):
        return self._clusters

    @property
    def cluster_bias(self):
        return self._cluster_bias

    def store_inter_exceedance_histogram(self, bin_width, bias='left'):
        """
        bin_width:
            type <int / float>

        bias:
            type <str>
        """
        inter_exceedance_times = self.events['time']
        w = int(bin_width)
        rightmost_edge = int(np.ceil(np.max(inter_exceedance_times)) / float(w)) * int(w)
        edges = np.arange(0, rightmost_edge +1, bin_width)
        H = Histogram(inter_exceedance_times, distribution_model=None)
        H.initialize_histogram(edges, criterion='edges', bias=bias, remove_empty_boundaries=True)
        self._histograms['inter-exceedance times'] = H

    def store_inverse_transform_sample_histogram(self):
        inter_exceedance_times = self.events['time']
        inter_exceedance_histogram = self._histograms['inter-exceedance times']
        x = np.random.uniform(size=inter_exceedance_times.size)
        y = np.log(1/x) * np.mean(inter_exceedance_times) / np.mean(np.log(1/x))
        H = Histogram(y, distribution_model=None)
        H.initialize_histogram(inter_exceedance_histogram.edges, criterion='edges', bias=inter_exceedance_histogram.bias)
        self._histograms['inverse-transform sample'] = H

    def store_moment_estimators(self, baseline=None):
        """
        baseline:
            type <int / float> or None
        """
        inter_exceedance_times = self.events['time']
        sub1, sub2 = inter_exceedance_times - 1, inter_exceedance_times - 2
        self._moment_estimators['first-order'] = 2 * np.square(np.sum(sub1)) / (inter_exceedance_times.size * np.sum(sub1 * sub2))
        self._moment_estimators['threshold'] = 2 * np.square(np.sum(inter_exceedance_times)) / (inter_exceedance_times.size * np.sum(np.square(inter_exceedance_times)))
        if baseline is not None:
            self._moment_estimators['baseline'] = baseline

    def store_time_thresholds(self):
        ts = np.sort(self.events['time'])
        n = ts.size + 1
        for bias, estimator in self.moment_estimators.items():
            index = int(estimator * n)
            try:
                self._time_thresholds[bias] = ts[index]
            except:
                self._time_thresholds[bias] = np.nan

    def store_temporal_clusters(self, bias='threshold'):
        """
        bias:
            type <str>
        """
        if bias not in list(self.time_thresholds.keys()):
            raise ValueError("invalid method bias: {}".format(bias))
        tc = self.time_thresholds[bias]
        parameters, inter_exceedance_times = self.events['parameter'], self.events['time']
        condition = (tc <= inter_exceedance_times)
        indices = np.where(condition)[0] + 1
        # S = EventSearcher({'dt' : inter_exceedance_times})
        # indices = S.search_events(parameters='dt', conditions='greater than or equal', values=tc)[1] +1
        result = {key : np.array(np.split(value, indices)) for key, value in parameters.items()}
        self._clusters.update(result)
        self._cluster_bias = bias

class TimeSeries(TimeSeriesViewer):

    def __init__(self, DB, timestep, event_type, ref_parameter, directory, ticksize=7, labelsize=8, textsize=9, titlesize=11, headersize=20, cellsize=15):
        """
        DB:
            type <custom class>

        timestep:
            type <str>

        event_type:
            type <str>

        ref_parameter:
            type <str>

        directory:
            type <str>

        ticksize:
            type <int / float>

        labelsize:
            type <int / float>

        textsize:
            type <int / float>

        titlesize:
            type <int / float>

        headersize:
            type <int / float>

        cellsize:
            type <int / float>
        """
        DB.verify_event_type(event_type)
        super().__init__(timestep, event_type, ref_parameter, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize)
        self.DB = DB
        self._events = []
        self._available.update(dict(DB.available))

    @property
    def events(self):
        return self._events

    def load_events(self, extreme_parameter, extrema='maximum', search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, nan_kwargs=None):
        """
        extreme_parameter:
            type <str>

        extrema:
            type <str>

        ...

        nan_kwargs:
            type <dict> or None
        """
        ## ORGANIZE SEARCH KWARGS
        search_kwargs = dict()
        search_kwargs['parameters'] = search_parameters
        search_kwargs['conditions'] = search_conditions
        search_kwargs['values'] = search_values
        search_kwargs['apply_to'] = apply_to
        search_kwargs['modifiers'] = modifiers
        ## ORGANIZE NAN KWARGS
        if nan_kwargs is None:
            nan_kwargs = dict(parameter=None, policy=None, repl=None)
        elif not isinstance(nan_kwargs, dict):
            raise ValueError("invalid type(nan_kwargs): {}".format(type(nan_kwargs)))
        ## ORGANIZE IDENTIFIERS
        identifiers = dict()
        identifiers['extreme parameter'] = extreme_parameter
        identifiers['search kwargs'] = search_kwargs
        identifiers['nan kwargs'] = nan_kwargs
        ## LOAD DATA (SEARCH IS OPTIONAL)
        data = getattr(self.DB, self.event_type)
        parameters = self.DB.filter_nans(data, **nan_kwargs)
        parameters = self.DB.reorganize_parameters(parameters, self.ref_parameter, extreme_parameter)
        S = EventSearcher(parameters)
        parameters = S.search_events(**search_kwargs)[0]
        ## GROUP EVENTS BY REGULARIZED TIME
        parameters['elapsed'] = self.DB.get_elapsed_time(parameters, self.timestep)
        indices = npi.group_by(parameters['elapsed']).argmax(parameters[self.ref_parameter])[1]
        parameters = {key : arr[indices] for key, arr in parameters.items()}
        tpad = self.DB.pad_datetimes(parameters, self.timestep)
        parameters = self.DB.group_regular_parameters(tpad, parameters, self.timestep, parameters['elapsed'])
        ## SAVE EVENTS
        result = {'parameter' : parameters, 'identifier' : identifiers, 'inter-exceedance' : dict()}
        self._events.append(result)

    def load_inter_exceedances(self, extreme_values=None, extreme_conditions='greater than or equal', apply_to='all', modifiers=None):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None

        extreme_conditions:
            type <str / tuple / list / array>

        ...
        """
        extreme_values = self.autocorrect_extreme_values(extreme_values)
        for events in self.events:
            parameters, identifiers = events['parameter'], events['identifier']
            S = EventSearcher(parameters)
            for extreme_value in extreme_values:
                if extreme_value not in list(events['inter-exceedance'].keys()):
                    result = dict()
                    if extreme_value is None:
                        inter_exceedance_times = np.diff(parameters['elapsed'])
                        result['parameter'] = parameters
                        result['search kwargs'] = None
                    else:
                        search_kwargs = dict()
                        search_kwargs['parameters'] = self.ref_parameter
                        search_kwargs['conditions'] = extreme_conditions
                        search_kwargs['values'] = extreme_value
                        search_kwargs['apply_to'] = apply_to
                        search_kwargs['modifiers'] = modifiers
                        _parameters = S.search_events(**search_kwargs)[0]
                        inter_exceedance_times = np.diff(_parameters['elapsed'])
                        result['parameter'] = _parameters
                        result['search kwargs'] = search_kwargs
                    result['time'] = inter_exceedance_times
                    events['inter-exceedance'][extreme_value] = result

    def load_inter_exceedance_distributions(self, extreme_values=None, bin_width=30, bias='left'):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None

        bin_width:
            type <int / float>

        bias:
            type <str>
        """
        extreme_values = self.autocorrect_extreme_values(extreme_values)
        for events in self.events:
            for extreme_value in extreme_values:
                if extreme_value not in list(events['inter-exceedance'].keys()):
                    raise ValueError("events of extreme value {} have not been loaded".format(extreme_value))
                if 'distribution' in list(events['inter-exceedance'][extreme_value].keys()):
                    IED = events['inter-exceedance'][extreme_value]['distribution']
                    IED.store_inter_exceedance_histogram(bin_width, bias)
                    IED.store_inverse_transform_sample_histogram()
                else:
                    IED = InterExceedanceDistribution(events['inter-exceedance'][extreme_value])
                    IED.store_inter_exceedance_histogram(bin_width, bias)
                    IED.store_inverse_transform_sample_histogram()
                    events['inter-exceedance'][extreme_value]['distribution'] = IED

    def load_unbiased_estimators(self, jndices, distribution_model='linear equation', walpha=0.05, wintercept=0.05, wtheta=0.1, alpha_model=None, intercept_model=None, theta_model=None, exclude_original_sample=False, ddof=0, bias='left', nresamples=1000, nshuffles=3, nan_repl=None, **kwargs):
        """
        jndices:
            type <array / slice>

        distribution_model:
            type <str>

        walpha:
            type <int / float>

        wintercept:
            type <int / float>

        wtheta:
            type <float>

        alpha_model:
            type <str> or None

        intercept_model:
            type <str> or None

        theta_model:
            type <str> or None

        exclude_original_sample:
            type <bool>

        ddof:
            type <int>

        bias:
            type <str>

        nresamples:
            type <int>

        nshuffles:
            type <int>

        nan_repl:
            type <int / float> or None
        """
        for events in self.events:
            parameters = events['parameter']
            UB = UnbiasedEstimators(parameters[self.ref_parameter], jndices, distribution_model, nresamples, nshuffles, nan_repl, exclude_original_sample)
            UB.store_max_spectra(ddof, **kwargs)
            UB.store_extremal_index()
            UB.store_point_estimators()
            UB.store_histograms(walpha, wintercept, wtheta, alpha_model, intercept_model, theta_model, bias)
            UB.store_statistics(ddof)
            UB.store_extreme_value_estimates()
            events['unbiased estimators'] = UB

    def load_moment_estimators_and_time_thresholds(self, extreme_values=None, baseline=None):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None

        baseline:
            type <int / float> or None
        """
        extreme_values = self.autocorrect_extreme_values(extreme_values)
        for events in self.events:
            for extreme_value in extreme_values:
                if extreme_value not in list(events['inter-exceedance'].keys()):
                    raise ValueError("events of extreme value {} have not been loaded".format(extreme_value))
                if 'distribution' in list(events['inter-exceedance'][extreme_value].keys()):
                    IED = events['inter-exceedance'][extreme_value]['distribution']
                    IED.store_moment_estimators(baseline)
                    IED.store_time_thresholds()
                else:
                    IED = InterExceedanceDistribution(events['inter-exceedance'][extreme_value])
                    IED.store_moment_estimators(baseline)
                    IED.store_time_thresholds()
                    events['inter-exceedance'][extreme_value]['distribution'] = IED

    def load_clusters(self, extreme_values=None, bias='threshold'):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None

        bias:
            type <str>
        """
        extreme_values = self.autocorrect_extreme_values(extreme_values)
        for events in self.events:
            for extreme_value in extreme_values:
                if extreme_value not in list(events['inter-exceedance'].keys()):
                    raise ValueError("events of extreme value {} have not been loaded".format(extreme_value))
                if 'distribution' not in list(events['inter-exceedance'][extreme_value].keys()):
                    raise ValueError("inter_exceedance distribution for extreme value {} has not been loaded".format(extreme_value))
                IED = events['inter-exceedance'][extreme_value]['distribution']
                IED.store_temporal_clusters(bias)

    def subview_inter_exceedance_histograms(self, ax, events, extreme_value, facecolors, xlim, yspace):
        """
        ax:
            type <matplotlib oject>

        events:
            type <dict>

        extreme_value:
            type <int / float> or None

        facecolors:
            type <tuple / list / array>

        xlim:
            type <tuple / list / array>

        yspace:
            type <int>
        """
        if extreme_value not in list(events['inter-exceedance'].keys()):
            raise ValueError("events of extreme value {} have not been loaded".format(extreme_value))
        inter_exceedances = events['inter-exceedance'][extreme_value]
        IED = inter_exceedances['distribution']
        search_kwargs = events['inter-exceedance'][extreme_value]['search kwargs']
        nevents_label = '${:,}$ {}s'.format(IED.events['time'].size +1, self.event_type_label)
        search_label = self.get_search_label(search_kwargs)
        bbox = {'facecolor': 'gray', 'alpha': 0.2, 'pad': 2}
        sc_label = self.get_solar_cycle_label(events)
        ex_label = self.get_extreme_parameter_label(events)
        title = '{}: {}'.format(sc_label, ex_label)
        legend_title = '{} {}'.format(self.event_type_label, search_label)
        largest_count = 0
        for (key, histogram), facecolor in zip(IED.histograms.items(), facecolors):
            ax.bar(histogram.midpoints, histogram.observed_counts, width=histogram.bin_widths, facecolor=facecolor, label=key.title(), alpha=0.5)
            xticks = histogram.edges
            yprime = np.max(histogram.observed_counts)
            if yprime > largest_count:
                largest_count = yprime
        ax.set_xticks(xticks[::2])
        ax.set_xticks(xticks[1::2], minor=True)
        ymax = self.round_up(largest_count * 1.3, nearest=yspace)
        # yticks = np.arange(0, ymax, yspace)
        # ax.set_yticks(yticks[::2])
        # ax.set_yticks(yticks[1::2], minor=True)
        ax.set_xlim(xlim)
        # ax.set_ylim([0, ymax])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.text(0.95, 0.95, nevents_label, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, fontsize=self.labelsize, bbox=bbox)
        ax.set_title(title, fontsize=self.titlesize)
        ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
        return ax, legend_title, ymax

    def view_inter_exceedance_histograms(self, extreme_values, layout, facecolors=('r', 'b'), xlim=(0, 300), yspace=25, save=False, **kwargs):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None

        layout:
            type <str>

        facecolors:
            type <tuple / list / array>

        xlim:
            type <tuple / list / array>

        yspace:
            type <int>

        save:
            type <bool>
        """
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < 2:
            raise ValueError("{} facecolors for 2 histograms".format(ncolors))
        if layout == 'overlay':
            raise ValueError("invalid layout '{}' for this figure".method(layout))
        vecdir = self.get_vecdir(layout)
        legend_kwargs = dict(ncol=2, loc='lower center', mode='expand', fontsize=self.labelsize)
        extreme_values = self.autocorrect_extreme_values(extreme_values)
        xlabel = 'Inter-Exceedance Times ({}s)'.format(self.timestep)
        ylabel = 'Observed Frequency'
        header = 'Distribution of Inter-Exceedance Times Against Exponential Decay'
        for extreme_value in extreme_values:
            ymax = []
            fig, axes = self.get_dim2_figure_and_axes(layout, n=len(self.events), **kwargs)
            try:
                for ax, events in zip(axes.ravel(), self.events):
                    ax, legend_title, _ymax = self.subview_inter_exceedance_histograms(ax, events, extreme_value, facecolors, xlim, yspace)
                    ymax.append(_ymax)
                handles, labels = axes.ravel()[0].get_legend_handles_labels()
                legend_kwargs['handles'] = [handles[0], handles[1]]
                legend_kwargs['labels'] = [labels[0], labels[1]]
            except:
                for events in self.events:
                    axes, legend_title, _ymax = self.subview_inter_exceedance_histograms(axes, events, extreme_value, facecolors, xlim, yspace)
                    ymax.append(_ymax)
            ymax = np.max(ymax)
            yticks = np.arange(0, ymax, yspace)
            try:
                for ax in axes.ravel():
                    ax.set_yticks(yticks[::2])
                    ax.set_yticks(yticks[1::2], minor=True)
                    ax.set_ylim([0, ymax])
            except:
                axes.set_yticks(yticks[::2])
                axes.set_yticks(yticks[1::2], minor=True)
                axes.set_ylim([0, ymax])
            fig.suptitle(header, fontsize=self.titlesize)
            fig, axes = self.add_shared_labels(fig, axes, xlabel, ylabel, vecdir)
            fig.align_ylabels()
            fig.subplots_adjust(bottom=0.15, wspace=0.2, hspace=0.2)
            leg = fig.legend(**legend_kwargs)
            leg.set_title(legend_title, prop={'size': self.textsize})
            leg._legend_box.align = "center"
            frame = leg.get_frame()
            frame.set_edgecolor('k')
            if save == True:
                savename = 'inter_exceedance_distribution_{}_{}'.format(extreme_value, layout)
            else:
                savename = None
            self.display_image(fig, savename=savename)

    def subview_unbiased_estimators_histogram(self, ax, events, attribute, sym, counts_type, alpha, statistics, facecolor):
        """
        ax:
            type <matplotlib object>

        events:
            type <dict>

        attribute:
            type <str>

        sym:
            type <str>

        counts_type:
            type <str>

        alpha:
            type <int / float>

        statistics:
            type <tuple / list / array> or None

        facecolor:
            type <str>
        """
        UB = events['unbiased estimators']
        data = getattr(UB, attribute)
        H, arr = data['histogram'], data['value']
        w = H.edges[1] - H.edges[0]
        y = getattr(H, '{}_counts'.format(counts_type))
        sc_label = self.get_solar_cycle_label(events)
        ex_label = self.get_extreme_parameter_label(events)
        title = '{}: {}'.format(sc_label, ex_label)
        legend_title = '{} Distribution via {} Resamples'.format(sym, UB.nresamples)
        label = title[:]
        if statistics is None:
            bottom = 0.15
        else:
            bottom = 0.05 * len(statistics) + 0.15
            for statistic in statistics:
                stat_value = data['statistics'][statistic]
                label += '\n{}$(${}$)$ $=$ ${:.2f}$'.format(statistic, sym, stat_value)
        ax.bar(H.midpoints, y, width=H.bin_widths, facecolor=facecolor, label=label, alpha=alpha)
        if attribute == 'theta':
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(w/2))
            ax.xaxis.set_major_locator(ticker.MultipleLocator(w))
        else:
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(w))
            ax.xaxis.set_major_locator(ticker.MultipleLocator(w*2))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        if attribute == 'theta':
            ax.set_xlim([0, 1])
        ax.tick_params(axis='both', which='both', labelsize=self.ticksize)
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_title(title, fontsize=self.titlesize)
        ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
        return ax, legend_title, bottom

    def view_unbiased_estimators_histograms(self, attribute, layout, counts_type='observed', statistics=None, facecolors=('steelblue', 'orange', 'green', 'salmon'), save=False, **kwargs):
        """
        attribute:
            type <str>

        layout:
            type <str>

        counts_type:
            type <str>

        statistics:
            type <str / tuple / list / array> or None

        facecolors:
            type <tuple / list / array>

        save:
            type <bool>
        """
        if attribute not in ('alpha', 'intercept', 'theta'):
            raise ValueError("invalid attribute: {}".format(attribute))
        nevents = len(self.events)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < nevents:
            raise ValueError("{} facecolors for {} events".format(ncolors, nevents))
        if statistics is not None:
            if isinstance(statistics, str):
                statistics = [statistics]
            elif not isinstance(statistics, (tuple, list, np.ndarray)):
                raise ValueError("invalid type(statistics): {}".format(type(statistics)))
        vecdir = self.get_vecdir(layout)
        legend_kwargs = dict(loc='lower center', mode='expand', fontsize=self.labelsize)
        legend_kwargs['ncol'] = nevents //2 if nevents > 3 else nevents
        sym = self.available['label']['unbiased estimators'][attribute]
        ylabel = '{} Frequency'.format(counts_type[:].title())
        header = 'Histogram of {} Distribution'.format(sym)
        alpha = 1 / nevents
        fig, axes = self.get_dim2_figure_and_axes(layout, nevents, **kwargs)
        if layout == 'overlay':
            for events, facecolor in zip(self.events, facecolors):
                axes, legend_title, bottom = self.subview_unbiased_estimators_histogram(axes, events, attribute, sym, counts_type, alpha, statistics, facecolor)
            bottom += 0.05
            axes.set_title('')
        else:
            try:
                for ax, events, facecolor in zip(axes.ravel(), self.events, facecolors):
                    ax, legend_title, bottom = self.subview_unbiased_estimators_histogram(ax, events, attribute, sym, counts_type, alpha, statistics, facecolor)
            except:
                for events, facecolor in zip(self.events, facecolors):
                    axes, legend_title, bottom = self.subview_unbiased_estimators_histogram(axes, events, attribute, sym, counts_type, alpha, statistics, facecolor)
        fig.suptitle(header, fontsize=self.titlesize)
        fig, axes = self.add_shared_labels(fig, axes, sym, ylabel, vecdir)
        fig.align_ylabels()
        fig.subplots_adjust(bottom=bottom)
        leg = fig.legend(**legend_kwargs)
        leg.set_title(legend_title, prop={'size': self.textsize})
        leg._legend_box.align = "center"
        frame = leg.get_frame()
        frame.set_edgecolor('k')
        if save == True:
            savename = 'unbiased_estimators_histogram_{}_{}_{}'.format(attribute, counts_type, layout)
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def subview_max_spectrum(self, ax, events, keys, ith_resample, facecolors):
        """
        ax:
            type <matplotlib object>

        events:
            type <dict>

        keys:
            type <tuple / list / array>

        ith_resample:
            type <int>

        facecolors:
            type <tuple / list / array>
        """
        identifiers = events['identifier']
        UB = events['unbiased estimators']
        max_spectrum = UB.max_spectra[ith_resample]
        ymin, ymax = [], []
        for key, facecolor in zip(keys, facecolors):
            if key == 'yj optimized':
                ax.plot(UB.xj, max_spectrum[key], color=facecolor, alpha=1, label='Power-Law')
                ymin.append(np.min(max_spectrum[key]))
                ymax.append(np.max(max_spectrum[key]))
            elif key == 'yj initial':
                ax.scatter(UB.js, max_spectrum[key], color=facecolor, label='Max Spectrum', marker='.')
                ymin.append(np.min(max_spectrum[key]))
                ymax.append(np.max(max_spectrum[key]))
            elif key in ('standard deviation', 'standard error'):
                ymin.append(np.min(max_spectrum['yj initial'] - max_spectrum[key]))
                ymax.append(np.max(max_spectrum['yj initial'] + max_spectrum[key]))
                ax.errorbar(UB.js, max_spectrum['yj initial'], yerr=max_spectrum[key], alpha=0.7, capsize=5, label=key.title(), ecolor=facecolor, fmt='none')
        ax.set_xticks(UB.js[1::2], minor=True)
        ax.set_xticks(UB.js[::2])
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='both', which='both', labelsize=self.ticksize)
        ax.set_xlim([UB.js[0] - 0.5, UB.js[-1] + 0.5])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        return ax, int(np.floor(np.min(ymin))), int(np.ceil(np.max(ymax)))

    def subview_max_spectrum_logscale(self, ax, ylim, xlabel=None, ylabel=None, ax_color='gray'):
        """
        ax:
            type <matplotlib object>

        ylim:
            type <array>

        xlabel:
            type <str> or None

        ylabel:
            type <str> or None

        ax_color:
            type <str>
        """
        mirror_ax = self.get_mirror_ax(ax, basex=2, basey=2, xscale='log', yscale='log', ax_color=ax_color)
        xlim = np.array(ax.get_xlim())
        try:
            mirror_ax.set_xlim(2**xlim)
        except:
            mirror_ax.set_xlim(2**np.array(xlim))
        try:
            mirror_ax.set_ylim(2**ylim)
        except:
            mirror_ax.set_ylim(2**np.array(ylim))
        if xlabel is None:
            mirror_ax.set_xticklabels([])
        else:
            mirror_ax.set_xlabel(xlabel, fontsize=self.labelsize, labelpad=10)
        if ylabel is None:
            mirror_ax.set_yticklabels([])
        else:
            mirror_ax.set_ylabel(ylabel, fontsize=self.labelsize, labelpad=12)
        return mirror_ax

    def view_max_spectrum(self, layout, show_errors=False, show_fit=False, show_logscale=False, ith_resample=0, yspace=0.5, facecolors=('darkorange', 'r'), save=False, **kwargs):
        """
        layout:
            type <str>

        show_errors:
            type <bool>

        show_fit:
            type <bool>

        show_logscale:
            type <bool>

        ith_resample:
            type <int>

        yspace:
            type <int / float>

        facecolors:
            type <tuple / list / array>

        save:
            type <bool>
        """
        savename_suffix = []
        keys = ['yj initial']
        if show_errors == True:
            keys += ['standard deviation', 'standard error']
            header = 'Max Spectrum Errors'
            savename_suffix.append('errors')
        if show_fit == True:
            keys += ['yj optimized']
            header = 'Max Spectrum Power-Law'
            savename_suffix.append('fit')
        if show_logscale == True:
            savename_suffix.append('log')
        nevents = len(self.events)
        minimum_colors = len(keys)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < minimum_colors:
            raise ValueError("{} facecolors provided but {} required".format(ncolors, minimum_colors))
        if layout == 'overlay':
            raise ValueError("invalid layout '{}' for this figure".method(layout))
        if layout != 'square':
            if 'vertical' not in layout:
                raise ValueError("invalid layout '{}' for this figure".format(layout))
        vecdir = self.get_vecdir(layout)
        sym = self.available['label'][self.event_type]['unit'][self.ref_parameter]
        x1label = r'$j$ $(log_2$ {}s$)$'.format(self.timestep)
        x2label = r'$2^j$ $({}s)$'.format(self.timestep)
        y1label = r'$Y(j)$ $(log_2$ {}$)$'.format(sym)
        y2label = r'$2^{Y(j)}$ ' + r'$(${}$)$'.format(sym)
        legend_kwargs = dict(loc='lower center', ncol=minimum_colors, mode='expand', fontsize=self.labelsize)
        ymin, ymax = [], []
        fig, axes = self.get_dim2_figure_and_axes(layout, nevents, **kwargs)
        try:
            for ax, events in zip(axes.ravel(), self.events):
                ax, _ymin, _ymax = self.subview_max_spectrum(ax, events, keys, ith_resample, facecolors)
                ymin.append(_ymin)
                ymax.append(_ymax)
            handles, labels = ax.get_legend_handles_labels()
            legend_kwargs['handles'] = handles
            legend_kwargs['labels'] = labels
            yticks = np.arange(np.min(ymin), np.max(ymax) +yspace, yspace)
            ylim = np.array([yticks[0], yticks[-1]])
            for ax in axes.ravel():
                ax.set_yticks(yticks[1::2], minor=True)
                ax.set_yticks(yticks[::2])
                ax.set_ylim(ylim)
                ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
                ax.grid(color='k', linestyle=':', alpha=0.3)
            if show_logscale == True:
                for ax, events in zip(axes.ravel(), self.events):
                    UB = events['unbiased estimators']
                    mean_value = UB.alpha['statistics']['mean']
                    sym = self.available['label']['unbiased estimators']['alpha']
                    alpha_label = r'mean$(${}$)$ $=$ ${:.2f}$'.format(sym, mean_value)
                    sc_label = self.get_solar_cycle_label(events)
                    ex_label = self.get_extreme_parameter_label(events)
                    title = '{}: {}'.format(sc_label, ex_label)
                    text_box = ax.text(0.05, 0.95, '{}\n{}'.format(title, alpha_label), fontsize=self.textsize, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
                    text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
                try:
                    for ax in axes[0, :-1].ravel():
                        mirror_ax = self.subview_max_spectrum_logscale(ax, ylim, xlabel=x2label, ylabel=None)
                        mirror_ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                    for ax in axes[1:, -1].ravel():
                        mirror_ax = self.subview_max_spectrum_logscale(ax, ylim, xlabel=None, ylabel=y2label)
                        mirror_ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                    mirror_ax = self.subview_max_spectrum_logscale(axes[0, -1], ylim, xlabel=x2label, ylabel=y2label)
                    mirror_ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                    mirror_ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                except:
                    mirror_ax = self.subview_max_spectrum_logscale(axes[0], ylim, xlabel=x2label, ylabel=y2label)
                    mirror_ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                    mirror_ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                    for ax in axes[1:].ravel():
                        mirror_ax = self.subview_max_spectrum_logscale(ax, ylim, xlabel=None, ylabel=y2label)
                        mirror_ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
            else:
                for ax, events in zip(axes.ravel(), self.events):
                    UB = events['unbiased estimators']
                    mean_value = UB.alpha['statistics']['mean']
                    sym = self.available['label']['unbiased estimators']['alpha']
                    alpha_label = r'mean$(${}$)$ $=$ ${:.2f}$'.format(sym, mean_value)
                    sc_label = self.get_solar_cycle_label(events)
                    ex_label = self.get_extreme_parameter_label(events)
                    title = '{}: {}'.format(sc_label, ex_label)
                    text_box = ax.text(0.05, 0.95, alpha_label, fontsize=self.textsize, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
                    text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
                    ax.set_title(title, fontsize=self.titlesize)
        except:
            for events in self.events:
                axes, _ymin, _ymax = self.subview_max_spectrum(axes, events, keys, ith_resample, facecolors)
                ymin.append(_ymin)
                ymax.append(_ymax)
            yticks = np.arange(np.min(ymin), np.max(ymax) +yspace, yspace)
            ylim = (yticks[0], yticks[-1])
            axes.set_yticks(yticks[1::2], minor=True)
            axes.set_yticks(yticks[::2])
            axes.set_ylim(ylim)
            axes.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
            axes.grid(color='k', linestyle=':', alpha=0.3)
            if show_logscale == True:
                UB = self.events[0]['unbiased estimators']
                mean_value = UB.alpha['statistics']['mean']
                sym = self.available['label']['unbiased estimators']['alpha']
                alpha_label = r'mean$(${}$)$ $=$ ${:.2f}$'.format(sym, mean_value)
                sc_label = self.get_solar_cycle_label(self.events[0])
                ex_label = self.get_extreme_parameter_label(self.events[0])
                title = '{}: {}'.format(sc_label, ex_label)
                text_box = axes.text(0.05, 0.95, '{}\n{}'.format(title, alpha_label), fontsize=self.textsize, horizontalalignment='left', verticalalignment='top', transform=axes.transAxes)
                text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
                mirror_ax = self.subview_max_spectrum_logscale(axes, ylim, xlabel=x2label, ylabel=y2label)
                mirror_ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
                mirror_ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
            else:
                UB = self.events[0]['unbiased estimators']
                mean_value = UB.alpha['statistics']['mean']
                sym = self.available['label']['unbiased estimators']['alpha']
                alpha_label = r'mean$(${}$)$ $=$ ${:.2f}$'.format(sym, mean_value)
                text_box = axes.text(0.05, 0.95, alpha_label, fontsize=self.textsize, horizontalalignment='left', verticalalignment='top', transform=axes.transAxes)
                text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
                sc_label = self.get_solar_cycle_label(self.events[0])
                ex_label = self.get_extreme_parameter_label(self.events[0])
                title = '{}: {}'.format(sc_label, ex_label)
                axes.set_title(title, fontsize=self.titlesize)
            handles, labels = axes.get_legend_handles_labels()
            if (show_fit == True) and (show_errors == True):
                legend_kwargs['handles'] = [handles[0], handles[1], handles[2], handles[3]]
                legend_kwargs['labels'] = [labels[0], labels[1], labels[2], labels[3]]
            elif (show_fit == True) and (show_errors == False):
                legend_kwargs['handles'] = [handles[0], handles[1]]
                legend_kwargs['labels'] = [labels[0], labels[1]]
            elif (show_fit == False) and (show_errors == True):
                legend_kwargs['handles'] = [handles[0], handles[1]]
                legend_kwargs['labels'] = [labels[0], labels[1]]
            else:
                legend_kwargs['handles'] = [handles[0]]
                legend_kwargs['labels'] = [labels[0]]
        fig.suptitle(header, fontsize=self.titlesize)
        fig, axes = self.add_shared_labels(fig, axes, x1label, y1label, vecdir='row')
        fig.align_ylabels()
        fig.subplots_adjust(hspace=0.3, wspace=0.3)
        fig.legend(**legend_kwargs)
        if save == True:
            savename = 'max_spectrum'
            for suffix in savename_suffix:
                savename = '{}_{}'.format(savename, suffix)
            savename = '{}_{}'.format(savename, layout)
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def subview_point_estimators(self, ax, events, facecolors, yticks):
        """
        ax:
            type <matplotlib object>

        events:
            type <dict>

        facecolors:
            type <tuple / list / array>

        yticks:
            type <tuple / list / array>
        """
        UB = events['unbiased estimators']
        sc_label = self.get_solar_cycle_label(events)
        ex_label = self.get_extreme_parameter_label(events)
        title = '{}: {}'.format(sc_label, ex_label)
        markers = ('o', '*', '_', '_')
        keys = ('mean', 'median', 'maximum', 'minimum')
        for key, facecolor, marker in zip(keys, facecolors, markers):
            if key in ('mean', 'median'):
                label = r'{}$(\theta_j)$'.format(key)
            elif key == 'maximum':
                label = r'bounds$(\theta_j)$'
            else:
                label = None
            ax.scatter(UB.xj, UB.point_estimators[key], marker=marker, color=facecolor, label=label)
        ax.set_xticks(UB.xj[1::2], minor=True)
        ax.set_xticks(UB.xj[::2])
        ax.set_yticks(yticks[1::2], minor=True)
        ax.set_yticks(yticks[::2])
        ax.tick_params(axis='both', which='both', labelsize=self.ticksize)
        ax.set_xlim([UB.xj[0] - 0.5, UB.xj[-1] + 0.5])
        ax.set_ylim([0, 1])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_title(title, fontsize=self.titlesize)
        ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
        return ax

    def view_point_estimators(self, layout, facecolors=('steelblue', 'darkorange', 'k'), save=False, **kwargs):
        """
        layout:
            type <str>

        facecolors:
            type <tuple / list / array>

        save:
            type <bool>
        """
        nevents = len(self.events)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        if not isinstance(facecolors, list):
            facecolors = list(facecolors)
        facecolors.append(facecolors[-1])
        ncolors = len(facecolors)
        if ncolors < 3:
            raise ValueError("{} facecolors provided but 3 required".format(ncolors))
        if layout == 'overlay':
            raise ValueError("invalid layout '{}' for this figure".method(layout))
        if layout != 'square':
            if 'vertical' not in layout:
                raise ValueError("invalid layout '{}' for this figure".format(layout))
        xlabel, ylabel = r'$j$', r'$\theta_j$'
        yticks = np.arange(0, 1.1, 0.1)
        legend_kwargs = dict(loc='lower center', ncol=3, mode='expand', fontsize=self.labelsize)
        fig, axes = self.get_dim2_figure_and_axes(layout, nevents, **kwargs)
        try:
            for ax, events in zip(axes.ravel(), self.events):
                ax = self.subview_point_estimators(ax, events, facecolors, yticks)
            handles, labels = ax.get_legend_handles_labels()
            legend_kwargs['handles'] = handles
            legend_kwargs['labels'] = labels
        except:
            for events in self.events:
                axes = self.subview_point_estimators(axes, events, facecolors, yticks)
            handles, labels = axes.get_legend_handles_labels()
            legend_kwargs['handles'] = [handles[0], handles[1], handles[2]]
            legend_kwargs['labels'] = [labels[0], labels[1], labels[2]]
        fig, axes = self.add_shared_labels(fig, axes, xlabel, ylabel, vecdir='row')
        fig.subplots_adjust(hspace=0.3, wspace=0.3)
        fig.suptitle('Point Estimators', fontsize=self.titlesize)
        fig.legend(**legend_kwargs)
        if save == True:
            savename = 'point_estimators_{}'.format(layout)
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def subview_moment_estimates_with_time_threshold(self, ax_top, ax_btm, events, facecolors, ncolors, theta_ticks, time_ticks):
        """
        ax_top:
            type <matplotlib object>

        ax_btm:
            type <matplotlib object>

        events:
            type <dict>

        facecolors:
            type <tuple / list / array>

        ncolors:
            type <int>

        theta_ticks:
            type <array>

        time_ticks:
            type <array>
        """
        biases = ('first-order', 'threshold', 'baseline')
        sc_label = self.get_solar_cycle_label(events)
        ex_label = self.get_extreme_parameter_label(events)
        title = '{}: {}'.format(sc_label, ex_label)
        xs = []
        thetas = {key : [] for key in biases}
        times = {key : [] for key in biases}
        for extreme_value, inter_exceedance in events['inter-exceedance'].items():
            xs.append(extreme_value)
            IED = events['inter-exceedance'][extreme_value]['distribution']
            for key in list(IED.moment_estimators.keys()):
                thetas[key].append(IED.moment_estimators[key])
                times[key].append(IED.time_thresholds[key])
        if len(thetas['baseline']) == 0:
            include_baseline = False
            del thetas['baseline']
            del times['baseline']
        else:
            include_baseline = True
            if ncolors < 3:
                raise ValueError("3 facecolors required")
        xs = np.array(xs)
        indices = np.argsort(xs)
        xs = xs[indices]
        for key in list(thetas.keys()):
            thetas[key] = np.array(thetas[key])[indices]
            times[key] = np.array(times[key])[indices]
        ax_btm.set_yscale('log', basey=5)
        for ax in (ax_btm, ax_top):
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.yaxis.set_minor_formatter(ticker.NullFormatter())
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        for key, facecolor in zip(list(thetas.keys()), facecolors):
            if key == 'baseline':
                label = 'baseline'
            else:
                label = '{} bias'.format(key)
            ax_top.plot(xs, thetas[key], color=facecolor, alpha=0.45, label=label)
            ax_btm.plot(xs, times[key], color=facecolor, alpha=0.45)
        ax_top.set_yticks(theta_ticks[1::2], minor=True)
        ax_top.set_yticks(theta_ticks[::2])
        ax_top.set_ylim([theta_ticks[0], theta_ticks[-1]])
        ax_btm.set_ylim([1, time_ticks[-1]])
        for ax in (ax_top, ax_btm):
            ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
            ax.set_xlim([xs[0], xs[-1]])
            ax.grid(color='k', linestyle=':', alpha=0.3)
        ax_top.set_title(title, fontsize=self.labelsize)
        return ax_top, ax_btm, include_baseline

    def subview_maximum_extreme_value(self, ax, events, facecolor, label=None):
        """
        ax:
            type <matplotlib object>

        events:
            type <dict>

        facecolor:
            type <str>

        label:
            type <str> or None
        """
        parameters = events['parameter']
        _p = parameters[self.ref_parameter]
        xmax = np.max(_p)
        ylim = ax.get_ylim()
        handle = ax.axvline(xmax, ymin=ylim[0], ymax=ylim[-1], color=facecolor, label=label, linestyle='-', marker=None)
        return ax, xmax, handle

    def view_moment_estimators_with_time_threshold(self, show_max_extreme_value=False, facecolors=('green', 'darkorange', 'steelblue'), save=False, **kwargs):
        """
        show_max_extreme_value:
            type <bool>

        facecolors:
            type <tuple / list / array>

        save:
            type <bool>
        """
        nevents = len(self.events)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < 2:
            raise ValueError("at least 2 facecolors required")
        theta_ticks = np.arange(0, 1.1, 0.1)
        time_ticks = np.geomspace(1, 3125, 6).astype(int)
        xlabel = '{} $(${}$)$'.format(self.ref_parameter.title(), self.unit_label)
        theta_label = 'Moment Estimator\n' + self.available['label']['unbiased estimators']['theta']
        time_label = 'Time Threshold\n' + r'$T_C$ ({}s)'.format(self.timestep)
        header = r'Moment Estimators & Time Thresholds'
        legend_kwargs = dict(loc='lower center', ncol=3, mode='expand', fontsize=self.labelsize)
        label_max = 'Maximum {}'.format(self.ref_parameter.title())
        fig, axes = self.get_dim2_figure_and_axes(layout='double-horizontal', n=nevents*2, **kwargs)
        try:
            for ax_top, ax_btm, events in zip(axes[0, :].ravel(), axes[1, :].ravel(), self.events):
                ax_top, ax_btm, include_baseline = self.subview_moment_estimates_with_time_threshold(ax_top, ax_btm, events, facecolors, ncolors, theta_ticks, time_ticks)
                if show_max_extreme_value == True:
                    for ax in (ax_top, ax_btm):
                        ax, _, _ = self.subview_maximum_extreme_value(ax, events, facecolor='k', label=label_max)
            for ax in axes[-1, :].ravel():
                ax.set_xlabel(xlabel, fontsize=self.labelsize)
            for ax in axes[:, 1:].ravel():
                ax.set_yticklabels([])
            axes[0, 0].set_ylabel(theta_label, fontsize=self.labelsize)
            axes[-1, 0].set_ylabel(time_label, fontsize=self.labelsize)
        except:
            (ax_top, ax_btm) = axes
            for events in self.events:
                ax_top, ax_btm, include_baseline = self.subview_moment_estimates_with_time_threshold(ax_top, ax_btm, events, facecolors, ncolors, theta_ticks, time_ticks)
                if show_max_extreme_value == True:
                    for ax in (ax_top, ax_btm):
                        ax, _, _ = self.subview_maximum_extreme_value(ax, events, facecolor='k', label=label_max)
                # ax_top.set_xticklabels([])
                ax_btm.set_xlabel(xlabel, fontsize=self.labelsize)
                ax_top.set_ylabel(theta_label, fontsize=self.labelsize)
                ax_btm.set_ylabel(time_label, fontsize=self.labelsize)
        fig.align_ylabels()
        handles, labels = ax_top.get_legend_handles_labels()
        if include_baseline == True:
            if show_max_extreme_value == True:
                legend_kwargs['handles'] = [handles[0], handles[1], handles[2], handles[3]]
                legend_kwargs['labels'] = [labels[0], labels[1], labels[2], labels[3]]
            else:
                legend_kwargs['handles'] = [handles[0], handles[1], handles[2]]
                legend_kwargs['labels'] = [labels[0], labels[1], labels[2]]
        else:
            if show_max_extreme_value == True:
                legend_kwargs['handles'] = [handles[0], handles[1], handles[2]]
                legend_kwargs['labels'] = [labels[0], labels[1], labels[2]]
            else:
                legend_kwargs['handles'] = [handles[0], handles[1]]
                legend_kwargs['labels'] = [labels[0], labels[1]]
        fig.suptitle(header, fontsize=self.titlesize)
        fig.subplots_adjust(bottom=0.2, hspace=0.3, wspace=0.3)
        fig.legend(**legend_kwargs)
        if save == True:
            savename = 'moment_estimators_with_time_thresholds'
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def subview_chronological_clusters(self, ax, events, extreme_value, parameter, dt_key, f, dx, facecolors, search_kwargs, major='year', minor='month', fmt="%Y-%m-%d", rotation=None):
        """
        ax:
            type <matplotlib object>

        events:
            type <dict>

        extreme_value:
            type <int / float> or None

        parameter:
            type <str>

        dt_key:
            type <str>

        f:
            type <function> or None

        dx:
            type <datetime object>

        facecolors:
            type <tuple / list / array>

        search_kwargs:
            type <dict>

        major:
            type <str>

        minor:
            type <str>

        fmt:
            type <str>

        rotation:
            type <int / float> or None
        """
        sc_label = self.get_solar_cycle_label(events)
        ex_label = self.get_extreme_parameter_label(events)
        title = '{}: {}'.format(sc_label, ex_label)
        IED = events['inter-exceedance'][extreme_value]['distribution']
        try:
            TC = TemporalClusters(IED.clusters, IED.cluster_bias, **search_kwargs)
        except:
            TC = TemporalClusters(IED.clusters, IED.cluster_bias)
        t_clusters, p_clusters = TC.clusters[dt_key], TC.clusters[parameter]
        nclusters = len(t_clusters)
        nevents = np.concatenate(t_clusters, axis=0).size
        nlabel = '${:,}$ {}s \nvia ${:,}$ Clusters'.format(nevents, self.event_type_label, nclusters)
        xs, ys, ws, cs = [], [], [], []
        for ith_cluster, (ts, ps) in enumerate(zip(t_clusters, p_clusters)):
            if len(ts) == 1:
                x = ts[0]
                w = date2num(x + dx) - date2num(x)
            else:
                x_left, x_right = ts[0], ts[1]
                x = x_left + (x_right - x_left) / 2
                w = date2num(x_right) - date2num(x_left)
            if ith_cluster % 2 == 0:
                _color = facecolors[0]
            else:
                _color = facecolors[1]
            if f is None:
                for p in ps:
                    ys.append(p)
                    xs.append(x)
                    ws.append(w)
                    cs.append(_color)
            else:
                xs.append(x)
                ys.append(f(ps))
                ws.append(w)
                cs.append(_color)
        ymax = max(ys) * 1.3
        ax.bar(xs, ys, width=ws, color=cs)
        ax = self.transform_x_as_datetime(ax, major, minor, fmt, rotation)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        text_box = ax.text(0.05, 0.95, nlabel, fontsize=self.labelsize, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_ylim([0, ymax])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_title(title, fontsize=self.titlesize)
        ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
        return ax

    def view_chronological_clusters(self, parameter, statistic=None, extreme_values=None, layout='vertical', facecolors=('r', 'b'), search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, major='year', minor='month', fmt="%Y-%m-%d", rotation=None, ddof=0, dt_key='datetime', save=False, **kwargs):
        """
        parameter:
            type <str>

        statistic:
            type <str> or None

        extreme_values:
            type <int / float / tuple / list / array> or None

        layout:
            type <str>

        facecolors:
            type <tuple / list / array>

        ...

        major:
            type <str>

        minor:
            type <str>

        fmt:
            type <str>

        rotation:
            type <int / float> or None

        ddof:
            type <int>

        dt_key:
            type <str>

        save:
            type <bool>
        """
        nevents = len(self.events)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < 2:
            raise ValueError("2 facecolors required")
        if layout == 'overlay':
            raise ValueError("invalid layout '{}' for this figure".format(layout))
        if layout != 'square':
            if 'vertical' not in layout:
                raise ValueError("invalid layout '{}' for this figure".format(layout))
        vecdir = self.get_vecdir(layout)
        header = 'Chronological Clusters'
        xlabel = 'Date'
        if statistic is None:
            try:
                unit_label = self.available['label'][self.event_type]['unit'][parameter]
                ylabel = '{} $(${}$)$'.format(parameter.title(), unit_label)
            except:
                ylabel = '{}'.format(parameter.title())
            f = None
        else:
            if statistic == 'standard deviation':
                f = lambda args : np.std(args, ddof=ddof)
            else:
                fkeys = ('mean', 'median', 'max', 'min', 'delta')
                fstats = (np.mean, np.median, np.max, np.min, lambda args : max(args) - min(args))
                fmap = dict(zip(fkeys, fstats))
                f = fmap[statistic]
            try:
                unit_label = self.available['label'][self.event_type]['unit'][parameter]
                ylabel = '{} {} $(${}$)$'.format(statistic.title(), parameter.title(), unit_label)
            except:
                ylabel = '{} {}'.format(statistic.title(), parameter.title())
        if (search_parameters is None) and (search_conditions is None) and (search_values is None):
            search_kwargs = None
        else:
            search_kwargs = dict(parameters=search_parameters, conditions=search_conditions, values=search_values, apply_to=apply_to, modifiers=modifiers)
        search_label = self.get_search_label(search_kwargs)
        xyi, xyf = (0, 0), (0.1, 0.1)
        handle_a = Rectangle(xyi, *xyf, fill=True, visible=True, facecolor=facecolors[0], edgecolor='none')
        handle_b = Rectangle(xyi, *xyf, fill=True, visible=True, facecolor=facecolors[1], edgecolor='none')
        handle_e = Rectangle(xyi, *xyf, alpha=0)
        legend_kwargs = dict(ncol=3, loc='lower center', mode='expand', fontsize=self.labelsize)
        legend_kwargs['handles'] = [handle_e, tuple((handle_a, handle_b)), handle_e]
        legend_kwargs['labels'] = [' ', 'Consecutive Clusters', ' ']
        legend_kwargs['handler_map'] = {tuple: HandlerTuple(None)}
        dx = datetime.timedelta(**{'{}s'.format(self.timestep) : 1})
        _kw = dict(**kwargs)
        if 'sharex' in list(_kw.keys()):
            x_is_same = True
        else:
            x_is_same = False
        for extreme_value in self.autocorrect_extreme_values(extreme_values):
            fig, axes = self.get_dim2_figure_and_axes(layout=layout, n=nevents, **kwargs)
            try:
                for ax, events in zip(axes.ravel(), self.events):
                    ax = self.subview_chronological_clusters(ax, events, extreme_value, parameter, dt_key, f, dx, facecolors, search_kwargs, major, minor, fmt, rotation)
                try:
                    for ax in axes[:, 0].ravel():
                        ax.set_ylabel(ylabel, fontsize=self.labelsize)
                    for ax in axes[:, 1:].ravel():
                        ax.set_yticklabels([])
                    for ax in axes[-1, :].ravel():
                        ax.set_xlabel(xlabel, fontsize=self.labelsize)
                    # if x_is_same == True:
                    #     for ax in axes[1:, :].ravel():
                    #         ax.set_xticklabels([])
                except:
                    for ax in axes.ravel():
                        ax.set_ylabel(ylabel, fontsize=self.labelsize)
                    # if x_is_same == True:
                    #     for ax in axes[:-1].ravel():
                    #         ax.set_xticklabels([])
                    axes[-1].set_xlabel(xlabel, fontsize=self.labelsize)
            except:
                for events in self.events:
                    axes = self.subview_chronological_clusters(axes, events, extreme_value, parameter, dt_key, f, dx, facecolors, search_kwargs, major, minor, fmt, rotation)
                    axes.set_ylabel(ylabel, fontsize=self.labelsize)
                    axes.set_xlabel(xlabel, fontsize=self.labelsize)
            stitle = fig.suptitle(header, fontsize=self.titlesize)
            fig, axes = self.add_shared_labels(fig, axes, xlabel, ylabel, vecdir)
            fig.align_ylabels()
            if x_is_same == True:
                fig.subplots_adjust(hspace=0.4, bottom=0.2)
            else:
                fig.subplots_adjust(hspace=0.85, bottom=0.2)
            leg = fig.legend(**legend_kwargs)
            if search_label is None:
                leg.set_title('All Clusters', prop={'size': self.labelsize})
            else:
                leg.set_title(search_label, prop={'size': self.labelsize})
            leg._legend_box.align = "center"
            frame = leg.get_frame()
            frame.set_edgecolor('k')
            if save == True:
                savename = 'chronological_{}_{}'.format(parameter, extreme_value)
                if statistic is not None:
                    savename = '{}_{}'.format(savename, statistic)
            else:
                savename = None
            bbox_extra_artists = [stitle, leg]
            self.display_image(fig, savename=savename)

    def subview_relative_cluster_statistics(self, ax_top, ax_mid, ax_btm, events, extreme_value, facecolors, ylabels, search_kwargs):
        """
        ax_top:
            type <matplotlib object>

        ax_btm:
            type <matplotlib object>

        events:
            type <dict>

        extreme_value:
            type <int / float> or None

        facecolors:
            type <tuple / list / array>

        ylabels:
            type <tuple / list / array>

        search_kwargs:
            type <dict> or None
        """
        sc_label = self.get_solar_cycle_label(events)
        ex_label = self.get_extreme_parameter_label(events)
        title = '{}: {}'.format(sc_label, ex_label)
        IED = events['inter-exceedance'][extreme_value]['distribution']
        try:
            TC = TemporalClusters(IED.clusters, IED.cluster_bias, **search_kwargs)
        except:
            TC = TemporalClusters(IED.clusters, IED.cluster_bias)
        TC.store_relative_statistics()
        x, ytop, ymid, ybtm = TC.relative_statistics['cluster size'], TC.relative_statistics['number events'], TC.relative_statistics['number clusters'], TC.relative_statistics['relative probability']
        xticks = np.arange(0, x[-1] +2, 1).astype(int)
        xmax = max(xticks)
        ymax = []
        for ax, y, facecolor, ylabel in zip((ax_top, ax_mid, ax_btm), (ytop, ymid, ybtm), facecolors, ylabels):
            if ylabel is None:
                _label = None
            else:
                _label = ylabel[:].replace('\n', ' ')
            ax.bar(x, y, width=1, facecolor=facecolor, label=_label)
            ymax.append(np.max(y))
            if xmax <= 10:
                ax.set_xticks(xticks[1::2], minor=True)
                ax.set_xticks(xticks[::2])
            elif 10 < xmax <= 45:
                ax.set_xticks(xticks, minor=True)
                ax.set_xticks(xticks[::5])
            else:
                ax.set_xticks(xticks, minor=True)
                ax.set_xticks(xticks[::10])
            ax.set_xlim([0, xticks[-1]])
            ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
        ax_top.set_title(title, fontsize=self.labelsize)
        return ax_top, ax_mid, ax_btm, ymax

    def view_relative_cluster_statistics(self, extreme_values=None, facecolors=('darkorange', 'green', 'purple'), search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, yspace=50, save=False, **kwargs):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None

        facecolors:
            type <tuple / list / array>

        ...

        yspace:
            type <int>

        save:
            type <bool>
        """
        nevents = len(self.events)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < 3:
            raise ValueError("at least 3 facecolors required, but only {} provided".format(ncolors))
        xlabel = 'Cluster Size'
        ylabels = ('Number of\nClusters', 'Number of\nExtreme Events', 'Relative\nProbability')
        none_labels = [None for ylabel in ylabels]
        header = 'Relative Cluster Statistics'
        if (search_parameters is None) and (search_conditions is None) and (search_values is None):
            search_kwargs = None
        else:
            search_kwargs = dict(parameters=search_parameters, conditions=search_conditions, values=search_values, apply_to=apply_to, modifiers=modifiers)
        search_label = self.get_search_label(search_kwargs)
        legend_kwargs = dict(ncol=3, loc='lower center', mode='expand', fontsize=self.labelsize)
        ybtm_ticks = np.arange(0, 1.1, 0.1)
        for extreme_value in self.autocorrect_extreme_values(extreme_values):
            ytop_max, ymid_max = 0, 0
            fig, axes = self.get_dim2_figure_and_axes(layout='triple-horizontal', n=nevents*3, **kwargs)
            try:
                for idx, (ax_top, ax_mid, ax_btm, events) in enumerate(zip(axes[0, :].ravel(), axes[1, :].ravel(), axes[2, :].ravel(), self.events)):
                    if idx == 0:
                        ax_top, ax_mid, ax_btm, (_ytop_max, _ymid_max, _) = self.subview_relative_cluster_statistics(ax_top, ax_mid, ax_btm, events, extreme_value, facecolors, ylabels, search_kwargs)
                    else:
                        ax_top, ax_mid, ax_btm, (_ytop_max, _ymid_max, _) = self.subview_relative_cluster_statistics(ax_top, ax_mid, ax_btm, events, extreme_value, facecolors, none_labels, search_kwargs)
                    if _ytop_max > ytop_max:
                        ytop_max = _ytop_max
                    if _ymid_max > ymid_max:
                        ymid_max = _ymid_max
            except:
                (ax_top, ax_mid, ax_btm) = axes
                for idx, events in enumerate(self.events):
                    if idx == 0:
                        ax_top, ax_mid, ax_btm, (_ytop_max, _ymid_max, _) = self.subview_relative_cluster_statistics(ax_top, ax_mid, ax_btm, events, extreme_value, facecolors, ylabels, search_kwargs)
                    else:
                        ax_top, ax_mid, ax_btm, (_ytop_max, _ymid_max, _) = self.subview_relative_cluster_statistics(ax_top, ax_mid, ax_btm, events, extreme_value, facecolors, none_labels, search_kwargs)
                    if _ytop_max > ytop_max:
                        ytop_max = _ytop_max
                    if _ymid_max > ymid_max:
                        ymid_max = _ymid_max
            ytop_max = self.round_up(ytop_max, nearest=yspace)
            ymid_max = self.round_up(ymid_max, nearest=yspace)
            ytop_ticks = np.arange(0, ytop_max + yspace + 1, yspace).astype(int)
            ymid_ticks = np.arange(0, ymid_max + yspace + 1, yspace).astype(int)
            try:
                for ax in axes[0, :].ravel():
                    ax.set_yticks(ytop_ticks[::2])
                    ax.set_yticks(ytop_ticks[1::2], minor=True)
                for ax in axes[1, :].ravel():
                    ax.set_yticks(ymid_ticks[::2])
                    ax.set_yticks(ymid_ticks[1::2], minor=True)
                for ax in axes[2, :].ravel():
                    ax.set_yticks(ybtm_ticks[::2])
                    ax.set_yticks(ybtm_ticks[1::2], minor=True)
                for ax in axes[:-1, :].ravel():
                    ax.set_xticklabels([])
                for ax in axes[:, 1:].ravel():
                    ax.set_yticklabels([])
                for ax in axes[-1, :].ravel():
                    ax.set_xlabel(xlabel, fontsize=self.labelsize)
                for ax, ylabel in zip(axes[:, 0].ravel(), ylabels):
                    ax.set_ylabel(ylabel, fontsize=self.labelsize)
            except:
                # (ax_top, ax_mid, ax_btm) = axes
                for ax, ticks, ylabel in zip(axes.ravel(), (ytop_ticks, ymid_ticks, ybtm_ticks), ylabels):
                    ax.set_yticks(ticks[::2])
                    ax.set_yticks(ticks[1::2], minor=True)
                    ax.set_xlabel(xlabel, fontsize=self.labelsize)
                    ax.set_ylabel(ylabel, fontsize=self.labelsize)
            for ax in axes.ravel():
                ax.xaxis.label.set_size(self.labelsize)
                ax.yaxis.label.set_size(self.labelsize)
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(self.ticksize)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(self.ticksize)
                ax.grid(color='k', alpha=0.3, linestyle=':')
            fig.subplots_adjust(hspace=0.25, wspace=0.25, bottom=0.2)
            fig.suptitle(header, fontsize=self.titlesize)
            leg = fig.legend(**legend_kwargs)
            if search_label is None:
                leg.set_title('All Clusters', prop={'size': self.labelsize})
            else:
                leg.set_title(search_label, prop={'size': self.labelsize})
            leg._legend_box.align = "center"
            frame = leg.get_frame()
            frame.set_edgecolor('k')
            if save == True:
                savename = 'relative_cluster_statistics_{}'.format(extreme_value)
            else:
                savename = None
            self.display_image(fig, savename=savename)

    def subview_cluster_duration_histograms(self, ax, attribute, TC, tc, facecolor, tc_color, title, attribute_label, tc_label, alpha, show_time_line=False):
        """

        """
        durs = getattr(TC, attribute)
        histogram = durs['histogram']
        handle_bar = ax.bar(histogram.midpoints, histogram.observed_counts, width=histogram.bin_widths, color=facecolor, label=attribute_label, alpha=alpha)
        handles = [handle_bar]
        # if 'intra' in attribute:
        ylim = ax.get_ylim()
        if show_time_line == True:
            handle_tc = ax.axvline(tc, ymin=0, ymax=ylim[-1], color=tc_color, linestyle='-', marker=None, label=tc_label)
            ax.legend(handles=[handle_tc], labels=[tc_label], bbox_to_anchor=(0.95, 0.95), loc='upper right', borderaxespad=0., fontsize=self.labelsize)
            handles.append(handle_tc)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
        xlim = ax.get_xlim()
        ax.grid(color='k', linestyle=':', alpha=0.3)
        if title is not None:
            ax.set_title(title, fontsize=self.labelsize)
        return ax, handles

    def view_cluster_duration_histograms(self, extreme_values, layout, zoom_in=False, show_intra_times=False, show_intra_durations=False, show_inter_durations=False, bias='left', intra_times_model=None, intra_durations_model=None, inter_durations_model=None, facecolors=('mediumpurple', 'darkgreen', 'steelblue', 'salmon', 'darkorange'), search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, save=False, **kwargs):
        """
        extreme_values:
            type <int / float / tuple / list / array> or None

        layout:
            type <str>

        show_intra_times:
            type <bool>

        show_intra_durations:
            type <bool>

        show_inter_durations:
            type <bool>

        facecolors:
            type <tuple / list / array>

        ...

        save:
            type <bool>
        """
        show_args = np.array([show_intra_times, show_intra_durations, show_inter_durations])
        if np.any(show_args == True) != True:
            raise ValueError("all show_ args are set to False")
        show_indices = np.where(show_args == True)[0]
        _attributes = np.array(['intra_times', 'intra_durations', 'inter_durations'])
        attributes = _attributes[show_indices]
        nshows = len(attributes)
        nevents = len(self.events)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < nshows:
            raise ValueError("not enough facecolors provided")
        tc_color = facecolors[-1]
        rev_colors = np.array(facecolors)[::-1]
        if (search_parameters is None) and (search_conditions is None) and (search_values is None):
            search_kwargs = None
        else:
            search_kwargs = dict(parameters=search_parameters, conditions=search_conditions, values=search_values, apply_to=apply_to, modifiers=modifiers)
        search_label = self.get_search_label(search_kwargs)
        xlabel = 'Time $(${}s$)$'.format(self.timestep)
        ylabel = 'Observed Frequency'
        tc_partial_label = r'Time Threshold $t_C$'
        if layout == 'overlay':
            alpha = 1 / nevents
            ncol = nshows +1 if nshows < 3 else (nshows +1)//2
        else:
            alpha = 1
            ncol = 2
        vecdir = self.get_vecdir(layout)
        legend_kwargs = dict(loc='lower center', mode='expand', fontsize=self.labelsize)
        legend_kwargs['ncol'] = ncol
        time_thresholds, temporal_clusters, titles = [], [], []
        for extreme_value in self.autocorrect_extreme_values(extreme_values):
            time_thresholds, temporal_clusters, titles, time_threshold_labels = [], [], [], []
            for events in self.events:
                sc_label = self.get_solar_cycle_label(events)
                ex_label = self.get_extreme_parameter_label(events)
                title = '{}: {}'.format(sc_label, ex_label)
                IED = events['inter-exceedance'][extreme_value]['distribution']
                tc = IED.time_thresholds[IED.cluster_bias]
                tc_label = r'$t_C$ $=$ ${}$ {}s'.format(tc, self.timestep)
                try:
                    TC = TemporalClusters(IED.clusters, IED.cluster_bias, **search_kwargs)
                except:
                    TC = TemporalClusters(IED.clusters, IED.cluster_bias)
                TC.store_cluster_times_and_durations()
                if show_intra_times == True:
                    TC.store_intra_times_histogram(bias, distribution_model=intra_times_model)
                if show_intra_durations == True:
                    TC.store_intra_durations_histogram(bias, distribution_model=intra_durations_model)
                if show_inter_durations == True:
                    TC.store_inter_durations_histogram(bias, distribution_model=inter_durations_model)
                time_thresholds.append(tc)
                time_threshold_labels.append(tc_label)
                temporal_clusters.append(TC)
                titles.append(title)
            mean_time_threshold = np.mean(time_thresholds)
            for ith_attribute, (attribute, _facecolor) in enumerate(zip(attributes, facecolors)):
                attribute_label = attribute.replace('_', ' ').title().replace(' ', '-')
                header = '{} of Clusters'.format(attribute_label)
                handle_bars, handle_times = [], []
                fig, axes = self.get_dim2_figure_and_axes(layout=layout, n=nevents, **kwargs)
                try:
                    show_time_line = True
                    for i, (ax, TC, tc, title, tc_label) in enumerate(zip(axes.ravel(), temporal_clusters, time_thresholds, titles, time_threshold_labels)):
                        if i == 0:
                            _attribute_label = attribute_label[:]
                        else:
                            _attribute_label = None
                        ax, _handles = self.subview_cluster_duration_histograms(ax, attribute, TC, tc, _facecolor, tc_color, title, _attribute_label, tc_label, alpha, show_time_line)
                    try:
                        handles = [_handles[0], _handles[1]]
                        labels = [attribute_label, tc_partial_label]
                    except:
                        handles = [_handles[0]]
                        labels = [attribute_label]
                    legend_kwargs['handles'] = handles
                    legend_kwargs['labels'] = labels
                except:
                    if layout == 'overlay':
                        show_time_line = False
                        for i, (TC, tc, title, facecolor, tc_color, tc_label) in enumerate(zip(temporal_clusters, time_thresholds, titles, facecolors, rev_colors, time_threshold_labels)):
                            _attribute_label = '{}\n{}'.format(title, tc_label)
                            _title = None
                            _tc_label = None
                            axes, _ = self.subview_cluster_duration_histograms(axes, attribute, TC, tc, facecolor, tc_color, _title, _attribute_label, _tc_label, alpha, show_time_line)
                    else:
                        show_time_line = True
                        for i, (TC, tc, title, tc_label) in enumerate(zip(temporal_clusters, time_thresholds, titles, time_threshold_labels)):
                            if i == 0:
                                _attribute_label = attribute_label[:]
                            else:
                                _attribute_label = None
                        axes, _handles = self.subview_cluster_duration_histograms(axes, attribute, TC, tc, _facecolor, tc_color, title, _attribute_label, tc_label, alpha, show_time_line)
                        try:
                            handles = [_handles[0], _handles[1]]
                            labels = [attribute_label, tc_partial_label]
                        except:
                            handles = [_handles[0]]
                            labels = [attribute_label]
                        legend_kwargs['handles'] = handles
                        legend_kwargs['labels'] = labels
                if zoom_in == True:
                    try:
                        xlim = axes.ravel()[0].get_xlim()
                        for ax in axes.ravel():
                            if attribute == 'intra_durations':
                                ax.set_xlim([0, 0.25 * xlim[-1]])
                            elif attribute == 'inter_durations':
                                ax.set_xlim([-1 * mean_time_threshold, 0.125 * xlim[-1]])
                    except:
                        xlim = axes.get_xlim()
                        if attribute == 'intra_durations':
                            axes.set_xlim([0, 0.25 * xlim[-1]])
                        elif attribute == 'inter_durations':
                            axes.set_xlim([-1 * mean_time_threshold, 0.125 * xlim[-1]])
                fig, axes = self.add_shared_labels(fig, axes, xlabel, ylabel, vecdir)
                fig.align_ylabels()
                fig.suptitle(header, fontsize=self.titlesize)
                fig.subplots_adjust(hspace=0.4, wspace=0.3, bottom=0.2)
                leg = fig.legend(**legend_kwargs)
                if search_label is None:
                    leg.set_title('All Clusters', prop={'size': self.labelsize})
                else:
                    leg.set_title(search_label, prop={'size': self.labelsize})
                leg._legend_box.align = "center"
                frame = leg.get_frame()
                frame.set_edgecolor('k')
                if save == True:
                    savename = 'cluster_{}_{}'.format(attribute, extreme_value)
                    if (zoom_in == True) and (attribute != 'intra_times'):
                        savename = '{}_zoom'.format(savename)
                    savename = '{}_{}'.format(savename, layout)
                else:
                    savename = None
                self.display_image(fig, savename=savename)

    def collect_tabular_event_specs(self, events, extreme_value, rowLabels, cellText, idx):
        """
        events:
            type <dict>

        ...

        rowLabels:
            type <list>

        cellText:
            type <list>

        idx:
            type <int>
        """
        parameters, identifiers = events['parameter'], events['identifier']
        UB = events['unbiased estimators']
        S = EventSearcher(parameters)
        indices = S.search_events(parameters=self.ref_parameter, conditions='not equal', values=UB.nan_repl, apply_to='all', modifiers=None)[1]
        nevents = np.sum(indices)
        nelapsed = parameters['elapsed'][-1] - parameters['elapsed'][0]
        if extreme_value is None:
            if idx == 0:
                rowLabels.append('Number of Events')
                rowLabels.append('Total Time Elapsed ({}s)'.format(self.timestep))
            cellText.append('{:,}'.format(nevents))
            cellText.append('{:,}'.format(nelapsed))
        else:
            specs = events['inter-exceedance'][extreme_value]
            search_kwargs = specs['search kwargs']
            IED = specs['distribution']
            nextreme = IED.events['time'].size + 1
            if idx == 0:
                search_label = self.get_search_label(search_kwargs)
                row_label = 'Number of Extreme Events\n{}'.format(search_label)
                rowLabels.append('Number of Events')
                rowLabels.append(row_label)
                rowLabels.append('Total Time Elapsed ({}s)'.format(self.timestep))
            cellText.append('{:,}'.format(nevents))
            cellText.append('{:,}'.format(nextreme))
            cellText.append('{:,}'.format(nelapsed))
        return rowLabels, cellText

    def collect_tabular_unbiased_estimators_analysis(self, events, rowLabels, cellText, idx):
        """
        events:
            type <dict>

        rowLabels:
            type <list>

        cellText:
            type <list>

        idx:
            type <int>
        """
        parameters, identifiers = events['parameter'], events['identifier']
        UB = events['unbiased estimators']
        attributes = ('alpha', 'intercept', 'theta')
        for attribute in attributes:
            label = self.available['label']['unbiased estimators'][attribute]
            if idx == 0:
                if 'theta' in label:
                    row_label = 'Extremal Index {}'.format(label)
                else:
                    row_label = 'Mean {}'.format(label)
                rowLabels.append(row_label)
            ub_value = getattr(UB, attribute)['statistics']['mean']
            rounded_ub_value = np.round(ub_value, decimals=2)
            cellText.append(r'${:.2f}$'.format(rounded_ub_value))
        for statistic, value in UB.extreme_value_estimates.items():
            if idx == 0:
                row_label = '{} Extreme {} Threshold Estimate $(${}$)$'.format(statistic.title(), self.ref_parameter.title(), self.unit_label)
                rowLabels.append(row_label)
            cellText.append(value)
        return rowLabels, cellText

    def collect_tabular_clusters_by_size_analysis(self, events, extreme_value, search_kwargs, ddof=0):
        """
        events:
            type <dict>

        extreme_value:
            type <int / float> or None

        search_kwargs:
            type <dict>

        ddof:
            type <int>

        ...
        """
        if 'inter-exceedance' not in list(events.keys()):
            raise ValueError("inter-exceedances have not yet been initialized")
        if extreme_value not in list(events['inter-exceedance'].keys()):
            raise ValueError("events of extreme value {} have not been loaded".format(extreme_value))
        if 'distribution' not in list(events['inter-exceedance'][extreme_value].keys()):
            raise ValueError("inter-exceedance distributions have not yet been initialized")
        IED = events['inter-exceedance'][extreme_value]['distribution']
        TC = TemporalClusters(IED.clusters, IED.cluster_bias, **search_kwargs)
        TC.store_relative_statistics()
        _search_kwargs = dict()
        _search_kwargs['parameters'] = 'cluster size'
        _search_kwargs['conditions'] = 'equal'
        cellText = []
        for cluster_size in TC.relative_statistics['cluster size']:
            _search_kwargs['values'] = cluster_size
            _TC = TemporalClusters(TC.clusters, TC.bias, **_search_kwargs)
            _TC.store_relative_statistics()
            _TC.store_cluster_times_and_durations()
            _nclusters = _TC.relative_statistics['number clusters'][0]
            _nevents = _TC.relative_statistics['number events'][0]
            # _nclusters = len(_elapsed_clusters)
            # _nevents = np.concatenate(_nclusters, axis=0).size
            # _elapsed_clusters = _TC.clusters['elapsed']
            f = lambda args, decimals=2 : r'${:.2f}$ $±$ ${:.2f}$'.format(np.round(np.mean(args), decimals=decimals), np.round(np.std(args, ddof=ddof), decimals=decimals))
            _intra_time_label = f(_TC.intra_times['value'])
            _intra_dur_label = f(_TC.intra_durations['value'])
            _inter_dur_label = f(_TC.inter_durations['value'])
            cellText.append(cluster_size)
            cellText.append(_nclusters)
            cellText.append(_nevents)
            cellText.append(_intra_time_label)
            cellText.append(_intra_dur_label)
            cellText.append(_inter_dur_label)
        return cellText

    def collect_tabular_cluster_analysis(self, events, extreme_value, rowLabels, cellText, idx, search_kwargs, ddof=0):
        """
        events:
            type <dict>

        extreme_value:
            type <int / float> or None

        rowLabels:
            type <list>

        cellText:
            type <list>

        idx:
            type <int>

        search_kwargs:
            type <dict>

        ddof:
            type <int>
        """
        if 'inter-exceedance' not in list(events.keys()):
            raise ValueError("inter-exceedances have not yet been initialized")
        if extreme_value not in list(events['inter-exceedance'].keys()):
            raise ValueError("events of extreme value {} have not been loaded".format(extreme_value))
        if 'distribution' not in list(events['inter-exceedance'][extreme_value].keys()):
            raise ValueError("inter-exceedance distributions have not yet been initialized")
        IED = events['inter-exceedance'][extreme_value]['distribution']
        TC = TemporalClusters(IED.clusters, IED.cluster_bias, **search_kwargs)
        TC.store_cluster_times_and_durations()
        elapsed_clusters = TC.clusters['elapsed']
        if idx == 0:
            identifiers, specs = events['identifier'], events['inter-exceedance'][extreme_value]
            rowLabels.append(r'Extremal Moment Estimator $\theta$')
            rowLabels.append(r'Time Threshold $T_C$ ({}s)'.format(self.timestep))
            if extreme_value is None:
                row_label = 'Number of Extreme Events'
            else:
                search_label = self.get_search_label(search_kwargs)
                row_label = 'Number of Extreme Events\n{}'.format(search_label)
            rowLabels.append(row_label)
            rowLabels.append('Number of Clusters')
            rowLabels.append('Mean Intra-Time ({}s)'.format(self.timestep))
            rowLabels.append('Mean Intra-Duration ({}s)'.format(self.timestep))
            rowLabels.append('Mean Inter-Duration ({}s)'.format(self.timestep))
            rowLabels.append('Mean Cluster Size')
            rowLabels.append(r'Expected Cluster Size $(\frac{1}{\theta})$')
        nextreme = np.concatenate(elapsed_clusters, axis=0).size
        nsizes = np.array([cluster.size for cluster in elapsed_clusters])
        nclusters = len(elapsed_clusters)
        theta = IED.moment_estimators[IED.cluster_bias]
        tc = IED.time_thresholds[IED.cluster_bias]
        f = lambda args : '${:.2f}$ ± ${:.2f}$'.format(np.round(np.mean(args), decimals=2), np.round(np.std(args, ddof=ddof), decimals=2))
        intra_time_label = f(np.concatenate(TC.intra_times['value'], axis=0))
        intra_duration_label = f(TC.intra_durations['value'])
        inter_duration_label = f(TC.inter_durations['value'])
        size_label = f(nsizes)
        # avg_dur = np.mean(TC.intra_durations['value'])
        # std_dur = np.std(TC.intra_durations['value'], ddof=ddof)
        # avg_size = np.mean(nsizes)
        # std_size = np.std(nsizes, ddof=ddof)
        cellText.append('${:.2f}$'.format(np.round(theta, decimals=2)))
        cellText.append(tc)
        cellText.append(nextreme)
        cellText.append(nclusters)
        cellText.append(intra_time_label)
        cellText.append(intra_duration_label)
        cellText.append(inter_duration_label)
        cellText.append(size_label)
        # cellText.append('${:.2f}$ ± ${:.2f}$'.format(np.round(avg_dur, decimals=2), np.round(std_dur, decimals=2)))
        # cellText.append('${:.2f}$ ± ${:.2f}$'.format(np.round(avg_size, decimals=2), np.round(std_size, decimals=2)))
        cellText.append('${:.2f}$'.format(np.round(1/theta, decimals=2)))
        return rowLabels, cellText

    def view_extreme_event_analysis_table(self, extreme_value=None, facecolors=('peachpuff', 'bisque'), colColours='orange', rowColours='orange', save=False, **kwargs):
        """
        ...

        facecolors:
            type <tuple / list / array>

        colColours:
            type <str>

        rowColours:
            type <str>

        save:
            type <bool>
        """
        ncols = len(self.events)
        colLabels, rowLabels, cellText = [], [], []
        for idx, events in enumerate(self.events):
            sc_label = self.get_solar_cycle_label(events)
            ex_label = self.get_extreme_parameter_label(events)
            title = '{}: {}'.format(sc_label, ex_label)
            colLabels.append(title)
            rowLabels, cellText = self.collect_tabular_event_specs(events, extreme_value, rowLabels, cellText, idx)
            rowLabels, cellText = self.collect_tabular_unbiased_estimators_analysis(events, rowLabels, cellText, idx)
        header = 'Extreme Event Analysis'
        fig, axes, table = self.subview_tabular_data(header, colLabels, rowLabels, cellText, facecolors, colColours, rowColours, loc='center', cellLoc='center', **kwargs)
        if save == True:
            savename = 'table_extreme_event_analysis_{}'.format(extreme_value)
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def view_cluster_analysis_table(self, extreme_value=None, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, facecolors=('peachpuff', 'bisque'), rowColours='orange', colColours='orange', ddof=0, save=False, **kwargs):
        """
        facecolors:
            type <tuple / list / array>

        rowColours:
            type <str>

        colColours:
            type <str>

        save:
            type <bool>
        """
        ncols = len(self.events)
        colLabels, rowLabels, cellText = [], [], []
        header = 'Cluster Analysis'
        if search_parameters is None:
            search_kwargs = dict()
        else:
            search_kwargs = dict(parameters=search_parameters, conditions=search_conditions, values=search_values, apply_to=apply_to, modifiers=modifiers)
            search_label = self.get_search_label(search_kwargs)
            header = '{}\n{}'.format(header, search_label)
        for idx, events in enumerate(self.events):
            sc_label = self.get_solar_cycle_label(events)
            ex_label = self.get_extreme_parameter_label(events)
            title = '{}: {}'.format(sc_label, ex_label)
            colLabels.append(title)
            rowLabels, cellText = self.collect_tabular_cluster_analysis(events, extreme_value, rowLabels, cellText, idx, search_kwargs, ddof)
        fig, axes, table = self.subview_tabular_data(header, colLabels, rowLabels, cellText, facecolors, colColours, rowColours, loc='center', cellLoc='center', **kwargs)
        if save == True:
            savename = 'table_cluster_analysis_{}'.format(extreme_value)
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def view_clusters_by_size_table(self, extreme_value=None, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, facecolors=('peachpuff', 'bisque'), colColours='orange', ddof=0, save=False, **kwargs):
        """
        facecolors:
            type <tuple / list / array>

        rowColours:
            type <str>

        colColours:
            type <str>

        save:
            type <bool>
        """
        ncols = len(self.events)
        colLabels, rowLabels, cellText = [], None, []
        rowColours = None
        header = 'Cluster Analysis'
        if search_parameters is None:
            search_kwargs = dict()
        else:
            search_kwargs = dict(parameters=search_parameters, conditions=search_conditions, values=search_values, apply_to=apply_to, modifiers=modifiers)
            search_label = self.get_search_label(search_kwargs)
            header = '{}\n{}'.format(header, search_label)
        colLabels = ['Cluster Size', 'Number of Clusters', 'Number of Extreme Events']
        colLabels += ['Mean Intra-Time ({}s)'.format(self.timestep)]
        colLabels += ['Mean Intra-Duration ({}s)'.format(self.timestep)]
        colLabels += ['Mean Inter-Duration ({}s)'.format(self.timestep)]
        for idx, events in enumerate(self.events):
            sc_label = self.get_solar_cycle_label(events)
            ex_label = self.get_extreme_parameter_label(events)
            title = '{}: {}'.format(sc_label, ex_label)
            cellText = self.collect_tabular_clusters_by_size_analysis(events, extreme_value, search_kwargs, ddof)
            fig, axes, table = self.subview_column_data(header, colLabels, cellText, facecolors, colColours, loc='center', cellLoc='center', **kwargs)
            # fig, axes, table = self.subview_tabular_data(header, colLabels, rowLabels, cellText, facecolors, colColours, rowColours, loc='center', cellLoc='center', **kwargs)
            if save == True:
                savename = 'table_cluster_analysis_sizes_{}_{}'.format(extreme_value, search_label.replace(' ', '_'))
            else:
                savename = None
            self.display_image(fig, savename=savename, dpi=300)

    def view_frechet_distribution(self, facecolors=('red', 'green', 'blue', 'orange', 'purple', 'black'), save=False, **kwargs):
        """
        facecolors:
            type <tuple / list / array

        save:
            type <bool>
        """
        x = np.arange(0.1, 5.006, .005)
        xticks = np.arange(0, x[-1], 0.5)
        yticks_pdf = np.arange(0, 1.3, 0.25)
        yticks_cdf = np.arange(0, 1.1, 0.25)
        alpha = np.array([1, 1, 2, 2, 3, 3])
        mu = np.zeros(alpha.size, dtype=int)
        sigma = np.array([1, 2, 1, 2, 1, 2])
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < sigma.size:
            raise ValueError("{} facecolors provided for {} varied parameters".format(ncolors, sigma.size))
        fig, axes = self.get_dim2_figure_and_axes(layout='double-horizontal', n=2, **kwargs)
        (ax_top, ax_btm) = axes
        for a, m, s, facecolor in zip(alpha, mu, sigma, facecolors):
            label = r'$\alpha = {}, \mu = {}, \sigma = {}$'.format(a, m, s)
            tmp = (x - m)/s
            ypdf = (a/s) * tmp**(-1 - a) * np.exp(-1 * (tmp ** -a))
            ycdf = np.exp(-1 * (tmp ** -a))
            ax_top.plot(x, ypdf, color=facecolor, label=label, alpha=0.5)
            ax_btm.plot(x, ycdf, color=facecolor, alpha=0.5)
            for ax in (ax_top, ax_btm):
                ax.set_xticks(xticks, minor=True)
                ax.set_xticks(xticks[::2])
                ax.set_xlim([0, 5])
            for ax, yticks in zip((ax_top, ax_btm), (yticks_pdf, yticks_cdf)):
                ax.set_yticks(yticks, minor=True)
                ax.set_yticks(yticks[::2])
                ax.set_ylim([0, yticks[-1]])
                ax.grid(color='k', alpha=0.3, linestyle=':')
        pdf_formula_label = r'$PDF(x | \alpha, \mu, \sigma) = \frac{\alpha}{\sigma} (\frac{x-\mu}{\sigma})^{-1-\alpha} e^{- (\frac{x-\mu}{\sigma})^{-\alpha}}$'
        cdf_formula_label = r'$CDF(x | \alpha, \mu, \sigma) = e^{- (\frac{x-\mu}{\sigma})^{-\alpha}}$'
        for ax, formula_label in zip(axes, (pdf_formula_label, cdf_formula_label)):
            ax.set_title(formula_label, fontsize=self.titlesize)
        ax_top.set_ylabel('PDF', fontsize=self.labelsize)
        ax_btm.set_xlabel('x', fontsize=self.labelsize)
        ax_btm.set_ylabel('CDF', fontsize=self.labelsize)
        # fig.tight_layout()
        fig.suptitle(r'Fr$\acute{e}$chet Distribution', fontsize=self.titlesize)
        leg = fig.legend(loc=7, fontsize=self.labelsize)
        fig.subplots_adjust(hspace=0.2, wspace=0.3, right=0.75)
        leg.set_title('Parameters', prop={'size': self.titlesize})
        leg._legend_box.align = "center"
        frame = leg.get_frame()
        frame.set_edgecolor('k')
        if save == True:
            savename = 'frechet_distribution'
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def view_cluster_example(self, facecolors=('darkorange', 'steelblue', 'purple'), time_unit='hour', save=False, **kwargs):
        """

        """
        time_threshold = 5
        clusters = np.array([np.array([2, 3, 5]), np.array([11, 13, 17]), np.array([22, 26])]) - 2
        ncluster = len(clusters)
        xticks = np.unique(np.sum([cluster.copy().tolist() for cluster in clusters]))
        header = 'Example: $3$ Clusters'
        intra_time, intra_duration, inter_duration = [], [], []
        intra_arrowprops = {'arrowstyle': '|-|', 'color' : 'k'}
        inter_arrowprops = {'arrowstyle': '<->', 'color' : 'gray'}
        fig, ax = self.get_dim2_figure_and_axes(layout='overlay', n=1, **kwargs)
        for ith_cluster, curr_cluster in enumerate(clusters):
            ti, tf = curr_cluster[0], curr_cluster[-1]
            intra_time = np.diff(curr_cluster)
            intra_duration = tf - ti
            ax.annotate('', xy=(tf, 0.95), xycoords='data', xytext=(ti, 0.95), textcoords='data', fontsize=self.labelsize, arrowprops=intra_arrowprops)
            ax.text(ti +1, 0.925, s=r'$\Delta T_{intra}$', fontsize=self.labelsize)
            if ith_cluster < ncluster - 1:
                next_cluster = clusters[ith_cluster+1]
                _ti, _tf = next_cluster[0], next_cluster[-1]
                inter_duration = _ti - tf
                ax.annotate('', xy=(_ti, 1.05), xycoords='data', xytext=(tf, 1.05), textcoords='data', fontsize=self.labelsize, arrowprops=inter_arrowprops)
                ax.text(_tf +2, 1.075, s=r'$\Delta T_{inter}$')
            ax.scatter(curr_cluster, np.ones(curr_cluster.size), color=facecolors[ith_cluster], label='Cluster #{}'.format(ith_cluster+1))
        ax.plot([np.nan], [np.nan], color='none', label='$T_C = {}$ {}s'.format(time_threshold, time_unit))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, fontsize=self.ticksize)
        ax.set_xlabel('Consecutive Differences of Elapsed {}s'.format(time_unit.title()), fontsize=self.labelsize)
        ax.set_yticks([0.85, 1.15])
        ax.set_yticklabels([])
        ax.set_xlim([-2, 30])
        ax.set_ylim([0.85, 1.15])
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.suptitle(header, fontsize=self.titlesize)
        fig.legend(ncol=4, bbox_transform=fig.transFigure, bbox_to_anchor=(0, 0.65, 1, 0.2), fancybox=True, loc='upper center', mode='expand', fontsize=self.labelsize)
        if save == True:
            savename = 'example_clusters'
        else:
            savename = None
        self.display_image(fig, savename=savename)





##
