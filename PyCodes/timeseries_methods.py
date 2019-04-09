import numpy as np
import datetime
from scipy.stats import sem, kurtosis, skew
from scipy.optimize import minimize, basinhopping
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from adaptive_histogram import *
from speed_optimization import *
from search_methods import *

class InterExceedance():

    """
    Purpose:
        (1) Select extreme events by specifying speed threshold.
        (2) Get time-deltas in-between consecutive extreme events.
    """

    def __init__(self, vthreshold):
        """
        vthreshold          :   type <int / float> or None
        """
        self.vthreshold = vthreshold
        self.result = {}

    def store_data(self, data):
        """
        data                :   type <dict>
        """
        if self.vthreshold is None:
            self.result['data'] = data
        elif isinstance(self.vthreshold, (float, int, np.float, np.int, np.int64)):
            Searcher = SearchEjecta(data)
            kwargs = {'search_conditions' : 'greater than', 'search_values' : self.vthreshold, 'search_parameters' : 'speed'}
            self.result['data'] = Searcher.search(**kwargs)
        else:
            raise ValueError("vthreshold = type <int / float> or None")

    @property
    def data(self):
        return self.result['data']

    @property
    def times(self):
        return np.diff(self.result['data']['elapsed'])


class InterExceedanceDecay():

    """
    Purpose:
        (1) Compute a histogram of inter-exceedance times between fast CMEs.
        (2) Compute a histogram of the inverse-transform sample data.
        (3) Compare the rate of decay between the two histograms.
    """

    def __init__(self, IE):
        """
        IE                  :   type <cls>
        """
        self.IE = IE
        self.result = {'Histogram' : {}}

    @property
    def data(self):
        return self.IE.data

    @property
    def times(self):
        return self.IE.times

    def generate_inverse_transform_sample(self):
        x = np.random.uniform(size=self.times.size)
        y_tmp = np.log(1/x)
        y = y_tmp * np.mean(self.times) / np.mean(y_tmp)
        self.result['x'] = x
        self.result['y'] = y

    @property
    def x(self):
        return self.result['x']

    @property
    def y(self):
        return self.result['y']

    def generate_histograms(self, bias='left', bin_width=15):
        """
        bias        :   type <str>
        """
        edges = np.arange(0, 301, bin_width, dtype=int)
        for key, data in zip(('inverse sample', 'time-delta'), (self.y, self.times)):
            Histogram = AdaptiveHistogram(data, threshold=None)
            Histogram.initialize_edges(n=edges, edges_by='custom')
            Histogram.initialize_observed_counts(bias)
            Histogram.initialize_midpoints()
            Histogram.initialize_bin_widths()
            Histogram.initialize_normalization()
            self.result['Histogram'][key] = Histogram

    def Histogram(self, htype):
        """
        htype       :   type <str>
        """
        available_htypes = ('inverse sample', 'time-delta')
        if htype not in available_htypes:
            raise ValueError("unknown htype: {}; available htypes: {}".format(htype, available_htypes))
        return self.result['Histogram'][htype]

class MaxSpectrum():

    """
    Purpose:
        (1) Compute D(j, k) for all j and for all k.
        (2) Compute initial Y(j) for all j.
        (3) Compute linear fit of Y(j) over sub-interval of j.
    """

    # def __init__(self, data, jinitial=4, jfinal=13):
    def __init__(self, data, jinitial=4, jfinal=12):
        """
        data        :   type <array>
        jinitial    :   type <int>
        jfinal      :   type <int>
        """
        self.data = data
        self.jndices = np.linspace(jinitial, jfinal, jfinal - jinitial + 1, dtype=int)
        # self.jndices = np.linspace(jinitial, jfinal, jfinal - jinitial + 1, dtype=int)[:-1]
        self.result = {}

    @property
    def jmax(self):
        return int(np.floor(np.log2(self.data.size)))

    @property
    def js(self):
        return np.linspace(1, self.jmax, self.jmax, dtype=int)

    @property
    def js_fit(self):
        return self.js[self.jndices]

    @staticmethod
    def get_jth_kth_index(j):
        """
        j           :   type <int>
        """
        jprime = 2**j
        return (-1, jprime)

    def get_shapeshifted_data(self, j):
        """
        j           :   type <int>
        """
        desired_size_factor = np.prod([n for n in j if n != -1])
        if -1 in j:  ## IMPLICIT ARRAY SIZE
            desired_size = self.data.size // desired_size_factor * desired_size_factor
        else:
            desired_size = desired_size_factor
        return self.data.flat[:desired_size].reshape(j)

    @staticmethod
    def get_data_maximums(vshift, pad_char=None):
        """
        vshift      :   type <array>
        pad_char    :   type <int / float> or None
        """
        res = np.max(vshift.T, axis=0)
        if pad_char is not None:
            res = res[res != pad_char]
        return res[res > 0]

    @staticmethod
    def get_error_bars(dlog, ddof=0):
        """
        dlog        :   type <array>
        ddof        :   type <int>
        """
        standard_deviation = np.std(dlog, ddof=ddof)
        standard_error = sem(dlog, ddof=ddof)
        return standard_deviation, standard_error

    def generate_djk(self, pad_char=1, ddof=0):
        """
        pad_char    :   type <int / float> or None
        ddof        :   type <int>
        """
        keys = ('dlog', 'npad', 'standard deviation', 'standard error')
        res = {key : [] for key in keys}
        for j in self.js:
            index = self.get_jth_kth_index(j)
            vshift = self.get_shapeshifted_data(index)
            vmax = self.get_data_maximums(vshift, pad_char)
            dlog = np.log2(vmax)
            npad = dlog.size
            st_dev, st_err = self.get_error_bars(dlog, ddof)
            args = (dlog, npad, st_dev, st_err)
            for key, value in zip(keys, args):
                res[key].append(value)
        for key in keys:
            res[key] = np.array(res[key])
        self.result.update(res)

    @property
    def number_of_nonpads(self):
        return self.result['npad']

    @property
    def dlog(self):
        return self.result['dlog']

    def D(self, j, k):
        """
        j           :   type <int> or None
        k           :   type <int> or None
        """
        if j is None:
            if k is None:
                return self.dlog
            elif isinstance(k, int):
                return self.dlog.T[k] ## NOT res[:, k] BECAUSE ARRAY IS RAGGED
            else:
                raise ValueError("k = type <int> or None")
        elif isinstance(j, int):
            if k is None:
                return self.dlog[j]
            elif isinstance(k, int):
                return self.dlog[j][k] ## NOT res[j, k] BECAUSE ARRAY IS RAGGED
            else:
                raise ValueError("k = type <int> or None")
        else:
            raise ValueError("j = type <int> or None")

    def generate_yj_initial(self):
        res = []
        for dlog, npad in zip(self.dlog, self.number_of_nonpads):
            yj_init = np.sum(dlog) / npad
            res.append(yj_init)
        self.result['yj initial'] = np.array(res)

    @property
    def yj_init(self):
        return self.result['yj initial']

    def generate_yj_fitted(self, method='Nelder-Mead', scale='local', **kwargs):
        """
        method      :   type <str>
        scale       :   type <str>
        **kwargs    :   type <dict>; pass to scipy minimize
        """
        GLS = GeneralizedLeastSquares(self.js_fit, self.yj_init[self.jndices])
        wts = self.number_of_nonpads[self.jndices]
        res = GLS.fit(wts=wts, method=method, scale=scale, **kwargs)
        prms = np.array([1/res.x[0], res.x[1]])
        # yj_fit = GLS.f_linear(prms)
        yj_fit = self.js_fit / prms[0] + prms[1]
        self.result['parameters'] = prms
        self.result['yj optimized'] = yj_fit

    @property
    def alpha(self):
        return self.result['parameters'][0]

    @property
    def intercept(self):
        return self.result['parameters'][1]

    @property
    def yj_fit(self):
        return self.result['yj optimized']

    def Y(self, j=None, fit=True):
        """
        j           :   type <int> or None
        fit         :   type <bool>
        """
        if fit is True:
            ydata = self.yj_fit
        elif fit is False:
            ydata = self.yj_init
        else:
            raise ValueError("fit = True or False")
        if j is None:
            return ydata
        elif isinstance(j, int):
            index = j - self.jndices[0]
            return ydata[index]
        else:
            raise ValueError("j = type <int> or None")

class UnbiasedEstimators():

    """
    Purpose:
        (1) Compute max spectrum of data.
        (2) Compute multiple resamples of the data and compute corresponding
            max spectra.
        (3) Collect statistics (such as alpha hat) from multiple resamples.
    """

    def __init__(self, data, nresamples=100, nshuffles=10, hat_func=np.mean):
        """
        data                :   type <array>
        nresamples          :   type <int>
        nshuffles           :   type <int>
        hat_func            :   type <function>
        """
        self.data = data
        self.nresamples = nresamples
        self.nshuffles = nshuffles
        self.indices = np.array(list(range(data.size)), dtype=int)
        self.f_alpha = {'mean' : np.mean, 'median' : np.median, 'hat' : hat_func}
        self.storage = {}
        self.view_keys = ('optimized fit', 'alpha histogram')

    @property
    def original_indices(self):
        return np.argsort(self.indices)

    @staticmethod
    def compute_max_spectrum(data, pad_char=1, ddof=0, method='Nelder-Mead', scale='local', **kwargs):
        """
        data                :   type <array>
        pad_char            :   type <int / float> or None
        ddof                :   type <int>
        method              :   type <str>
        scale               :   type <str>
        **kwargs            :   type <dict>; pass to scipy minimize
        """
        MS = MaxSpectrum(data)
        MS.generate_djk(pad_char, ddof)
        MS.generate_yj_initial()
        MS.generate_yj_fitted(method, scale, **kwargs)
        return MS

    def resample_indices(self):
        if self.nshuffles > 0:
            for n in range(self.nshuffles):
                np.random.shuffle(self.indices)
        np.random.shuffle(self.indices)

    def store_max_spectra(self):
        MS_orig = self.compute_max_spectrum(self.data)
        res = []
        for n in range(self.nresamples):
            self.resample_indices()
            MS = self.compute_max_spectrum(self.data[self.indices])
            res.append(MS)
        res.insert(0, MS_orig)
        self.storage['max spectra'] = np.array(res)

    @property
    def max_spectra(self):
        return self.storage['max spectra']

    def select_ith_resample(self, i=0):
        """
        i                   :   type <int>
        """
        if not isinstance(i, int):
            raise ValueError("i = type <int>")
        MS = self.max_spectra[i]
        return MS

    def D(self, j, k, i=0):
        """
        j                   :   type <int> or None
        k                   :   type <int> or None
        i                   :   type <int>
        """
        MS = self.select_ith_resample(i)
        return MS.D(j, k)

    def Y(self, j=None, fit=True, i=0):
        """
        j                   :   type <int> or None
        fit                 :   type <bool>
        i                   :   type <int>
        """
        MS = self.select_ith_resample(i)
        return MS.Y(j, fit)

    def number_of_nonpads(self, i=0):
        """
        i                   :   type <int>
        """
        MS = self.select_ith_resample(i)
        return MS.number_of_nonpads

    def error(self, err, i=0):
        """
        err                 :   type <str>
        i                   :   type <int>
        """
        available_err = ('standard deviation', 'standard error')
        if err not in available_err:
            raise ValueError("unknown err: {}; available err: {}".format(err, available_err))
        MS = self.select_ith_resample(i)
        return MS.result[err]

    @staticmethod
    def select_initial_index(include_original):
        """
        include_original    :   type <bool>
        """
        if not isinstance(include_original, bool):
            raise ValueError("include_original = True or False")
        if include_original is True:
            lowest_index = 0
        else:
            lowest_index = 1
        return lowest_index

    def alpha(self, i='hat', include_original=True):
        """
        i                   :   type <int> or type <str> or None
        include_original    :   type <bool>
        """
        if i is None:
            res = []
            lowest_index = self.select_initial_index(include_original)
            for i in range(lowest_index, self.nresamples+1):
                MS = self.select_ith_resample(i)
                res.append(MS.alpha)
            return np.array(res)
        elif isinstance(i, int):
            MS = self.select_ith_resample(i)
            return MS.alpha
        elif isinstance(i, str):
            if i in list(self.f_alpha.keys()):
                f = self.f_alpha[i]
            elif i == 'skew':
                f = lambda arg : skew(arg)
            elif i == 'kurtosis':
                f = lambda arg : kurtosis(arg)
            else:
                available_i = tuple(list(self.f_alpha.keys()) + ['skew', 'kurtosis'])
                raise ValueError("unknown i: {}; available i: {}".format(i, available_i))
            # if i not in list(self.f_alpha.keys()):
            #     raise ValueError("unknown i: {}; available i: {}".format(i, tuple(list(self.f_alpha.keys()))))
            # f = self.f_alpha[i]
            lowest_index = self.select_initial_index(include_original)
            res = []
            for i in range(lowest_index, self.nresamples+1):
                MS = self.select_ith_resample(i)
                res.append(MS.alpha)
            return f(np.array(res))
        else:
            raise ValueError("i = 'type <int> or type <str> or None")

    def intercept(self, i='hat', include_original=True):
        """
        i                   :   type <int> or type <str> or None
        include_original    :   type <bool>
        """
        if i is None:
            res = []
            lowest_index = self.select_initial_index(include_original)
            for i in range(lowest_index, self.nresamples+1):
                MS = self.select_ith_resample(i)
                res.append(MS.intercept)
            return np.array(res)
        elif isinstance(i, int):
            MS = self.select_ith_resample(i)
            return MS.intercept
        elif isinstance(i, str):
            if i in list(self.f_alpha.keys()):
                f = self.f_alpha[i]
            elif i == 'skew':
                f = lambda arg : skew(arg)
            elif i == 'kurtosis':
                f = lambda arg : kurtosis(arg)
            else:
                available_i = tuple(list(self.f_alpha.keys()) + ['skew', 'kurtosis'])
                raise ValueError("unknown i: {}; available i: {}".format(i, available_i))
            # if i not in list(self.f_alpha.keys()):
            #     raise ValueError("unknown i: {}; available i: {}".format(i, tuple(list(self.f_alpha.keys()))))
            # f = self.f_alpha[i]
            lowest_index = self.select_initial_index(include_original)
            res = []
            for i in range(lowest_index, self.nresamples+1):
                MS = self.select_ith_resample(i)
                res.append(MS.intercept)
            return f(np.array(res))
        else:
            raise ValueError("i = 'type <int> or type <str> or None")

    def store_alpha_histogram(self, include_original=True, bias='left'):
        """
        n                   :   type <int / float> or type <array>
        edges_by            :   type <str>
        include_original    :   type <bool>
        bias                :   type <str>
        """
        edges = np.arange(3.2, 4.01, 0.02)
        data = self.alpha(i=None, include_original=include_original)
        Histogram = AdaptiveHistogram(data, threshold=None)
        Histogram.initialize_edges(n=edges, edges_by='custom')
        Histogram.initialize_observed_counts(bias)
        Histogram.initialize_midpoints()
        Histogram.initialize_bin_widths()
        Histogram.initialize_normalization()
        self.storage['Histogram'] = Histogram

    @property
    def Histogram(self):
        return self.storage['Histogram']

    def store_speed_thresholds(self, include_original=True):
        """ """
        if include_original is True:
            lower_index = 0
        elif include_original is False:
            lower_index = 1
        else:
            raise ValueError("include_original = True or False")
        tmp = []
        for i in range(lower_index, self.nresamples):
            MS = self.select_ith_resample(i)
            MS.yj_fit
            tmp.append(MS.yj_fit)
        tmp = np.array(tmp)
        ymean = np.mean(tmp, axis=0)
        ymedian = np.median(tmp, axis=0)
        f = lambda arg : int(round(2**min(arg)))
        self.storage['speed threshold'] = {'mean' : f(ymean), 'median' : f(ymedian)}

    @property
    def speed_thresholds(self):
        """ """
        return self.storage['speed threshold']

class ExtremalIndex():

    """

    """

    def __init__(self, UB):
        """ """
        self.UB = UB
        self.storage = {}

    @property
    def js(self):
        """ """
        MS = self.UB.select_ith_resample(i=0)
        return MS.js_fit

    def store_theta_hat(self):
        """ """
        MS_orig = self.UB.select_ith_resample(i=0)
        res = []
        for i in range(1, self.UB.nresamples):
            MS = self.UB.select_ith_resample(i)
            power = MS.alpha * (MS_orig.yj_fit - MS.yj_fit).T
            res.append(power)
        self.storage['theta hat'] = 2**np.array(res)

    @property
    def theta_hat(self):
        """ """
        return self.storage['theta hat']

    def store_point_estimators(self):
        """ """
        keys = ('mean', 'median', 'minimum', 'maximum')
        f_stats = (np.mean, np.median, np.min, np.max)
        res = {key : f(self.theta_hat, axis=0) for key, f in zip(keys, f_stats)}
        self.storage['point estimators'] = res

    @property
    def point_estimators(self):
        """ """
        return self.storage['point estimators']

    def store_theta_histogram(self, n, edges_by, bias='left'):
        """ """
        data = self.theta_hat.copy().reshape(-1)
        Histogram = AdaptiveHistogram(data, threshold=None)
        Histogram.initialize_edges(n, edges_by)
        Histogram.initialize_observed_counts(bias)
        Histogram.initialize_midpoints()
        Histogram.initialize_bin_widths()
        Histogram.initialize_normalization()
        self.storage['Histogram'] = Histogram

    @property
    def Histogram(self):
        """ """
        return self.storage['Histogram']

class ClusterDurationHistograms():

    def __init__(self, clusters, time_threshold):
        self.clusters = clusters
        self.time_threshold = time_threshold
        self.available_htypes = ('intra time', 'intra duration', 'inter duration')

    @staticmethod
    def flatten_data(data):
        """
        data                :   type <array>
        """
        return np.concatenate(data, axis=0)

    @property
    def intra_time(self):
        diffs = np.array([np.diff(cluster).tolist() for cluster in self.clusters['elapsed']])
        return self.flatten_data(diffs)

    @property
    def intra_duration(self):
        data = self.clusters['elapsed']
        res = np.array([subarr[-1] - subarr[0] for subarr in data])
        return res

    @property
    def inter_duration(self):
        data = self.clusters['elapsed']
        return np.array([data[idx][0] - data[idx-1][-1] for idx in range(1, len(data))])

    @staticmethod
    def get_intra_bins(delta):
        """
        delta                   :   type <array>
        """
        return np.linspace(0, delta, delta+1, dtype=int)

    @staticmethod
    def get_inter_bins(delta, method='freedman-diaconis'):
        """
        delta                   :   type <array>
        method                  :   type <str>
        """
        if method == 'freedman-diaconis':
            iqr_bounds = np.percentile(delta, np.array([25, 75]))
            iqr_delta = max(iqr_bounds) - min(iqr_bounds)
            wbin = (iqr_delta * 2) // (len(delta))**(1/3)
            lbin = np.floor(np.amin(delta))
            rbin = np.ceil(np.amax(delta))
            nbin = 2*int((rbin - lbin) // wbin)
            res = np.linspace(lbin, rbin, nbin, dtype=int)
        else:
            raise ValueError("not yet implemented")
        return res

    @staticmethod
    def histogram_subroutine(data, edges, bias='left'):
        """ """
        Histogram = AdaptiveHistogram(data, threshold=None)
        Histogram.initialize_edges(edges, edges_by='custom')
        Histogram.initialize_observed_counts(bias)
        Histogram.initialize_midpoints()
        Histogram.initialize_bin_widths()
        Histogram.initialize_normalization()
        return Histogram

    @property
    def intra_time_histogram(self):
        """ """
        peak = np.max(self.intra_time)
        edges = self.get_intra_bins(peak)
        Histogram = self.histogram_subroutine(self.intra_time, edges)
        return Histogram

    @property
    def intra_duration_histogram(self):
        """ """
        peak = np.max(self.intra_duration)
        edges = self.get_intra_bins(peak)
        Histogram = self.histogram_subroutine(self.intra_duration, edges)
        return Histogram

    @property
    def inter_duration_histogram(self):
        """ """
        edges = self.get_inter_bins(self.inter_duration, method='freedman-diaconis')
        Histogram = self.histogram_subroutine(self.inter_duration, edges)
        return Histogram

class TemporalClusters():

    """
    Purpose:

        (1) calculate the moment estimator of the extremal index
        (2) calculate the time threshold that separates intra-
            clusters from inter-clusters
        (3) organize CMEs into CME clusters
        (4) search CME clusters by condition(s)
    """

    def __init__(self, IE, bias='threshold'):
        """
        IE                          :   type <cls>
        bias                        :   type <str>
        """
        self.IE = IE
        self.bias = bias
        self.result = {}
        self.histograms = {}
        self.relative_statistics = {}

    @property
    def data(self):
        return self.IE.data

    @property
    def interexceedance_times(self):
        return self.IE.times

    def initialize_extremal_moment_estimator(self, initial_theta=None):
        """
        initial_theta               :   type <int / float> or None
        """
        if initial_theta is None:
            if isinstance(self.bias, str):
                bias_keys = ('threshold', 'first order')
                if self.bias not in bias_keys:
                    raise ValueError("unknown bias: {}; available bias: {}".format(self.bias, bias_keys))
                if self.bias == 'threshold':
                    self.result['extremal moment estimator'] = 2 * np.sum(self.interexceedance_times)**2 / (self.interexceedance_times.size * np.sum(self.interexceedance_times**2))
                else:
                    self.result['extremal moment estimator'] = 2 * np.sum(self.interexceedance_times - 1)**2 / (self.interexceedance_times.size * np.sum((self.interexceedance_times - 1) * (self.interexceedance_times - 2)))
            else:
                raise ValueError("bias = type <str>")
        else:
            if isinstance(initial_theta, (int, float)):
                self.result['extremal moment estimator'] = initial_theta
            else:
                raise ValueError("initial_theta = type <int / float> or None")

    @property
    def extremal_moment_estimator(self):
        return self.result['extremal moment estimator']

    def initialize_time_threshold(self, initial_tc=None):
        """
        initial_tc              :   type <int / float> or None
        """
        if initial_tc is None:
            values = np.sort(self.interexceedance_times)
            size = values.size + 1
            index = int(np.floor(self.extremal_moment_estimator * size))
            if index == values.size:
                index -= 1
            # index = int(np.round(self.extremal_moment_estimator * size))
            self.result['time threshold'] = values[index]
        else:
            self.result['time threshold'] = initial_tc

    @property
    def time_threshold(self):
        return self.result['time threshold']

    @property
    def cluster_indices(self):
        condition = (self.time_threshold < self.interexceedance_times)
        return np.where(condition)[0] + 1

    def initialize_clusters(self):
        indices = self.cluster_indices
        self.result['clusters'] = {key : np.array(np.split(value, indices)) for key, value in self.data.items()}

    @property
    def clusters(self):
        return self.result['clusters']

    def search(self, search_parameters, search_conditions, search_values, fkeys=None, load_keys=('ith cluster', 'cluster size'), apply_to='all', view_keys=None, use_abs=False, use_halo=None):
        """
        search_parameters       :   type <str> or type <tuple / list / array>
        search_conditions       :   type <str> or type <tuple / list / array>
        search_values           :   type <int / float / str> or type <tuple / list / array>
        fkeys                   :   None or type <tuple / list / array> or type <str>
        load_keys               :   None or type <tuple / list / array>
        apply_to                :   type <str>
        view_keys               :   None or type <str> or type <tuple / list / array>
        use_abs                 :   type <tuple / list / array> or type <bool>
        """
        SC = SearchClusters(self.clusters, self.cluster_indices)
        res = SC.search(search_parameters, search_conditions, search_values, fkeys, load_keys, apply_to, view_keys, use_abs, use_halo)
        return res

    def get_cluster_histograms(self, cluster_size='all', search_conditions='exact match', **search_kwargs):
        """
        cluster_size            :   'all' or type <int>
        search_conditions       :   type <str>
        search_kwargs           :   type <dict>
        """
        load_keys = ('ith cluster', 'cluster size')
        if cluster_size == 'all':
            # data = self.clusters
            data = self.search(search_parameters='cluster size', search_conditions='greater than', search_values=1, load_keys=load_keys)
        else:
            data = self.search(search_parameters='cluster size', search_conditions=search_conditions, search_values=cluster_size, load_keys=load_keys, **search_kwargs)
        return ClusterDurationHistograms(data, self.time_threshold)

    def initialize_relative_statistics(self, skip_lone_ejecta=True):
        """
        skip_lone_ejecta        :   type <bool>
        """
        if skip_lone_ejecta is True:
            smin = 2
        elif skip_lone_ejecta is False:
            smin = 1
        else:
            raise ValueError("skip_lone_ejecta = True or False")
        key = list(self.clusters.keys())[0]
        rep_sizes = np.array([subarr.size for subarr in self.clusters[key] if subarr.size >= smin])
        unq_sizes, unq_counts = np.unique(rep_sizes, return_counts=True)
        res = {'cluster size' : unq_sizes, 'number of clusters' : unq_counts}
        res['number of ejecta'] = unq_sizes * unq_counts
        f_prob = lambda arg : arg / np.sum(arg)
        res['probability'] = f_prob(res['number of ejecta'])
        self.relative_statistics.update(res)

class SolarCorrelations():

    """
    Purpose:
        (1) Store temporal sunspot frequencies.
        (2) Compare temporal ejection frequencies by condition.
        (3) Compare temporal cluster frequencies by condition.
    """

    def __init__(self, path):
        """
        path                    :   type <str>
        """
        self.path = path
        self.base_keys = ('datetime number', 'count', 'year', 'month', 'day')
        self.result = {}

    @property
    def ymapping(self):
        """ """
        return {'daily' : dict(delta=25, nearest=100), 'monthly' : dict(delta=500, nearest=10**3), 'yearly' : dict(delta=5*10**3, nearest=10**4)}

    @staticmethod
    def get_key_from_counts_by(counts_by):
        """ """
        if counts_by in ('monthly', 'yearly'):
            key = counts_by[:].replace('ly', '')
        elif counts_by == 'daily':
            key = 'day'
        else:
            raise ValueError("counts_by = 'monthly', 'yearly', 'daily'")
        return key

    @staticmethod
    def get_datetime_number_from_string(dt):
        """ """
        return date2num(datetime.datetime.strptime(dt, "%Y/%m/%d"))

    def to_daily(self, data, counts_by):
        """ """
        res = {key : [] for key in self.base_keys}
        indices = np.where((np.diff(data['day']) != 0) | (np.diff(data['month']) > 0))[0] + 1
        mod_data = {key : np.split(value, indices) for key, value in data.items() if key in self.base_keys}
        for ith_subarray in range(len(mod_data['year'])):
            arr_yy, arr_mm, arr_dd = mod_data['year'][ith_subarray], mod_data['month'][ith_subarray], mod_data['day'][ith_subarray]
            yy, mm, dd = arr_yy[0], arr_mm[0], arr_dd[0]
            dt = self.get_datetime_number_from_string('{}/{}/{}'.format(yy, mm, dd))
            cc = arr_yy.size
            args = (dt, cc, yy, mm, dd)
            for key, arg in zip(self.base_keys, args):
                res[key].append(arg)
        return {key : np.array(value) for key, value in res.items()}

    def from_daily(self, data, counts_by):
        """ """
        res = {key : [] for key in self.base_keys}
        key = self.get_key_from_counts_by(counts_by)
        indices = np.where(np.diff(data[key]) != 0)[0] + 1
        mod_data = {key : np.split(value, indices) for key, value in data.items() if key in self.base_keys}
        for ith_subarray in range(len(mod_data['year'])):
            arr_yy, arr_mm, arr_dd = mod_data['year'][ith_subarray], mod_data['month'][ith_subarray], mod_data['day'][ith_subarray]
            yy, mm, dd = arr_yy[0], arr_mm[0], arr_dd[0]
            dt = self.get_datetime_number_from_string('{}/{}/{}'.format(yy, mm, dd))
            cc = np.sum(mod_data['count'][ith_subarray])
            args = (dt, cc, yy, mm, dd)
            for key, arg in zip(self.base_keys, args):
                res[key].append(arg)
        return {key : np.array(value) for key, value in res.items()}

    def initialize_daily_sunspots(self):
        """ """
        year, month, day, dcount = np.loadtxt(self.path, skiprows=0, delimiter=None, unpack=True, usecols=(0, 1, 2, 4), dtype=int)
        dcount[dcount < 0] = 0
        daily = {'year' : year, 'month' : month, 'day' : day, 'count' : dcount}
        daily['datetime number'] = np.array([self.get_datetime_number_from_string('{}/{}/{}'.format(yy, mm, dd)) for yy, mm, dd in zip(year, month, day)])
        self.result.update({'daily' : daily})

    def initialize_monthly_sunspots(self):
        """ """
        monthly = self.from_daily(self.result['daily'], counts_by='monthly')
        self.result.update({'monthly' : monthly})

    def initialize_yearly_sunspots(self):
        """ """
        yearly = self.from_daily(self.result['daily'], counts_by='yearly')
        self.result.update({'yearly' : yearly})

    def subselect_sunspots(self, counts_by, TimeSeries=None):
        """ """
        if TimeSeries is None:
            return self.result[counts_by]
        else:
            yi, yf = TimeSeries.sc_data['original data']['year'][0], TimeSeries.sc_data['original data']['year'][-1]
            year = self.result[counts_by]['year']
            indices = np.where((year >= yi) & (year <= yf))[0]
            return {key : value[indices] for key, value in self.result[counts_by].items()}

    def subselect_ejecta(self, counts_by, TimeSeries, search_parameters=None, search_conditions=None, search_values=None, **search_kwargs):
        """ """
        Searcher = SearchEjecta(TimeSeries.sc_data['original data'])
        if search_parameters is None:
            data = Searcher.data
        else:
            data = Searcher.search(search_parameters, search_conditions, search_values, **search_kwargs)
        daily = self.to_daily(data, counts_by)
        if counts_by == 'daily':
            return daily
        elif counts_by in ('monthly', 'yearly'):
            return self.from_daily(daily, counts_by)
        else:
            raise ValueError("counts_by = 'monthly', 'yearly', or 'daily'")

class TimeSeriesAnalysis():

    """
    Purpose:
        (1) can run any/all of the following methods:
            *-  inter-exceedance decay comparison
            *-  unbiased estimators
            *-  extremal index (*)
            *-  temporal clustering
        (2) search CMEs
        (3) search CME clusters
    """

    def __init__(self, sc_data, keys):
        self.sc_data = sc_data
        self.keys = keys
        self.storage = {}
        self.result = {}

    @property
    def identity(self):
        """ """
        return self.sc_data['identity']

    @property
    def available_keys(self):
        """ """
        return ('inter-exceedance decay', 'unbiased estimators', 'extremal index', 'temporal clustering', 'cluster thresholds', 'relative cluster statistics', 'sunspot correlation', 'speed distribution fit')

    def store_sunspot_correlator(self, path):
        """ """
        Correlator = SolarCorrelations(path)
        Correlator.initialize_daily_sunspots()
        Correlator.initialize_monthly_sunspots()
        Correlator.initialize_yearly_sunspots()
        self.storage['SunSpot'] = Correlator

    @property
    def SunSpotCorrelator(self):
        """ """
        return self.storage['SunSpot']

    def store_speed_distribution_fit(self, n, edges_by, bin_threshold=5, vmin=20, method='Nelder-Mead', scale='local', minimizer_kwargs=None, htype='threshold', bias='left', **kwargs):
        """ """
        errs = ('maximum likelihood estimation', 'minimum gtest estimation') #, 'minimum chi square estimation')
        vs = self.sc_data['original data']['speed']
        if vmin is not None:
            vs = vs[vs >= vmin]
        SD = SpeedDistribution(vs)
        SD.initialize_preliminary_histogram(n, edges_by, bias)
        SD.initialize_thresholded_histogram(n, edges_by, bin_threshold, bias)
        SD.fit(errs, method, scale, minimizer_kwargs, htype, **kwargs)
        self.storage['Speed Distribution'] = SD
        self.storage['Speed Error-Space'] = {}
        for err in errs:
            EDM = ExtraDimensionalMapper(SD)
            # EDM.initialize_dim2_grid(err, x=None, y=None, xfrac=0.05, xn=21, xspace_by='number', yfrac=0.1, yn=21, yspace_by='number') ## ELLIPTICAL CONTOURS
            EDM.initialize_dim2_grid(err, x=None, y=None, xfrac=0.05, xn=21, xspace_by='number', yfrac=0.2, yn=21, yspace_by='number')
            EDM.initialize_dim3_grid(err, htype)
            self.storage['Speed Error-Space'][err] = EDM

    @property
    def SpeedDistribution(self):
        """ """
        return self.storage['Speed Distribution']

    def store_interexceedance_instance(self, vthreshold):
        """
        vthreshold          :   type <int / float>
        """
        IE = InterExceedance(vthreshold)
        IE.store_data(self.sc_data['padded data'])
        self.storage['InterExceedance'] = IE

    @property
    def InterExceedance(self):
        return self.storage['InterExceedance']

    def store_interexceedance_decays(self):
        IED = InterExceedanceDecay(self.InterExceedance)
        IED.generate_inverse_transform_sample()
        IED.generate_histograms()
        self.result['inter-exceedance decay'] = IED

    @property
    def InterExceedanceDecay(self):
        return self.result['inter-exceedance decay']

    def store_unbiased_estimators(self, nresamples):
        """
        nresamples          :   type <int>
        """
        UB = UnbiasedEstimators(self.sc_data['padded data']['speed'], nresamples=nresamples)
        UB.store_max_spectra()
        UB.store_alpha_histogram()
        UB.store_speed_thresholds()
        self.result['unbiased estimators'] = UB

    @property
    def UnbiasedEstimators(self):
        return self.result['unbiased estimators']

    def store_extremal_index(self, n, edges_by, bias='left'):
        """ """
        if 'unbiased estimators' not in list(self.result.keys()):
            raise ValueError("must run 'unbiased estimators' before running 'extremal index'")
        EI = ExtremalIndex(self.result['unbiased estimators'])
        EI.store_theta_hat()
        EI.store_point_estimators()
        EI.store_theta_histogram(n, edges_by, bias)
        self.result['extremal index'] = EI

    @property
    def ExtremalIndex(self):
        """ """
        return self.result['extremal index']

    def store_temporal_clusters(self, initial_theta=None, initial_tc=None):
        """
        initial_theta       :   None or type <int / float>
        initial_tc          :   None or type <int / float>
        """
        TC = TemporalClusters(self.InterExceedance, bias='threshold')
        TC.initialize_extremal_moment_estimator(initial_theta)
        TC.initialize_time_threshold(initial_tc)
        TC.initialize_clusters()
        self.result['temporal clustering'] = TC

    @property
    def TemporalClustering(self):
        return self.result['temporal clustering']

    @property
    def clusters(self):
        return self.TemporalClustering.clusters

    @property
    def declustering_time(self):
        return self.TemporalClustering.time_threshold

    @property
    def extremal_moment(self):
        return self.TemporalClustering.extremal_moment_estimator

    def store_cluster_thresholds(self, vmin=500, vmax=1500, vn=10):
        """ """
        vthresholds = np.arange(vmin, vmax+vn-1, vn, dtype=int)
        res = {}
        for bias in ('threshold', 'baseline', 'first order'):
            res[bias] = {'moment estimator' : [], 'time threshold' : []}
            for vthreshold in vthresholds:
                IE = InterExceedance(vthreshold)
                IE.store_data(self.sc_data['padded data'])
                if bias == 'baseline':
                    initial_theta = 0.5
                else:
                    initial_theta = None
                TemporalData = TemporalClusters(IE, bias)
                TemporalData.initialize_extremal_moment_estimator(initial_theta)
                TemporalData.initialize_time_threshold(initial_tc=None)
                res[bias]['moment estimator'].append(TemporalData.extremal_moment_estimator)
                res[bias]['time threshold'].append(TemporalData.time_threshold)
            for key, value in res[bias].items():
                res[bias][key] = np.array(value)
            res[bias]['vthreshold'] = vthresholds
        self.storage['estimates by speed'] = res

    @property
    def cluster_thresholds_by_speed(self):
        """ """
        return self.storage['estimates by speed']

    def store_relative_cluster_statistics(self, skip_lone_ejecta=True):
        """ """
        self.TemporalClustering.initialize_relative_statistics(skip_lone_ejecta)
        self.storage['relative cluster statistics'] = self.TemporalClustering.relative_statistics

    @property
    def relative_cluster_statistics(self):
        """ """
        return self.storage['relative cluster statistics']

    # def search_cmes(self, search_parameters, search_conditions, search_values):
    #     """
    #     search_parameters   :   type <tuple / list / array> or type <str>
    #     search_conditions   :   type <tuple / list / array> or type <str>
    #     search_values       :   type <tuple / list / array> or type <int / float>
    #     """
    #     raise ValueError("not yet implemented")

    def search_cme_clusters(self, clusters, search_parameters, search_conditions, search_values, fkeys=None, load_keys=('ith cluster', 'cluster size'), apply_to='all', view_keys=None, use_abs=False, use_halo=None):
        """
        search_parameters   :   type <tuple / list / array> or type <str>
        search_conditions   :   type <tuple / list / array> or type <str>
        search_values       :   type <tuple / list / array> or type <int / float>
        """
        Searcher = SearchClusters(clusters, np.cumsum([arr.size for arr in clusters['speed']]))
        return Searcher.search(search_parameters, search_conditions, search_values, fkeys, load_keys, apply_to, view_keys, use_abs, use_halo)

    def get_cluster_duration_histogram_instance(self, clusters):
        """ """
        CDH = ClusterDurationHistograms(clusters, self.declustering_time)
        return CDH

    def autocorrect_keys(self):
        if isinstance(self.keys, str):
            self.keys = [self.keys]
        elif not isinstance(self.keys, (tuple, list, np.ndarray)):
            raise ValueError("keys = type <tuple / list / array> or type <str>")

    def run(self, vthreshold, unbiased_estimator_kw={}, extremal_index_kw={}, cluster_statistic_kw={}, speed_distribution_kw={}, sunspot_kw={}):
        """ """
        self.autocorrect_keys()
        self.storage['vthreshold'] = vthreshold
        self.store_interexceedance_instance(vthreshold)
        for key in self.keys:
            if key not in self.available_keys:
                raise ValueError("unknown key: {}; available keys: {}".format(key, self.available_keys))
            if key == 'inter-exceedance decay':
                self.store_interexceedance_decays()
            elif key == 'unbiased estimators':
                self.store_unbiased_estimators(**unbiased_estimator_kw)
            elif key == 'extremal index':
                self.store_extremal_index(**extremal_index_kw)
            elif key == 'temporal clustering':
                self.store_temporal_clusters()
            elif key == 'cluster thresholds':
                self.store_cluster_thresholds()
            elif key == 'relative cluster statistics':
                self.store_relative_cluster_statistics(**cluster_statistic_kw)
            elif key == 'speed distribution fit':
                self.store_speed_distribution_fit(**speed_distribution_kw)
            elif key == 'sunspot correlation':
                self.store_sunspot_correlator(**sunspot_kw)


























##
