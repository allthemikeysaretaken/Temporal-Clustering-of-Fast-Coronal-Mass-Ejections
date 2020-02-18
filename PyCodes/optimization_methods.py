import numpy as np
from scipy.stats import kurtosis, skew
from scipy.integrate import quad
from scipy.optimize import minimize, basinhopping, OptimizeResult

class NormalDistribution():

    def __init__(self):
        self.nparams = 2

    def get_statistics(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        result = dict()
        result['mean'] = prms[0]
        result['standard deviation'] = prms[1]
        return result

    @staticmethod
    def pdf(prms, x):
        """
        prms:
            type <tuple / list / array>

        x:
            type <array>
        """
        return np.exp(- np.square((x - prms[0])/prms[1]) / 2) / (prms[1] * np.sqrt(2 * np.pi))

    @staticmethod
    def get_parameter_guess(x, **kwargs):
        """
        x:
            type <array>
        """
        mu, sigma = np.mean(x), np.std(x, **kwargs)
        return (mu, sigma)

    def get_kwargs(self):
        result = {}
        result['get_statistics'] = self.get_statistics
        result['pdf'] = self.pdf
        result['get_parameter_guess'] = self.get_parameter_guess
        result['nparams'] = self.nparams
        return result

class LogNormalDistribution():

    def __init__(self):
        self.nparams = 2

    @staticmethod
    def from_normal_to_lognormal(prms):
        """
        prms:
            type <tuple / list / array>
        """
        mu, sigma = prms
        mu_mod = np.exp(mu + sigma**2 /2)
        sigma_mod = np.sqrt(np.exp(2*mu + sigma**2) * (np.exp(sigma**2) - 1))
        return (mu_mod, sigma_mod)

    @staticmethod
    def from_lognormal_to_normal(prms):
        """
        prms:
            type <tuple / list / array>
        """
        mu, sigma = prms
        variance = sigma**2
        tmp = variance / mu**2 + 1
        mu_mod = np.log(mu/np.sqrt(tmp))
        sigma_mod = np.sqrt(np.log(tmp))
        return (mu_mod, sigma_mod)

    def get_statistics(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        (mu, sigma) = self.from_normal_to_lognormal(prms)
        result = dict()
        result['mean'] = mu
        result['standard deviation'] = sigma
        result['median'] = np.sqrt((np.exp(prms[1]**2 -1)) * (np.exp(2*prms[0] + prms[1]**2)))
        result['mode'] = np.exp(prms[0] - prms[1]**2)
        return result

    @staticmethod
    def pdf(prms, x):

        """
        prms:
            type <tuple / list / array>

        x:
            type <array>
        """
        return np.exp(- np.square((np.log(x) - prms[0]) / prms[1]) / 2) / (x * prms[1] * np.sqrt(2 * np.pi))

    def get_parameter_guess(self, x, **kwargs):
        """
        x:
            type <array>
        """
        mu, sigma = np.mean(x), np.std(x, **kwargs)
        return self.from_lognormal_to_normal([mu, sigma])

    def get_kwargs(self):
        result = {}
        result['get_statistics'] = self.get_statistics
        result['pdf'] = self.pdf
        result['get_parameter_guess'] = self.get_parameter_guess
        result['from_lognormal_to_normal'] = self.from_lognormal_to_normal
        result['from_normal_to_lognormal'] = self.from_normal_to_lognormal
        result['nparams'] = self.nparams
        return result

class LinearEquation():

    def __init__(self):
        self.nparams = 2

    @staticmethod
    def f(prms, x):
        """
        prms:
            type <tuple / list / array>

        x:
            type <array>
        """
        return prms[0] * x + prms[1]

    @staticmethod
    def get_parameter_guess(x, y):
        """
        x:
            type <array>

        y:
            type <array>
        """
        dy = y[-1] - y[0]
        dx = x[-1] - x[0]
        m = dy / dx
        b = np.mean([y[idx] - m * x[idx] for idx in (0, -1)])
        return (m, b)

    def get_kwargs(self):
        result = {}
        result['f'] = self.f
        result['get_parameter_guess'] = self.get_parameter_guess
        result['nparams'] = self.nparams
        return result

class DistributionModel():

    def __init__(self, distribution_model):
        """
        distribution_model:
            type <str>
        """
        self._distribution_model = None
        self.initialize_model(distribution_model)

    @property
    def available_distribution_models(self):
        result = {}
        result['linear equation'] = LinearEquation()
        result['normal distribution'] = NormalDistribution()
        result['lognormal distribution'] = LogNormalDistribution()
        return result

    def initialize_model(self, distribution_model):
        """
        distribution_model:
            type <str>
        """
        if distribution_model not in list(self.available_distribution_models.keys()):
            raise ValueError("invalid distribution_model: {}".format(distribution_model))
        self._distribution_model = distribution_model
        model = self.available_distribution_models[distribution_model]
        for key, attribute in model.get_kwargs().items():
            setattr(self, key, attribute)
        try:
            f = lambda x, prms : self.pdf(prms, x)
            setattr(self, 'integrable_pdf', f)
            setattr(self, 'is_probability_density', True)
        except:
            setattr(self, 'is_probability_density', False)

class Extremizer():

    @staticmethod
    def search_local_extremum(f, parameter_guess, args=(), method=None, jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None):
        return minimize(f, parameter_guess, args, method, jac, hess, hessp, bounds, constraints, tol, callback, options)

    @staticmethod
    def search_global_extremum(f, parameter_guess, niter=100, T=1.0, stepsize=0.5, minimizer_kwargs=None, take_step=None, accept_test=None, callback=None, interval=50, disp=False, niter_success=None, seed=None):
        return basinhopping(f, parameter_guess, niter, T, stepsize, minimizer_kwargs, take_step, accept_test, callback, interval, disp, niter_success, seed)

    def search_extremum(self, f, parameter_guess=None, scale='local', minimizer_kwargs=None, **kwargs):
        """
        f:
            type <function>

        parameter_guess:
            type <tuple / list / array>

        scale:
            type <str>

        minimizer_kwargs:
            type <dict> or None
        """
        result = {}
        if scale == 'local':
            if minimizer_kwargs is None:
                minimizer_kwargs = {}
            minimizer_kwargs.update(**kwargs)
            res = self.search_local_extremum(f, parameter_guess, **minimizer_kwargs)
        elif scale == 'global':
            res = self.search_global_extremum(f, parameter_guess, minimizer_kwargs=minimizer_kwargs, **kwargs).lowest_optimization_result
        else:
            raise ValueError("invalid scale: {}".format(scale))
        if res.success != True:
            raise ValueError("optimization was not successful:\n{}\n".format(res))
        result['parameters'] = res.x
        result['fun'] = res.fun
        result['success'] = res.success
        return result









class GeneralizedLeastSquares(Extremizer):

    def __init__(self, data, model):
        """
        data:
            type <array>

        model:
            type <custom class>
        """
        super().__init__()
        self.x = data[0]
        self.y = data[1]
        self.model = model
        self._err_func = None

    @property
    def err_func(self):
        return self._err_func

    def get_unweighted_residuals(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        return np.square(self.y - self.model.f(prms, self.x))

    def get_unweighted_residual_sum(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        residuals = self.get_unweighted_residuals(prms)
        return np.sum(residuals)

    def get_residual_weights(residuals, tol=1e-10):
        """
        residuals:
            type <array>

        tol:
            type <float>
        """
        err = np.sum(residuals)
        if total < tol:
            return np.ones(err.size)
        else:
            return residuals / total

    def get_weighted_residual_sum(self, prms, wts):
        """
        prms:
            type <tuple / list / array>

        wts:
            type <array>
        """
        residuals = self.get_unweighted_residuals(prms)
        return np.sum(wts * residuals)

    def get_residually_weighted_residual_sum(self, prms, tol=1e-10):
        """
        prms:
            type <tuple / list / array>

        tol:
            type <float>
        """
        residuals = self.get_unweighted_residuals(prms)
        wts = self.get_residual_weights(residuals, tol)
        return np.sum(wts * residuals)

    def set_err_func(self, weights, tol):
        """

        """
        if weights is None:
            f = lambda prms : self.get_unweighted_residual_sum(prms)
        elif isinstance(weights, np.ndarray):
            f = lambda prms : self.get_weighted_residual_sum(prms, weights)
        elif isinstance(weights, str):
            if weights == 'residual':
                f = lambda prms : self.get_residually_weighted_residual_sum(prms, tol)
            else:
                raise ValueError("invalid weights: {}".format(weights))
        else:
            raise ValueError("invalid type(weights): {}".format(type(weights)))
        self._err_func = f

    def fit(self, parameter_guess=None, scale='local', minimizer_kwargs=None, weights=None, tol=1e-10, **kwargs):
        """
        prms:
            type <tuple / list / array>

        scale:
            type <str>

        minimizer_kwargs:
            type <dict> or None

        weights:
            type <str / array> or None

        tol:
            type <float>
        """
        self.set_err_func(weights, tol)
        if parameter_guess is None:
            parameter_guess = self.model.get_parameter_guess(self.x, self.y)
            # indices = np.argsort(self.x)
            # xi, yi = np.copy(self.x[indices]), np.copy(self.y[indices])
            # parameter_guess = self.model.get_parameter_guess(xi, yi)
        return self.search_extremum(self.err_func, parameter_guess, scale, minimizer_kwargs, **kwargs)

class MaximumLikelihoodEstimation(Extremizer):

    def __init__(self, data, model):
        """
        data:
            type <array>

        model:
            type <custom class>
        """
        super().__init__()
        self.data = data
        self.model = model
        self._err_func = None

    @property
    def err_func(self):
        return self._err_func

    def log_likelihood(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        return np.log(self.model.pdf(prms, self.data))

    def set_err_func(self, weights):
        """
        prms:
            type <tuple / list / array>

        weights:
            ...
        """
        if weights is None:
            f = lambda prms : -1 * np.sum(self.log_likelihood(prms))
        else:
            if not isinstance(weights, np.ndarray):
                raise ValueError("invalid type(weights): {}".format(type(weights)))
            f = lambda prms : -1 * np.sum(self.log_likelihood(prms) * weights)
        self._err_func = f

    def fit(self, parameter_guess=None, scale='local', minimizer_kwargs=None, weights=None, **kwargs):
        """
        parameter_guess:
            type <tuple / list / array> or None

        scale:
            type <str>

        minimizer_kwargs:
            type <dict> or None

        weights:
            ...
        """
        self.set_err_func(weights)
        if parameter_guess is None:
            parameter_guess = self.model.get_parameter_guess(np.sort(self.data))
        return self.search_extremum(self.err_func, parameter_guess, scale, minimizer_kwargs, **kwargs)

class BinStatistics(Extremizer):

    def __init__(self, histogram):
        """
        histogram:
            type <custom class>
        """
        super().__init__()
        self.histogram = histogram
        self._err_func = None

    @property
    def err_func(self):
        return self._err_func

    def reduce_statistic(self, err):
        """
        err:
            type <int / float>
        """
        return err / (self.histogram.bin_widths.size - self.histogram.model.nparams - 1)

    def get_expected_values(self, args=(), **kwargs):
        """
        args:
            type <tuple / list / array>
        """
        result = [quad(self.histogram.model.integrable_pdf, lbound, ubound, args, **kwargs)[0] * self.histogram.data.size for lbound, ubound in zip(self.histogram.edges[:-1], self.histogram.edges[1:])]
        return np.array(result)

    def get_chi_square(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        expected_counts = self.get_expected_values(prms)
        csq = (self.histogram.observed_counts - expected_counts)**2 / expected_counts
        return np.sum(csq)

    def get_reduced_chi_square(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        csq = self.get_chi_square(prms)
        return self.reduce_statistic(csq)

    def get_gstatistic(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        expected_counts = self.get_expected_values(prms)
        gte = 2 * self.histogram.observed_counts * np.log(self.histogram.observed_counts / expected_counts)
        return np.sum(gte)

    def get_reduced_gstatistic(self, prms):
        """
        prms:
            type <tuple / list / array>
        """
        err = self.get_gstatistic(prms)
        return self.reduce_statistic(err)

    def set_err_func(self, error_metric, reduce_statistic):
        """

        """
        if error_metric == 'chi square':
            if reduce_statistic == True:
                f = lambda prms : self.get_reduced_chi_square(prms)
            else:
                f = lambda prms : self.get_chi_square(prms)
        elif error_metric == 'g-test':
            if reduce_statistic == True:
                f = lambda prms : self.get_reduced_gstatistic(prms)
            else:
                f = lambda prms : self.get_gstatistic(prms)
        else:
            raise ValueError("invalid error_metric: {}".format(error_metric))
        self._err_func = f

    def fit(self, error_metric, parameter_guess=None, scale='local', reduce_statistic=False, minimizer_kwargs=None, **kwargs):
        """
        error_metric:
            type <str>

        parameter_guess:
            type <tuple / list / array> or None

        scale:
            type <str>

        reduce_statistic:
            type <bool>

        minimizer_kwargs:
            type <dict> or None
        """
        self.set_err_func(error_metric, reduce_statistic)
        if parameter_guess is None:
            parameter_guess = self.histogram.model.get_parameter_guess(np.sort(self.histogram.data))
        return self.search_extremum(self.err_func, parameter_guess, scale, minimizer_kwargs, **kwargs)

class Histogram():

    """
    This class contains methods to construct a histogram with
    fitting methods (such as chi square and g-statistic) and
    adaptable binning.
    """

    def __init__(self, data, distribution_model=None):
        """
        data:
            type <array>

        distribution_model:
            type <str> or None
        """
        self.data = data
        self.x = np.sort(data)
        self._bias = None
        self._threshold = None
        self._threshold_condition = None
        self._edges = None
        self._observed_counts = None
        self._optimizations = {}
        if distribution_model is None:
            self.model = None
        else:
            self.model = DistributionModel(distribution_model)

    def __str__(self):
        bias = '\n .. BIAS:\n{}\n'.format(self.bias)
        threshold = '\n .. THRESHOLD:\n{}\n'.format(self.threshold)
        normalization_constant = '\n .. NORMALIZATION CONSTANT:\n{}\n'.format(self.normalization_constant, self.normalization_constant)
        bin_widths = '\n .. {} BIN WIDTHS:\n{}\n'.format(self.bin_widths.size, self.bin_widths)
        edges = '\n .. {} EDGES:\n{}\n'.format(self.edges.size, self.edges)
        midpoints = '\n .. {} MIDPOINTS:\n{}\n'.format(self.midpoints.size, self.midpoints)
        observed_counts = '\n .. {} OBSERVED COUNTS:\n{}\n'.format(self.observed_counts.size, self.observed_counts)
        return '{}{}{}{}{}{}{}'.format(bias, threshold, normalization_constant, bin_widths, edges, midpoints, observed_counts)

    @property
    def optimizations(self):
        return self._optimizations

    @property
    def bias(self):
        return self._bias

    @property
    def threshold(self):
        return self._threshold

    @property
    def threshold_condition(self):
        return self._threshold_condition

    @property
    def edges(self):
        return self._edges

    @property
    def midpoints(self):
        return (self.edges[1:] + self.edges[:-1])/2

    @property
    def bin_widths(self):
        return np.diff(self.edges)

    @property
    def observed_counts(self):
        return self._observed_counts

    @property
    def normalization_constant(self):
        return np.sum(self.bin_widths * self.observed_counts)

    @property
    def normalized_counts(self):
        return self.observed_counts / self.normalization_constant

    def initialize_edges(self, value, criterion='number of edges'):
        """
        value:
            type <int / float>

        criterion:
            type <str>
        """
        if criterion == 'number of edges':
            self._edges = np.linspace(min(self.data), max(self.data), value)
        elif criterion == 'number of bins':
            self.initialize_edges(value +1, 'number of edges')
        elif criterion == 'bin width':
            self._edges = np.arange(min(self.data), max(self.data) + 0.9 * value, value)
        elif criterion == 'edges':
            self._edges = value
        else:
            raise ValueError("invalid criterion: {}".format(criterion))

    def initialize_observed_counts(self, bias='left'):
        """
        bias:
            type <str>
        """
        if bias == 'left':
            counts, edges = np.histogram(self.data, self.edges)
            self._edges = edges
        elif bias == 'right':
            counts = np.zeros(len(self.edges) - 1, dtype=int)
            for idx, val in zip(*np.unique(np.searchsorted(self.edges, self.data, side='left'), return_counts=True)):
                counts[idx - 1] = val
        else:
            raise ValueError("invalid bias: {}".format(self.bias))
        self._bias = bias
        self._observed_counts = counts

    def apply_bin_threshold(self, threshold=None, threshold_condition='greater than or equal'):
        """
        threshold:
            type <int> or None

        threshold_condition:
            type <str>
        """
        if threshold is not None:
            if not isinstance(threshold, int):
                raise ValueError("invalid threshold: {}".format(threshold))
            if threshold < 0:
                raise ValueError("threshold must be non-negative")
            loc = np.argmax(self.observed_counts)
            idx = 0
            left_edges = [self.edges[0]]
            left_counts = []
            while idx <= loc:
                cs, es = self.observed_counts[idx], self.edges[idx +1]
                while cs < threshold:
                    idx += 1
                    cs += self.observed_counts[idx]
                    es = self.edges[idx +1]
                left_edges.append(es)
                left_counts.append(cs)
                idx += 1
            idx = len(self.observed_counts) -1
            right_edges = []
            right_counts = []
            while idx > loc:
                cs, es = self.observed_counts[idx], self.edges[idx +1]
                while cs < threshold:
                    idx -= 1
                    cs += self.observed_counts[idx]
                    es = self.edges[idx +1]
                right_edges.append(es)
                right_counts.append(cs)
                idx -= 1
            right_edges[0] = self.edges[-1]
            self._edges = np.array(left_edges + right_edges[::-1])
            self._observed_counts = np.array(left_counts + right_counts[::-1])
            self._threshold = threshold
            self._threshold_condition = threshold_condition

    def erase_empty_bins_at_boundaries(self):
        condition = (self._observed_counts > 0)
        loc = np.where(condition)[0]
        xmin, xmax = np.min(loc), np.max(loc)
        edges = np.copy(self._edges[xmin : xmax +1])
        counts = np.copy(self._observed_counts[xmin : xmax])
        self._edges = edges
        self._observed_counts = counts

    def initialize_histogram(self, value, criterion='number of edges', threshold=None, condition='greater than or equal', bias='left', remove_empty_boundaries=False):
        """
        value:
            type <int / float>

        criterion:
            type <str>

        threshold:
            type <int> or None

        condition:
            type <str>

        bias:
            type <str>

        remove_empty_boundaries:
            type <bool>
        """
        self.initialize_edges(value, criterion)
        self.initialize_observed_counts(bias)
        if remove_empty_boundaries == True:
            self.erase_empty_bins_at_boundaries()
        self.apply_bin_threshold(threshold, condition)

    def initialize_optimization(self, error_metric, parameter_guess=None, scale='local', reduce_statistic=False, minimizer_kwargs=None, **kwargs):
        """
        error_metric:
            type <str>

        parameter_guess:
            type <tuple / list / array> or None

        scale:
            type <str>

        reduce_statistic:
            type <bool>
        """
        if self.model is None:
            raise ValueError("distribution_model was not initialized")
        if self.model.is_probability_density != True:
            raise ValueError("not yet implemented")
        BS = BinStatistics(self)
        result = BS.fit(error_metric, parameter_guess, scale, reduce_statistic, minimizer_kwargs, **kwargs)
        if error_metric in list(self._optimizations.keys()):
            if result['fun'] < self._optimizations[error_metric]['fun']:
                self._optimizations[error_metric] = result
        else:
            self._optimizations[error_metric] = result

class ErrorDimensionalization():

    def __init__(self, Error, prms):
        self.Error = Error
        self.prms = prms
        self._x = None
        self._y = None
        self._z = None
        self._X = None
        self._Y = None
        self._Z = None

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def X(self):
        return self._X

    @property
    def Y(self):
        return self._Y

    @property
    def Z(self):
        return self._Z

    def get_parameter_space(self, prm, frac=1, n=1, space_by='number'):
        """

        """
        delta = prm * frac
        pmin, pmax = prm - delta, prm + delta
        if pmin < 0:
            pmin = 0
        if space_by == 'number':
            pspace = np.linspace(pmin, pmax, n)
        elif space_by == 'width':
            pspace = np.arange(pmin, pmax + n, n)
        else:
            raise ValueError("invalid space_by: {}".format(space_by))
        return pspace

    def load_parameter_space(self, x=None, y=None, xfrac=0.5, xn=7, xspace_by='number', yfrac=0.5, yn=7, yspace_by='number'):
        """

        """
        if x is None:
            x = self.get_parameter_space(self.prms[0], xfrac, xn, xspace_by)
        if not isinstance(x, np.ndarray):
            raise ValueError("invalid type(x): {}".format(type(x)))
        if y is None:
            y = self.get_parameter_space(self.prms[1], yfrac, yn, yspace_by)
        if not isinstance(y, np.ndarray):
            raise ValueError("invalid type(y): {}".format(type(y)))
        X, Y = np.meshgrid(x, y)
        self._x = x
        self._y = y
        self._X = X
        self._Y = Y
        z = np.array([self.Error.err_func([xi, yi]) for xi, yi in zip(self.X.reshape(-1), self.Y.reshape(-1))])
        self._z = z
        self._Z = z.reshape(X.shape)







































# prms = (50, 10)
# data = np.random.normal(loc=50, scale=10, size=1000)
# x = np.sort(data)
# MLE = MaximumLikelihoodEstimation(x, DistributionModel('normal distribution'))
# result = MLE.fit((30, 20), 'local', None, **{'method' : 'Nelder-Mead'})
#
# _prms = (4.5, 0.7)
# prms = LogNormalDistribution().from_lognormal_to_normal(_prms)
# print(prms)
# data = np.random.lognormal(*prms, size=1000)
# x = np.sort(data)
# MLE = MaximumLikelihoodEstimation(x, DistributionModel('lognormal distribution'))
# result = MLE.fit((2, 0.2), 'local', None, **{'method' : 'Nelder-Mead'})
# print(result)


# _prms = (4.5, 0.7)
# prms = LogNormalDistribution().from_lognormal_to_normal(_prms)
# print(prms)
# data = np.random.lognormal(*prms, size=1000)
# x = np.sort(data)
# edges = np.arange(0, 10.1, 0.2)
# H = Histogram(x, 'lognormal distribution')
# H.initialize_histogram(edges, 'edges', remove_empty_boundaries=True)
# H.initialize_optimization('chi square', parameter_guess=(1.5, 0.15), reduce_statistic=True, scale='local')
# H.initialize_optimization('g-test', parameter_guess=(1.5, 0.15), reduce_statistic=True, scale='local')
# print(H.optimizations)




# prms = (50, 10)
# data = np.random.normal(loc=50, scale=10, size=1000)
# x = np.sort(data)
# edges = np.arange(0, 101, 10).astype(int)
# H = Histogram(x, 'normal distribution')
# H.initialize_histogram(edges, 'edges', remove_empty_boundaries=True)
# H.initialize_optimization('chi square', parameter_guess=(30, 15), reduce_statistic=True, scale='local')
# H.initialize_optimization('g-test', parameter_guess=(30, 15), reduce_statistic=True, scale='local')
# print(H.optimizations)

# _prms = (4, 0.5)
# data = np.random.lognormal(mean=_prms[0], sigma=_prms[1], size=1000)
# prms = LogNormalDistribution().from_lognormal_to_normal([np.mean(data), np.std(data)])
# # prms = LogNormalDistribution().from_normal_to_lognormal(_prms)
# print(prms)
# x = np.sort(data)
# # print(x)
# edges = np.arange(0, 10, 1).astype(int)
# H = Histogram(x, 'normal distribution')
# H.initialize_histogram(edges, 'edges', remove_empty_boundaries=True)
# H.initialize_optimization('chi square', parameter_guess=(3, 0.4), reduce_statistic=True, scale='local')
# H.initialize_optimization('g-test', parameter_guess=(3, 0.4), reduce_statistic=True, scale='local')
# print(H.optimizations)





# np.random.seed(327)
# prms = (50, 10)
# data = np.random.normal(loc=prms[0], scale=prms[1], size=1000).astype(int)
# x = np.sort(data)
#
# edges = np.arange(0, 101, 10).astype(int)
# for scale in ('local', 'global'):
#     for error_metric in ('chi square', 'g-test'):
#         H = Histogram(x, 'normal distribution')
#         H.initialize_histogram(edges, 'edges', remove_empty_boundaries=True)
#         H.initialize_optimization(error_metric, parameter_guess=(30, 15), reduce_statistic=True, scale=scale)
#         print("\n{}\n".format(H.optimizations))
#

# model = DistributionModel('normal distribution')
# MLE = MaximumLikelihoodEstimation(data, model)
# mle_fit = MLE.fit(parameter_guess=(44, 15))
# print(mle_fit)

# x = np.arange(20)
# noise = np.random.uniform(low=-0.5, high=0.5, size=x.size)
# m, b = np.exp(1), -np.pi
# parameter_guess = (3, -5)
# y = m * x + b + noise
# data = [x, y]
# minimizer_kwargs = {'method' : 'Nelder-Mead'}
# model = DistributionModel('linear equation')
# GLS = GeneralizedLeastSquares(data, model)
# result = GLS.fit(parameter_guess, minimizer_kwargs=minimizer_kwargs)
# print(result)








##
