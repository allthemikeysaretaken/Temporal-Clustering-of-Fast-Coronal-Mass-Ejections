import numpy as np
from scipy import special
from scipy.optimize import minimize, basinhopping
from scipy.stats import sem, chisquare, mode

from visual_configurations import *

class DistributionMap(VisualEditor):

    def __init__(self):
        super().__init__()
        self._model = None
        self._f = None
        self._parametrize = None

    @property
    def model(self):
        return self._model

    @property
    def f(self):
        return self._f

    @property
    def parametrize(self):
        return self._parametrize

    @staticmethod
    def from_lognormal_to_normal(prms):
        """ """
        mu, sigma = prms
        variance = prms[1]**2
        tmp = variance / prms[0]**2 + 1
        mu_mod = np.log(prms[0]/np.sqrt(tmp))
        sigma_mod = np.sqrt(np.log(tmp))
        return (mu_mod, sigma_mod)

    @staticmethod
    def from_normal_to_lognormal(prms):
        """ """
        mu_mod = np.exp(prms[0] + prms[1]**2 /2)
        sigma_mod = np.sqrt(np.exp(2*prms[0] + prms[1]**2) * (np.exp(prms[1]**2) - 1))
        return (mu_mod, sigma_mod)

    @property
    def linear(self):
        eqn = lambda prms, x : x * prms[0] + prms[1]
        get_initial_parameter_guess = lambda x, y : np.array([(y[-1] - y[0]) / (x[-1] - x[0]), np.mean([y[idx] - x[idx] * (y[-1] - y[0]) / (x[-1] - x[0]) for idx in (0, -1)])])
        return {'eqn' : eqn, 'parametrize' : get_initial_parameter_guess}

    @property
    def normal(self):
        eqn = lambda prms, x : np.exp(- (x - prms[0])**2 / (2* prms[1]**2)) / (prms[1] * np.sqrt(2*np.pi))
        get_initial_parameter_guess = lambda x : (np.mean(x), np.std(x))
        get_stats = lambda x : {'mean' : np.mean(x), 'median' : np.median(x), 'mode' : mode(x, axis=None)[0][0], 'standard deviation' : np.std(x)}
        return {'pdf' : eqn, 'parametrize' : get_initial_parameter_guess, 'statistize' : get_stats}

        (mu, sigma) = SD.from_normal_to_lognormal(parameters)
        median = np.sqrt((np.exp(parameters[1]**2 -1)) * (np.exp(2*parameters[0] + parameters[1]**2)))
        mode = np.exp(parameters[0] - parameters[1]**2)

    def get_lognormal_stats(self, prms):
        """ """
        (mu, sigma) = self.from_normal_to_lognormal(prms)
        median = np.sqrt((np.exp(prms[1]**2 -1)) * (np.exp(2*prms[0] + prms[1]**2)))
        md = np.exp(prms[0] - prms[1]**2)
        return {'mean' : mu, 'median' : median, 'mode' : md, 'standard deviation' : sigma}

    @property
    def lognormal(self):
        eqn = lambda prms, x : np.exp(- (np.log(x) - prms[0])**2 / (2* prms[1]**2)) / (x * prms[1] * np.sqrt(2*np.pi))
        get_initial_parameter_guess = lambda x : self.from_lognormal_to_normal((np.mean(x), np.std(x)))
        get_stats = lambda prms : self.get_lognormal_stats(prms)
        return {'pdf' : eqn, 'parametrize' : get_initial_parameter_guess, 'statistize' : get_stats}

    @property
    def models(self):
        res = {}
        res['linear'] = self.linear
        res['normal'] = self.normal
        res['lognormal'] = self.lognormal
        return res

    def configure_model(self, model):
        """ """
        available_models = list(self.models.keys())
        if model in available_models:
            dmap = self.models[model]
            try:
                self._f = dmap['pdf']
            except:
                self._f = dmap['eqn']
            self._parametrize = dmap['parametrize']
        else:
            raise ValueError("invalid model: {}; available_models: {}".format(model, available_models))

class ExtraDimensionalMapper(DistributionMap):

    def __init__(self):
        super().__init__()
        self.dim2_space = {}
        self.dim3_space = {}

    @staticmethod
    def get_parameter_space(prm, frac, n, space_by):
        """
        prm         :   type <float / int>
        frac        :   type <float>
        n           :   type <int>
        space_by    :   type <str>
        """
        delta = prm * frac
        if space_by == 'number':
            pmin = prm - delta
            if pmin < 0:
                pmin = 0
            pspace = np.linspace(pmin, prm + delta, n)
        elif space_by == 'width':
            pmin = prm - n * delta
            pmax = prm + n * delta
            if pmin < 0:
                # pmin, pmax = 0, pmax + pmin
                pmin = 0
            pspace = np.arange(pmin, pmax, delta)
        else:
            raise ValueError("space_by = 'number' or 'width'")
        return pspace

    def initialize_dim2_space(self, prms, estimator, x=None, y=None, xfrac=0.5, xn=7, xspace_by='number', yfrac=0.5, yn=7, yspace_by='number'):
        """
        Inputs 'x' and/or 'y' will supersede other inputs if not None.
        x           :   type <array> or None
        y           :   type <array> or None
        xfrac       :   type <float>
        xn          :   type <int>
        xspace_by   :   type <str>
        yfrac       :   type <float>
        yn          :   type <int>
        yspace_by   :   type <str>
        """
        self.dim2_space[estimator] = {}
        if x is None:
            self.dim2_space[estimator]['x'] = self.get_parameter_space(prms[0], xfrac, xn, xspace_by)
        elif isinstance(x, (tuple, list, np.ndarray)):
            self.dim2_space[estimator]['x'] = np.array(x)
        else:
            raise ValueError("x = type <array> or None")
        if y is None:
            self.dim2_space[estimator]['y'] = self.get_parameter_space(prms[1], yfrac, yn, yspace_by)
        elif isinstance(y, (tuple, list, np.ndarray)):
            self.dim2_space[estimator]['y'] = np.array(y)
        else:
            raise ValueError("y = type <array> or None")

    def initialize_dim3_space(self, prms, estimator, f, args=()):
        """ """
        if len(self.dim2_space.keys()) == 0:
            raise ValueError("initialize dim2_space before initializing dim3_space")
        self.dim3_space[estimator] = {}
        X, Y = np.meshgrid(self.dim2_space[estimator]['x'], self.dim2_space[estimator]['y'])
        xx = X.copy().reshape(-1)
        yy = Y.copy().reshape(-1)
        z = np.array([f(np.array([xi, yi]), *args) for xi, yi in zip(xx, yy)]).reshape(-1)
        if estimator in ('maximum likelihood',):
            z = -z
        Z = z.reshape(X.shape)
        self.dim3_space[estimator]['x'] = xx
        self.dim3_space[estimator]['y'] = yy
        self.dim3_space[estimator]['z'] = z
        self.dim3_space[estimator]['X'] = X
        self.dim3_space[estimator]['Y'] = Y
        self.dim3_space[estimator]['Z'] = Z

    # @property
    # def estimator_symbols(self):
    #     res = {}
    #     res['maximum likelihood'] = r'$- ln(L)$'
    #     res['gtest'] = r'$G_{min}$'
    #     res['reduced gtest'] = r'$G_{min, red}$'
    #     res['chi square'] = r'$\chi_{min}^2$'
    #     res['reduced chi square'] = r'$\chi_{min, red}^2$'
    #     return res

class Extremizer(ExtraDimensionalMapper):

    def __init__(self):
        super().__init__()

    @property
    def available_local_methods(self):
        res = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP']
        res += ['trust-constr', 'dogleg', 'trust-ncd', 'trust-exact', 'trust-krylov']
        return tuple(res)

    @property
    def available_global_methods(self):
        res = ('basin-hop',) # 'brute-force', 'simulated annealing')
        return res

    @staticmethod
    def autocorrect_bounds(bounds, local_method):
        """ """
        if bounds is not None:
            if local_method not in ('L-BFGS-B', 'TNC', 'SLSQP'):
                raise ValueError("bounds = None for this local_method")
        return bounds

    @staticmethod
    def autocorrect_constraints(constraints, local_method):
        """ """
        if constraints is None:
            constraints = ()
        else:
            if local_method not in ('COBYLA', 'SLSQP'):
                raise ValueError("constraints = None for this local_method")
        return constraints

    def search_local_extremum(self, f_err, initial_parameter_guess, args=(), bounds=None, constraints=None, local_method='Nelder-Mead', **kwargs):
        """ """
        if local_method not in self.available_local_methods:
            raise ValueError("invalid method: {}; available local methods: {}".format(local_method, self.available_local_methods))
        bounds = self.autocorrect_bounds(bounds, local_method)
        constraints = self.autocorrect_constraints(constraints, local_method)
        kwargs = dict(**kwargs)
        res = minimize(f_err, x0=initial_parameter_guess, args=args, bounds=bounds, constraints=constraints, method=local_method, **kwargs)
        return res

    def search_global_extremum(self, f_err, initial_parameter_guess, args=(), bounds=None, constraints=None, local_method='Nelder-Mead', global_method='basin-hop', **kwargs):
        """ """
        if global_method not in self.available_global_methods:
            raise ValueError("invalid method: {}; available global methods: {}".format(global_method, self.available_global_methods))
        if global_method == 'basin-hop':
            if local_method in self.available_local_methods:
                bounds = self.autocorrect_bounds(bounds, local_method)
                constraints = self.autocorrect_constraints(constraints, local_method)
                minimizer_kwargs = {'args' : args, 'bounds' : bounds, 'constraints' : constraints}
                res = basinhopping(f_err, x0=initial_parameter_guess, minimizer_kwargs=minimizer_kwargs, **kwargs)
            else:
                raise ValueError("invalid method: {}; available local methods: {}".format(local_method, self.available_local_methods))
        else:
            raise ValueError("not yet implemented")
        return res

    def search_for_extremum(self, f_err, initial_parameter_guess, args=(), bounds=None, constraints=None, local_method='Nelder-Mead', global_method=None, **kwargs):
        """ """
        if global_method is None:
            res = self.search_local_extremum(f_err, initial_parameter_guess, args, bounds, constraints, local_method, **kwargs)
        else:
            res = self.search_global_extremum(f_err, initial_parameter_guess, args, bounds, constraints, local_method, global_method, **kwargs)
        return res

class GeneralizedLeastSquares(Extremizer):

    def __init__(self, x, y):
        """ """
        super().__init__()
        self.configure_model('linear')
        self.x = x
        self.y = y
        self._f_err = None

    @property
    def f_err(self):
        return self._f_err

    def get_residuals(self, prms):
        """ """
        return (self.y - self.f(prms, self.x))**2

    @staticmethod
    def get_inverse_residual_weights(residuals, tolerance=10**-10):
        """ """
        inv = np.abs(residuals)
        tot = np.sum(inv)
        if tot < tolerance:
            res = np.ones(inv.size)
        else:
            residuals = residuals / tot
        return residuals

    def ordinary(self, prms):
        """ """
        residuals = self.get_residuals(prms)
        return np.sum(residuals)

    def custom_weighted(self, prms, weights):
        """ """
        residuals = self.get_residuals(prms)
        return np.sum(weights * residuals)

    def residually_weighted(self, prms, tolerance=10**-10):
        """ """
        residuals = self.get_residuals(prms)
        weights = self.get_inverse_residual_weights(residuals, tolerance)
        return np.sum(weights * residuals)

    def fit(self, initial_parameter_guess=None, weights=None, bounds=None, constraints=None, local_method='Nelder-Mead', global_method=None, **kwargs):
        """ """
        if initial_parameter_guess is None:
            initial_parameter_guess = self.parametrize(self.x, self.y)
        if weights is None:
            f_err = self.ordinary
            args = ()
        elif isinstance(weights, (tuple, list, np.ndarray)):
            f_err = self.custom_weighted
            args = (weights,)
        elif weights == 'residual':
            f_err = self.residually_weighted
            args = ()
        self._f_err = f_err
        res = self.search_for_extremum(f_err, initial_parameter_guess, args, bounds, constraints, local_method, global_method, **kwargs)
        return res

class MaximumLikelihoodEstimation(Extremizer):

    def __init__(self, x, model):
        """ """
        super().__init__()
        self.configure_model(model)
        self.x = x
        self._f_err = None

    @property
    def f_err(self):
        return self._f_err

    def get_log_likelihood(self, prms):
        """ """
        return np.log(self.f(prms, self.x))

    # def f_err(self, prms):
    #     """ """
    #     logl = self.get_log_likelihood(prms)
    #     return -1 * np.sum(logl)

    def fit(self, initial_parameter_guess=None, bounds=None, constraints=None, local_method='Nelder-Mead', global_method=None, show_initial_parameter_guess=False, **kwargs):
        """ """
        if initial_parameter_guess is None:
            initial_parameter_guess = self.parametrize(self.x)
        if show_initial_parameter_guess == True:
            print("\n .. initial parameter guess:\n\t{}\n".format(initial_parameter_guess))
        f_err = lambda prms : -1 * np.sum(self.get_log_likelihood(prms))
        self._f_err = f_err
        args = ()
        res = self.search_for_extremum(f_err, initial_parameter_guess, args, bounds, constraints, local_method, global_method, **kwargs)
        return res

class BinnedStatisticEstimation(Extremizer):

    def __init__(self, x, histogram, model):
        """ """
        super().__init__()
        self.configure_model(model)
        self.x = x
        self.histogram = histogram
        self.integrable_pdf = lambda x, prms : self.f(prms, x)
        self.nprms = len(self.parametrize(self.x))
        self._f_err = None

    @property
    def f_err(self):
        return self._f_err

    def reduce_statistic(self, err):
        """ """
        return err / (self.histogram.bin_widths.size - self.nprms - 1)

    def get_chi_square(self, prms):
        """ """
        expected_counts = self.histogram.get_expectation_values(f=self.integrable_pdf, args=prms)
        csq = (self.histogram.observed_counts - expected_counts)**2 / expected_counts
        return np.sum(csq)

    # def pvalue(self, csq):
    #     """ """
    #     return 1 - stats.chi2.cdf(csq, self.nprms)

    def get_gtest(self, prms):
        """ """
        expected_counts = self.histogram.get_expectation_values(f=self.integrable_pdf, args=prms)
        gte = 2 * self.histogram.observed_counts * np.log(self.histogram.observed_counts / expected_counts)
        return np.sum(gte)

    @property
    def error_map(self):
        res = {}
        res['chi square'] = lambda prms : self.get_chi_square(prms)
        res['reduced chi square'] = lambda prms : self.reduce_statistic(self.get_chi_square(prms))
        res['gtest'] = lambda prms : self.get_gtest(prms)
        res['reduced gtest'] = lambda prms : self.reduce_statistic(self.get_gtest(prms))
        return res

    def fit(self, error_metric, initial_parameter_guess=None, bounds=None, constraints=None, local_method='Nelder-Mead', global_method=None, show_initial_parameter_guess=False, **kwargs):
        """ """
        available_error_metrics = list(self.error_map.keys())
        if error_metric in available_error_metrics:
            f_err = self.error_map[error_metric]
        else:
            raise ValueError("invalid error_metric: {}; available error_metrics: {}".format(error_metric, available_error_metrics))
        self._f_err = f_err
        if initial_parameter_guess is None:
            initial_parameter_guess = self.parametrize(self.x)
        if show_initial_parameter_guess == True:
            print("\n .. initial parameter guess:\n\t{}\n".format(initial_parameter_guess))
        args = ()
        res = self.search_for_extremum(f_err, initial_parameter_guess, args, bounds, constraints, local_method, global_method, **kwargs)
        return res
