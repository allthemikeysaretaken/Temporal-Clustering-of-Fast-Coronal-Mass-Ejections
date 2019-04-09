import numpy as np
from scipy import special
from scipy.optimize import minimize, basinhopping
from scipy.stats import chisquare, chi2

from adaptive_histogram import *

class GeneralizedLeastSquares():

    """
    Purpose:
        Perform ordinary least squares or weighted least squares. For the
        weighted routines, one can input a weights vector or can use the
        residuals to weight each point.
    """

    def __init__(self, x, y):
        """
        x           :   type <array>
        y           :   type <array>
        """
        self.x = x
        self.y = y
        self.result = {}

    @staticmethod
    def get_weights(resid, tolerance=10**-10):
        """
        resid       :   type <array>
        tolerance   :   type <float>
        """
        inv = np.abs(resid)
        tot = np.sum(inv)
        if tot < tolerance:
            res = np.ones(inv.size)
            print('')
            print("tolerance condition violated; weights equalized")
            print('')
        else:
            resid = resid / tot
        return resid

    def f_linear(self, prms):
        """
        prms        :   type <tuple / list / array>
        """
        return self.x * prms[0] + prms[1]

    def get_residuals(self, prms):
        """
        prms        :   type <tuple / list / array>
        """
        return (self.y - self.f_linear(prms))**2

    def residually_weighted(self, prms):
        """
        prms        :   type <tuple / list / array>
        """
        resid = self.get_residuals(prms)
        wts = self.get_weights(resid)
        err = wts * resid
        return np.sum(err)

    def constant_weighted(self, prms, wts):
        """
        prms        :   type <tuple / list / array>
        wts         :   type <array>
        """
        resid = self.get_residuals(prms)
        err = wts * resid
        return np.sum(err)

    def unweighted(self, prms):
        """
        prms        :   type <tuple / list / array>
        """
        resid = self.get_residuals(prms)
        return np.sum(resid)

    def dispatch_f_err(self, wts):
        """
        wts         :   type <array>
        """
        if wts is None:
            f = lambda prms : self.unweighted(prms)
        elif isinstance(wts, (tuple, list, np.ndarray)):
            f = lambda prms : self.constant_weighted(prms, np.array(wts))
        elif wts == 'residual':
            f = lambda prms : self.residually_weighted(prms)
        else:
            raise ValueError("unknown wts: {}; available wts: (None, type <tuple / list / array>, or 'residual')".format(wts))
        return f

    def get_initial_parameter_guess(self):
        init_slope = (self.y[-1] - self.y[0]) / (self.x[-1] - self.x[0])
        init_intercept = np.mean([self.y[idx] - self.x[idx] * init_slope for idx in (0, -1)])
        return np.array([init_slope, init_intercept])

    def fit(self, prms=None, wts=None, method='Nelder-Mead', scale='local', **kwargs):
        """
        prms        :   type <tuple / list / array> or None
        wts         :   type <array>
        method      :   type <str>
        scale       :   type <str>
        **kwargs    :   type <dict>; pass to scipy minimize
        """
        if prms is None:
            prms = self.get_initial_parameter_guess()
        f_err = self.dispatch_f_err(wts)
        if scale == 'local':
            res = minimize(f_err, x0=prms, method=method, **kwargs)
        elif scale == 'global':
            minimizer_kwargs = dict(**kwargs)
            minimizer_kwargs['method'] = method
            res = basinhopping(f_err, x0=prms, minimizer_kwargs=minimizer_kwargs)
        else:
            raise ValueError("scale = 'local' or 'global'")
        return res

class SpeedDistribution():

    def __init__(self, data):
        """ """
        self.data = data
        self.x = np.sort(data)
        self.n = data.size
        self.p = 2
        self.storage = {}
        self.histograms = {}
        self.result = {}


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

    def initialize_parameter_guess(self, init_prms=None):
        """ """
        if init_prms is None:
            mu, sigma = np.mean(self.data), np.std(self.data)
            prms = (mu, sigma)
            m, s = self.from_lognormal_to_normal(prms)
            self.storage['initial parameter guess'] = (m, s)
        elif len(init_prms) == 2:
            self.storage['initial parameter guess'] = init_prms
        else:
            raise ValueError("init_prms = None or type <tuple / list / array> of size 2")

    @staticmethod
    def lognormal_pdf(prms, x):
        """ """
        return np.exp(- (np.log(x) - prms[0])**2 / (2* prms[1]**2)) / (x * prms[1] * np.sqrt(2*np.pi))

    def integrable_lognormal_pdf(self, x, prms):
        """ """
        return self.lognormal_pdf(prms, x)

    def initialize_preliminary_histogram(self, n, edges_by, bias='left'):
        """ """
        H = AdaptiveHistogram(self.data, None)
        H.initialize_edges(n, edges_by)
        H.initialize_observed_counts(bias)
        H.initialize_midpoints()
        H.initialize_bin_widths()
        H.initialize_normalization()
        self.histograms['preliminary'] = H

    def initialize_thresholded_histogram(self, n, edges_by, bin_threshold, bias='left'):
        """ """
        H = AdaptiveHistogram(self.data, bin_threshold)
        H.initialize_edges(n, edges_by)
        H.initialize_observed_counts(bias)
        H.apply_bin_threshold()
        H.initialize_midpoints()
        H.initialize_bin_widths()
        H.initialize_normalization()
        self.histograms['threshold'] = H

    def reduce_error(self, error, htype='threshold'):
        """
        err         :   type <float>
        """
        return error / (self.histograms[htype].bin_widths.size - self.p - 1)

    def chi_square(self, prms, htype):
        """ """
        expected_counts = self.histograms[htype].get_expectation_values(self.integrable_lognormal_pdf, prms)
        csq = (self.histograms[htype].observed_counts - expected_counts)**2 / expected_counts
        return self.reduce_error(np.sum(csq))

    def pvalue(self, csq):
        """
        csq         :   type <float>
        """
        return 1 - chi2.cdf(csq, self.p)

    def gtest(self, prms, htype='threshold'):
        """
        prms        :   type <tuple / list / array>
        """
        expected_counts = self.histograms[htype].get_expectation_values(self.integrable_lognormal_pdf, prms)
        res = 2 * np.sum(self.histograms[htype].observed_counts * np.log(self.histograms[htype].observed_counts / expected_counts))
        return self.reduce_error(res)

    def log_likelihood(self, prms):
        """ """
        return np.sum(np.log(self.lognormal_pdf(prms, self.x)))

    @property
    def errors_to_colors(self):
        """ """
        return {'maximum likelihood estimation' : 'r', 'minimum gtest estimation' : 'b', 'minimum chi square estimation' : 'darkorange', 'maximum pvalue estimation' : 'purple'}

    def dispatch_objective_error(self, key, htype='threshold'):
        """ """
        if key == 'maximum likelihood estimation':
            f = lambda prms : -1 * self.log_likelihood(prms)
        elif key == 'minimum chi square estimation':
            f = lambda prms : self.chi_square(prms, htype)
        elif key == 'minimum gtest estimation':
            f = lambda prms : self.gtest(prms, htype)
        elif key == 'maximum pvalue estimation':
            f = lambda prms : -1 * self.pvalue(self.chi_square(prms, htype))
        else:
            raise ValueError("unknown key: {}; available keys: ('maximum  likelihood estimation', 'minimnum chi square estimation', 'minimum gtest estimation')")
        return f

    def search_local_extrema(self, key, method='Nelder-Mead', htype='threshold', **kwargs):
        """ """
        f = self.dispatch_objective_error(key, htype)
        res = minimize(f, x0=self.storage['initial parameter guess'], method=method, **kwargs)
        if res.success is True:
            return {'parameters' : res.x, 'objective' : res.fun}
        else:
            raise ValueError("optimization was not successful")

    def search_global_extrema(self, key, method='Nelder-Mead', minimizer_kwargs=None, htype='threshold', **kwargs):
        """ """
        f = self.dispatch_objective_error(key, htype)
        if minimizer_kwargs is None:
            res = basinhopping(f, x0=self.storage['initial parameter guess'], **kwargs)
        elif isinstance(minimizer_kwargs, dict):
            res = basinhopping(f, x0=self.storage['initial parameter guess'], minimizer_kwargs=minimizer_kwargs, **kwargs)
        else:
            raise ValueError("see SCIPY docs")
        return {'parameters' : res.x, 'objective' : res.fun}

    def fit(self, errs, method='Nelder-Mead', scale='local', minimizer_kwargs=None, htype='threshold', **kwargs):
        """ """
        if 'initial parameter guess' not in list(self.storage.keys()):
            self.initialize_parameter_guess()
        for err in errs:
            if scale == 'local':
                res = self.search_local_extrema(err, method, htype, **kwargs)
            elif scale == 'global':
                res = self.search_local_extrema(err, method, minimizer_kwargs, htype, **kwargs)
            else:
                raise ValueError("scale = 'local' or 'global'")
            res['y'] = self.lognormal_pdf(res['parameters'], self.x)
            self.result[err] = res

    def get_y(self, err, normalized=False, htype='threshold'):
        """ """
        y = self.result[err]['y']
        if normalized is True:
            return y
        elif normalized is False:
            H = self.histograms[htype]
            return y * H.normalization_constant
        else:
            raise ValueError("normalized = True or False")


class ExtraDimensionalMapper():

    def __init__(self, SD):
        """ """
        self.SD = SD
        self.result = {}

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
                pmin, pmax = 0, pmax + pmin
            pspace = np.arange(pmin, pmax, delta)
        else:
            raise ValueError("space_by = 'number' or 'width'")
        return pspace

    def initialize_dim2_grid(self, err, x=None, y=None, xfrac=0.5, xn=7, xspace_by='number', yfrac=0.5, yn=7, yspace_by='number'):
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
        if x is None:
            self.result['x'] = self.get_parameter_space(self.SD.result[err]['parameters'][0], xfrac, xn, xspace_by)
        elif isinstance(x, (tuple, list, np.ndarray)):
            self.result['x'] = np.array(x)
        else:
            raise ValueError("x = type <array> or None")
        if y is None:
            self.result['y'] = self.get_parameter_space(self.SD.result[err]['parameters'][1], yfrac, yn, yspace_by)
        elif isinstance(y, (tuple, list, np.ndarray)):
            self.result['y'] = np.array(y)
        else:
            raise ValueError("y = type <array> or None")

    @property
    def x(self):
        return self.result['x']

    @property
    def y(self):
        return self.result['y']

    @staticmethod
    def check_invert_z(err):
        """ """
        yes = ('maximum likelihood estimation', 'maximum pvalue estimation')
        no = ('minimum chi square estimation', 'minimum gtest estimation')
        if err in yes:
            invert_z = True
        elif err in no:
            invert_z = False
        else:
            raise ValueError("unknown err: {}; available err: {} or {}".format(err, yes, no))
        return invert_z

    def initialize_dim3_grid(self, err, htype='threshold'):
        """ """
        f = self.SD.dispatch_objective_error(err, htype)
        X, Y = np.meshgrid(self.x, self.y)
        xx = X.copy().reshape(-1)
        yy = Y.copy().reshape(-1)
        z = np.array([f(np.array([xi, yi])) for xi, yi in zip(xx, yy)])
        invert_z = self.check_invert_z(err)
        if invert_z is True:
            z = - z
        elif invert_z is not False:
            raise ValueError("invert_z = True or False")
        self.result['z'] = z
        Z = z.reshape(-1).reshape(X.shape)
        self.result['X'] = X
        self.result['Y'] = Y
        self.result['Z'] = Z

    @property
    def z(self):
        return self.result['z']

    @property
    def X(self):
        return self.result['X']

    @property
    def Y(self):
        return self.result['Y']

    @property
    def Z(self):
        return self.result['Z']

    def get_zlabel(self, err):
        """ """
        if err == 'maximum likelihood estimation':
            zlabel = r'- $ln(L)$ ' + r'$= - {:.2f}$'.format(self.SD.result[err]['objective'])
        elif err == 'minimum gtest estimation':
            zlabel = r'$G_{min}$ ' + r'$= {:.2f}$'.format(self.SD.result[err]['objective'])
        else:
            zlabel = err[:].replace(' ', '\n')
        return zlabel











##
