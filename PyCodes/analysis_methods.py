from histogram_methods import *
from optimization_methods import *
from data_processing import *
# from search_methods import *
from cluster_optimizers import *
from visual_configurations import *

class FrechetDistribution(VisualEditor):

    def __init__(self):
        """ """
        super().__init__()
        self.x = np.arange(0.1, 5.006, .005)
        self.xticks = np.arange(0, self.x[-1], 0.5)
        self.alpha = np.array([1, 1, 2, 2, 3, 3])
        self.mu = np.zeros(self.alpha.size, dtype=int)
        self.sigma = np.array([1, 2, 1, 2, 1, 2])
        self.facecolors = ('red', 'green', 'blue', 'orange', 'purple', 'black')

    def subview_probability_density(self, ax, ticksize=7, labelsize=8, textsize=8):
        """ """
        yticks = np.arange(0, 1.3, 0.25)
        y, labels = [], []
        for a, m, s in zip(self.alpha, self.mu, self.sigma):
            tmp = (self.x - m)/s
            yi = (a/s) * tmp**(-1 - a) * np.exp(-1 * (tmp ** -a))
            label = r'$\alpha = {}, \mu = {}, \sigma = {}$'.format(a, m, s)
            y.append(yi)
            labels.append(label)
        for yi, facecolor, label in zip(y, self.facecolors, labels):
            ax.plot(self.x, yi, color=facecolor, alpha=0.5, label=label)
        ax.set_xticks(self.xticks, minor=True)
        ax.set_xticks(self.xticks[::2])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlim([0, 5])
        ax.set_ylim([0, 1.3])
        formula_label = r'$PDF(x | \alpha, \mu, \sigma) = \frac{\alpha}{\sigma} (\frac{x-\mu}{\sigma})^{-1-\alpha} e^{- (\frac{x-\mu}{\sigma})^{-\alpha}}$'
        text_box = ax.text(0.95, 0.95, formula_label, fontsize=textsize, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlabel('x')
        ax.set_ylabel('Probability Density')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_cumulative_density(self, ax, ticksize=7, labelsize=8, textsize=8):
        """ """
        yticks = np.arange(0, 1.1, 0.25)
        y, labels = [], []
        for a, m, s in zip(self.alpha, self.mu, self.sigma):
            tmp = (self.x - m)/s
            yi = np.exp(-1 * (tmp ** -a))
            label = r'$\alpha = {}, \mu = {}, \sigma = {}$'.format(a, m, s)
            y.append(yi)
            labels.append(label)
        for yi, facecolor, label in zip(y, self.facecolors, labels):
            ax.plot(self.x, yi, color=facecolor, alpha=0.5, label=label)
        ax.set_xticks(self.xticks, minor=True)
        ax.set_xticks(self.xticks[::2])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlim([0, 5])
        ax.set_ylim([0, 1])
        formula_label = r'$CDF(x | \alpha, \mu, \sigma) = e^{- (\frac{x-\mu}{\sigma})^{-\alpha}}$'
        text_box = ax.text(0.95, 0.05, formula_label, fontsize=textsize, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlabel('x')
        ax.set_ylabel('Cumulative Density')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

class EventFrequencyComparison(VisualEditor):

    def __init__(self, sunspot_data, cme_data):
        """ """
        super().__init__()
        self.sunspot_data = sunspot_data
        self.Searcher = SearchEvents(cme_data)
        self.available_timescales = ('daily', 'monthly', 'yearly')
        self.cme_facecolors = iter(['steelblue', 'green', 'purple', 'r', 'k'])
        self._converter = None

    @property
    def converter(self):
        return self._converter

    def initialize_converter(self):
        TSC = TimeSeriesConversion()
        if self.sunspot_data is not None:
            TSC.update_daily_sunspots(self.sunspot_data)
        self._converter = TSC

    def subview_sunspot_data(self, ax, desired_timescale='daily', ticksize=7, labelsize=8, textsize=8):
        """ """
        if desired_timescale not in self.available_timescales:
            raise ValueError("invalid desired_timescale: {}; available desired_timescale: {}".format(desired_timescale, self.available_timescales))
        if desired_timescale in ('monthly', 'yearly'):
            sunspot_data = self._converter.count_events_subroutine(self.converter.daily_sunspots, desired_timescale)
            if desired_timescale == 'monthly':
                yticks = np.arange(0, int(10e3), int(10e2))
                yspace = 2
            else:
                yticks = np.arange(0, int(10e4) + 1, int(10e3))
                yspace = 2
        else: ## daily
            sunspot_data = deepcopy(self.converter.daily_sunspots)
            yticks = np.arange(0, 401, 25)
            yspace = 4
        label = '${:,}$ Sunspots'.format(np.sum(sunspot_data['count']))
        facecolor = 'darkorange'
        if desired_timescale == 'yearly':
            handle = ax.scatter(sunspot_data['datetime object'], sunspot_data['count'], color=facecolor, marker='x')
            base = 5000 # int(10e4)
        else:
            if desired_timescale == 'monthly':
                base = 1000
            else:
                base = 50 # 100
            handle = ax.plot(sunspot_data['datetime object'], sunspot_data['count'], color=facecolor, linestyle='-', linewidth=1)
        ax = self.transform_x_as_datetime(ax)
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::yspace])
        ax.get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.grid(color='k', alpha=0.3, linestyle=':')
        # ax.set_xlim([date2num(sunspot_data['datetime object'][0]), date2num(sunspot_data['datetime object'][-1])])
        ax.set_ylim([0, yticks[-1]])
        ax.set_xlabel('Date')
        ax.set_ylabel('{}\nSunspots'.format(desired_timescale.title()))
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        ax.tick_params(axis='x', colors=facecolor, which='both')
        ax.tick_params(axis='y', colors=facecolor, which='both')
        ax.xaxis.label.set_color(facecolor)
        ax.yaxis.label.set_color(facecolor)
        return ax, handle, label

    def subview_cme_data(self, ax, yticks, yspace, search_parameters='speed', search_conditions='greater than', search_values=1, apply_to='all', desired_timescale='daily', ticksize=7, labelsize=8, textsize=8):
        """ """
        if desired_timescale not in self.available_timescales:
            raise ValueError("invalid desired_timescale: {}; available desired_timescale: {}".format(desired_timescale, self.available_timescales))
        search_inputs = np.array([search_parameters, search_conditions])
        if np.all(search_inputs == None):
            tmp_data = deepcopy(self.Searcher.events)
            ylabel = '{} CMEs'.format(desired_timescale.title())
        else:
            tmp_data = self.Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            S = Symbolizer()
            search_title = '{} {} ${}$ {}'.format(S.parameter_to_symbol_map[search_parameters], S.condition_to_symbol_map[search_conditions], search_values, S.unit_to_symbol_map[search_parameters])
            ylabel = '{} CMEs\n{}'.format(desired_timescale.title(), search_title)
        self._converter.count_daily_cmes(tmp_data)
        if desired_timescale in ('monthly', 'yearly'):
            cme_data = self._converter.count_events_subroutine(self.converter.daily_cmes, desired_timescale)
        else: ## daily
            cme_data = deepcopy(self.converter.daily_cmes)
        label = '${:,}$ CMES'.format(np.sum(cme_data['count']))
        facecolor = next(self.cme_facecolors)
        if desired_timescale == 'yearly':
            handle = ax.scatter(cme_data['datetime object'], cme_data['count'], color=facecolor, marker='x')
        else:
            handle = ax.plot(cme_data['datetime object'], cme_data['count'], color=facecolor, linestyle='-', linewidth=1)
        ax = self.transform_x_as_datetime(ax)
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::yspace])
        ax.get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.grid(color='k', alpha=0.3, linestyle=':')
        # ax.set_xlim([date2num(cme_data['datetime object'][0]), date2num(cme_data['datetime object'][-1])])
        ax.set_ylim([0, yticks[-1]])
        ax.set_xlabel('Date')
        ax.set_ylabel(ylabel)
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        ax.tick_params(axis='x', colors=facecolor, which='both')
        ax.tick_params(axis='y', colors=facecolor, which='both')
        ax.xaxis.label.set_color(facecolor)
        ax.yaxis.label.set_color(facecolor)
        return ax, handle, label

class SpeedDistribution(VisualEditor):

    def __init__(self, speeds, identifier=None):
        """ """
        super().__init__()
        self.speeds = speeds[speeds > 20]
        self.x = np.sort(self.speeds)
        self.identifier = identifier
        self._original_histogram = None
        self._thresholded_histogram = None
        self.optimizations = {}

    @property
    def original_histogram(self):
        return self._original_histogram

    @property
    def thresholded_histogram(self):
        return self._thresholded_histogram

    def initialize_histograms(self, n=None, w=None, edges=None, bin_threshold=5, bias='left'):
        """ """
        AH = AdaptiveHistogram()
        AH.initialize_histogram(self.speeds, n, w, edges, bin_threshold=None, bias=bias) #, exclude_empty_boundaries=True)
        self._original_histogram = AH
        if bin_threshold is None:
            self._thresholded_histogram = AH
        else:
            BH = AdaptiveHistogram()
            BH.initialize_histogram(self.speeds, n, w, edges, bin_threshold=bin_threshold, bias=bias) #, exclude_empty_boundaries=True)
            self._thresholded_histogram = BH

    def initialize_binned_optimization(self, initial_parameter_guess=None, estimator='reduced gtest', show_optimized_result=False):
        """ """
        BSE = BinnedStatisticEstimation(self.speeds, self.thresholded_histogram, model='lognormal')
        if (('chi square' in estimator) or ('gtest' in estimator)):
            optim_result = BSE.fit(estimator, initial_parameter_guess=initial_parameter_guess)
        else:
            raise ValueError("invalid bin_estimator: {}".format(estimator))
        if show_optimized_result == True:
            print("\n ..BIN OPTIMIZATION RESULT:\n{}\n".format(optim_result))
        self.optimizations[estimator] = {'optimized parameters' : optim_result.x, 'cls' : BSE, 'objective' : optim_result.fun}

    def initialize_alternate_optimization(self, initial_parameter_guess=None, estimator='maximum likelihood', show_optimized_result=False):
        """ """
        if estimator == 'maximum likelihood':
            OPT = MaximumLikelihoodEstimation(self.speeds, model='lognormal')
        else:
            raise ValueError("not yet implemented")
        optim_result = OPT.fit(initial_parameter_guess=initial_parameter_guess)
        if show_optimized_result == True:
            print("\n .. ALT OPTIMIZATION RESULT:\n{}\n".format(optim_result))
        self.optimizations[estimator] = {'optimized parameters' : optim_result.x, 'cls' : OPT, 'objective' : optim_result.fun}

    @property
    def estimator_labels(self):
        res = {}
        res['maximum likelihood'] = 'Maximum Likelihood Estimation'
        res['reduced gtest'] = 'Minimum Reduced G-Test Estimation'
        res['gtest'] = 'Minimum G-Test Estimation'
        res['chi square'] = 'Minimum Chi Square Estimation'
        res['reduced chi square'] = 'Minimum Reduced Chi Square Estimation'
        return res

    @property
    def estimator_colors(self):
        res = {}
        res['maximum likelihood'] = 'darkorange'
        res['reduced gtest'] = 'darkgreen'
        res['gtest'] = 'darkgreen'
        res['chi square'] = 'r'
        res['reduced chi square'] = 'r'
        return res

    @property
    def estimator_symbols(self):
        res = {}
        res['maximum likelihood'] = r'$ln(L_{max})$'
        res['reduced gtest'] = r'$G_{min, red}$'
        res['gtest'] = r'$G_{min}$'
        res['chi square'] = r'$\chi_{min}^2$'
        res['reduced chi square'] = r'$\chi_{min, red}^2$'
        return res

    @property
    def estimator_levels(self):
        res = {}
        res['maximum likelihood'] = np.array([0, 25, 50, 100, 150, 225, 300, 450, 600, 750, 900, 1150, 1300, 1450, 1600, 2000, 3000, 5000])[::-1] * -1000
        res['reduced gtest'] = (0, 100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000, 30000, 40000, 50000, 60000, 75000, 100000, 150000, 200000, 300000)
        res['gtest'] = 21
        res['chi square'] = 21
        res['reduced chi square'] = 21
        return res

    @staticmethod
    def get_statistical_parameters(optimizer, prms):
        """ """
        (mu, sigma) = optimizer.from_normal_to_lognormal(prms)
        median = np.sqrt((np.exp(prms[1]**2 -1)) * (np.exp(2*prms[0] + prms[1]**2)))
        mode = np.exp(prms[0] - prms[1]**2)
        stat_args = (mu, sigma, median, mode)
        stat_keys = ('Mean', 'Standard Deviation', 'Median', 'Mode')
        stat_colors = ('darkred', None, 'mediumblue', 'darkorchid')
        return stat_args, stat_keys, stat_colors

    @staticmethod
    def subview_statistical_parameters(ax, ymax, stat_args, stat_keys, stat_colors, textsize):
        """ """
        stat_handles, stat_labels = [], []
        for key, arg, facecolor in zip(stat_keys, stat_args, stat_colors):
            label = '{}: ${:.0f}$ '.format(key, np.round(arg, decimals=0)) + r'$\frac{km}{s}$'
            if facecolor is None:
                handle = ax.axvline(0, ymin=0, ymax=0, color='none', linestyle=':', marker=None, label=label)
            else:
                handle = ax.axvline(arg, ymin=0, ymax=ymax, color=facecolor, linestyle=':', marker=None, label=label, alpha=0.5)
            stat_handles.append(handle)
            stat_labels.append(label)
        ax.legend(handles=stat_handles, labels=stat_labels, loc='upper right', ncol=1, fontsize=textsize)
        return ax

    def subview_highlighted_tail(self, ax, y):
        """ """
        label = r'Tail ($V_{CME}$ $≥$ $800$ $\frac{km}{s}$)'
        handle = ax.fill_between(self.x, 0, y, where=(self.x >= 800), color='red', alpha=0.875, label=label)
        return ax, handle, label

    def subview_confidence_interval(self, ax, y, stat_args, fill_color='b'):
        """ """
        (mu, sigma, median, mode) = stat_args
        label_dev1 = r'$\mu \pm \sigma$'
        label_dev2 = r'$\mu \pm 2 \sigma$'
        label_dev3 = r'$\mu \pm 3 \sigma$'
        handle_dev1 = ax.fill_between(self.x, 0, y, where=((self.x > (mu - sigma)) & (self.x < (mu + sigma))), color=fill_color, alpha=0.1, label=label_dev1)
        ax.fill_between(self.x, 0, y, where=((self.x > (mu - 2*sigma)) & (self.x < (mu - sigma))), color=fill_color, alpha=0.3)
        handle_dev2 = ax.fill_between(self.x, 0, y, where=((self.x > (mu + sigma)) & (self.x < (mu + 2*sigma))), color=fill_color, alpha=0.3, label=label_dev2)
        ax.fill_between(self.x, 0, y, where=((self.x > (mu - 3*sigma)) & (self.x < (mu - 2*sigma))), color=fill_color, alpha=0.6)
        handle_dev3 = ax.fill_between(self.x, 0, y, where=((self.x > (mu + 2*sigma)) & (self.x < (mu + 3*sigma))), color=fill_color, alpha=0.6, label=label_dev3)
        ax.fill_between(self.x, 0, y, where=(self.x > (mu + 3*sigma)), color=fill_color, alpha=1)
        handles = [handle_dev1, handle_dev2, handle_dev3]
        labels = [label_dev1, label_dev2, label_dev3]
        return ax, handles, labels

    @staticmethod
    def subview_arrow_to_tail(ax, textsize, normalized=False):
        """ """
        arrowprops = {'arrowstyle': '->', 'color' : 'k'}
        if normalized == True:
            y1, y2 = 0.0004, 0.00065
        else:
            y1, y2 = 600, 825
        ax.annotate(r'Tail ($V_{CME}$ $≥$ $800$ $\frac{km}{s}$)', xy=(800, y1), xycoords='data', xytext=(1000, y2), textcoords='data', fontsize=textsize, arrowprops=arrowprops)
        return ax

    def subview_tail_histogram(self, ax, ticksize=7, labelsize=8, textsize=8, titlesize=10):
        """ """
        counts_type = 'observed frequency'
        ax.bar(self.original_histogram.midpoints, self.original_histogram.observed_counts, width=self.original_histogram.bin_widths, align='center', facecolor='k', alpha=0.375)
        ax.step(self.original_histogram.midpoints, self.original_histogram.observed_counts, color='darkorange', linewidth=0.5, where='mid')
        ax.set_xscale('log', basex=5)
        ax.set_yscale('log', basey=5)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.axvline(800, ymin=0, ymax=3125, color='r', linestyle='--', label=r'$800$ $\frac{km}{s}$')
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlabel(r'$V_{CME}$ ($\frac{km}{s}$)')
        ax.set_ylabel(counts_type.title())
        ax.set_xlim([625, 15625])
        ax.set_ylim([0.5, 3125])
        Searcher = SearchEvents({'x' : self.x})
        devents = Searcher.search(search_parameters='x', search_conditions='greater than or equal', search_values=800)
        label_top = r'${:,}$ Extreme Events '.format(devents['x'].size) + r'($V_{CME}$ $≥$ ' + r'${}$'.format(800) + r'$\frac{km}{s}$)'
        label_btm = r'max($V_{CME}$)' + '$ = {:,}$ '.format(int(self.x[-1])) + r'$\frac{km}{s}$'
        text_box = ax.text(0.95, 0.95, '{}\n{}'.format(label_top, label_btm), fontsize=textsize, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        if self.identifier is not None:
            ax.set_title(self.identifier, fontsize=titlesize)
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_probability_density_fitted_histogram(self, ax, estimator, ticksize=7, labelsize=8, textsize=8, titlesize=10, show_confidence_interval=False, highlight_tail=False, point_to_tail=False, show_statistics=False):
        """ """
        available_estimators = list(self.optimizations.keys())
        if estimator not in available_estimators:
            raise ValueError("invalid estimator: {}; available estimator: {}".format(estimator, available_estimators))
        optims = self.optimizations[estimator]
        optimizer, prms = optims['cls'], optims['optimized parameters']
        y = optimizer.f(prms, self.x)
        stat_args, stat_keys, stat_colors = self.get_statistical_parameters(optimizer, prms)
        counts_type = 'probability density'
        hist_label = 'Normalized Histogram'
        # hist_handle, = ax.step(self.original_histogram.midpoints, self.original_histogram.normalized_counts, color='k', linewidth=0.5, where='mid', label=hist_label)
        hist_handle, = ax.step(self.thresholded_histogram.midpoints, self.thresholded_histogram.normalized_counts, color='k', linewidth=0.5, where='mid', label=hist_label)
        fit_label = 'Probability Density Fit via {}'.format(self.estimator_labels[estimator])
        fit_handle, = ax.plot(self.x, y, color=self.estimator_colors[estimator], linewidth=1, linestyle='-', label=fit_label)
        handles = [hist_handle, fit_handle]
        labels = [hist_label, fit_label]
        xticks = np.arange(0, int(self.x[-1] + 2*self.original_histogram.bin_widths[0] +1), 50)
        yticks = np.arange(0, 0.0031, 0.00025)
        ymax = yticks[-1]
        if show_confidence_interval == True:
            ax, con_handles, con_labels = self.subview_confidence_interval(ax, y, stat_args, fill_color='b')
            handles.extend(con_handles)
            labels.extend(con_labels)
        if show_statistics == True:
            ax = self.subview_statistical_parameters(ax, ymax=10, stat_args=stat_args, stat_keys=stat_keys, stat_colors=stat_colors, textsize=textsize)
        if highlight_tail == True:
            ax, fill_handle, fill_label = self.subview_highlighted_tail(ax, y)
            handles.append(fill_handle)
            labels.append(fill_label)
        if point_to_tail == True:
            ax = self.subview_arrow_to_tail(ax, textsize, normalized=True)
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::5])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::4])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlim([0, 2000])
        ax.set_ylim([0, yticks[-1]])
        ax.set_xlabel(r'$V_{CME}$ ($\frac{km}{s}$)')
        ax.set_ylabel(counts_type.title())
        if self.identifier is not None:
            ax.set_title(self.identifier, fontsize=titlesize)
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax, handles, labels

    def subview_empirically_fitted_histogram(self, ax, estimator, ticksize=7, labelsize=8, textsize=8, titlesize=10, show_confidence_interval=False, highlight_tail=False, point_to_tail=False, show_statistics=False):
        """ """
        available_estimators = list(self.optimizations.keys())
        if estimator not in available_estimators:
            raise ValueError("invalid estimator: {}; available estimator: {}".format(estimator, available_estimators))
        optims = self.optimizations[estimator]
        optimizer, prms = optims['cls'], optims['optimized parameters']
        y = optimizer.f(prms, self.x) * self.thresholded_histogram.normalization_constant
        stat_args, stat_keys, stat_colors = self.get_statistical_parameters(optimizer, prms)
        xticks = np.arange(0, int(self.x[-1] + 2*self.original_histogram.bin_widths[0] +1), 50)
        if self.original_histogram.bin_widths[0] > 70:
            yticks = np.arange(0, 5001, 100)
        else:
            yticks = np.arange(0, 3001, 50)
        ymax = yticks[-1]
        counts_type = 'observed frequency'
        hist_label = 'Histogram'
        # hist_handle = ax.bar(self.original_histogram.midpoints, self.original_histogram.observed_counts, width=self.original_histogram.bin_widths, align='center', facecolor='k', edgecolor='none', linewidth=0, alpha=0.6, label=hist_label)
        hist_handle = ax.bar(self.thresholded_histogram.midpoints, self.thresholded_histogram.observed_counts, width=self.thresholded_histogram.bin_widths, align='center', facecolor='k', edgecolor='none', linewidth=0, alpha=0.6, label=hist_label)
        handles = [hist_handle]
        labels = [hist_label]
        fit_label = 'Empirical Fit via {}'.format(self.estimator_labels[estimator])
        fit_handle, = ax.plot(self.x, y, color=self.estimator_colors[estimator], linewidth=1, linestyle='-', label=fit_label)
        handles.append(fit_handle)
        labels.append(fit_label)
        if show_confidence_interval == True:
            ax, con_handles, con_labels = self.subview_confidence_interval(ax, y, stat_args, fill_color='b')
            handles.extend(con_handles)
            labels.extend(con_labels)
        if show_statistics == True:
            ax = self.subview_statistical_parameters(ax, ymax, stat_args, stat_keys, stat_colors, textsize)
        if highlight_tail == True:
            ax, fill_handle, fill_label = self.subview_highlighted_tail(ax, y)
            handles.append(fill_handle)
            labels.append(fill_label)
        if point_to_tail == True:
            ax = self.subview_arrow_to_tail(ax, textsize, normalized=False)
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::5])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::5])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlim([0, 2000])
        ax.set_ylim([0, yticks[-1]])
        ax.set_xlabel(r'$V_{CME}$ ($\frac{km}{s}$)')
        ax.set_ylabel(counts_type.title())
        if self.identifier is not None:
            ax.set_title(self.identifier, fontsize=titlesize)
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax, handles, labels

    def subview_fit_comparisons(self, ax, ticksize=7, labelsize=8, textsize=8, titlesize=10, highlight_tail=False, point_to_tail=False, counts_type='observed frequency'):
        """ """
        available_estimators = list(self.optimizations.keys())
        xticks = np.arange(0, int(self.x[-1] + 2*self.original_histogram.bin_widths[0] +1), 50)
        for estimator in available_estimators:
            optims = self.optimizations[estimator]
            optimizer, prms = optims['cls'], optims['optimized parameters']
            if counts_type == 'probability density':
                y = optimizer.f(prms, self.x)
                minor_yticks = np.arange(0, 0.0031, 0.00025)
                major_yticks = minor_yticks[::4]
                hist_label = 'Normalized Histogram'
            elif counts_type == 'observed frequency':
                y = optimizer.f(prms, self.x) * self.thresholded_histogram.normalization_constant
                hist_label = 'Histogram'
                if self.original_histogram.bin_widths[0] > 70:
                    minor_yticks = np.arange(0, 5001, 100)
                else:
                    minor_yticks = np.arange(0, 3001, 50)
                major_yticks = minor_yticks[::5]
            else:
                raise ValueError("invalid counts_type: {}; available counts_type = 'probability density' or 'observed frequency'".format(counts_type))
            ax.plot(self.x, y, color=self.estimator_colors[estimator], linewidth=1, linestyle='-', label=self.estimator_labels[estimator], alpha=1/len(available_estimators))
            if highlight_tail == True:
                ax, fill_handle, fill_label = self.subview_highlighted_tail(ax, y)
            if point_to_tail == True:
                ax = self.subview_arrow_to_tail(ax, textsize, normalized=normalized)
        ax.bar(self.thresholded_histogram.midpoints, self.thresholded_histogram.observed_counts, width=self.thresholded_histogram.bin_widths, align='center', facecolor='k', edgecolor='none', linewidth=0, alpha=0.6, label=hist_label)
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::5])
        ax.set_yticks(minor_yticks, minor=True)
        ax.set_yticks(major_yticks)
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlim([0, 2000])
        ax.set_ylim([0, minor_yticks[-1]])
        ax.set_xlabel(r'$V_{CME}$ ($\frac{km}{s}$)')
        ax.set_ylabel(counts_type.title())
        if self.identifier is not None:
            ax.set_title(self.identifier, fontsize=titlesize)
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_fit_surface(self, fig, ax, estimator, ticksize, labelsize, textsize, titlesize, cmap='plasma', elev=None, azim=None, antialiased=False, shade=False, show_contours=False):
        """ """
        optims = self.optimizations[estimator]
        optimizer, prms, objective = optims['cls'], optims['optimized parameters'], optims['objective']
        trans_prms = optimizer.from_normal_to_lognormal(prms)
        xlabel = r'$\mu_N$'
        ylabel = r'$\sigma_N$'
        unit_label = r'$\frac{km}{s}$'
        normal_label = r'{} $=$ ${:.2f}$ '.format(xlabel, prms[0]) + '{} '.format(unit_label) + r', {} $=$ ${:.2f}$ '.format(ylabel, prms[1]) + '{} '.format(unit_label)
        lognormal_label = r'$\mu$ $=$ ${:.2f}$ '.format(trans_prms[0]) + '{} '.format(unit_label) + r', $\sigma$ $=$ ${:.2f}$ '.format(trans_prms[1]) + '{} '.format(unit_label)
        if estimator in ('maximum likelihood',):
            central_color = 'k'
            if elev == 'auto':
                elev = -60 # 30
            if azim == 'auto':
                azim = -30 # 206.25
            objective *= -1
            sub_label = "{0:,.2f}".format(objective)
            zlabel = '$-${}'.format(self.estimator_symbols[estimator])
            objective_label = '{} $=$ ${}$'.format(zlabel, sub_label)
            zloc = -1
            zbound = -1 * self.round_up_to(-1 * objective, base=100)
        else:
            central_color = 'white'
            if elev == 'auto':
                elev = 30 # None # -60 # 30
            if azim == 'auto':
                azim = 0 # None # 30 # 345
            sub_label = "{0:,.2f}".format(objective)
            zlabel = self.estimator_symbols[estimator]
            objective_label = '{} $=$ ${}$'.format(zlabel, sub_label)
            zloc = 0
            zbound = self.round_up_to(objective, base=10)
        prms_label = '{}\n{}\n{}'.format(normal_label, lognormal_label, objective_label)
        X, Y, Z = optimizer.dim3_space[estimator]['X'], optimizer.dim3_space[estimator]['Y'], optimizer.dim3_space[estimator]['Z']
        z = optimizer.dim3_space[estimator]['z']
        norm = self.get_normalized_colormap(np.min(z), np.max(z), estimator=estimator)
        if show_contours == True:
            ax.contour(X, Y, Z, self.estimator_levels[estimator], norm=norm, linewidth=3, cmap=cmap, linestyles='solid', offset=-1)
            ax.contour(X, Y, Z, self.estimator_levels[estimator], norm=norm, linewidth=3, colors='k', linestyles='solid')
            ax.scatter(prms[0], prms[1], objective, color=central_color, marker='*', linewidth=0, s=50)
            opacity = 0.5
        else:
            opacity = 1
        surf_handle = ax.plot_surface(X, Y, Z, cmap=cmap, norm=norm, antialiased=antialiased, shade=shade, alpha=opacity)
        for ticklabel in ax.get_xticklabels():
            ticklabel.set_ha("center")
            ticklabel.set_va("center")
        for ticklabel in ax.get_yticklabels():
            ticklabel.set_ha("center")
            ticklabel.set_va("center")
        zticks = ax.get_zticks()
        mod_zticks = zticks.copy()
        mod_zticks[zloc] = int(zbound)
        fmt = FuncFormatter(lambda x, p: format(int(x), ','))
        ax.w_zaxis.set_major_formatter(fmt)
        # cbar = fig.colorbar(surf_handle, ax=ax, orientation='horizontal', shrink=0.75, pad=0.1, norm=norm, format=fmt, ticks=mod_zticks, extend='both')
        cbar = fig.colorbar(surf_handle, ax=ax, orientation='horizontal', shrink=0.75, pad=0.1, norm=norm, format=fmt, ticks=zticks, extend='both')
        cbar.ax.tick_params(labelsize=ticksize)
        # cbar.ax.set_title(objective_label, fontsize=labelsize)
        cbar.ax.set_title(prms_label, fontsize=labelsize)
        for ax_ijk in (ax.xaxis, ax.yaxis, ax.zaxis):
            ax_ijk.set_rotate_label(False)
            ax_ijk._axinfo['label']['space_factor'] = 3.0
            ax_ijk.set_rotate_label(False)
        # ax.xaxis._axinfo['label']['space_factor'] = 0.
        # ax.yaxis._axinfo['label']['space_factor'] = 0.
        # ax.zaxis._axinfo['label']['space_factor'] = 0.
        # ax.xaxis._axinfo['tick']['inward_factor'] = 0.6
        # ax.xaxis._axinfo['tick']['outward_factor'] = 0.2
        # ax.yaxis._axinfo['tick']['inward_factor'] = 0.6
        # ax.yaxis._axinfo['tick']['outward_factor'] = 0.2
        # ax.zaxis._axinfo['tick']['inward_factor'] = 0.2
        # ax.zaxis._axinfo['tick']['outward_factor'] = 1.
        ax.set_xlabel(xlabel, rotation=0, labelpad=0) # labelpad=-5
        ax.set_ylabel(ylabel, rotation=0, labelpad=0) # labelpad=-5
        # ax.set_zlabel(zlabel, rotation=0)
        if ((elev is not None) or (azim is not None)):
            ax.view_init(elev, azim)
        if self.identifier is not None:
            if estimator in ('maximum likelihood',):
                if ((elev == -60) and (azim == -30)):
                    ytit = 1.15
                else:
                    ytit = 1.1
            else:
                ytit = 1
            ax.set_title(self.identifier, fontsize=titlesize, y=ytit)
        ax = self.set_tick_size(ax, ticksize, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize, labelsize)
        return fig, ax, cbar

    def subview_fit_contours(self, ax, estimator, norm, xticks, yticks, ticksize, labelsize, textsize, titlesize, cmap='plasma'):
        """ """
        optims = self.optimizations[estimator]
        optimizer, prms, objective = optims['cls'], optims['optimized parameters'], optims['objective']
        trans_prms = optimizer.from_normal_to_lognormal(prms)
        unit_label = r'$\frac{km}{s}$'
        xlabel = r'$\mu_N$'
        ylabel = r'$\sigma_N$'
        normal_label = r'{} $=$ ${:.2f}$ '.format(xlabel, prms[0]) + '{} '.format(unit_label) + r', {} $=$ ${:.2f}$ '.format(ylabel, prms[1]) + '{} '.format(unit_label)
        lognormal_label = r'$\mu$ $=$ ${:.2f}$ '.format(trans_prms[0]) + '{} '.format(unit_label) + r', $\sigma$ $=$ ${:.2f}$ '.format(trans_prms[1]) + '{} '.format(unit_label)
        if estimator in ('maximum likelihood',):
            central_color = 'k'
            sub_label = "{0:,.2f}".format(objective * -1)
            objective_label = '$-${} $=$ ${}$'.format(self.estimator_symbols[estimator], sub_label)
        else:
            central_color = 'white'
            sub_label = "{0:,.2f}".format(objective)
            objective_label = '{} $=$ ${}$'.format(self.estimator_symbols[estimator], sub_label)
        prms_label = '{}\n{}\n{}'.format(normal_label, lognormal_label, objective_label)
        X, Y, Z = optimizer.dim3_space[estimator]['X'], optimizer.dim3_space[estimator]['Y'], optimizer.dim3_space[estimator]['Z']
        fill_handle = ax.contourf(X, Y, Z, self.estimator_levels[estimator], cmap=cmap, norm=norm)
        line_handle = ax.contour(X, Y, Z, self.estimator_levels[estimator], colors=central_color, linewidths=0.5)
        prms_handle = ax.scatter(prms[0], prms[1], color=central_color, marker='*', linewidth=0, s=50, label=prms_label)
        z_posmin = np.min(np.abs(Z.reshape(-1)))
        if z_posmin > 10:
            ax.clabel(line_handle, inline=True, fontsize=ticksize, fmt=FuncFormatter(lambda x, p: format(int(x), ',')))
        else:
            if z_posmin > 1:
                ax.clabel(line_handle, inline=True, fontsize=ticksize, fmt='%.0f')
            else:
                ax.clabel(line_handle, inline=True, fontsize=ticksize, fmt='%.1f')
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::2])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlabel(xlabel) #, labelpad=20)
        ax.set_ylabel(ylabel) #, labelpad=20)
        ax.set_xlim([xticks[0], xticks[-1]])
        ax.set_ylim([yticks[0], yticks[-1]])
        leg = ax.legend(handles=[prms_handle], labels=[prms_label], bbox_to_anchor=(0.5, -0.15), loc='upper center', ncol=1, fontsize=textsize)
        leg._legend_box.align = "center"
        if central_color in ('w', 'white'):
            frame = leg.get_frame()
            frame.set_facecolor('darkorange')
            frame.set_edgecolor('k')
        if self.identifier is not None:
            ax.set_title(self.identifier, fontsize=titlesize)
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

class InterExceedanceDistribution(VisualEditor):

    def __init__(self, inter_exceedance_times, w=15, max_edge=301):
        """ """
        super().__init__()
        self.inter_exceedance_times = inter_exceedance_times
        self.edges = np.arange(0, max_edge, w, dtype=int)
        self._time_delta_histogram = None
        self._inverse_sample_histogram = None

    def initialize_time_delta_histogram(self, bias='left'):
        """ """
        AH = AdaptiveHistogram()
        AH.initialize_histogram(self.inter_exceedance_times, edges=self.edges, bin_threshold=None, bias=bias)
        self._time_delta_histogram = AH

    @property
    def time_delta_histogram(self):
        return self._time_delta_histogram

    def initialize_inverse_sample_histogram(self, bias='left'):
        """ """
        x = np.random.uniform(size=self.inter_exceedance_times.size)
        intermediate_term = np.log(1/x)
        y = intermediate_term * np.mean(self.inter_exceedance_times) / np.mean(intermediate_term)
        AH = AdaptiveHistogram()
        AH.initialize_histogram(y, edges=self.edges, bin_threshold=None, bias=bias)
        self._inverse_sample_histogram = AH

    @property
    def inverse_sample_histogram(self):
        return self._inverse_sample_histogram

    def subview_histogram_comparison(self, ax, ticksize=7, labelsize=8, textsize=8, counts_type='observed frequency'):
        """ """
        histograms = (self.inverse_sample_histogram, self.time_delta_histogram)
        if counts_type == 'observed frequency':
            get_y = lambda histogram : histogram.observed_counts
        elif counts_type == 'probability density':
            get_y = lambda histogram : histogram.normalized_counts
        else:
            raise ValueError("counts_type = 'observed frequency' or 'probability density'")
        yticks = np.arange(0, 501, 25).astype(int)
        facecolors = ('r', 'b')
        labels = ('Inverse-Transform Sample', 'Time-Deltas')
        for histogram, facecolor, label in zip(histograms, facecolors, labels):
            ax.bar(histogram.midpoints, get_y(histogram), width=histogram.bin_widths, facecolor=facecolor, label=label, alpha=0.5, align='center')
        ax.set_xticks(self.time_delta_histogram.edges[1::2], minor=True)
        ax.set_xticks(self.time_delta_histogram.edges[::2])
        # ax.set_yticks(yticks[1::2], minor=True)
        # ax.set_yticks(yticks[::2])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::3])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        text_box = ax.text(0.95, 0.95, '${:,}$ Extreme Events'.format(self.inter_exceedance_times.size +1), fontsize=textsize, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlim([self.time_delta_histogram.edges[0], self.time_delta_histogram.edges[-1]])
        ax.set_xlabel('Elapsed Hours in-between Extreme Events')
        ax.set_ylabel(counts_type.title())
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

class MaxSpectrum(VisualEditor):

    def __init__(self):
        super().__init__()

    @staticmethod
    def get_jth_kth_index(j):
        """
        j           :   type <int>
        """
        jprime = 2**j
        return (-1, jprime)

    @staticmethod
    def get_shapeshifted_data(data, jprime):
        """
        j           :   type <int>
        """
        desired_size_factor = np.prod([n for n in jprime if n != -1])
        if -1 in jprime:  ## IMPLICIT ARRAY SIZE
            desired_size = data.size // desired_size_factor * desired_size_factor
        else:
            desired_size = desired_size_factor
        return data.flat[:desired_size].reshape(jprime)

    @staticmethod
    def get_data_maximums(data, pad_char):
        """
        vshift      :   type <array>
        pad_char    :   type <int / float> or None
        """
        res = np.max(data.T, axis=0)
        if pad_char is not None:
            res = res[res != pad_char]
        return res[res > 0]

    def get_max_spectrum(self, js, data, pad_char=1, ddof=0):
        """ """
        keys = ('dlog', 'npad', 'standard deviation', 'standard error', 'initial Y(j)')
        res = {key : [] for key in keys}
        for j in js:
            jprime = self.get_jth_kth_index(j)
            shapeshifted_data = self.get_shapeshifted_data(data, jprime)
            partial_spectrum = self.get_data_maximums(shapeshifted_data, pad_char)
            dlog = np.log2(partial_spectrum)
            st_dev = np.std(dlog, ddof=ddof)
            st_err = sem(dlog, ddof=ddof)
            yj_init = np.sum(dlog) / dlog.size # np.mean(dlog)
            args = (dlog, dlog.size, st_dev, st_err, yj_init)
            for key, arg in zip(keys, args):
                res[key].append(arg)
        return {key : np.array(value) for key, value in res.items()}

    @staticmethod
    def get_optimized_estimators(js, max_spectrum, js_subinterval, local_method='Nelder-Mead', global_method=None, **kwargs):
        """ """
        xi = js[js_subinterval]
        yi = max_spectrum['initial Y(j)'][js_subinterval]
        weights = max_spectrum['npad'][js_subinterval]
        GLS = GeneralizedLeastSquares(xi, yi)
        # fit_results = GLS.fit(weights='residual', local_method=local_method, global_method=global_method, **kwargs)
        fit_results = GLS.fit(weights=weights, local_method=local_method, global_method=global_method, **kwargs)
        yj_fit = GLS.f(fit_results.x, xi)
        return {'fitted Y(j)' : yj_fit, 'alpha' : 1/fit_results.x[0], 'intercept' : fit_results.x[1]}

class UnbiasedEstimators(MaxSpectrum):

    def __init__(self, speeds):
        """ """
        super().__init__()
        self.speeds = speeds
        self._nresamples = None
        self.max_spectra = []
        self._alphas = []
        self._intercepts = []
        self._js_subinterval = None
        self._js_prime = None
        self._alpha_hat_histogram = None
        self._thetas = None
        self._theta_hat_histogram = None
        self.point_estimators = {}
        self.speed_thresholds = {}

    @property
    def nresamples(self):
        return self._nresamples

    @property
    def alphas(self):
        return np.array(self._alphas)

    @property
    def intercepts(self):
        return np.array(self._intercepts)

    @property
    def jmax(self):
        return int(np.floor(np.log2(self.speeds.size)))

    @property
    def js(self):
        return np.linspace(1, self.jmax, self.jmax, dtype=int)

    @staticmethod
    def resample_indices(indices, nshuffles=3):
        """ """
        if nshuffles > 0:
            for n in range(nshuffles):
                np.random.shuffle(indices)
        np.random.shuffle(indices)
        return indices

    def store_max_spectra(self, js_subinterval, nshuffles=3, nresamples=1000, pad_char=1, ddof=0, local_method='Nelder-Mead', global_method=None, **kwargs):
        """ """
        self._js_subinterval = js_subinterval
        self._js_prime = self.js[js_subinterval]
        max_spectrum = self.get_max_spectrum(self.js, self.speeds, pad_char, ddof)
        opt_spectrum = self.get_optimized_estimators(self.js, max_spectrum, js_subinterval, local_method, global_method, **kwargs)
        max_spectrum.update(opt_spectrum)
        self.max_spectra.append(max_spectrum)
        self._alphas.append(max_spectrum['alpha'])
        self._intercepts.append(max_spectrum['intercept'])
        indices = np.arange(self.speeds.size)
        for ith_resample in range(nresamples):
            indices = self.resample_indices(indices, nshuffles)
            max_spectrum = self.get_max_spectrum(self.js, self.speeds[indices], pad_char, ddof)
            opt_spectrum = self.get_optimized_estimators(self.js, max_spectrum, js_subinterval, local_method, global_method, **kwargs)
            max_spectrum.update(opt_spectrum)
            self.max_spectra.append(max_spectrum)
            self._alphas.append(max_spectrum['alpha'])
            self._intercepts.append(max_spectrum['intercept'])
        self._nresamples = nresamples

    @property
    def js_subinterval(self):
        return self._js_subinterval

    @property
    def js_prime(self):
        return self._js_prime

    def subview_power_law_errors(self, ax, ticksize=7, labelsize=8, textsize=8, ith_resample=0):
        """ """
        xlim = np.array([0, 17])
        ylim = np.array([7, 12])
        yticks = np.arange(ylim[0], ylim[-1]+0.1, 0.5)
        max_spectrum = self.max_spectra[ith_resample]
        for key, facecolor in zip(('standard deviation', 'standard error'), ('red', 'blue')):
            ax.errorbar(self.js, max_spectrum['initial Y(j)'], yerr=max_spectrum[key], alpha=0.7, capsize=5, label=key.title(), ecolor=facecolor, fmt='none')
        ax.set_xticks(self.js[1::2], minor=True)
        ax.set_xticks(self.js[::2])
        ax.set_yticks(yticks[1::2], minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(r'$j$ ($log_2$ hours)')
        ax.set_ylabel(r'$Y(j)$ ($log_2$ $\frac{km}{s}$)')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_linear_power_law(self, ax, ticksize=7, labelsize=8, textsize=8, ith_resample=0):
        """ """
        xlim = np.array([0.1, 17])
        ylim = np.array([7, 12])
        yticks = np.arange(ylim[0], ylim[-1]+0.1, 0.5)
        max_spectrum = self.max_spectra[ith_resample]
        ax.scatter(self.js, max_spectrum['initial Y(j)'], color='gray', marker='.', label='Initial Max Spectrum')
        ax.plot(self.js_prime, max_spectrum['fitted Y(j)'], color='red', label='Power-Law')
        ax.set_xticks(self.js[1::2], minor=True)
        ax.set_xticks(self.js[::2])
        ax.set_yticks(yticks[1::2], minor=True)
        ax.set_yticks(yticks[::2])
        # ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(r'$j$ ($log_2$ hours)')
        ax.set_ylabel(r'$Y(j)$ ($log_2$ $\frac{km}{s}$)')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_logarithmic_power_law(self, ax, ticksize=7, labelsize=8, textsize=8, ith_resample=0, identifier=None):
        """ """
        xlim = 2 ** np.array([0.1, 17])
        ylim = 2 ** np.array([8, 12])
        # yticks = 2 ** np.arange(ylim[0], ylim[-1]+0.1, 0.5)
        max_spectrum = self.max_spectra[ith_resample]
        ax.set_xscale('log', basex=2)
        ax.set_yscale('log', basey=2)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.plot(2 ** self.js_prime, max_spectrum['fitted Y(j)'], color='red')
        ax.scatter(2 ** self.js, max_spectrum['initial Y(j)'], color='gray', marker='.')
        ax.set_xticks(2 ** self.js[::3])
        # ax.set_yticks(2 ** yticks[::3])
        ax.tick_params(axis='x', colors='gray')
        ax.tick_params(axis='y', colors='gray')
        ax.xaxis.label.set_color('gray')
        ax.yaxis.label.set_color('gray')
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.grid(color='k', alpha=0.3, linestyle=':', which='both')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        alpha_label = r'mean($\hat\alpha) = {0:.3f}$'.format(np.mean(self.alphas[1:]))
        # alpha_label = r'$\hat\alpha = {0:.2f}$'.format(self.alphas[ith_resample])
        if identifier is None:
            tlabel = alpha_label
        else:
            tlabel = '{}:\n{}'.format(identifier, alpha_label)
        text_box = ax.text(0.05, 0.95, tlabel, fontsize=textsize, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlabel(r'$2^j$ (hours)')
        ax.set_ylabel(r'$2^{Y(j)}$ ($\frac{km}{s}$)')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def initialize_alpha_hat_histogram(self, edges=None, n=None, w=None, bias='left'):
        """ """
        AH = AdaptiveHistogram()
        AH.initialize_histogram(self.alphas[1:], edges=edges, n=n, w=w, bin_threshold=None, bias=bias, exclude_empty_boundaries=False)
        self._alpha_hat_histogram = AH

    @property
    def alpha_hat_histogram(self):
        return self._alpha_hat_histogram

    def subview_alpha_hat_histogram(self, ax, ticksize=7, labelsize=8, textsize=8, counts_type='observed frequency', transparency=1, identifier=None):
        """ """
        histogram = self.alpha_hat_histogram
        if counts_type == 'observed frequency':
            y = histogram.observed_counts
            maj_mult = 20
            min_mult = 10
        elif counts_type == 'probability density':
            y = histogram.normalized_counts
            maj_mult = 2
            min_mult = 1
        else:
            raise ValueError("counts_type = 'observed frequency' or 'probability density'")
        skew_label = r'skew($\hat\alpha) = {0:.2f}$'.format(histogram.skew)
        kurt_label = r'kurtosis($\hat\alpha) = {0:.2f}$'.format(histogram.kurtosis)
        if identifier is None:
            legend_label = '{}\n{}'.format(skew_label, kurt_label)
        else:
            legend_label = '{}\n{}\n{}'.format(identifier, skew_label, kurt_label)
        ax.bar(histogram.midpoints, y, width=histogram.bin_widths, align='center', label=legend_label, alpha=transparency)
        ax.set_xticks(histogram.edges, minor=True)
        ax.set_xticks(histogram.edges[::5])
        ax.yaxis.set_major_locator(MultipleLocator(maj_mult))
        ax.yaxis.set_minor_locator(MultipleLocator(min_mult))
        ax.set_xlim([histogram.edges[0], histogram.edges[-1]])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlabel(r'$\hat\alpha$')
        ax.set_ylabel(counts_type.title())
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def initialize_thetas(self):
        exponents = []
        original_max_spectrum = self.max_spectra[0]
        for max_spectrum in self.max_spectra[1:]:
            delta = original_max_spectrum['fitted Y(j)'] - max_spectrum['fitted Y(j)']
            exponent = max_spectrum['alpha'] * delta.T
            exponents.append(exponent)
        self._thetas = 2 ** np.array(exponents)

    @property
    def thetas(self):
        return self._thetas

    def initialize_point_estimators(self):
        keys = ('mean', 'median', 'maximum', 'minimum')
        fs = (np.mean, np.median, np.max, np.min)
        point_estimators = {key : f(self.thetas, axis=0) for key, f in zip(keys, fs)}
        self.point_estimators.update(point_estimators)

    def subview_point_estimators(self, ax, ticksize=7, labelsize=8, textsize=8):
        """ """
        yticks = np.arange(0, 1.01, 0.05) # np.arange(0.2, 0.81, 0.05)
        for key, label in zip(('maximum', 'minimum'), ('max/min', None)):
            ax.scatter(self.js_prime, self.point_estimators[key], color='k', label=label, marker='_')
        for key, facecolor, marker in zip(('mean', 'median'), ('steelblue', 'darkorange'), ('o', '*')):
            ax.scatter(self.js_prime, self.point_estimators[key], color=facecolor, alpha=0.5, label=r'{} $\theta(j)$'.format(key.title()), marker=marker)
        ax.set_xticks(self.js_prime[1::2], minor=True)
        ax.set_xticks(self.js_prime[::2])
        ax.set_yticks(yticks[1::2], minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlabel(r'$j$ ($log_2$ hours)')
        ax.set_ylabel(r'$\theta(j)$')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def initialize_theta_hat_histogram(self, edges=None, n=None, w=None, bias='left'):
        """ """
        AH = AdaptiveHistogram()
        AH.initialize_histogram(self.thetas.copy().reshape(-1), edges=edges, n=n, w=w, bin_threshold=None, bias=bias)
        self._theta_hat_histogram = AH

    @property
    def theta_hat_histogram(self):
        return self._theta_hat_histogram

    def subview_theta_hat_histogram(self, ax, ticksize=7, labelsize=8, textsize=8, counts_type='observed frequency', transparency=1, identifier=None):
        """ """
        histogram = self.theta_hat_histogram
        if counts_type == 'observed frequency':
            y = histogram.observed_counts
            yticks = np.arange(0, 4001, 250).astype(int)
        elif counts_type == 'probability density':
            y = histogram.normalized_counts
            yticks = np.arange(0, 21, 2.5).astype(int)
        else:
            raise ValueError("counts_type = 'observed frequency' or 'probability density'")
        mean_label = r'mean($\hat\theta$) = ${0:.2f}$'.format(np.mean(histogram.data))
        median_label = r'median($\hat\theta$) = ${0:.2f}$'.format(np.median(histogram.data))
        if identifier is None:
            legend_label = '{}\n{}'.format(mean_label, median_label)
        else:
            legend_label = '{}\n{}\n{}'.format(identifier, mean_label, median_label)
        ax.bar(histogram.midpoints, y, width=histogram.bin_widths, align='center', label=legend_label, alpha=transparency)
        ax.set_xticks(histogram.edges[1::2], minor=True)
        ax.set_xticks(histogram.edges[::2])
        ax.set_yticks(yticks[1::2], minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlim([0, 1])
        ax.set_xlabel(r'$\hat\theta$')
        ax.set_ylabel(counts_type.title())
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def initialize_speed_thresholds(self):
        """ """
        yj = np.array([max_spectrum['fitted Y(j)'] for max_spectrum in self.max_spectra])
        keys = ('mean', 'median')
        fs = (np.mean, np.median)
        for key, f in zip(keys, fs):
            vthreshold = np.min(2 ** f(yj, axis=0))
            self.speed_thresholds[key] = int(round(vthreshold))

class ClusterParametrization(VisualEditor):

    def __init__(self, inter_exceedance_times):
        """ """
        super().__init__()
        self.inter_exceedance_times = inter_exceedance_times
        self.available_bias = ('first-order', 'threshold', 'baseline')
        self.moment_estimators = {}
        self.time_thresholds = {}
        self._density_optimizer = None
        self._agglomerative_optimizer = None
        self._clusters = None
        self._bias = None

    @property
    def bias(self):
        return self._bias

    def initialize_bias(self, bias='threshold'):
        self._bias = bias

    def initialize_moment_estimators(self, baseline=None):
        """ """
        self.moment_estimators['first-order'] = 2 * np.sum(self.inter_exceedance_times - 1)**2 / (self.inter_exceedance_times.size * np.sum((self.inter_exceedance_times - 1) * (self.inter_exceedance_times - 2)))
        self.moment_estimators['threshold'] = 2 * np.sum(self.inter_exceedance_times)**2 / (self.inter_exceedance_times.size * np.sum(self.inter_exceedance_times**2))
        self.moment_estimators['baseline'] = baseline

    def initialize_time_thresholds(self, baseline=None):
        """ """
        values = np.sort(self.inter_exceedance_times)
        size = self.inter_exceedance_times.size + 1
        for bias in list(self.moment_estimators.keys()):
            if self.moment_estimators[bias] is None:
                self.time_thresholds[bias] = baseline
            else:
                index = int(np.floor(self.moment_estimators[bias] * size))
                try:
                    self.time_thresholds[bias] = values[index]
                    # self.time_thresholds[bias] = int(np.round(self.extremal_moment_estimator[bias] * size))
                except:
                    self.time_thresholds[bias] = np.nan

    @property
    def bias_colors(self):
        return {bias : facecolor for bias, facecolor in zip(self.available_bias, ('darkgreen', 'darkorange', 'steelblue'))}

    def subview_moment_estimators(self, ax, speed_thresholds, moment_estimators, ticksize=7, labelsize=8):
        """ """
        yticks = np.arange(0, 1.01, 0.1)
        biases = list(moment_estimators.keys())
        for bias in biases:
            if bias == 'baseline':
                label = bias.title()
            else:
                label = '{} Bias'.format(bias.title())
            ax.plot(speed_thresholds, moment_estimators[bias], color=self.bias_colors[bias], linestyle='-', alpha=1/len(biases), label=label)
        ax.set_xticks(speed_thresholds[::10], minor=True)
        ax.set_xticks(speed_thresholds[::50])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::2])
        ax.set_xlim([speed_thresholds[0], speed_thresholds[-1]])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlabel(r'Speed Threshold $(\frac{km}{s})$')
        ax.set_ylabel(r'Moment Estimator $\hat\theta$')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_time_thresholds(self, ax, speed_thresholds, time_thresholds, ticksize=7, labelsize=8):
        """ """
        yticks = np.arange(0, 5001, 25)
        biases = list(time_thresholds.keys())
        for bias in biases:
            if bias == 'baseline':
                label = bias.title()
            else:
                label = '{} Bias'.format(bias.title())
            ax.plot(speed_thresholds, time_thresholds[bias], color=self.bias_colors[bias], linestyle='-', alpha=1/len(biases), label=label)
        ax.set_xticks(speed_thresholds[::10], minor=True)
        ax.set_xticks(speed_thresholds[::50])
        ax.set_yscale('log', basey=5)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        # ax.set_yticks(yticks, minor=True)
        # ax.set_yticks(yticks[::6])
        ax.set_xlim([speed_thresholds[0], speed_thresholds[-1]])
        ax.set_ylim([yticks[0], yticks[-1]])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlabel(r'Speed Threshold $(\frac{km}{s})$')
        ax.set_ylabel(r'Time Threshold (hours)')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def get_clusters(self, data, bias):
        """ """
        # condition = (self.time_thresholds[bias] < self.inter_exceedance_times)
        condition = (self.time_thresholds[bias] <= self.inter_exceedance_times)
        indices = np.where(condition)[0] + 1
        return {key : np.array(np.split(value, indices)) for key, value in data.items()}

    @staticmethod
    def select_optimized_groups(groups, keys, include_levels=False):
        """ """
        _x = [groups[key][0] for key in keys]
        x = np.concatenate(_x, axis=0)
        if include_levels == True:
            levels = np.concatenate([[jdx for idx in range(len(_x[jdx]))] for jdx in range(len(_x))], axis=0)
            return x, levels
        else:
            return x

    def initialize_density_optimization(self, bias, radius=None, minimum_npoints=None, method='dbscan', identifiers='auto', distance_metric='euclidean'):
        """ """
        if radius is None:
            radius = self.time_thresholds[bias]
        optimizer = DensityMethods(np.array([self.inter_exceedance_times]))
        optimizer.dispatch(method, identifiers=identifiers, distance_metric=distance_metric, radius=radius, minimum_npoints=minimum_npoints)
        self._density_optimizer = optimizer

    @property
    def density_optimizer(self):
        return self._density_optimizer

    def subview_density_optimized_clusters(self, ax, cmap='plasma', show_clusters=False, show_noise=False, label_clusters=False, label_noise=False, zoom_clusters=False, ticksize=7, labelsize=8, textsize=8):
        """ """
        cluster_keys = np.arange(self.density_optimizer.k).astype(str)
        search_inputs = np.array([show_clusters, show_noise])
        if ((self.density_optimizer.identification_labels is None) or (isinstance(self.density_optimizer.identification_labels, str))):
            identification_labels = np.array([self.density_optimizer.identification_labels for key in cluster_keys])
        else:
            identification_labels = np.array(self.density_optimizer.identification_labels)
        nlabels = len(identification_labels)
        if nlabels < self.density_optimizer.k:
            raise ValueError("{} labels for {} clusters".format(nlabels, self.k))
        x_clusters = None
        if np.any(search_inputs == True):
            if show_clusters == True:
                x_clusters, x_levels = self.select_optimized_groups(self.density_optimizer.clusters, self.density_optimizer.cluster_keys, include_levels=True)
                ax.scatter(x_clusters, np.zeros(x_clusters.shape).astype(int), cmap=cmap, c=x_levels, alpha=0.8, label='clusters', marker='.')
                if ((label_clusters == True) and (np.all(identification_labels is not None))):
                    facecolors = self.get_discrete_facecolors_from_colormap(cmap, density_optimizer.k)
                    for key, c, label in zip(self.density_optimizer.cluster_keys, facecolors, identification_labels):
                        xi = self.density_optimizer.clusters[key][0]
                        ax.text(np.median(xi), 2.5, label, color=c, fontsize=textsize)
            if show_noise == True:
                x_noise = self.select_optimized_groups(self.density_optimizer.noise, self.density_optimizer.noise_keys, include_levels=False)
                ax.scatter(x_noise, np.zeros(x_noise.shape).astype(int), facecolor='k', alpha=0.3, label='noise', marker='.')
                if label_noise == True:
                    for key in self.density_optimizer.noise_keys:
                        xi = self.density_optimizer.noise[key][0]
                        ax.text(np.median(xi), 2.5, 'noise', color=c, fontsize=textsize)
        else:
            ax.scatter(self.inter_exceedance_times, np.zeros(self.inter_exceedance_times.shape).astype(int), facecolor='k', label='data', marker='.')
        ax.tick_params(axis='x', which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=False, right=False, labelleft=False, labelright=False)
        ax.set_yticks([])
        if zoom_clusters == True:
            if x_clusters is None:
                x_clusters = self.select_optimized_groups(self.density_optimizer.clusters, self.density_optimizer.cluster_keys)
            xmax = np.max(x_clusters)
            _xmin = np.min(x_clusters)
            xmin = self.round_down_to(_xmin, 10)
            if xmin == 1:
                xmin = 0
            xticks = np.arange(xmin, xmax +50, 50).astype(int)
            ax.set_xticks(xticks, minor=True)
            ax.set_xticks(xticks[::4])
            ax.set_xlim([xmin, xmax +50])
        else:
            xmax = np.max(self.inter_exceedance_times)
            xticks = np.arange(0, xmax +100, 10).astype(int)
            if len(xticks) <= 10:
                ax.set_xticks(xticks[::10], minor=True)
                ax.set_xticks(xticks[::50])
            else:
                ax.set_xticks(xticks[::25], minor=True)
                ax.set_xticks(xticks[::150])
            ax.set_xlim([-10, xmax +100])
        ax.set_ylim([-10, 10])
        text_box = ax.text(0.95, 0.95, '${}$ Clusters via ${}$ Events'.format(self.density_optimizer.k, self.inter_exceedance_times.size +1), fontsize=textsize, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlabel('Inter-Exceedance Times (hours)')
        ax.set_ylabel('')
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def initialize_agglomerative_optimization(self, method='single', identifiers='auto', distance_metric='euclidean'):
        """ """
        optimizer = AgglomerativeHierarchicalMethods(np.array([self.inter_exceedance_times]))
        optimizer.dispatch(method, identifiers, distance_metric=distance_metric)
        self._agglomerative_optimizer = optimizer

    @property
    def agglomerative_optimizer(self):
        return self._agglomerative_optimizer

    def subview_agglomerative_optimized_dendrogram(self, ax, facecolors, label_dendrogram=False, tick_connectors=False, ticksize=7, labelsize=8):
        """ """
        ymax = np.max(self.agglomerative_optimizer.levels)
        if ymax <= 1000:
            yticks = np.arange(0, ymax +50, 50)
        elif ymax <= 3000:
            yticks = np.arange(0, ymax +125, 125)
        else:
            yticks = np.arange(0, ymax +250, 250)
        ax = self.agglomerative_optimizer.subview_dendrogram(ax, facecolors, alpha=0.8, label_dendrogram=label_dendrogram, tick_connectors=tick_connectors, ticksize=ticksize, labelsize=labelsize)
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::4])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_ylabel('Levels (hours)')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax


class ClusteringStatistics(VisualEditor):

    def __init__(self, clusters):
        """ """
        super().__init__()
        self.clusters = clusters
        self._intra_times = []
        self._intra_durations = []
        self._inter_durations = []
        self._intra_times_histogram = None
        self._intra_durations_histogram = None
        self._inter_durations_histogram = None
        self.relative_statistics = {}
        self.mean_duration = {}
        self.median_duration = {}
        self.stdev_duration = {}
        self.max_duration = {}
        self.min_duration = {}
        self._duration_labels = []

    def initialize_times_and_durations(self):
        for ith_cluster, cluster in enumerate(self.clusters['elapsed hour']):
            intra_times = np.diff(cluster)
            intra_duration = cluster[-1] - cluster[0]
            try:
                next_cluster = self.clusters['elapsed hour'][ith_cluster +1]
                inter_duration = next_cluster[0] - cluster[-1]
                self._inter_durations.append(inter_duration)
            except:
                pass
            self._intra_times.append(intra_times)
            self._intra_durations.append(intra_duration)

    @property
    def intra_times(self):
        return np.array(self._intra_times)

    @property
    def intra_durations(self):
        return np.array(self._intra_durations)

    @property
    def inter_durations(self):
        return np.array(self._inter_durations)

    def initialize_intra_times_histogram(self, bias='left'):
        """ """
        flat_times = np.concatenate(self.intra_times, axis=0)
        # edges = np.arange(0, max(flat_times) + 1)
        edges = np.arange(-0.5, max(flat_times) + 1.5)
        AH = AdaptiveHistogram()
        AH.initialize_histogram(flat_times, edges=edges, n=None, w=None, bin_threshold=None, bias=bias)
        self._intra_times_histogram = AH

    @property
    def intra_times_histogram(self):
        return self._intra_times_histogram

    def initialize_intra_durations_histogram(self, bias='left'):
        """ """
        # edges = np.arange(0, max(self.intra_durations) + 1)
        edges = np.arange(-0.5, max(self.intra_durations) + 1.5)
        AH = AdaptiveHistogram()
        AH.initialize_histogram(self.intra_durations, edges=edges, n=None, w=None, bin_threshold=None, bias=bias)
        self._intra_durations_histogram = AH

    @property
    def intra_durations_histogram(self):
        return self._intra_durations_histogram

    def initialize_inter_durations_histogram(self, bias='left'):
        """ """
        iqr_bounds = np.percentile(self.inter_durations, np.array([25, 75]))
        iqr_delta = max(iqr_bounds) - min(iqr_bounds)
        wbin = (iqr_delta * 2) // self.inter_durations.size ** (1/3)
        lbin = np.floor(np.min(self.inter_durations))
        rbin = np.ceil(np.max(self.inter_durations))
        nbin = 2 * int((rbin - lbin) // wbin)
        edges = np.linspace(lbin, rbin, nbin, dtype=int)
        AH = AdaptiveHistogram()
        AH.initialize_histogram(self.inter_durations, edges=edges, n=None, w=None, bin_threshold=None, bias=bias)
        self._inter_durations_histogram = AH

    @property
    def inter_durations_histogram(self):
        return self._inter_durations_histogram

    def subview_intra_time_histogram(self, ax, ticksize=7, labelsize=8, textsize=8, counts_type='observed frequency', time_threshold=None):
        """ """
        if counts_type == 'observed frequency':
            y = self.intra_times_histogram.observed_counts
            yticks = np.arange(0, max(y) * 1.5, 5)
        elif counts_type == 'probability density':
            y = self.intra_times_histogram.normalized_counts
            yticks = np.arange(0, max(y) * 1.5, 1)
        else:
            raise ValueError("counts_type = 'observed frequency' or 'probability density'")
        ax.bar(self.intra_times_histogram.midpoints, y, width=self.intra_times_histogram.bin_widths, facecolor='mediumpurple', alpha=0.5, align='center')
        if time_threshold is not None:
            ax.axvline(time_threshold, ymin=0, ymax=yticks[-1], color='darkorange', linestyle='-', marker=None, label='$T_C = {}$ hours'.format(time_threshold))
        ax.set_xticks(self.intra_times_histogram.midpoints, minor=True)
        ax.set_xticks(self.intra_times_histogram.midpoints[1::2])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlim([0, self.intra_times_histogram.edges[-1] +1])
        ax.legend(loc='upper left', fontsize=textsize)
        # text_box = ax.text(0.05, 0.95, '$T_C = {}$ hours'.format(time_threshold), fontsize=textsize, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
        # text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlabel('Intra-Event Times (hours)')
        ax.set_ylabel(counts_type.title())
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_intra_duration_histogram(self, ax, ticksize=7, labelsize=8, textsize=8, counts_type='observed frequency', time_threshold=None):
        """ """
        if counts_type == 'observed frequency':
            y = self.intra_durations_histogram.observed_counts
            ymax = max(y)
            if ymax <= 10:
                yticks = np.arange(0, ymax * 1.3, 2)
            else:
                yticks = np.arange(0, ymax * 1.3, 5)
        elif counts_type == 'probability density':
            y = self.intra_durations_histogram.normalized_counts
            yticks = np.arange(0, max(y) * 1.3, 1)
        else:
            raise ValueError("counts_type = 'observed frequency' or 'probability density'")
        ax.bar(self.intra_durations_histogram.midpoints, y, width=self.intra_durations_histogram.bin_widths, facecolor='darkgreen', alpha=0.5, align='center')
        if time_threshold is not None:
            ax.axvline(time_threshold, ymin=0, ymax=yticks[-1], color='darkorange', linestyle='-', marker=None, label='$T_C = {}$ hours'.format(time_threshold))
        ax.set_xticks(self.intra_durations_histogram.midpoints, minor=True)
        if self.intra_durations_histogram.edges[-1] <= 50:
            ax.set_xticks(self.intra_durations_histogram.midpoints[::5])
        elif self.intra_durations_histogram.edges[-1] <= 100:
            ax.set_xticks(self.intra_durations_histogram.midpoints[::10])
        else:
            ax.set_xticks(self.intra_durations_histogram.midpoints[::20])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        ax.set_xlim([0, self.intra_durations_histogram.edges[-1] +1])
        ax.legend(loc='upper right', fontsize=textsize)
        # text_box = ax.text(0.05, 0.95, '$T_C = {}$ hours'.format(time_threshold), fontsize=textsize, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
        # text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlabel('Intra-Cluster Durations (hours)')
        ax.set_ylabel(counts_type.title())
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_inter_duration_histogram(self, ax, ticksize=7, labelsize=8, textsize=8, counts_type='observed frequency', time_threshold=None):
        """ """
        if counts_type == 'observed frequency':
            y = self.inter_durations_histogram.observed_counts
            yticks = np.arange(0, max(y) * 1.3, 5)
        elif counts_type == 'probability density':
            y = self.inter_durations_histogram.normalized_counts
            yticks = np.arange(0, max(y) * 1.3, 1)
        else:
            raise ValueError("counts_type = 'observed frequency' or 'probability density'")
        ax.bar(self.inter_durations_histogram.midpoints, y, width=self.inter_durations_histogram.bin_widths, facecolor='steelblue', alpha=0.5, align='center')
        if time_threshold is not None:
            ax.axvline(time_threshold, ymin=0, ymax=yticks[-1], color='darkorange', linestyle='-', marker=None, label='$T_C = {}$ hours'.format(time_threshold))
        xticks = np.arange(0, max(self.inter_durations_histogram.edges) + 100, 100)
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::3]) # xticks[::2]
        # ax.set_xticks(self.inter_durations_histogram.midpoints, minor=True)
        # ax.set_xticks(self.inter_durations_histogram.midpoints[::20])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', linestyle=':', alpha=0.3)
        # xlim = [-100, 0.25 * self.inter_durations_histogram.edges[-1]]
        loc = np.where(self.inter_durations_histogram.normalized_counts >= 2.5 * 10**-4)[0]
        index = loc[-1]
        if len(loc) == 0:
            xlim = [-100, self.inter_durations_histogram.edges[-1]]
        else:
            xmax = self.round_up_to(self.inter_durations_histogram.edges[index], base=100)
            xlim = [-100, xmax]
        ax.set_xlim(xlim)
        ax.legend(loc='upper right', fontsize=textsize)
        # text_box = ax.text(0.05, 0.95, '$T_C = {}$ hours'.format(time_threshold), fontsize=textsize, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
        # text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlabel('Inter-Cluster Durations (hours)')
        ax.set_ylabel(counts_type.title())
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def initialize_relative_statistics(self):
        initial_cluster_sizes = np.array([cluster.size for cluster in self.clusters['elapsed hour']])
        cluster_sizes, size_counts = np.unique(initial_cluster_sizes, return_counts=True)
        nejecta = cluster_sizes * size_counts
        self.relative_statistics['cluster size'] = cluster_sizes
        self.relative_statistics['number of clusters'] = size_counts
        self.relative_statistics['number of ejecta'] = nejecta
        self.relative_statistics['probability'] = nejecta / np.sum(nejecta)
        # self.relative_statistics['mean duration'] = ...

    def initialize_duration_statistics(self):
        """ """
        self.mean_duration[None] = np.mean(self._intra_durations)
        self.median_duration[None] = np.median(self._intra_durations)
        self.stdev_duration[None] = np.std(self._intra_durations)
        self.max_duration[None] = np.max(self._intra_durations)
        self.min_duration[None] = np.min(self._intra_durations)
        for cluster_size in self.relative_statistics['cluster size']:
            Searcher = SearchClusters(self.clusters)
            Searcher.add_scale_parameters()
            searched_clusters = Searcher.search(search_parameters='cluster size', search_conditions='exact match', search_values=cluster_size)
            intra_durs = np.array([cluster[-1] - cluster[0] for cluster in searched_clusters['elapsed hour']])
            mcentral = np.mean(intra_durs)
            stdev = np.std(intra_durs)
            dur_label = '{0:.2f} ± {1:.2f}'.format(np.round(mcentral, decimals=2), np.round(stdev, decimals=2))
            self.mean_duration[cluster_size] = mcentral
            self.stdev_duration[cluster_size] = stdev
            self.median_duration[cluster_size] = np.median(intra_durs)
            self.max_duration[cluster_size] = np.max(intra_durs)
            self.min_duration[cluster_size] = np.min(intra_durs)
            self._duration_labels.append(dur_label)

    @property
    def duration_labels(self):
        return np.array(self._duration_labels)

    def subview_relative_cluster_frequency(self, ax, ticksize=7, labelsize=8, textsize=8, compare=False):
        """ """
        x = self.relative_statistics['cluster size']
        y = self.relative_statistics['number of clusters']
        xticks = np.arange(0, x[-1] +2)
        if compare == True:
            yticks = np.arange(0, 151, 5)
        else:
            yticks = np.arange(0, max(y)*1.3, 5)
        ax.bar(x, y, width=1, facecolor='darkorange', edgecolor='b')
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::2])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::5])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        # text_box = ax.text(0.95, 0.95, '${:,}$ Clusters'.format(self.clusters['elapsed hour'].size), fontsize=textsize, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        # text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlim([0, xticks[-1]])
        ax.set_xlabel('Cluster Size')
        ax.set_ylabel('Number of\nClusters')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_relative_ejecta_frequency(self, ax, ticksize=7, labelsize=8, textsize=8, compare=False):
        """ """
        x = self.relative_statistics['cluster size']
        y = self.relative_statistics['number of ejecta']
        xticks = np.arange(0, x[-1] +2)
        if compare == True:
            yticks = np.arange(0, 301, 5)
        else:
            yticks = np.arange(0, max(y)*1.3, 5)
        ax.bar(x, y, width=1, facecolor='darkgreen', edgecolor='r')
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::2])
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[::10])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        # nevents = np.concatenate(self.clusters['elapsed hour'], axis=0).size
        # text_box = ax.text(0.95, 0.95, '${:,}$ Extreme Events'.format(nevents), fontsize=textsize, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        # text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlim([0, xticks[-1]])
        ax.set_xlabel('Cluster Size')
        ax.set_ylabel('Number of\nExtreme Events')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_relative_probability(self, ax, ticksize=7, labelsize=8, textsize=8, compare=False):
        """ """
        x = self.relative_statistics['cluster size']
        y = self.relative_statistics['probability']
        xticks = np.arange(0, x[-1] +2)
        if compare == True:
            yticks = np.arange(0, 1.01, 0.1)
        else:
            yticks = np.arange(0, max(y) + 0.11, 0.1)
        ax.bar(x, y, width=1, facecolor='r', edgecolor='k')
        ax.set_xticks(xticks, minor=True)
        ax.set_xticks(xticks[::2])
        ax.set_yticks(yticks[1::2], minor=True)
        ax.set_yticks(yticks[::2])
        ax.grid(color='k', alpha=0.3, linestyle=':')
        ax.set_xlim([0, xticks[-1]])
        ax.set_xlabel('Cluster Size')
        ax.set_ylabel('Relative\nProbability')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax

    def subview_chronological_cluster_sizes(self, ax, ticksize=7, labelsize=8, textsize=8, ydelta=0):
        """ """
        alternating_facecolors = ('r', 'b')
        x, y, w, c = [], [], [], []
        f = lambda dt : date2num(datetime.datetime.strptime(dt.split(' ')[-1], '%H:%M:%S'))
        lone_width = f('1999/10/27 01:00:00') - f('1999/10/27 00:00:00')
        for ith_cluster, elapsed_cluster in enumerate(self.clusters['elapsed hour']):
            dt_cluster = self.clusters['datetime object'][ith_cluster]
            dt_initial, dt_final = date2num(dt_cluster[0]), date2num(dt_cluster[-1])
            facecolor = np.roll(alternating_facecolors, ith_cluster)[0]
            c.append(facecolor)
            y.append(elapsed_cluster.size)
            if elapsed_cluster.size == 1:
                x.append(dt_initial)
                w.append(lonewidth)
            else:
                tmp = [dt_initial, dt_final]
                x.append(np.mean(tmp))
                w.append(np.diff(tmp)[0])
        nclusters = self.clusters['elapsed hour'].size
        nevents = np.sum(y)
        yticks = np.arange(max(y) + 3 + ydelta)
        ax.bar(x, y, width=w, color=c, align='center')
        ax.set_yticks(yticks, minor=True)
        ax.set_yticks(yticks[1::2])
        ax = self.transform_x_as_datetime(ax)
        ax.grid(color='k', alpha=0.3, linestyle=':')
        # ax.legend(loc='upper left', fontsize=textsize)
        text_box = ax.text(0.95, 0.95, '${:,}$ Extreme Events\n via ${:,}$ Clusters'.format(nevents, nclusters), fontsize=textsize, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        text_box.set_bbox(dict(facecolor='gray', alpha=0.25, edgecolor='k'))
        ax.set_xlabel('Date')
        ax.set_ylabel('Cluster Size')
        ax = self.set_tick_size(ax, ticksize, ticksize)
        ax = self.set_label_size(ax, labelsize, labelsize)
        return ax
