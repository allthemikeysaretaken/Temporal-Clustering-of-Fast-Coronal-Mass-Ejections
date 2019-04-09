import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.legend_handler
from matplotlib.ticker import ScalarFormatter, LogFormatter, NullFormatter
import datetime
from matplotlib.colors import Normalize
from matplotlib.legend import Legend


class TicksLabelsMethods():

    @staticmethod
    def round_to_nearest(number, direction, nearest=10):
        """ """
        if direction == 'up':
            return int(np.ceil(number + nearest) / nearest) * nearest
        elif direction == 'down':
            return int(np.floor(number) / nearest) * nearest
        else:
            raise ValueError("direction = 'up' or 'down'")

    @staticmethod
    def get_every_other(ticks):
        """ """
        return ticks[::2], ticks[1::2]

    @staticmethod
    def get_every_n(ticks, n, shift=0):
        """ """
        majors, minors = [], []
        for idx in range(len(ticks)):
            if (idx + shift) % n == 0:
                majors.append(ticks[idx])
            else:
                minors.append(ticks[idx])
        return np.array(majors).astype(int), np.array(minors).astype(int)

    @property
    def step_mapping(self):
        """ """
        res = {}
        res['every other'] = {'f' : self.get_every_other, 'kwargs' : dict()}
        res['every three'] = {'f' : self.get_every_n, 'kwargs' : dict(n=3, shift=0)}
        res['every ten'] = {'f' : self.get_every_n, 'kwargs' : dict(n=10, shift=0)}
        return res

    @staticmethod
    def configure_xaxis(ax, mirror_ticks=False, mirror_labels=False, erase_ticks=False, erase_labels=False):
        """ """
        if mirror_ticks is True:
            ax.xaxis.tick_top()
        if mirror_labels is True:
            ax.xaxis.set_label_position('top')
        if erase_ticks is True:
            ax.set_xticklabels([])
        if erase_labels is True:
            ax.set_xlabel('')

    @staticmethod
    def configure_yaxis(ax, mirror_ticks=False, mirror_labels=False, erase_ticks=False, erase_labels=False):
        """ """
        if mirror_ticks is True:
            ax.yaxis.tick_right()
        if mirror_labels is True:
            ax.yaxis.set_label_position('right')
        if erase_ticks is True:
            ax.set_yticklabels([])
        if erase_labels is True:
            ax.set_ylabel('')

    def standard_configuration(self, axes):
        """ """
        iter_axes = axes.ravel()
        n = len(iter_axes)
        ndim = len(axes.shape)
        if n == 1:
            xerase_axes, xmirror_axes, yerase_axes, ymirror_axes = None, None, None, None
        elif n > 1:
            if ndim == 1:
                xmirror_axes = None
                xerase_axes = None
                ymirror_axes = np.array([axes[-1]])
                if n > 2:
                    yerase_axes = np.array([iter_axes[idx] for idx in range(n) if idx not in (0, n-1)])
                else:
                    yerase_axes = None
            else:
                xmirror_axes = None
                xerase_axes = axes[:-1, :]
                ymirror_axes = axes[:, -1]
                if axes.shape[-1] > 2:
                    yerase_axes = axes[:, 1:-1]
                else:
                    yerase_axes = None
        return xerase_axes, xmirror_axes, yerase_axes, ymirror_axes

    def edge_column_configuration(self, axes):
        """ """
        nrows, ncols = axes.shape
        if nrows > 1:
            if nrows > 2:
                xerase_axes = axes[1:-1, :]
                xmirror_axes = axes[0, :]
            else:
                xerase_axes = None
                xmirror_axes = None
        if ncols > 2:
            ymirror_axes = axes[:, ncols//2:]
            # yerase_axes = axes[:, 1:-1]
            yerase_axes = None
        else:
            ymirror_axes = axes[:, -1]
            yerase_axes = None
        return xerase_axes, xmirror_axes, yerase_axes, ymirror_axes


    @property
    def configuration_mapping(self):
        """ """
        res = {}
        res['standard'] = dict(f=self.standard_configuration, kwargs={})
        res['edge column'] = dict(f=self.edge_column_configuration, kwargs={})
        # res['linear'] = dict(f=self.linear_configuration, kwargs={})
        # res['logarithmic'] = dict(f=self.logarithmic_configuration, kwargs={})
        return res

    def autoconfigure_axes(self, axes, key='standard'):
        """ """
        skip = False
        if not isinstance(axes, np.ndarray):
            axes = np.array(axes)
        config_map = self.configuration_mapping[key]
        f, kwargs = config_map['f'], config_map['kwargs']
        xerase_axes, xmirror_axes, yerase_axes, ymirror_axes = f(axes, **kwargs)
        if xmirror_axes is not None:
            for ax in xmirror_axes.ravel():
                TicksLabelsMethods.configure_xaxis(ax, mirror_ticks=True, mirror_labels=True)
        if ymirror_axes is not None:
            for ax in ymirror_axes.ravel():
                TicksLabelsMethods.configure_yaxis(ax, mirror_ticks=True, mirror_labels=True)
        if xerase_axes is not None:
            for ax in xerase_axes.ravel():
                TicksLabelsMethods.configure_xaxis(ax, erase_ticks=True, erase_labels=True)
        if yerase_axes is not None:
            for ax in yerase_axes.ravel():
                TicksLabelsMethods.configure_yaxis(ax, erase_ticks=True, erase_labels=True)

    @staticmethod
    def make_doublesided(axes_dom, axes_sub, nrows, ncols=None):
        """ """
        if nrows > 1:
            axes_lin = np.array(axes_dom).reshape((nrows, ncols))
            axes_log = np.array(axes_sub).reshape((nrows, ncols))
            for ax in axes_lin[:, 1:].ravel():
                TicksLabelsMethods.configure_yaxis(ax, erase_ticks=True, erase_labels=True)
            for ax in axes_log[:, :-1].ravel():
                TicksLabelsMethods.configure_yaxis(ax, erase_ticks=True, erase_labels=True)
            for ax in axes_lin[:-1, :].ravel():
                TicksLabelsMethods.configure_xaxis(ax, erase_ticks=True, erase_labels=True)
            for ax in axes_log[1:, :].ravel():
                TicksLabelsMethods.configure_xaxis(ax, erase_ticks=True, erase_labels=True)
        else:
            if not isinstance(axes_dom, np.ndarray):
                axes_dom = np.array(axes_dom)
            if not isinstance(axes_sub, np.ndarray):
                axes_sub = np.array(axes_sub)
            for ax in axes_dom[1:].ravel():
                TicksLabelsMethods.configure_yaxis(ax, erase_ticks=True, erase_labels=True)
            for ax in axes_sub[:-1].ravel():
                TicksLabelsMethods.configure_yaxis(ax, erase_ticks=True, erase_labels=True)
        return axes_dom, axes_sub

    @staticmethod
    def get_intra_time_ticks(H):
        """ """
        if H.edges[-1] > 20:
            if H.edges[-1] >= 50:
                xticks = H.edges[::5]
            else:
                xticks = H.edges[::2]
        else:
            xticks = H.edges
        if max(H.observed_counts) > 40:
            ynearest = 20
            ydelta = 10
        else:
            ynearest = 10
            ydelta = 5
        ymax = TicksLabelsMethods().round_to_nearest(max(H.observed_counts), direction='up', nearest=ynearest)
        if ydelta is None:
            yticks = np.unique(np.linspace(0, ymax*1.5, 20).astype(int))
        else:
            yticks = np.arange(0, ymax + ydelta - 1, ydelta//2, dtype=int)
        xlim = (xticks[0], xticks[-1])
        ylim = (yticks[0], yticks[-1])
        return xticks, yticks, xlim, ylim

    @staticmethod
    def get_intra_duration_ticks(H):
        """ """
        if H.edges[-1] >= 80:
            xticks = H.edges[::20]
        elif H.edges[-1] >= 50:
            xticks = H.edges[::5]
        else:
            xticks = H.edges[::2]
        yprime = max(H.observed_counts)
        if yprime >= 100:
            ynearest = 60
            ydelta = 30
        elif yprime >= 50:
            ynearest = 40
            ydelta = 20
        elif yprime > 40:
            ynearest = 20
            ydelta = 10
        else:
            ynearest = 10
            ydelta = 5
        ymax = TicksLabelsMethods().round_to_nearest(yprime, direction='up', nearest=ynearest)
        yticks = np.arange(0, ymax + ydelta - 1, ydelta//2, dtype=int)
        xlim = (xticks[0], xticks[-1] + xticks[1] - xticks[0])
        ylim = (yticks[0], yticks[-1])
        return xticks, yticks, xlim, ylim

    @staticmethod
    def get_inter_duration_ticks(H):
        """ """
        loc = np.where(H.normalized_counts >= 2.5e-4)[0]
        if len(loc) == 0:
            loc = np.where(H.normalized_counts > 0)[0]
        xmax, xmin = H.edges[max(loc)], -100
        xnearest, ynearest = 10, 50
        xdelta = TicksLabelsMethods().round_to_nearest(int((xmax - xmin) / xnearest), direction='up', nearest=xnearest)
        xticks = np.arange(xmin, xmax, int(xdelta), dtype=int)
        ydelta = 20
        ymax = TicksLabelsMethods().round_to_nearest(max(H.observed_counts), direction='up', nearest=ynearest)
        if ydelta is None:
            yticks = np.unique(np.linspace(0, ymax*1.5, 20).astype(int))
        else:
            yticks = np.arange(0, ymax + ydelta - 1, ydelta//2, dtype=int)
        xlim = (xticks[0], xticks[-1] + xticks[1] - xticks[0])
        ylim = (yticks[0], yticks[-1])
        return xticks, yticks, xlim, ylim

class RowColumnOperations():

    def __init__(self, MultipleTimeSeries):
        """ """
        if isinstance(MultipleTimeSeries, (tuple, list, np.ndarray)):
            self.MultipleTimeSeries = MultipleTimeSeries
        else:
            self.MultipleTimeSeries = [MultipleTimeSeries]
        self.n = len(self.MultipleTimeSeries)
        self.result = {}

    def get_rows_columns(self, nrows=None, ncols=None, swap_rc=False, apply_overlay=False):
        """ """
        if apply_overlay is True:
            nrows, ncols = 1, 1
        else:
            condition = ((nrows is None) or (ncols is None))
            if condition is False:
                ntotal = nrows * ncols
                if ((ntotal > self.n) or (ntotal < 1)):
                    raise ValueError("cannot figure {} rows by {} columns for {} plots".format(nrows, ncols, self.n))
            if self.n == 1:
                if condition is True:
                    nrows, ncols = 1, 1
            elif self.n == 2:
                if condition is True:
                    nrows, ncols = 1, 2
            elif self.n == 3:
                if condition is True:
                    nrows, ncols = 3, 1
            elif self.n >= 4:
                if condition is True:
                    nrows = self.n // 2
                    ncols = int(self.n / nrows)
            if swap_rc is True:
                nrows, ncols = ncols, nrows
        return nrows, ncols

    @property
    def key_mapping(self):
        """ """
        res = {}
        res['inter-exceedance distribution decay'] = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        res['alpha-hat histogram'] = dict(nrows=1, ncols=1, swap_rc=False, apply_overlay=True)
        res['max spectrum power-law'] = None
        res['unbiased estimator errors'] = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        res['point estimators'] = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        res['extremal index histogram'] = dict(nrows=1, ncols=1, swap_rc=False, apply_overlay=True)
        res['moment estimates'] = None
        res['speed tail'] = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        res['speed distribution'] = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        res['sunspot correlations'] = None
        res['chronological cluster size'] = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        res['relative cluster statistics'] = None
        res['durational histograms'] = None
        res['dim3 speed distribution'] = None
        res['frechet distribution'] = None
        res['illustrative cluster example'] = dict(nrows=1, ncols=1, swap_rc=False, apply_overlay=True)
        return res

    def initialize(self, key, sharex=False, sharey=False, **kwargs):
        """ """
        if key is None:
            return None, None
        else:
            init_kwargs = self.key_mapping[key]
            if init_kwargs is None:
                fig = plt.figure(**kwargs)
                axes = None
            else:
                nrows, ncols = self.get_rows_columns(**init_kwargs)
                fig, axes = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey, **kwargs)
                if isinstance(axes, (tuple, list)):
                    axes = np.array(axes)
                elif not isinstance(axes, np.ndarray):
                    axes = np.array([axes])
            return fig, axes

class AxesMethods(RowColumnOperations):

    def __init__(self, MultipleTimeSeries, ticksize, labelsize, titlesize, saveloc, extension):
        """ """
        super().__init__(MultipleTimeSeries)
        self.ticksize = ticksize
        self.labelsize = labelsize
        self.titlesize = titlesize
        self.saveloc = saveloc
        self.extension = extension

    @staticmethod
    def get_partial_title(TimeSeries, include_speed_type=False):
        """ """
        partial_title = r'SC ${}$'.format(TimeSeries.identity['solar cycle'])
        if include_speed_type is True:
            if TimeSeries.identity['speed type'] == 'linear speed':
                symbol = r'$V_{linear}$'
            elif TimeSeries.identity['speed type'] == 'second order initial speed':
                symbol = r'$V_{20 R_{\odot}, i}$'
            elif TimeSeries.identity['speed type'] == 'second order final speed':
                symbol = r'$V_{20 R_{\odot}, f}$'
            elif TimeSeries.identity['speed type'] == 'second order 20R speed':
                symbol = r'$V_{20 R_{\odot}}$'
            else:
                raise ValueError("unknown speed type: {}".format(TimeSeries.identity['speed type']))
            partial_title += ' ({})'.format(symbol)
        return partial_title

    @staticmethod
    def get_rectangular_handle(fill=False, edgecolor='none', visible=False, xyi=(0, 0), xyf=(0.1, 0.1), facecolor='none'):
        """ """
        return matplotlib.patches.Rectangle(xyi, *xyf, fill=fill, edgecolor=edgecolor, visible=visible, facecolor=facecolor)

    @staticmethod
    def get_axis_ticks(steps, ticks):
        """ """
        step_mapping = TicksLabelsMethods().step_mapping[steps]
        f, kwargs = step_mapping['f'], step_mapping['kwargs']
        res = f(ticks, **kwargs)
        return res

    def develop_x_axis(self, ax, ticks=None, label=None, steps='every other', bounds=None):
        """ """
        if ticks is not None:
            if steps is None:
                ax.set_xticks(ticks)
                ax.set_xticklabels(ticks, fontsize=self.ticksize)
            else:
                xticks = self.get_axis_ticks(steps, ticks)
                ax.set_xticks(xticks[0])
                ax.set_xticks(xticks[1], minor=True)
        ax.tick_params(axis='x', labelsize=self.ticksize)
        if label is not None:
            ax.set_xlabel(label, fontsize=self.labelsize)
        if bounds is not None:
            ax.set_xlim(bounds)
        return ax

    def develop_y_axis(self, ax, ticks=None, label=None, steps='every other', bounds=None):
        """ """
        if ticks is not None:
            if steps is None:
                ax.set_yticks(ticks)
                ax.set_yticklabels(ticks, fontsize=self.ticksize)
            else:
                yticks = self.get_axis_ticks(steps, ticks)
                ax.set_yticks(yticks[0])
                ax.set_yticks(yticks[1], minor=True)
        ax.tick_params(axis='y', labelsize=self.ticksize)
        if label is not None:
            ax.set_ylabel(label, fontsize=self.labelsize)
        if bounds is not None:
            ax.set_ylim(bounds)
        return ax

    def autoposition_standard_legend(self, fig, axes, loc='auto', ncol='auto', **kwargs):
        """ """
        fig.tight_layout()
        if self.n == 1:
            if loc == 'auto':
                loc = 'upper right'
            if ncol == 'auto':
                ncol = 1
            try:
                leg = axes[0].legend(loc=loc, ncol=ncol, fontsize=self.labelsize)
            except:
                leg = axes[0, 0].legend(loc=loc, ncol=ncol, fontsize=self.labelsize)
        else:
            try:
                handles, labels = axes[0].get_legend_handles_labels()
            except:
                handles, labels = axes[0, 0].get_legend_handles_labels()
            if ncol == 'auto':
                ncol = len(labels)
            if loc == 'auto':
                fig.subplots_adjust(bottom=0.125)
                if len(labels) > 1:
                    leg = fig.legend(handles, labels, loc='lower center', ncol=ncol, fontsize=self.labelsize, mode='expand', fancybox=True, borderaxespad=0.2, **kwargs)
                else:
                    leg = fig.legend(handles, labels, loc='lower center', ncol=ncol, fontsize=self.labelsize, fancybox=True, borderaxespad=0.2, **kwargs)
            else:
                leg = fig.legend(handles, labels, loc=loc, ncol=ncol, fontsize=self.labelsize, fancybox=True, **kwargs)
        return leg

    def autoposition_overlay_legend(self, ax, handle, label, idx, facecolor, borderaxespad=0.2, frameon=True, shift=0.4, loc='upper right', include_speed_type=False):
        """ """
        if loc == 'upper right':
            bboxs = [(0.98, 0.98 - (i * 0.225)) for i in range(self.n)]
        elif loc == 'upper left':
            bboxs = [(0.02, 0.98 - (i * 0.225)) for i in range(self.n)]
        else:
            raise ValueError("not yet implemented")
        if isinstance(loc, str):
            locations = [loc for i in range(len(bboxs))]
        else:
            locations = loc
        bbox_to_anchor = bboxs[idx]
        loc = locations[idx]
        leg = Legend(ax, [handle], [label], loc=locations[idx], facecolor=facecolor, fontsize=self.labelsize, bbox_to_anchor=bbox_to_anchor, borderaxespad=borderaxespad, frameon=frameon)
        partial_title = self.get_partial_title(self.MultipleTimeSeries[idx], include_speed_type)
        leg.set_title(partial_title, prop={'size' : self.labelsize})
        for t in leg.get_texts():
            t.set_ha('left')
            t.set_position((shift, 0))
        ax.add_artist(leg)
        return leg

    def apply_title(self, fig, ax, header=None, partial_title=None, apply_overlay=False):
        """ """
        if self.n == 1:
            if header is not None:
                fig.suptitle(header, fontsize=self.titlesize, y=1.05)
        else:
            if apply_overlay is True:
                ax.set_title(header, fontsize=self.titlesize)
            else:
                if partial_title is not None:
                    ax.set_title(partial_title, fontsize=self.labelsize)

    def transform_as_subliminal(self, ax, xcolor='gray', ycolor='gray'):
        """ """
        ax.tick_params(axis='x', colors=xcolor, labelsize=self.ticksize)
        ax.tick_params(axis='y', colors=ycolor, labelsize=self.ticksize)
        ax.xaxis.label.set_color(xcolor)
        ax.yaxis.label.set_color(ycolor)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        return ax

    def transform_as_datetime(self, ax, axis='x'):
        """ """
        if axis == 'x':
            ax.xaxis.set_major_locator(mdates.YearLocator())
            ax.xaxis.set_minor_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            ax.tick_params(axis='x', labelsize=self.ticksize, rotation=15)
            ax.tick_params(axis='y', labelsize=self.ticksize)
        elif axis == 'y':
            ax.yaxis.set_major_locator(mdates.YearLocator())
            ax.yaxis.set_minor_locator(mdates.MonthLocator())
            ax.yaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            ax.tick_params(axis='x', labelsize=self.ticksize)
            ax.tick_params(axis='y', labelsize=self.ticksize, rotation=15)
        # ax.tick_params(axis=axis, labelsize=self.ticksize, rotation=15)
        return ax

    @staticmethod
    def select_specific_axis(axes, idx, apply_overlay):
        """ """
        if apply_overlay is True:
            try:
                ax = axes[0, 0]
            except:
                try:
                    ax = axes[0]
                except:
                    ax = axes
        else:
            ax = axes.ravel()[idx]
        return ax

    def subview_interexceedance_distribution_against_exponential_decay(self, fig, axes, facecolors=('r', 'b'), include_speed_type=False):
        """ """
        header = 'Comparison of the Distribution of Inter-Exceedance Times\n& the Corresponding Exponential Distribution'
        for idx, ax in enumerate(axes.ravel()):
            TimeSeries = self.MultipleTimeSeries[idx]
            htypes = ('inverse sample', 'time-delta')
            labels = ('Exponential Distribution', 'Distribution of Inter-Exceedance Times')
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            partial_title += r': ${}$ CMEs'.format(TimeSeries.InterExceedanceDecay.Histogram('time-delta').size)
            ys = []
            for htype, facecolor, label in zip(htypes, facecolors, labels):
                H = TimeSeries.InterExceedanceDecay.Histogram(htype)
                ys.append(H.observed_counts)
                ax.bar(H.midpoints, H.observed_counts, width=H.bin_widths, label=label, facecolor=facecolor, edgecolor='none', linewidth=0, alpha=0.5)
            ymax = TicksLabelsMethods.round_to_nearest(np.max(np.array(ys).reshape(-1)), 'up', 50)
            yticks = np.arange(0, ymax+9, 25, dtype=int)
            ax = self.develop_x_axis(ax, ticks=H.edges, label=r'Inter-Exceedance Times (hours)', steps='every other', bounds=[0, H.edges[-1]])
            ax = self.develop_y_axis(ax, ticks=yticks, label=r'Observed Frequency', steps='every other', bounds=None)
            ax.grid(color='k', alpha=0.3, linestyle=':')
            ax.set_title(partial_title, fontsize=self.labelsize)
        leg = self.autoposition_standard_legend(fig, axes, loc='auto', ncol='auto')
        fig.suptitle(header, fontsize=self.titlesize, y=1.05)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def subview_alpha_hat_histogram(self, fig, axes, density='observed', facecolor=('steelblue', 'darkorange', 'green', 'red'), apply_overlay=True, include_speed_type=False):
        """ """
        available_density = ('probability', 'observed')
        if density not in available_density:
            raise ValueError("unknown density: {}; available density: {}".format(density, available_density))
        if density == 'probability':
            f = lambda cls : cls.normalized_counts
            ylabel = r'Probability Density'
        elif density == 'observed':
            f = lambda cls : cls.observed_counts
            ylabel = r'Observed Frequency'
        else:
            raise ValueError("density = 'probability' or 'observed'")
        header = r'Histogram of $\hat\alpha$' + r' via ${}$ resamples'.format(self.MultipleTimeSeries[0].UnbiasedEstimators.nresamples)
        ys = []
        for idx, TimeSeries in enumerate(self.MultipleTimeSeries):
            xticks = np.round(TimeSeries.UnbiasedEstimators.Histogram.edges, decimals=2)[::3]
            y = f(TimeSeries.UnbiasedEstimators.Histogram)
            ymax = TicksLabelsMethods.round_to_nearest(max(y), 'up', 10)
            ys.append(ymax)
            yticks = np.arange(0, max(ys)+9, 10, dtype=int)
            alpha_skew = TimeSeries.UnbiasedEstimators.alpha('skew', include_original=False)
            alpha_kurt = TimeSeries.UnbiasedEstimators.alpha('kurtosis', include_original=False)
            skew_label = r'skew($\alpha$) $ = {:.3f}$'.format(alpha_skew)
            kurt_label = r'kurtosis($\alpha$) $ = {:.3f}$'.format(alpha_kurt)
            alpha_label = r'$\hat\alpha = {0:.3}$'.format(TimeSeries.UnbiasedEstimators.alpha('hat'))
            alpha_hat = np.round(TimeSeries.UnbiasedEstimators.alpha('hat'), decimals=2)
            label = '{}\n{}\n{}'.format(alpha_label, skew_label, kurt_label)
            handle = self.get_rectangular_handle(fill=False, edgecolor='none', visible=False, xyi=(0, 0), xyf=(0.1, 0.1))
            ax = self.select_specific_axis(axes, idx, apply_overlay)
            ax.bar(TimeSeries.UnbiasedEstimators.Histogram.midpoints, y, width=TimeSeries.UnbiasedEstimators.Histogram.bin_widths, label=None, facecolor=facecolor[idx], edgecolor='none', linewidth=0, alpha=1/self.n)
            ax = self.develop_x_axis(ax, ticks=xticks, label=r'$\hat\alpha$', steps='every other', bounds=[TimeSeries.UnbiasedEstimators.Histogram.edges[0], TimeSeries.UnbiasedEstimators.Histogram.edges[-1]])
            ax = self.develop_y_axis(ax, ticks=yticks, label=ylabel, steps='every other', bounds=None)
            ax.grid(color='k', alpha=0.3, linestyle=':')
            if apply_overlay is True:
                self.autoposition_overlay_legend(ax, handle, label, idx, facecolor[idx], loc='upper right', include_speed_type=include_speed_type)
            else:
                alpha_hat = np.round(TimeSeries.UnbiasedEstimators.alpha('hat'), decimals=2)
                string_alpha = r'$\hat\alpha = {0:.2}$'.format(alpha_hat)
                partial_title = self.get_partial_title(TimeSeries, include_speed_type)
                partial_title += ': {}'.format(string_alpha)
                ax.set_title(partial_title, fontsize=self.labelsize)
                ax.legend(handles=[handle], labels=[label], loc='upper right', fontsize=self.labelsize)
        self.fig.tight_layout()
        fig.suptitle(header, fontsize=self.titlesize, y=1.05)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def subview_unbiased_estimator_errors(self, fig, axes, facecolors=('r', 'b'), ith_resample=0, include_speed_type=False):
        """ """
        header = r'Unbiased Estimator Errors'
        error_types = ('standard deviation', 'standard error')
        symbols = (r'Standard Deviation $\sigma(j)$', r'Standard Error $\tilde\sigma_{\mu}(j)$')
        yticks = np.arange(7, 14, 1, dtype=int)
        for idx, ax in enumerate(axes.ravel()):
            TimeSeries = self.MultipleTimeSeries[idx]
            MS = TimeSeries.UnbiasedEstimators.select_ith_resample(ith_resample)
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            for err, facecolor, sym in zip(error_types, facecolors, symbols):
                ax.errorbar(MS.js, MS.Y(j=None, fit=False), yerr=MS.result[err], alpha=0.7, capsize=5, label=sym, ecolor=facecolor, fmt='none')
            ax = self.develop_x_axis(ax, ticks=MS.js, label=r'$j$ ($log_2$ hours)', steps='every other', bounds=[MS.js[0]-0.5, MS.js[-1]+0.5])
            ax = self.develop_y_axis(ax, ticks=yticks, label=r'$Y(j)$ ($log_2$ $\frac{km}{s}$)', steps='every other', bounds=[yticks[0], yticks[-1]])
            ax.set_title(partial_title, fontsize=self.labelsize)
            ax.grid(color='k', alpha=0.3, linestyle=':', which='both')
        leg = self.autoposition_standard_legend(fig, axes, loc='auto', ncol='auto')
        fig.suptitle(header, fontsize=self.titlesize, y=1.05)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def subview_power_law(self, fig, axes=None, basex=2, basey=2, ith_resample=0, include_speed_type=False):
        """ """
        if axes is None:
            axes = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        elif not isinstance(axes, dict):
            raise ValueError("axes = None or type <dict>")
        xlim = np.array([0.1, 17])
        ylim = np.array([8, 12])
        yticks = np.arange(ylim[0], ylim[-1]+0.1, 0.5)
        nrows, ncols = self.get_rows_columns(**axes)
        header = r'Power-Law via Weighted Fit of Max Spectrum'
        axes_lin, axes_log = [], []
        for idx in range(self.n):
            TimeSeries = self.MultipleTimeSeries[idx]
            ax_lin = fig.add_subplot(nrows, ncols, idx+1)
            ax_log = fig.add_subplot(nrows, ncols, idx+1, frame_on=False)
            ax_log.set_xscale('log', basex=basex)
            ax_log.set_yscale('log', basey=basey)
            ax_log.xaxis.set_major_formatter(ScalarFormatter())
            ax_log.yaxis.set_major_formatter(ScalarFormatter())
            MS = TimeSeries.UnbiasedEstimators.select_ith_resample(i=ith_resample)
            ax_lin.scatter(MS.js, MS.yj_init, color='gray', marker='.', label='initial Max Spectrum')
            ax_lin.plot(MS.js_fit, MS.yj_fit, color='r', label='Power-Law')
            ax_log.scatter(basex**MS.js, basey**MS.yj_init, color='gray', marker='.')
            ax_log.plot(basex**MS.js_fit, basey**MS.yj_fit, color='r')
            ax_lin = self.develop_x_axis(ax_lin, ticks=MS.js, label=r'$j$ ($log_2$ hours)', steps='every other', bounds=xlim)
            ax_lin = self.develop_y_axis(ax_lin, ticks=yticks, label=r'$Y(j)$ ($log_2$ $\frac{km}{s}$)', steps='every other', bounds=ylim)
            ax_log = self.develop_x_axis(ax_log, ticks=basex**MS.js, label=r'$2^j$ (hours)', steps='every three', bounds=basex**xlim)
            ax_log = self.develop_y_axis(ax_log, ticks=basey**yticks, label=r'$2^{Y(j)}$ ($\frac{km}{s}$)', steps='every three', bounds=basey**ylim)
            ax_log = self.transform_as_subliminal(ax_log, xcolor='gray', ycolor='gray')
            TicksLabelsMethods.configure_xaxis(ax_log, mirror_ticks=True, mirror_labels=True)
            TicksLabelsMethods.configure_yaxis(ax_log, mirror_ticks=True, mirror_labels=True)
            ax_log.grid(color='k', alpha=0.3, linestyle=':')
            labels = self.get_partial_title(self.MultipleTimeSeries[idx], include_speed_type)
            ax_lin.legend(loc='upper left', handles=[self.get_rectangular_handle()], labels=[labels], fontsize=self.labelsize)
            axes_lin.append(ax_lin)
            axes_log.append(ax_log)
        axes_lin, axes_log = TicksLabelsMethods.make_doublesided(axes_lin, axes_log, nrows, ncols)
        if self.n > 1:
            leg = self.autoposition_standard_legend(fig, axes_lin, loc='auto', ncol='auto')
        else:
            leg = self.autoposition_standard_legend(fig, axes_lin, loc='upper left', ncol='auto')
        fig.tight_layout()
        fig.suptitle(header, y=1.05, fontsize=self.titlesize)

    def subview_speed_tail(self, fig, axes, vthreshold=800, facecolor='darkorange', include_speed_type=False):
        """ """
        htype = 'preliminary'
        yt = [1, 10, 100, 1000]
        yticks = [0.001] + yt
        yticklabels = [0] + yt
        xticks = np.unique([vthreshold, 1000, 2000, 3000, 4000])
        ylabel = 'Observed\nFrequency'
        header = r'Comparison of the Tails' + '\n' + r'of the Speed Distributions'
        for idx, ax in enumerate(axes.ravel()):
            TimeSeries = self.MultipleTimeSeries[idx]
            SD = TimeSeries.SpeedDistribution
            H = SD.histograms[htype]
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            ax.bar(H.midpoints, H.observed_counts, width=H.bin_widths, label='CME Speeds', facecolor=facecolor, edgecolor='none', linewidth=0)
            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax = self.develop_x_axis(ax, ticks=xticks, label=r'CME Speed ($\frac{km}{s}$)', steps=None, bounds=[vthreshold, 4000])
            ax = self.develop_y_axis(ax, ticks=yticks, label='Observed\nFrequency', steps=None, bounds=[0.1, 10**3])
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.grid(color='k', alpha=0.3, linestyle=':')
            ax.set_title(partial_title, fontsize=self.labelsize)
        leg = self.autoposition_standard_legend(fig, axes, loc='auto', ncol='auto')
        fig.tight_layout()
        fig.suptitle(header, y=1.05, fontsize=self.titlesize)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def fill_speed_confidence_interval(self, ax, SD, err, y, facecolor):
        """ """
        parameters = SD.result[err]['parameters']
        (mu, sigma) = SD.from_normal_to_lognormal(parameters)
        median = np.sqrt((np.exp(parameters[1]**2 -1)) * (np.exp(2*parameters[0] + parameters[1]**2)))
        mode = np.exp(parameters[0] - parameters[1]**2)
        args = (mu, sigma, median, mode)
        ax.fill_between(SD.x, 0, y, where=((SD.x > (mu - sigma)) & (SD.x < (mu + sigma))), color=facecolor, alpha=0.1, label=r'$\mu \pm \sigma$')
        ax.fill_between(SD.x, 0, y, where=((SD.x > (mu - 2*sigma)) & (SD.x < (mu - sigma))), color=facecolor, alpha=0.3)
        ax.fill_between(SD.x, 0, y, where=((SD.x > (mu + sigma)) & (SD.x < (mu + 2*sigma))), color=facecolor, alpha=0.3, label=r'$\mu \pm 2 \sigma$')
        ax.fill_between(SD.x, 0, y, where=((SD.x > (mu - 3*sigma)) & (SD.x < (mu - 2*sigma))), color=facecolor, alpha=0.6)
        ax.fill_between(SD.x, 0, y, where=((SD.x > (mu + 2*sigma)) & (SD.x < (mu + 3*sigma))), color=facecolor, alpha=0.6, label=r'$\mu \pm 3 \sigma$')
        ax.fill_between(SD.x, 0, y, where=(SD.x > (mu + 3*sigma)), color=facecolor, alpha=1)
        return ax, args

    def subview_speed_distribution(self, fig, axes, err='minimum gtest estimation', fit_color='darkorange', con_color='b', hst_color='k', include_speed_type=False):
        """ """
        htype = 'preliminary'
        normalized = True
        xdelta = 250
        yticks = np.arange(0, 0.0036, 0.00025)
        header = r'Comparison of the Speed Distributions'
        for idx, ax in enumerate(axes.ravel()):
            TimeSeries = self.MultipleTimeSeries[idx]
            SD = TimeSeries.SpeedDistribution
            H = SD.histograms[htype]
            y = SD.get_y(err, normalized, htype)
            xticks = np.arange(0, H.edges[-1] - xdelta + 1, xdelta)
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            partial_title += r': ${}$ CMEs'.format(SD.x.size)
            ax.plot(SD.x, y, label=r'Optimized Fit', color=fit_color)
            if con_color is None:
                ax.bar(H.midpoints, H.normalized_counts, width=H.bin_widths, label='Normalized Histogram of CME Speeds', facecolor=hst_color, edgecolor='none', linewidth=0)
            else:
                ax.step(H.midpoints, H.normalized_counts, label='Normalized Histogram of CME Speeds', color=hst_color, linewidth=0.5, where='mid')
                ax, args = self.fill_speed_confidence_interval(ax, SD, err, y, facecolor=con_color)
                labels = ['Mean', 'Standard Deviation', 'Median', 'Mode']
                labels = ['{}: ${:.0f}$ '.format(label, np.round(arg, decimals=0)) + r'$\frac{km}{s}$' for label, arg in zip(labels, args)]
                handles = [self.get_rectangular_handle(fill=False, edgecolor='none', visible=False, xyi=(0, 0), xyf=(0.1, 0.1))]*len(labels)
                ax.legend(handles=handles, labels=labels, loc='upper right', fontsize=self.ticksize)
            ax.set_title(partial_title, fontsize=self.labelsize)
            ax = self.develop_x_axis(ax, ticks=xticks, label=r'CME Speed ($\frac{km}{s}$)', steps='every other', bounds=[0, H.edges[-1]])
            ax = self.develop_y_axis(ax, ticks=yticks, label='Probability\nDensity', steps='every other', bounds=[yticks[0], yticks[-1]])
            ax.grid(color='k', alpha=0.3, linestyle=':')
        # leg = Legend(ax, [handle], [label], loc=locations[idx], facecolor=facecolor, fontsize=self.labelsize, bbox_to_anchor=bbox_to_anchor, borderaxespad=borderaxespad, frameon=frameon)
        leg = self.autoposition_standard_legend(fig, axes, loc='auto', ncol='auto')
        fig.tight_layout()
        fig.suptitle(header, y=1.05, fontsize=self.titlesize)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def subview_point_estimators(self, fig, axes, facecolors=('steelblue', 'darkorange'), ith_resample=0, include_speed_type=False):
        """ """
        header = r'Point Estimators'
        keys = ('mean', 'median')
        labels = (r'Mean $\mu_{\theta}(j)$', r'Median $\tilde\mu_{\theta}(j)$')
        symbols = ('o', '*')
        yticks = np.arange(0, 1.1, 0.2)
        for idx, ax in enumerate(axes.ravel()):
            TimeSeries = self.MultipleTimeSeries[idx]
            xmin, xmax = TimeSeries.ExtremalIndex.js[0] - 1, TimeSeries.ExtremalIndex.js[-1] + 1
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            ax.scatter(TimeSeries.ExtremalIndex.js, TimeSeries.ExtremalIndex.point_estimators['minimum'], color='k', marker='_')
            ax.scatter(TimeSeries.ExtremalIndex.js, TimeSeries.ExtremalIndex.point_estimators['maximum'], color='k', marker='_')
            for key, fcolor, label, symbol in zip(keys, facecolors, labels, symbols):
                ax.scatter(TimeSeries.ExtremalIndex.js, TimeSeries.ExtremalIndex.point_estimators[key], color=fcolor, alpha=0.5, label=label, marker=symbol)
            ax = self.develop_x_axis(ax, ticks=TimeSeries.ExtremalIndex.js, label=r'$j$ ($log_2 hours$)', steps='every other', bounds=[xmin, xmax])
            ax = self.develop_y_axis(ax, ticks=yticks, label=r'$\theta(j)$', steps='every other', bounds=[yticks[0], yticks[-1]])
            ax.set_title(partial_title, fontsize=self.labelsize)
            ax.grid(color='k', alpha=0.3, linestyle=':')
        leg = self.autoposition_standard_legend(fig, axes, loc='auto', ncol='auto')
        fig.suptitle(header, fontsize=self.titlesize, y=1.05)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def subview_extremal_index_histogram(self, fig, axes, density='observed', facecolor=('steelblue', 'darkorange', 'green', 'red'), apply_overlay=True, include_speed_type=False):
        """ """
        available_density = ('probability', 'observed')
        xticks = np.arange(0, 1.05, 0.1)
        if density not in available_density:
            raise ValueError("unknown density: {}; available density: {}".format(density, available_density))
        if density == 'probability':
            f = lambda cls : cls.normalized_counts
            ylabel = r'Probability Density'
            yticks = None
            ylim = None
        elif density == 'observed':
            f = lambda cls : cls.observed_counts
            ylabel = r'Observed Frequency'
            yticks = np.arange(0, 1501, 50).astype(int)
            ylim = [0, 1500]
        else:
            raise ValueError("density = 'probability' or 'observed'")
        header = r'Histogram of Extremal Index $\hat\theta$' + '\n' + r' via ${}$ resamples'.format(self.MultipleTimeSeries[0].UnbiasedEstimators.nresamples)
        ys = []
        for idx, TimeSeries in enumerate(self.MultipleTimeSeries):
            mean_label = r'Mean $\mu_{\theta} =$ ' + r'${:.2f}$'.format(np.mean(TimeSeries.ExtremalIndex.Histogram.data))
            median_label = r'Median $\tilde\mu_{\theta} =$' + r'${:.2f}$'.format(np.median(TimeSeries.ExtremalIndex.Histogram.data))
            y = f(TimeSeries.ExtremalIndex.Histogram)
            ymax = TicksLabelsMethods.round_to_nearest(max(y), 'up', 10)
            ys.append(ymax)
            if yticks is None:
                yticks = np.arange(0, max(ys)+9, 10, dtype=int)
            label = '{}\n{}'.format(mean_label, median_label)
            ax = self.select_specific_axis(axes, idx, apply_overlay)
            ax.bar(TimeSeries.ExtremalIndex.Histogram.midpoints, y, width=TimeSeries.ExtremalIndex.Histogram.bin_widths, facecolor=facecolor[idx], alpha=1/self.n, edgecolor=None, label=label, align='center')
            ax = self.develop_x_axis(ax, ticks=xticks, label=r'$\theta$', steps='every other', bounds=[0, 1])
            ax = self.develop_y_axis(ax, ticks=yticks, label=ylabel, steps='every other', bounds=ylim)
            ax.grid(color='k', alpha=0.3, linestyle=':')
            if apply_overlay is True:
                handle = self.get_rectangular_handle(fill=False, edgecolor='none', visible=False, xyi=(0, 0), xyf=(0.1, 0.1))
                self.autoposition_overlay_legend(ax, handle, label, idx, facecolor[idx], loc='upper left', include_speed_type=include_speed_type)
        self.fig.tight_layout()
        fig.suptitle(header, fontsize=self.titlesize, y=1.05)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def subview_moment_estimates(self, fig, axes, facecolors=('darkgreen', 'darkorange', 'steelblue'), include_speed_type=False):
        """ """
        if axes is None:
            axes = []
        else:
            raise ValueError("this method initializes axes; axes = None")
        biases = ('threshold', 'first order', 'baseline')
        cmap = dict(zip(biases, facecolors))
        nrows, ncols = 2, self.n
        header = r'Moment Estimator $\theta$ & Declustering Time Threshold $T_C$' + '\n' + r'as a function of Speed Threshold $V_{threshold}$'
        for idx in range(2*self.n):
            ax = fig.add_subplot(nrows, ncols, idx+1)
            axes.append(ax)
        axes = np.array(axes).reshape((nrows, ncols))
        f_label = lambda bias : '{}'.format(bias.title()) if bias == 'baseline' else '{} Bias'.format(bias.title())
        yticks_theta = np.round(np.arange(0, 1.01, 0.1), decimals=1)
        for idx, TimeSeries in enumerate(self.MultipleTimeSeries):
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            ax_top = axes[0, idx]
            ax_btm = axes[1, idx]
            for bias in biases:
                x = TimeSeries.storage['estimates by speed'][bias]['vthreshold']
                y = TimeSeries.storage['estimates by speed'][bias]['moment estimator']
                z = TimeSeries.storage['estimates by speed'][bias]['time threshold']
                label = f_label(bias)
                ax_top.plot(x, y, color=cmap[bias], label=label, marker='.', alpha=1/3)
                ax_top = self.develop_x_axis(ax_top, ticks=x[::2], label=None, steps='every ten', bounds=[x[0], x[-1]])
                ax_top = self.develop_y_axis(ax_top, ticks=yticks_theta, label=r'$\theta$', steps='every other', bounds=None)
                ax_top.grid(color='k', alpha=0.3, linestyle=':')
                ax_top.set_title(partial_title, fontsize=self.labelsize)
                ax_btm.plot(x, z, color=cmap[bias], label=label, marker='.', alpha=1/3)
                ax_btm = self.develop_x_axis(ax_btm, ticks=x[::2], label=r'$V_{threshold}$ ($\frac{km}{s}$)', steps='every ten', bounds=[x[0], x[-1]])
                ax_btm = self.develop_y_axis(ax_btm, ticks=None, label=r'$T_C$ (hours)', steps='every other', bounds=None)
                ax_btm.grid(color='k', alpha=0.3, linestyle=':')
        TicksLabelsMethods().autoconfigure_axes(axes, key='edge column')
        for ax in axes[0, :].ravel():
            ax.set_xticklabels([])
        if self.n > 1:
            leg = self.autoposition_standard_legend(fig, axes, loc='auto', ncol='auto')
        else:
            fig.subplots_adjust(bottom=0.135)
            try:
                handles, labels = axes[0].get_legend_handles_labels()
            except:
                handles, labels = axes[0, 0].get_legend_handles_labels()
                leg = fig.legend(handles, labels, loc='lower center', ncol=len(labels), fontsize=self.labelsize, fancybox=True, borderaxespad=0.2)
        fig.suptitle(header, fontsize=self.titlesize, y=1.05)

    def subview_sunspot_correlations(self, fig, axes, vthresholds=(0, 800, 1000), search_parameters=['speed']*3, search_conditions=['greater than or equal']*3, counts_by='monthly', facecolor=('darkorange', 'steelblue', 'green', 'purple'), equalize_limits=True, search_kwargs=None, include_speed_type=False):
        """ """
        if axes is None:
            axes = []
        else:
            raise ValueError("this method initializes axes; axes = None")
        header = 'Correlation between\nSunspots & Extreme Events'
        nrows = len(vthresholds) + 1
        ylimits = {'daily' : (200, 30, 10, 6), 'monthly' : (8000, 300, 80, 60), 'yearly' : (30*10**3, 2500, 300, 150)}
        yspaces = {'daily' : (10, 10, 10, 5), 'monthly' : (400, 20, 10, 5), 'yearly' : (10, 10, 10, 5)}
        if search_kwargs is None:
            search_kwargs = [dict() for row in range(nrows-1)]
        for idx in range(self.n * nrows):
            ax = fig.add_subplot(nrows, self.n, idx+1)
            axes.append(ax)
        axes = np.array(axes).reshape((nrows, self.n))
        for ith_row in range(nrows):
            for ith_col, ax in enumerate(axes[ith_row, :]):
                TimeSeries = self.MultipleTimeSeries[ith_col]
                Correlator = TimeSeries.storage['SunSpot']
                if ith_row == 0:
                    data = Correlator.subselect_sunspots(counts_by, TimeSeries)
                    label = '${:,}$ Sunspots'.format(np.sum(data['count']))
                    ylabel = '{}\nSunspots'.format(counts_by.title())
                    partial_title = self.get_partial_title(TimeSeries, include_speed_type)
                    ax.set_title(partial_title, fontsize=self.labelsize)
                else:
                    data = Correlator.subselect_ejecta(counts_by, TimeSeries, search_parameters[ith_row-1], search_conditions[ith_row-1], search_values=vthresholds[ith_row-1], **search_kwargs[ith_row-1])
                    label = '${:,}$ CMEs'.format(np.sum(data['count']))
                    ylabel = '{}\nEjecta\n({} $≥ {}$'.format(counts_by.title(), search_parameters[ith_row-1].title(), vthresholds[ith_row-1]) + r' $\frac{km}{s}$)'
                ax.plot(data['datetime number'], data['count'], label=label, color=facecolor[ith_row], linestyle='-', linewidth=1)
                ax = self.transform_as_datetime(ax, axis='x')
                if ith_row == nrows - 1:
                    ax.set_xlabel('Date', fontsize=self.labelsize)
                else:
                    ax.set_xticklabels([])
                if equalize_limits is True:
                    ymax = ylimits[counts_by][ith_row]
                    yticks = np.arange(0, ymax, yspaces[counts_by][ith_row])
                    ax = self.develop_y_axis(ax, ticks=(np.ceil(yticks / 10) * 10).astype(int), label=None, steps='every other', bounds=(0, ymax))
                if ith_col in (0, len(axes[ith_row, :]) - 1):
                    ax.set_ylabel(ylabel, fontsize=self.labelsize)
                    if ith_col == len(axes[ith_row, :]) - 1:
                        ax.yaxis.tick_right()
                        ax.yaxis.set_label_position('right')
                else:
                    ax.set_yticklabels([])
                ax.grid(color='k', alpha=0.3, linestyle=':')
                ax.legend(loc='upper right', fontsize=self.ticksize)
        fig.align_ylabels()
        fig.suptitle(header, fontsize=self.titlesize, y=1.05)

    def subview_chronological_cluster_sizes(self, fig, axes, include_lone_ejecta=False, facecolors=('r', 'b'), include_speed_type=False):
        """ """
        if include_lone_ejecta is True:
            search_conditions = 'greater than or equal'
        elif include_lone_ejecta is False:
            search_conditions = 'greater than'
        else:
            raise ValueError("include_lone_ejecta = True or False")
        load_keys = ('ith cluster', 'cluster size')
        f = lambda args : datetime.datetime.strptime(args.split(' ')[0], "%Y/%m/%d")
        g = lambda dtime : mdates.date2num(f(dtime))
        lone_width = g('1999/10/27 01:00:00') - g('1999/10/27 00:00:00')
        for idx, ax in enumerate(axes.ravel()):
            TimeSeries = self.MultipleTimeSeries[idx]
            clusters = TimeSeries.TemporalClustering.search(search_parameters='cluster size', search_conditions=search_conditions, search_values=1, load_keys=load_keys)
            x, w, y = [], [], []
            icolors, handles = [], []
            for ith_cluster, data in enumerate(clusters['elapsed']):
                ith_color = facecolors[0] if ith_cluster % 2 == 0 else facecolors[1]
                icolors.append(ith_color)
                dt = clusters['datetime'][ith_cluster]
                dt_init = dt[0]
                dt_fin = dt[-1]
                if data.size == 1:
                    xi = g(dt_init)
                    wi = lone_width
                else:
                    tmp = np.array([g(dt_init), g(dt_fin)])
                    xi = np.mean(tmp)
                    wi = np.diff(tmp)
                x.append(xi)
                w.append(wi)
                y.append(data.size)
            ncluster = clusters['elapsed'].size
            nejecta = np.sum(y)
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            partial_title += ':\n${}$ Clusters from ${}$ Extreme Events'.format(ncluster, nejecta)
            ymax = max(y)
            yticks = np.linspace(0, ymax+1, ymax+2).astype(int)
            ax.bar(x, y, width=np.array(w).reshape(-1), color=icolors, align='center')
            ax = self.transform_as_datetime(ax, axis='x')
            ax = self.develop_y_axis(ax, ticks=yticks, label='Cluster Size', steps='every other', bounds=None)
            ax.grid(color='k', alpha=0.3, linestyle=':')
            ax.set_title(partial_title, fontsize=self.labelsize)
        handles = [tuple((self.get_rectangular_handle(fill=True, visible=True, facecolor=facecolors[ith_color]) for ith_color in range(len(facecolors))))]
        kwargs = {'handler_map' : {tuple: matplotlib.legend_handler.HandlerTuple(None)}}
        fig.subplots_adjust(bottom=0.125)
        fig.legend(handles=handles, labels=['Consecutive Clusters'], loc='lower center', fontsize=self.labelsize, **kwargs)
        fig.suptitle('Chronological Cluster Sizes', fontsize=self.titlesize, y=1.05)
        if self.n != 1:
            TicksLabelsMethods().autoconfigure_axes(axes, key='standard')

    def subview_relative_cluster_statistics(self, fig, axes, stats=('number of ejecta', 'number of clusters', 'probability'), facecolors=('r', 'darkorange', 'darkgreen'), edgecolors=('k', 'b', 'r'), include_speed_type=False):
        """ """
        if axes is None:
            axes = []
        else:
            raise ValueError("this method initializes axes; axes = None")
        nrows = len(stats)
        for idx in range(self.n * nrows):
            ax = fig.add_subplot(nrows, self.n, idx+1)
            axes.append(ax)
        axes = np.array(axes).reshape((nrows, self.n))
        handles, labels = [], []
        for ith_row, ith_statistic, ith_color, ith_edgecolor in zip(list(range(nrows)), stats, facecolors, edgecolors):
            ith_handle = self.get_rectangular_handle(fill=True, visible=True, facecolor=ith_color, edgecolor=ith_edgecolor)
            ith_label = ith_statistic.title()
            handles.append(ith_handle)
            labels.append(ith_label)
            for ith_col, ax in enumerate(axes[ith_row, :]):
                TimeSeries = self.MultipleTimeSeries[ith_col]
                x = TimeSeries.storage['relative cluster statistics']['cluster size']
                y = TimeSeries.storage['relative cluster statistics'][ith_statistic]
                xmin, xmax = min(x), max(x)
                xticks = np.linspace(xmin, xmax, xmax - xmin + 1).astype(int)
                if ith_row == 0:
                    partial_title = self.get_partial_title(TimeSeries, include_speed_type)
                    ax.set_title(partial_title, fontsize=self.labelsize)
                if ith_statistic == 'probability':
                    yticks = np.round(np.arange(0, 1.01, 0.1), decimals=1)
                else:
                    ymax = TicksLabelsMethods.round_to_nearest(max(y), 'up', 10)
                    ydelta = 10
                    yticks = np.arange(0, ymax + ydelta - 1, ydelta//2).astype(int)[::2]
                    if ith_statistic == 'number of ejecta':
                        yticks = yticks[::2]
                        # yticks = np.append(yticks, yticks[-1] + 2*ydelta)
                ax.bar(x, y, width=1, facecolor=ith_color, edgecolor=ith_edgecolor)
                ax = self.develop_x_axis(ax, ticks=xticks, label='Cluster Size', steps='every other', bounds=None)
                ax = self.develop_y_axis(ax, ticks=yticks, label=r'{}'.format(ith_statistic.title()), steps='every other', bounds=None)
                ax.grid(color='k', alpha=0.3, linestyle=':')
        # TicksLabelsMethods().autoconfigure_axes(axes, key='standard')
        TicksLabelsMethods().autoconfigure_axes(axes, key='edge column')
        for ax in axes[0, :].ravel():
            ax.set_xlabel('')
            ax.set_xticklabels([])
            ax.xaxis.tick_bottom()
        fig.align_ylabels()
        fig.suptitle(r'Relative Cluster Statistics', fontsize=self.titlesize, y=1.03)
        fig.subplots_adjust(bottom=0.13)
        fig.legend(handles=handles, labels=labels, loc='lower center', fontsize=self.labelsize, ncol=len(labels))

    def subview_durational_histograms(self, fig, axes, cluster_size='all', search_conditions='exact match', keys=('intra-time histogram', 'intra-duration histogram', 'inter-duration histogram'), facecolors=('r', 'darkgreen', 'b'), edgecolors=('none', 'none', 'none'), tc_color='darkorange', include_speed_type=False):
        """ """
        if axes is None:
            axes = []
        else:
            raise ValueError("this method initializes axes; axes = None")
        nrows = len(keys)
        for idx in range(self.n * nrows):
            ax = fig.add_subplot(nrows, self.n, idx+1)
            axes.append(ax)
        axes = np.array(axes).reshape((nrows, self.n))
        if cluster_size == 'all':
            header = r'Durational Histograms of All Clusters'
        else:
            if search_conditions in ('exact match', 'equal', 'equality'):
                header = r'Durational Histograms of Size-${}$ Clusters'.format(cluster_size)
            else:
                header = 'Durational Histograms\nfor Clusters of Size {} ${}$'.format(search_conditions.title(), cluster_size)
        handles, labels = [], []
        tc_booleans = []
        for ith_row, ith_key, ith_color, ith_edgecolor in zip(list(range(nrows)), keys, facecolors, edgecolors):
            htype = ith_key[:].replace('-', ' ').replace('histogram', '')[:-1]
            ith_handle = self.get_rectangular_handle(fill=True, visible=True, facecolor=ith_color, edgecolor=ith_edgecolor)
            ith_label = ith_key.title()
            handles.append(ith_handle)
            labels.append(ith_label)
            for ith_col, ax in enumerate(axes[ith_row, :]):
                TimeSeries = self.MultipleTimeSeries[ith_col]
                ClusterHistograms = TimeSeries.TemporalClustering.get_cluster_histograms(cluster_size, search_conditions)
                if ith_key == 'intra-time histogram':
                    H = ClusterHistograms.intra_time_histogram
                    xticks, yticks, xlim, ylim = TicksLabelsMethods().get_intra_time_ticks(H)
                    show_tc = False
                elif ith_key == 'intra-duration histogram':
                    H = ClusterHistograms.intra_duration_histogram
                    xticks, yticks, xlim, ylim = TicksLabelsMethods().get_intra_duration_ticks(H)
                    show_tc = True
                elif ith_key == 'inter-duration histogram':
                    H = ClusterHistograms.inter_duration_histogram
                    xticks, yticks, xlim, ylim = TicksLabelsMethods().get_inter_duration_ticks(H)
                    show_tc = True
                ax.bar(H.midpoints-0.5, H.observed_counts, width=H.bin_widths, facecolor=ith_color, edgecolor=ith_edgecolor)
                ax = self.develop_x_axis(ax, ticks=xticks, label='Elapsed Time (Hours)', steps='every other', bounds=xlim)
                ax = self.develop_y_axis(ax, ticks=yticks, label='Frequency', steps='every other', bounds=ylim)
                if ith_row == 0:
                    ncluster = len(ClusterHistograms.clusters['elapsed'])
                    nejecta = np.concatenate(ClusterHistograms.clusters['elapsed'], axis=0).size
                    partial_title = self.get_partial_title(TimeSeries, include_speed_type)
                    partial_title += '\n\n{} Clusters from {} Extreme Events\n$T_C = {}$ Hours'.format(ncluster, nejecta, TimeSeries.declustering_time)
                    ax.set_title(partial_title, fontsize=self.labelsize)
                if show_tc is True:
                    ax.axvline(TimeSeries.TemporalClustering.time_threshold, ymin=0, ymax=ylim[-1], color=tc_color, linestyle='-', marker=None)
                tc_booleans.append(show_tc)
                ax.grid(color='k', alpha=0.3, linestyle=':')
                if ith_col == len(axes[ith_row, :]) - 1:
                    ax.yaxis.tick_right()
                    ax.yaxis.set_label_position('right')
        # TicksLabelsMethods().autoconfigure_axes(axes, key='standard')
        # TicksLabelsMethods().autoconfigure_axes(axes, key='edge column')
        if np.any(np.array(tc_booleans) == True):
            tc_handle = self.get_rectangular_handle(fill=True, visible=True, facecolor=tc_color, edgecolor='none')
            tc_label = 'Time Threshold $T_C$'
            handles.append(tc_handle)
            labels.append(tc_label)
        fig.align_ylabels()
        fig.suptitle(header, fontsize=self.titlesize, y=1.065)
        fig.subplots_adjust(bottom=0.15, hspace=0.2, wspace=0.2)
        fig.legend(handles=handles, labels=labels, loc='lower center', fontsize=self.labelsize, ncol=len(labels))

    def subview_dim3_speed_distribution(self, fig, axes, err='minimum gtest estimation', elev=30, azim=345, ncontours=15, cmap='plasma', alpha=0.8, antialiased=True, shade=True, include_colorbar=True, include_speed_type=False):
        """ """
        if axes is None:
            axes = dict(nrows=None, ncols=None, swap_rc=False, apply_overlay=False)
        elif not isinstance(axes, dict):
            raise ValueError("axes = None or type <dict>")
        available_err = ('maximum likelihood estimation', 'minimum gtest estimation') ## LINE 978 OF timeseries_methods.py
        if err not in available_err:
            raise ValueError("unknown err: {}; available err: {}".format(err, available_err))
        nrows, ncols = self.get_rows_columns(**axes)
        for idx, TimeSeries in enumerate(self.MultipleTimeSeries):
            EDM = TimeSeries.storage['Speed Error-Space'][err]
            prms = EDM.SD.result[err]['parameters']
            zmin = EDM.SD.result[err]['objective']
            original_zlabel = EDM.get_zlabel(err)
            zlabel = '{} at ({:.2f}, {:.2f})'.format(original_zlabel, *prms)
            invert_z = EDM.check_invert_z(err)
            if err == 'maximum likelihood estimation':
                partial_zlabel = r'$-ln(L)$'
            elif err == 'minimum gtest estimation':
                partial_zlabel = r'$G$'
            if invert_z is True:
                zmin = - zmin
            if zmin >= 0:
                central_color = 'k'
            else:
                central_color = 'white'
            partial_title = self.get_partial_title(TimeSeries, include_speed_type)
            ax = fig.add_subplot(nrows, ncols, idx+1, projection='3d')
            norm = Normalize(np.amin(EDM.Z), np.amax(EDM.Z))
            ax.plot_surface(EDM.X, EDM.Y, EDM.Z, cmap=cmap, norm=norm, alpha=alpha, antialiased=antialiased, shade=shade)
            if ((elev is not None) or (azim is not None)):
                ax.view_init(elev, azim)
            ax.contour(EDM.X, EDM.Y, EDM.Z, ncontours, linewidth=3, cmap=cmap, linestyles='solid', offset=-1)
            ax.contour(EDM.X, EDM.Y, EDM.Z, ncontours, linewidth=3, colors='k', linestyles='solid')
            # ax.contour(EDM.X, EDM.Y, EDM.Z, ncontours, linewidth=3, colors=central_color, linestyles='solid')
            # ax.scatter(prms[0], prms[1], zmin, color=central_color, marker='*', linewidth=0, s=20, label='{} at ({:.2f}, {:.2f})'.format(zlabel, *prms))
            ax.grid(color='k', alpha=0.3, linestyle=':')
            if include_colorbar is True:
                zticks = ax.get_zticks()
                if err == 'maximum likelihood estimation':
                    zticklabels = ['{:,}'.format(zticks[idx]) if idx % 2 == 0 else '' for idx in range(len(zticks))]
                elif err == 'minimum gtest estimation':
                    zticklabels = ['{:,}'.format(zt) for zt in zticks]
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm._A = []
                cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', shrink=0.75, pad=0.1)
                # cbar = fig.colorbar(sm, ax=ax, orientation='vertical', shrink=0.75, pad=0.1)
                cbar.set_ticks(zticks)
                cbar.set_ticklabels(zticklabels)
                cbar.ax.tick_params(labelsize=self.ticksize)
                cbar.ax.set_title(zlabel, fontsize=self.labelsize)
            ax.set_title(partial_title, fontsize=self.labelsize)
            ax.xaxis.set_rotate_label(False)
            ax.yaxis.set_rotate_label(False)
            ax.zaxis.set_rotate_label(False)
            ax.tick_params(axis='both', labelsize=self.ticksize)
            ax.set_xlabel(r'$\mu$', fontsize=self.labelsize) #, labelpad=20, rotation=0)
            ax.set_ylabel(r'$\sigma$', fontsize=self.labelsize) #, labelpad=20, rotation=0)
            ax.set_zlabel(partial_zlabel, fontsize=self.labelsize)
            # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fontsize=self.ticksize)
        fig.subplots_adjust(wspace=0.2, hspace=0.35)
        fig.tight_layout()
        fig.suptitle('Speed Distribution Error-Space\nvia {}'.format(err.title()), fontsize=self.titlesize, y=1.05)
        # fig.subplots_adjust(bottom=0.3)
        # fig.legend(loc='lower center', fontsize=self.labelsize, mode='expand', ncol=)

    def subview_frechet_distribution(self, fig, axes, facecolors=('red', 'green', 'blue', 'orange', 'purple', 'black')):
        """ """
        if axes is not None:
            raise ValueError("axes = None")
        axes = []
        for idx in range(2):
            ax = fig.add_subplot(2, 1, idx+1)
            axes.append(ax)
        density = ('probability', 'cumulative')
        x = np.linspace(0.1, 5, 1001)
        alpha = np.array([1, 1, 2, 2, 3, 3])
        mu = np.zeros(alpha.size, dtype=int)
        sigma = np.array([1, 2, 1, 2, 1, 2])
        xticks = np.linspace(0, x[-1], x[-1]+1).astype(int)
        for key, ax in zip(density, axes):
            yn, labels = [], []
            if key == 'probability':
                yticks = np.arange(0, 1.3, 0.25)
                ylim = (0, 1.25)
                for a, m, s in zip(alpha, mu, sigma):
                    tmp = (x - m)/s
                    res = (a/s) * tmp**(-1 - a) * np.exp(-1 * (tmp ** -a))
                    label = r'$\alpha = {}, \mu = {}, \sigma = {}$'.format(a, m, s)
                    yn.append(res)
                    labels.append(label)
            else:
                ylim = (0, 1)
                yticks = np.arange(0, 1.1, 0.25)
                for a, m, s in zip(alpha, mu, sigma):
                    tmp = (x - m)/s
                    res = np.exp(-1 * (tmp ** -a))
                    label = r'$\alpha = {}, \mu = {}, \sigma = {}$'.format(a, m, s)
                    yn.append(res)
                    labels.append(label)
            for y, facecolor, label in zip(yn, facecolors, labels):
                ax.plot(x, y, color=facecolor, alpha=0.5, label=label)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks, fontsize=self.ticksize)
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks, fontsize=self.ticksize)
            ax.set_xlim([0, x[-1]])
            ax.set_ylim(ylim)
            ax.set_xlabel('x', fontsize=self.labelsize, labelpad=-7.5)
            ax.set_ylabel(r'{} Density'.format(key.title()), fontsize=self.labelsize)
            ax.grid(color='k', alpha=0.3, linestyle=':')
        axes[0].set_title(r'Fr$\acute{e}$chet Distribution', fontsize=self.titlesize)
        handles, labels = axes[0].get_legend_handles_labels()
        axes[0].set_xlabel('')
        axes[0].set_xticklabels([])
        axes[-1].set_title('')
        fig.legend(handles=handles, labels=labels, ncol=len(labels)//2, fancybox=True, loc='lower center', mode='expand', fontsize=self.labelsize)
        fig.subplots_adjust(bottom=0.175, hspace=0.3)
        fig.align_ylabels()
        fig.subplots_adjust(hspace=0.2, wspace=0.2, bottom=0.2)

    def subview_cluster_example(self, fig, axes, facecolors=('darkorange', 'steelblue', 'purple')):
        """ """
        ax = axes[0]
        time_threshold = 5
        clusters = np.array([np.array([2, 3, 5]), np.array([11, 13, 17]), np.array([24, 27])]) - 2
        ncluster = len(clusters)
        xticks = np.unique(np.sum([cluster.copy().tolist() for cluster in clusters]))
        empty_handle = self.get_rectangular_handle(fill=False, edgecolor='none', visible=False, xyi=(0, 0), xyf=(0.1, 0.1))
        intra_time, intra_duration, inter_duration = [], [], []
        intra_arrowprops = {'arrowstyle': '|-|', 'color' : 'k'}
        inter_arrowprops = {'arrowstyle': '<->', 'color' : 'gray'}
        for ith_cluster in range(len(clusters)):
            curr_cluster = clusters[ith_cluster]
            intra_time = np.diff(curr_cluster).tolist()
            intra_duration = curr_cluster[-1] - curr_cluster[0]
            ax.annotate('', xy=(curr_cluster[-1], 0.95), xycoords='data', xytext=(curr_cluster[0], 0.95), textcoords='data', arrowprops=intra_arrowprops)
            ax.text(curr_cluster[0]+1, 0.925, s=r'$\Delta T_{intra}$')
            if ith_cluster < ncluster - 1:
                next_cluster = clusters[ith_cluster+1]
                inter_duration = next_cluster[0] - curr_cluster[-1]
                ax.annotate('', xy=(next_cluster[0], 1.05), xycoords='data', xytext=(curr_cluster[-1], 1.05), textcoords='data', arrowprops=inter_arrowprops)
                ax.text(curr_cluster[-1]+2, 1.075, s=r'$\Delta T_{inter}$')
            ax.scatter(curr_cluster, np.ones(curr_cluster.size), color=facecolors[ith_cluster], label='Cluster #{}'.format(ith_cluster+1))
        ax.plot([np.nan], [np.nan], color='none', label='$T_C = {}$ hours'.format(time_threshold))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, fontsize=self.ticksize)
        ax.set_xlabel('Consecutive Differences of Elapsed Hours', fontsize=self.labelsize)
        ax.set_yticks([0.85, 1.15])
        ax.set_yticklabels([])
        ax.set_xlim([-2, 30])
        ax.set_ylim([0.85, 1.15])
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.legend(ncol=4, fancybox=True, loc='upper center', mode='expand', fontsize=self.labelsize)
        plt.subplots_adjust(top=0.8)

class FigureVisualizer(AxesMethods):

    def __init__(self, MultipleTimeSeries, key=None, saveloc=None, ticksize=7, labelsize=8, titlesize=10, extension='.png', sharex=False, sharey=False, include_speed_type=False, **kwargs):
        """ """
        super().__init__(MultipleTimeSeries, ticksize, labelsize, titlesize, saveloc, extension)
        self.key = key
        fig, axes = self.initialize(key, sharex=sharex, sharey=sharey, **kwargs)
        self.fig = fig
        self.axes = axes
        self.include_speed_type = include_speed_type

    @property
    def function_mapping(self):
        """ """
        res = {}
        res['inter-exceedance distribution decay'] = self.subview_interexceedance_distribution_against_exponential_decay
        res['alpha-hat histogram'] = self.subview_alpha_hat_histogram
        res['unbiased estimator errors'] = self.subview_unbiased_estimator_errors
        res['max spectrum power-law'] = self.subview_power_law
        res['point estimators'] = self.subview_point_estimators
        res['extremal index histogram'] = self.subview_extremal_index_histogram
        res['moment estimates'] = self.subview_moment_estimates
        res['speed tail'] = self.subview_speed_tail
        res['speed distribution'] = self.subview_speed_distribution
        res['sunspot correlations'] = self.subview_sunspot_correlations
        res['chronological cluster size'] = self.subview_chronological_cluster_sizes
        res['relative cluster statistics'] = self.subview_relative_cluster_statistics
        res['durational histograms'] = self.subview_durational_histograms
        res['dim3 speed distribution'] = self.subview_dim3_speed_distribution
        res['frechet distribution'] = self.subview_frechet_distribution
        res['illustrative cluster example'] = self.subview_cluster_example
        return res

    @property
    def available_keys(self):
        """ """
        return list(self.function_mapping.keys())

    @property
    def savepath(self):
        """ """
        if self.saveloc is None:
            return None
        elif isinstance(self.saveloc, str):
            fname = self.key[:].replace(' ', '_')
            return '{}{}{}'.format(self.saveloc, fname, self.extension)

    def finalize(self, dpi, bbox_inches, pad_inches):
        """
        path                :   type <str> or None
        dpi                 :   type <int>
        bbox_inches         :   type <str>
        pad_inches          :   type <int / float>
        """
        if self.savepath is None:
            plt.show()
            plt.close(self.fig)
        elif isinstance(self.saveloc, str):
            self.fig.savefig(self.savepath, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches)
            plt.close(self.fig)
        else:
            raise ValueError("path = None or type <str>")

    def dispatch(self, dpi=800, bbox_inches='tight', pad_inches=0.1, **kwargs):
        """ """
        if self.key not in self.available_keys:
            raise ValueError("unknown key: {}; available keys: {}".format(self.key, self.available_keys))
        f = self.function_mapping[self.key]
        if self.key not in ('frechet distribution', 'illustrative cluster example'):
            kwargs = dict(**kwargs)
            kwargs['include_speed_type'] = self.include_speed_type
        f(self.fig, self.axes, **kwargs)
        self.finalize(dpi, bbox_inches, pad_inches)

    @staticmethod
    def get_identifier(TimeSeries):
        """ """
        solar_cycle = TimeSeries.identity['solar cycle']
        cycle_type = TimeSeries.identity['cycle type']
        speed_type = TimeSeries.identity['speed type']
        vthreshold = TimeSeries.storage['vthreshold']
        identifier = 'SC-{}_subinterval-{}_{}_vthreshold-{}'.format(solar_cycle, cycle_type, speed_type, vthreshold)
        return solar_cycle, cycle_type, speed_type, vthreshold, identifier

    def get_text_path(self, partial_fname, identifier):
        """ """
        path = '{}{}{}'.format(self.saveloc, partial_fname, identifier)
        return path

    def save_analysis_results(self):
        """ """
        if self.saveloc is None:
            raise ValueError("specify saveloc as directory to save .txt file")
        partial_fname = 'analysis_results__'
        header = '\n** ANALYSIS RESULTS **\n'
        for TimeSeries in self.MultipleTimeSeries:
            solar_cycle, cycle_type, speed_type, vthreshold, identifier = self.get_identifier(TimeSeries)
            path = self.get_text_path(partial_fname, identifier)
            extremal_index = np.mean(TimeSeries.ExtremalIndex.theta_hat)
            clusters = TimeSeries.TemporalClustering.search('cluster size', 'greater than', 1, load_keys=('ith cluster', 'cluster size'))
            ncluster = len(clusters['speed'])
            nejecta = np.concatenate(clusters['speed'], axis=0).size
            avg_cluster_size = np.mean([arr.size for arr in clusters['speed']])
            string = '\n .. ALPHA HAT:\n{}\n\n .. INTERCEPT HAT:\n{}\n\n .. EXTREMAL INDEX:\n{}\n'.format(TimeSeries.UnbiasedEstimators.alpha('hat'), TimeSeries.UnbiasedEstimators.intercept('hat'), extremal_index)
            string += '\n .. EXTREMAL MOMENT:\n{}\n\n .. TIME THRESHOLD:\n{}\n'.format(TimeSeries.extremal_moment, TimeSeries.declustering_time)
            string += '\n .. NUMBER OF EJECTA:\n{}\n\n .. NUMBER OF CLUSTERS:\n{}\n'.format(nejecta, ncluster)
            string += '\n .. AVERAGE CLUSTER SIZE:\n{}\n\n .. EXPECTED AVG SIZE:\n{}\n\n'.format(avg_cluster_size, 1/extremal_index)
            with open(path, "w+") as text_file:
                text_file.write('{}\n'.format(header))
                text_file.write('{}'.format(string))
                for key in list(TimeSeries.UnbiasedEstimators.speed_thresholds.keys()):
                    text_file.write(" .. APPROXIMATE SPEED THRESHOLD ({}):\t{} km/s\n".format(key, TimeSeries.UnbiasedEstimators.speed_thresholds[key]))

    def save_clusters(self, search_parameters, search_conditions, search_values, fkeys=None, load_keys=('ith cluster', 'cluster size'), apply_to='all', view_keys=None, use_abs=False, use_halo=None):
        """ """
        if self.saveloc is None:
            raise ValueError("specify saveloc as directory to save .txt file")
        partial_fname = 'cluster_data__'
        header = '\n** CLUSTER DATA **\n'
        for TimeSeries in self.MultipleTimeSeries:
            solar_cycle, cycle_type, speed_type, vthreshold, identifier = self.get_identifier(TimeSeries)
            path = self.get_text_path(partial_fname, identifier)
            clusters = TimeSeries.search_cme_clusters(TimeSeries.clusters, search_parameters, search_conditions, search_values, fkeys, load_keys, apply_to, view_keys, use_abs, use_halo)
            with open(path, "w+") as text_file:
                text_file.write('{}\n'.format(header))
                for ith_cluster, (dt, hr, vi) in enumerate(zip(clusters['datetime'], clusters['elapsed'], clusters['speed'])):
                    string = '\nCLUSTER #{}:\n\n .. DATETIMES:\n{}\n\n .. ELAPSED HOURS:\n{}\n\n .. SPEEDS:\n{}\n\n'.format(ith_cluster+1, dt, hr, vi)
                    text_file.write('{}'.format(string))

    def save_cluster_statistics(self, search_parameters, search_conditions, search_values, fkeys=None, load_keys=('ith cluster', 'cluster size'), apply_to='all', view_keys=None, use_abs=False, use_halo=None, skip_lone_ejecta=True):
        """ """
        if self.saveloc is None:
            raise ValueError("specify saveloc as directory to save .txt file")
        partial_fname = 'cluster_stats__'
        header = '\n** CLUSTER STATISTICS **\n'
        if skip_lone_ejecta is True:
            condition = 'greater than or equal'
            value = 2
        elif skip_lone_ejecta is False:
            condition = 'greater than or equal'
            value = 1
        for TimeSeries in self.MultipleTimeSeries:
            solar_cycle, cycle_type, speed_type, vthreshold, identifier = self.get_identifier(TimeSeries)
            path = self.get_text_path(partial_fname, identifier)
            # base = self.TimeSeries.clusters
            base_clusters = TimeSeries.search_cme_clusters(TimeSeries.clusters, 'cluster size', condition, value, fkeys=None, load_keys=load_keys, apply_to='all', use_halo=use_halo)
            clusters = TimeSeries.search_cme_clusters(base_clusters, search_parameters, search_conditions, search_values, fkeys, load_keys, apply_to, view_keys, use_abs, use_halo)
            total_clusters, nclusters = len(base_clusters['speed']), len(clusters['speed'])
            total_ejecta, nejecta = np.concatenate(base_clusters['speed'], axis=0).size, np.concatenate(clusters['speed'], axis=0).size
            relative_probability = nejecta/total_ejecta
            CDH = TimeSeries.get_cluster_duration_histogram_instance(clusters)
            mean_duration = np.mean(CDH.intra_duration)
            stdev_duration = np.std(CDH.intra_duration)
            string = "\n .. search condition:\n\t{}\n".format(search_conditions)
            string += "\n .. search parameter:\n\t{}\n".format(search_parameters)
            string += "\n .. search value:\n\t{}\n".format(search_values)
            string += "\n .. number of clusters:\n\t{} of {}\n".format(nclusters, total_clusters)
            string += "\n .. number of ejecta:\n\t{} of {}\n".format(nejecta, total_ejecta)
            string += "\n .. relative probability:\n\t{}\n".format(relative_probability)
            string += "\n .. mean duration:\n\t{} ± {}\n".format(mean_duration, stdev_duration)
            with open(path, "w+") as text_file:
                text_file.write('{}\n'.format(header))
                text_file.write(string)
