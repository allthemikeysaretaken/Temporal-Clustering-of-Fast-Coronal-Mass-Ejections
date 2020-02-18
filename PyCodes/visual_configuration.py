import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, YearLocator, MonthLocator, DayLocator, WeekdayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.ticker as ticker
from matplotlib.legend_handler import HandlerTuple, HandlerRegularPolyCollection
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.mplot3d import Axes3D

class ScatterHandler(HandlerRegularPolyCollection):

    def update_prop(self, legend_handle, orig_handle, legend):
        """ """
        legend._set_artist_props(legend_handle)
        legend_handle.set_clip_box(None)
        legend_handle.set_clip_path(None)

    def create_collection(self, orig_handle, sizes, offsets, transOffset):
        """ """
        p = type(orig_handle)([orig_handle.get_paths()[0]], sizes=sizes, offsets=offsets, transOffset=transOffset, cmap=orig_handle.get_cmap(), norm=orig_handle.norm)
        a = orig_handle.get_array()
        if type(a) != type(None):
            p.set_array(np.linspace(a.min(), a.max(), len(offsets)))
        else:
            self._update_prop(p, orig_handle)
        return p

class Symbolizer():

    """
    This class contains mappings from strings to symbols
    for the purpose of plotting figures.
    """

    def __init__(self):
        super().__init__()
        labels = dict()
        labels['comparison'] = self.load_comparison_labels()
        labels['cme'] = self.load_cme_labels()
        labels['optimization'] = self.load_optimization_labels()
        # labels['flare'] = dict()
        labels['unbiased estimators'] = self.load_unbiased_estimators_labels()
        self._available = dict(label=labels)

    @property
    def available(self):
        return self._available

    @staticmethod
    def load_comparison_labels():
        keys = ('equal', 'greater than', 'greater than or equal', 'less than', 'less than or equal', 'not equal')
        labels = ('=', '>', '≥', '<', '≤', '≠')
        return dict(zip(keys, labels))

    @staticmethod
    def load_cme_labels():
        parameters = ('speed', 'mass', 'acceleration')
        units = (r'$\frac{km}{s}$', r'kg', r'$\frac{km}{s^2}$')
        speed_types = ('linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed')
        speed_labels = (r'$V_{CME, linear}$', r'$V_{CME, 20 R_{\odot}, i}$', r'$V_{CME, 20 R_{\odot}, f}$', r'$V_{CME, 20 R_{\odot}}$')
        return {'speed type' : dict(zip(speed_types, speed_labels)), 'unit' : dict(zip(parameters, units))}

    @staticmethod
    def load_optimization_labels():
        _available = dict()
        _available['maximum likelihood estimation'] = r'$ln$ $L_{max}$'
        _available['chi square'] = dict(reduced=r'$\chi_{min, red}^2$', ordinary=r'$\chi_{min}^2$')
        _available['g-test'] = dict(reduced=r'$G_{min, red}$', ordinary=r'$G_{min}$')
        return _available

    @staticmethod
    def load_unbiased_estimators_labels():
        keys = ('alpha', 'intercept', 'theta')
        labels = (r'$\hat\alpha$', r'$\hatC$', r'$\hat\theta$')
        return dict(zip(keys, labels))

    @staticmethod
    def round_up(values, nearest):
        """
        values:
            type <int / float / tuple / list / array>

        nearest:
            type <int>
        """
        if not isinstance(nearest, int):
            raise ValueError("invalid type(nearest): {}".format(type(nearest)))
        try:
            return int(np.ceil(max(values) / float(nearest)) * int(nearest))
        except:
            return int(np.ceil(values / float(nearest)) * int(nearest))

    @staticmethod
    def get_acronym(words):
        """

        """
        try:
            abbreviations = [phrase[0] for phrase in words.split(' ')]
        except:
            abbreviations = [phrase[0] for phrase in words.split('-')]
        return "".join(abbreviations)

class VisualConfiguration(Symbolizer):

    def __init__(self, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize):
        """
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
        super().__init__()
        self.directory = directory
        self.ticksize = ticksize
        self.labelsize = labelsize
        self.textsize = textsize
        self.titlesize = titlesize
        self.headersize = headersize
        self.cellsize = cellsize
        layouts = ('overlay', 'horizontal', 'vertical', 'square', 'double-vertical', 'double-horizontal')
        time_types = ('year', 'month', 'day', 'weekday', 'hour', 'minute', 'second')
        time_locators = (YearLocator, MonthLocator, DayLocator, WeekdayLocator, HourLocator, MinuteLocator, SecondLocator)
        _available = {'layout' : layouts, 'datetime locator' : dict(zip(time_types, time_locators))}
        self._available.update(_available)

    @staticmethod
    def get_vecdir(layout):
        """
        layout:
            type <str>
        """
        if 'horizontal' in layout:
            vecdir = 'column'
        else:
            vecdir = 'row'
        return vecdir

    @staticmethod
    def get_diagonal_table_colors(facecolors, nrows, ncols):
        """
        facecolors:
            type <tuple / list / array>

        nrows:
            type <int>

        ncols:
            type <int>
        """
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors != 2:
            raise ValueError("2 facecolors required but only {} provided".format(ncolors))
        cellColours = []
        for ith_row in range(nrows):
            for ith_col in range(ncols):
                if ith_row % 2 == 0:
                    if ith_col % 2 == 0:
                        cellColours.append(facecolors[0])
                    else:
                        cellColours.append(facecolors[1])
                else:
                    if ith_col % 2 == 0:
                        cellColours.append(facecolors[1])
                    else:
                        cellColours.append(facecolors[0])
        return np.array(cellColours).reshape((nrows, ncols))

    @staticmethod
    def get_row_column_numbers(layout, n):
        """
        layout:
            type <str>

        n:
            type <int>
        """
        if not isinstance(n, int):
            raise ValueError("invalid type(n): {}".format(type(n)))
        if n < 1:
            raise ValueError("n must be positive")
        if layout == 'overlay':
            nrows, ncols = 1, 1
        elif layout == 'horizontal':
            nrows, ncols = 1, n
        elif layout == 'vertical':
            nrows, ncols = n, 1
        elif layout == 'double-vertical':
            if n % 2 != 0:
                raise ValueError("cannot place rows evenly by 2 columns given odd number of subplots")
            nrows, ncols = n//2, 2
        elif layout == 'double-horizontal':
            if n % 2 != 0:
                raise ValueError("cannot place columns evenly by 2 rows given odd number of subplots")
            nrows, ncols = 2, n//2
        elif layout == 'triple-vertical':
            if n % 3 != 0:
                raise ValueError("cannot place rows evenly by 3 columns given odd number of subplots")
            nrows, ncols = n//3, 3
        elif layout == 'triple-horizontal':
            if n % 3 != 0:
                raise ValueError("cannot place columns evenly by 3 rows given odd number of subplots")
            nrows, ncols = 3, n//3
        elif layout == 'square':
            nrows = int(np.sqrt(n))
            if nrows * nrows == n:
                ncols = nrows
            else:
                raise ValueError("{} axes could not be arranged in equal numbers of rows and columns of equal spacing".format(n))
        else:
            raise ValueError("invalid layout: {}".format(layout))
        return nrows, ncols

    def get_dim2_figure_and_axes(self, layout, n, **kwargs):
        """
        layout:
            type <str>

        n:
            type <int>
        """
        nrows, ncols = self.get_row_column_numbers(layout, n)
        fig, ax = plt.subplots(nrows, ncols, **kwargs)
        return fig, ax

    def get_dim3_figure_and_axes(self, layout, n, **kwargs):
        """

        """
        fig = plt.figure()
        nrows, ncols = self.get_row_column_numbers(layout, n)
        rows, cols = list(range(nrows)), list(range(ncols))
        axes = []
        ith_sub = 0
        for ith_row, row in enumerate(rows):
            for ith_col, col in enumerate(cols):
                ith_sub += 1
                coord = (nrows, ncols, ith_sub)
                ax = fig.add_subplot(*coord, projection='3d')
                axes.append(ax)
        return fig, np.array(ax)

    # def get_dim3_figure_and_axes(self, layout, n, **kwargs):
    #     """
    #
    #     """
    #     fig = plt.figure()
    #     nrows, ncols = self.get_row_column_numbers(layout, n)
    #     axes = []
    #     _nrows, _ncols = list(range(nrows)), list(range(ncols))
    #     for _i, _row in enumerate(_nrows):
    #         i = _i + 1
    #         row = _row + 1
    #         for _j, _col in enumerate(_ncols):
    #             if _j != 0:
    #                 i += i + 1
    #             col = _col + 1
    #             coord = (row, col, i)
    #             print("\n .. coord:\t{}".format(coord))
    #             ax = fig.add_subplot(*coord, projection='3d')
    #             axes.append(ax)
    #     return fig, np.array(ax)

    def get_mirror_ax(self, ax, basex=None, basey=None, xscale='log', yscale='log', ax_color='gray'):
        """
        ax:
            type <matplotlib object>

        basex:
            type <int / float> or None

        basey:
            type <int / float> or None

        xscale:
            type <str> or None

        yscale:
            type <str> or None

        ax_color:
            type <str>
        """
        # newax = ax._make_twin_axes(frameon=False)
        newax = ax.figure.add_subplot(ax.get_subplotspec(), frameon=False)
        newax.xaxis.set(label_position='top')
        newax.yaxis.set(label_position='right', offset_position='right')
        newax.yaxis.get_label().set_rotation(-90)
        newax.yaxis.tick_right()
        newax.xaxis.tick_top()
        if basex is not None:
            newax.set_xscale(xscale, basex=basex)
            newax.xaxis.set_major_locator(ticker.LogLocator(base=2))
            newax.yaxis.set_minor_formatter(ticker.NullFormatter())
            newax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        if basey is not None:
            newax.set_yscale(yscale, basey=basey)
            newax.yaxis.set_major_locator(ticker.LogLocator(base=2))
            newax.yaxis.set_minor_formatter(ticker.NullFormatter())
            newax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        newax.xaxis.label.set_color(ax_color)
        newax.yaxis.label.set_color(ax_color)
        newax.tick_params(axis='both', which='both', labelsize=self.ticksize, colors=ax_color)
        return newax

    def transform_x_as_datetime(self, ax, major='year', minor='month', fmt="%Y-%m-%d", rotation=None):
        """
        ax:
            type <matplotlib object>

        major:
            type <str>

        minor:
            type <str>

        fmt:
            type <str>

        rotation:
            type <int / float> or None
        """
        maj_locator = self.available['datetime locator'][major]
        min_locator = self.available['datetime locator'][minor]
        ax.xaxis.set_major_locator(maj_locator())
        ax.xaxis.set_minor_locator(min_locator())
        ax.xaxis.set_major_formatter(DateFormatter(fmt))
        if rotation is not None:
            ax.tick_params(axis='x', rotation=rotation)
        return ax

    def add_shared_labels(self, fig, axes, xlabel, ylabel, vecdir='row'):
        """
        fig:
            type <matplotlib object>

        axes:
            type <array>

        xlabel:
            type <str>

        ylabel:
            type <str>

        vecdir:
            type <str>
        """
        try:
            try:
                for ax in axes[-1, :].ravel():
                    ax.set_xlabel(xlabel, fontsize=self.labelsize)
                for ax in axes[:, 0].ravel():
                    ax.set_ylabel(ylabel, fontsize=self.labelsize)
            except:
                if vecdir == 'row':
                    axes[-1].set_xlabel(xlabel, fontsize=self.labelsize)
                    for ax in axes.ravel():
                        ax.set_ylabel(ylabel, fontsize=self.labelsize)
                elif vecdir == 'column':
                    axes[0].set_ylabel(ylabel, fontsize=self.labelsize)
                    for ax in axes.ravel():
                        ax.set_xlabel(xlabel, fontsize=self.labelsize)
                else:
                    raise ValueError("invalid vecdir: {}".format(vecdir))
            fig.align_ylabels()
        except:
            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)
        return fig, axes

    def subview_column_data(self, header, colLabels, cellText, facecolors, colColours, loc='center', cellLoc='center', **kwargs):
        """
        header:
            type <str>

        colLabels:
            type <list>

        cellText:
            type <list>

        facecolors:
            type <tuple / list / array>

        colColours:
            type <str>

        loc:
            type <str>

        cellLoc:
            type <str>
        """
        cell_size = len(cellText)
        ncols = len(colLabels)
        nrows = cell_size // ncols
        cellText = np.array(cellText).reshape((nrows, ncols))
        cellColours = self.get_diagonal_table_colors(facecolors, nrows, ncols)
        colColours = [colColours for col in range(ncols)]
        fig, axes = self.get_dim2_figure_and_axes(layout='overlay', n=1, **kwargs)
        table = axes.table(colLabels=colLabels, cellText=cellText, colColours=colColours, cellColours=cellColours, loc=loc, cellLoc=cellLoc, **kwargs)
        axes.axis('off')
        axes.axis('tight')
        table.auto_set_font_size(False)
        table.set_fontsize(self.cellsize)
        for key, cell in table.get_celld().items():
            row, col = key
            if row == 0:
                cell.set_text_props(fontproperties=FontProperties(variant='small-caps', size=self.headersize)) # weight='semibold'
        xscale, yscale = (ncols, nrows)
        table.scale(xscale, yscale)
        return fig, axes, table

    def subview_tabular_data(self, header, colLabels, rowLabels, cellText, facecolors, colColours, rowColours, loc='center', cellLoc='center', **kwargs):
        """
        header:
            type <str>

        colLabels:
            type <list>

        rowLabels:
            type <list>

        cellText:
            type <list>

        facecolors:
            type <tuple / list / array>

        colColours:
            type <str>

        rowColours:
            type <str>

        loc:
            type <str>

        cellLoc:
            type <str>
        """
        cell_size = len(cellText)
        nrows, ncols = len(rowLabels), len(colLabels)
        cellText = np.array(cellText).reshape((ncols, nrows)).T
        cellColours = self.get_diagonal_table_colors(facecolors, nrows, ncols)
        colColours = [colColours for col in range(ncols)]
        rowColours = [rowColours for row in range(nrows)]
        fig, axes = self.get_dim2_figure_and_axes(layout='overlay', n=1, **kwargs)
        table = axes.table(colLabels=colLabels, rowLabels=rowLabels, cellText=cellText, colColours=colColours, rowColours=rowColours, cellColours=cellColours, loc=loc, cellLoc=cellLoc, **kwargs)
        axes.axis('off')
        axes.axis('tight')
        table.auto_set_font_size(False)
        table.set_fontsize(self.cellsize)
        for key, cell in table.get_celld().items():
            row, col = key
            if row == 0:
                cell.set_text_props(fontproperties=FontProperties(variant='small-caps', size=self.headersize)) # weight='semibold'
        xscale, yscale = (ncols, nrows)
        table.scale(xscale, yscale)
        return fig, axes, table

    def display_image(self, fig, savename=None, dpi=800, bbox_inches='tight', pad_inches=0.1, extension='.png', **kwargs):
        """
        fig:
            type <matplotlib object>

        savename:
            type <str> or None

        dpi:
            type <int>

        bbox_inches:
            type <str> or None

        pad_inches:
            type <float>

        extension:
            type <str>
        """
        if savename is None:
            plt.show()
        elif isinstance(savename, str):
            savepath = '{}{}{}'.format(self.directory, savename, extension)
            fig.savefig(savepath, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches, **kwargs)
        else:
            raise ValueError("invalid type(savename): {}".format(type(savename)))
        plt.close(fig)

    def convert_to_eps(self, extension='.png'):
        """ """
        n = len(extension)
        for dirpath, dirnames, filenames in os.walk(self.directory):
            for filename in [f for f in filenames if f.endswith(extension)]:
                img_path = os.path.join(dirpath, filename)
                save_path = '{}.eps'.format(img_path[:-n])
                img = Image.open(img_path)
                img.load()
                res = Image.new('RGB', img.size, (255, 255, 255))
                res.paste(img, mask=img.split()[3])
                res.save(save_path, 'EPS')
##
