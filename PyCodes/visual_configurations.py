import os
from PIL import Image
from copy import deepcopy
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, YearLocator, MonthLocator, DateFormatter
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties
from matplotlib.legend_handler import HandlerTuple, HandlerRegularPolyCollection
from matplotlib.ticker import FuncFormatter, ScalarFormatter, LogFormatter, NullFormatter, MultipleLocator, AutoLocator
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

class VisualEditor():

    @staticmethod
    def get_normalized_colormap(vmin, vmax, estimator=None, autorescale=False):
        """ """
        if estimator is None:
            norm = Normalize(vmin=vmin, vmax=vmax)
        else:
            if estimator in ('maximum likelihood',):
                if autorescale == True:
                    norm = Normalize(vmin=vmin /1000, vmax=0)
                else:
                    norm = Normalize(vmin=vmin, vmax=0)
            else:
                if autorescale == True:
                    norm = Normalize(vmin=0, vmax=vmax /1000)
                else:
                    norm = Normalize(vmin=0, vmax=vmax)
        return norm

    @staticmethod
    def set_tick_size(ax, xtick_size, ytick_size, ztick_size=None):
        """ """
        ax.tick_params(axis='x', labelsize=xtick_size)
        ax.tick_params(axis='y', labelsize=ytick_size)
        if ztick_size is not None:
            ax.tick_params(axis='z', labelsize=ztick_size)
        return ax

    @staticmethod
    def set_label_size(ax, xlabel_size, ylabel_size, zlabel_size=None):
        """ """
        ax.xaxis.label.set_size(ylabel_size)
        ax.yaxis.label.set_size(xlabel_size)
        if zlabel_size is not None:
            ax.zaxis.label.set_size(zlabel_size)
        return ax

    @staticmethod
    def round_up_to(value, base):
        """ """
        return (value + base - 1) // base * base

    @staticmethod
    def round_down_to(value, base):
        """ """
        return int(np.floor(value / base) * base)

    @staticmethod
    def transform_x_as_datetime(ax):
        """ """
        ax.xaxis.set_major_locator(YearLocator())
        ax.xaxis.set_minor_locator(MonthLocator())
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d"))
        ax.tick_params(axis='x', rotation=15)
        return ax

    @staticmethod
    def get_discrete_colors_from_colormap(cmap, args):
        """ """
        cmap = plt.cm.get_cmap(cmap)
        if isinstance(args, int):
            facecolors = [cmap(idx/args) for idx in range(args)]
        return facecolors

class FigureOptions(VisualEditor):

    def __init__(self, saveloc):
        """ """
        super().__init__()
        self.saveloc = saveloc

    def get_savepath(self, savename, extension):
        """ """
        if self.saveloc is None:
            return None
        else:
            return '{}{}{}'.format(self.saveloc, savename, extension)

    @staticmethod
    def display_figure(fig, savepath=None, dpi=800, bbox_inches='tight', pad_inches=0.1):
        """
        fig:
            *   type <matplotlib object>
        savepath:
            *   type <str> or None
        dpi:
            *   type <int>
        bbox_inches:
            *   type <str>
        pad_inches:
            *   type <float>
        """
        if savepath is None:
            plt.show()
        elif isinstance(savepath, str):
            fig.savefig(savepath, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches)
        else:
            raise ValueError("unknown type: {}; savepath = type <str> or None")
        plt.close(fig)

    @staticmethod
    def display_text(string, savepath=None):
        """
        string:
            *   type <str>
        savepath:
            *   type <str> or None
        """
        if savepath is None:
            print(string)
        else:
            if isinstance(string, str):
                np.savetxt(savepath, [string], fmt='%s')
            else:
                np.savetxt(savepath, string, fmt='%s')

    @staticmethod
    def get_rectangular_handle(facecolor='none', edgecolor='none', visible=False, fill=False, xyi=(0,0), xyf=(0.1, 0.1)):
        """ """
        return Rectangle(xyi, *xyf, facecolor=facecolor, edgecolor=edgecolor, visible=visible, fill=fill)

    @staticmethod
    def get_overlay_legend(fig, ax, n, bottom=0.215, textsize=8):
        """ """
        fig.subplots_adjust(bottom=bottom)
        if n == 1:
            handles, labels = ax.get_legend_handles_labels()
            empty_handle = Rectangle((0, 0), 1, 1, alpha=0)
            handles = [empty_handle] + handles + [empty_handle]
            labels = [' '] + labels + [' ']
            fig.legend(handles=handles, labels=labels, loc='lower center', fontsize=textsize, mode='expand', ncol=3)
        else:
            fig.legend(loc='lower center', fontsize=textsize, mode='expand', ncol=n)
        return fig, ax

    @staticmethod
    def configure_symmetrical_axes(nrows, ncols, axes, apply_labels=False, apply_ticklabels=False, autoconfigure=False, tops_and_bottoms=False, modify_columns=False, modify_singular=False):
        """ """
        if nrows > 1:
            try:
                for ax in axes[:-1, :].ravel():
                    if apply_labels == True:
                        ax.set_xlabel('')
                    if apply_ticklabels == True:
                        ax.set_xticklabels([])
            except:
                pass
        if ncols > 1:
            if nrows == 1:
                if apply_labels == True:
                    axes[-1].yaxis.set_label_position('right')
                if apply_ticklabels == True:
                    axes[-1].yaxis.tick_right()
            else:
                for ax in axes[:, -1].ravel():
                    if apply_labels == True:
                        ax.yaxis.set_label_position('right')
                    if apply_ticklabels == True:
                        ax.yaxis.tick_right()
                if ncols > 2:
                    for ax in axes[:, 1:-1].ravel():
                        if apply_labels == True:
                            ax.set_ylabel('')
                        if apply_ticklabels == True:
                            ax.set_yticklabels([])
        if autoconfigure == True:
            if nrows * ncols > 2:
                for ax in axes[0, 1:-1].ravel():
                    ax.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=False)
                    ax.set_ylabel('')
                for ax in axes[1:-1, 1:-1].ravel():
                    if tops_and_bottoms == True:
                        ax.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=False)
                    else:
                        ax.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=False)
                    ax.set_ylabel('')
                if tops_and_bottoms == True:
                    for ax in axes[1:-1, 0].ravel():
                        ax.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=False, left=True, right=True, labelleft=True, labelright=False)
                    for ax in axes[1:-1, -1].ravel():
                        ax.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=True)
                else:
                    for ax in axes[1:-1, 0].ravel():
                        ax.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=True, labelright=False)
                    for ax in axes[1:-1, -1].ravel():
                        ax.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=True)
                for ax in axes[-1, 1:-1].ravel():
                    ax.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=True, right=True, labelleft=False, labelright=False)
                    ax.set_ylabel('')
                if ncols > 2:
                    if tops_and_bottoms == True:
                        pass
                    else:
                        for ax in axes[1:, 1:-1].ravel():
                            ax.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=False)
                axes[0, 0].tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=True, labelright=False)
                axes[0, -1].tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=True)
                if tops_and_bottoms == True:
                    axes[-1, 0].tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=True, right=True, labelleft=True, labelright=False)
                    axes[-1, -1].tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=True, right=True, labelleft=False, labelright=True)
                else:
                    axes[-1, 0].tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=True, left=True, right=True, labelleft=True, labelright=False)
                    axes[-1, -1].tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=True, left=True, right=True, labelleft=False, labelright=True)
        if ((nrows == 1) and (ncols == 2)):
            axes[-1].set_ylabel('')
        if modify_columns == True:
            if ((nrows == 1) and (ncols == 2)):
                (ax_left, ax_right) = axes
                ax_left.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=True, left=True, right=True, labelleft=True, labelright=False)
                ax_right.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=True, left=True, right=True, labelleft=False, labelright=True)
        if modify_singular == True:
            if ((nrows == 1) and (ncols == 1)):
                axes[0].tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=True, right=True, labelleft=True, labelright=False)
        return axes

    def convert_to_eps(self, extension='.png'):
        """ """
        n = len(extension)
        for dirpath, dirnames, filenames in os.walk(self.saveloc):
            for filename in [f for f in filenames if f.endswith(extension)]:
                img_path = os.path.join(dirpath, filename)
                save_path = '{}.eps'.format(img_path[:-n])
                img = Image.open(img_path)
                img.load()
                res = Image.new('RGB', img.size, (255, 255, 255))
                res.paste(img, mask=img.split()[3])
                res.save(save_path, 'EPS')
