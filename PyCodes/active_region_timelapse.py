import os
import numpy as np
import aplpy
import moviepy.editor as mpy

class MichelsonDopplerImager():

    def __init__(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, rdirectory='', sdirectory=''):
        """
        ticksize:
            type <int / float>

        labelsize:
            type <int / float>

        textsize:
            type <int / float>

        titlesize:
            type <int / float>

        rdirectory:
            type <str>

        sdirectory:
            type <str>
        """
        self.ticksize = ticksize
        self.labelsize = labelsize
        self.textsize = textsize
        self.titlesize = titlesize
        self.rdirectory = rdirectory
        self.sdirectory = sdirectory
        self.titlesize = titlesize
        self.vmin = -60
        self.vmax = 60
        self.extension = '.fits'

    @staticmethod
    def search_directory(directory, extension):
        """
        directory:
            type <str>

        extension:
            type <str>
        """
        paths = ['{}{}'.format(directory, filename) for filename in os.listdir(directory) if filename.endswith(extension)]
        return sorted(paths)

    @staticmethod
    def get_filename_suffix(cmap):
        """
        cmap:
            type <str> or None
        """
        return '__{}'.format(cmap)

    @staticmethod
    def strip_filename(path):
        return path.rsplit('/', 1)[-1]

    def configure_carrington_map(self, fig, ax):
        """
        fig:
            type <matplotlib / apl object>

        ax:
            type <matplotlib / apl object>
        """
        xmin = 60
        xmax = 360
        xn = 60
        xscale = 360.
        xtext = 'Carrington Longitude'
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xticklabels = np.arange(xmin, xmax+1, xn)
        xticks = xticklabels / xscale * xlim[1] + xlim[0]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=self.ticksize)
        fig.axis_labels.set_xtext(xtext)
        ctime = fig._wcs.wcs_pix2world(xlim, ylim, origin=1)[0]
        title = r'Carrington Time: ${}$ - ${}$'.format(ctime[1], ctime[0])
        ax.set_title(title, fontsize=self.titlesize)
        return fig, ax

    def apply_color_transform(self, fig, cmap=None):
        """
        fig:
            type <matplotlib / apl object>

        cmap:
            type <str> or None
        """
        if cmap is None:
            fig.show_grayscale(vmin=self.vmin, vmax=self.vmax)
        else:
            fig.show_colorscale(vmin=self.vmin, vmax=self.vmax, cmap=cmap)
        return fig

    def save_image(self, path, cmap=None, extension='.png', dpi=800):
        """
        path:
            type <str>

        cmap:
            type <str> or None

        extension:
            type <str>

        dpi:
            type <int>
        """
        path_without_ext = path[:-len(self.extension)]
        filename = self.strip_filename(path_without_ext)
        suffix = self.get_filename_suffix(cmap)
        savepath = '{}{}{}{}'.format(self.sdirectory, filename, suffix, extension)
        if savepath not in self.search_directory(self.sdirectory, extension):
            fig = aplpy.FITSFigure(path)
            ax = fig._ax1
            fig, ax = self.configure_carrington_map(fig, ax)
            fig = self.apply_color_transform(fig, cmap)
            fig.savefig(savepath, dpi=dpi, format=None, transparent=False)

    def convert_from_fits(self, cmap=None, extension='.png', dpi=800):
        """
        cmap:
            type <str> or None

        extension:
            type <str>

        dpi:
            type <int>
        """
        for path in self.search_directory(self.rdirectory, '.fits'):
            self.save_image(path, cmap, extension, dpi)

    def save_timelapse(self, fps, cmap=None, extension='.mkv', codec='mpeg4', search_extension='.png'):
        """
        fps:
            type <int / float>

        cmap:
            type <str> or None

        extension:
            type <str>

        codec:
            type <str>

        search_extension:
            type <str>
        """
        suffix = self.get_filename_suffix(cmap)
        paths = [path for path in self.search_directory(self.sdirectory, search_extension) if suffix in path]
        n = len(paths)
        savepath = '{}_CarringtonRotations{}_fps{}{}'.format(self.sdirectory, suffix, fps, extension)
        clip = mpy.ImageSequenceClip(paths, fps=fps)
        clip.write_videofile(savepath, fps=fps, codec=codec, audio=False)

##
