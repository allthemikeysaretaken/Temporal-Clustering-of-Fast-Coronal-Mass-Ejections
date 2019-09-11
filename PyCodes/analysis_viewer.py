import aplpy
import moviepy.editor as mpy

from meta_analysis_methods import *

class FitsInterpreter():

    def __init__(self, fitsloc, saveloc=None):
        self.fitsloc = fitsloc
        if saveloc is None:
            self.saveloc = fitsloc
        else:
            self.saveloc = saveloc
        self.vmin = -60
        self.vmax = 60
        self.img_paths = []

    @staticmethod
    def get_savepath(savename, extension):
        """ """
        return '{}{}'.format(savename, extension)

    @staticmethod
    def get_filepaths(loc, extension=None):
        """ """
        if extension is None:
            filepaths = ['{}{}'.format(loc, filename) for filename in os.listdir(loc)]
        else:
            filepaths = ['{}{}'.format(loc, filename) for filename in os.listdir(loc) if filename.endswith(extension)]
        return filepaths

    @staticmethod
    def sort_filepaths(paths, sort_by='filename'):
        """ """
        if sort_by == 'filename':
            return sorted(paths)
        elif sort_by == 'extension':
            return paths.sort(key=lambda f: os.path.splitext(f)[1])
        else:
            raise ValueError("invalid sort_by: {}; available sort_by = 'filename' or 'extension'")

    @staticmethod
    def configure_carrington_map(fig, ax):
        """ """
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
        ax.set_xticklabels(xticklabels)
        fig.axis_labels.set_xtext(xtext)
        ctime = fig._wcs.wcs_pix2world(xlim, ylim, origin=1)[0]
        ax.set_title(r'Carrington Time: ${}$ - ${}$'.format(ctime[0], ctime[1]))
        return fig, ax

    def apply_color_transform(self, fig, cmap=None):
        """ """
        if cmap is None: ## grayscale
            fig.show_grayscale(vmin=self.vmin, vmax=self.vmax)
        else: ## color-mapped
            fig.show_colorscale(vmin=self.vmin, vmax=self.vmax, cmap=cmap)
        return fig

    def to_image(self, fitspath, savename, extension='.png', apply_false_color=False, apply_grayscale=False, dpi=800):
        """ """
        fig = aplpy.FITSFigure(fitspath)
        ax = fig._ax1
        fig, ax = self.configure_carrington_map(fig, ax)
        if apply_grayscale == True:
            fig = self.apply_color_transform(fig, cmap=None)
            savepath = self.get_savepath('{}__grayscale'.format(savename), extension)
            fig.savefig(savepath, dpi=dpi, format=None, transparent=False)
            self.img_paths.append(savepath)
        if apply_false_color == True:
            fig = self.apply_color_transform(fig, cmap='hot')
            savepath = self.get_savepath('{}__falsecolor'.format(savename), extension)
            fig.savefig(savepath, dpi=dpi, format=None, transparent=False)
            self.img_paths.append(savepath)
        plt.close()

    def save_images(self, extension='.png', apply_false_color=False, apply_grayscale=False, dpi=800):
        """ """
        n = len('.fits')
        for fitspath in self.get_filepaths(self.fitsloc, extension='.fits'):
            self.to_image(fitspath, savename=fitspath[:-n], extension=extension, apply_false_color=apply_false_color, apply_grayscale=apply_grayscale)

    def save_timelapse_animation(self, fps, savename, extension='.mkv', codec='mpeg4', include_false_color=False, include_grayscale=False):
        """ """
        conditional_paths = []
        if include_false_color == True:
            conditional_paths += [path for path in self.img_paths if '__falsecolor' in path]
        if include_grayscale == True:
            conditional_paths += [path for path in self.img_paths if '__grayscale' in path]
        if len(conditional_paths) == 0:
            raise ValueError("both inputs are False; set include_false_color = True and/or include_grayscale = True")
        img_paths = self.sort_filepaths(conditional_paths, sort_by='filename')
        savepath = '{}{}{}'.format(self.saveloc, savename, extension)
        clip = mpy.ImageSequenceClip(img_paths, fps=fps)
        clip.write_videofile(savepath, fps=fps, codec=codec, audio=False)

class CorrelatedSunSpotViewer(FigureOptions):

    def __init__(self, saveloc=None):
        super().__init__(saveloc)

    def view_cme_sunspot_frequency_comparison(self, sunspot_datas, cme_datas, desired_timescale='daily', denote_cycles=None, ticksize=7, labelsize=8, textsize=8, titlesize=10, add_legend=False, **kwargs):
        """ """
        if not isinstance(sunspot_datas, (tuple, list, np.ndarray)):
            sunspot_datas = [sunspot_datas]
        if not isinstance(cme_datas, (tuple, list, np.ndarray)):
            sunspot_datas = [cme_datas]
        n1 = len(sunspot_datas)
        n2 = len(cme_datas)
        if n1 != n2:
            raise ValueError("{} sunspot_datas with {} cme_datas".format(n1, n2))
        legend_labels, legend_handles = [], []
        fig, axes = plt.subplots(nrows=4, ncols=1, **kwargs)
        for sunspot_data, cme_data in zip(sunspot_datas, cme_datas):
            EFC = EventFrequencyComparison(sunspot_data, cme_data['data']['original'])
            EFC.initialize_converter()
            ax_top, handle_top, label_top = EFC.subview_sunspot_data(axes[0], desired_timescale=desired_timescale, ticksize=ticksize, labelsize=labelsize, textsize=textsize)
            legend_labels.append(label_top)
            legend_handles.append(handle_top)
            if desired_timescale == 'yearly':
                y1_ticks, y1_space = np.arange(0, 3001, 250), 4
                y2_ticks, y2_space = np.arange(0, 401, 25), 4
                y3_ticks, y3_space = np.arange(0, 201, 25), 2
            elif desired_timescale == 'monthly':
                y1_ticks, y1_space = np.arange(0, 301, 25), 4
                y2_ticks, y2_space = np.arange(0, 51, 1), 10
                y3_ticks, y3_space = np.arange(0, 31, 1), 5
            elif desired_timescale == 'daily':
                y1_ticks, y1_space = np.arange(0, 26, 1), 5
                y2_ticks, y2_space = np.arange(0, 16, 1), 5
                y3_ticks, y3_space = np.arange(0, 7, 1), 2
            else:
                raise ValueError("invalid desired_timescale: {}; available timescales = 'daily', 'monthly', or 'yearly'".format(desired_timescale))
            ax_umid, handle_umid, label_umid = EFC.subview_cme_data(axes[1], y1_ticks, y1_space, search_parameters='speed', search_conditions='greater than', search_values=0, apply_to='all', desired_timescale=desired_timescale, ticksize=ticksize, labelsize=labelsize, textsize=textsize)
            ax_bmid, handle_bmid, label_bmid = EFC.subview_cme_data(axes[2], y2_ticks, y2_space, search_parameters='speed', search_conditions='greater than', search_values=800, apply_to='all', desired_timescale=desired_timescale, ticksize=ticksize, labelsize=labelsize, textsize=textsize)
            ax_btm, handle_btm, label_btm = EFC.subview_cme_data(axes[3], y3_ticks, y3_space, search_parameters='speed', search_conditions='greater than', search_values=1000, apply_to='all', desired_timescale=desired_timescale, ticksize=ticksize, labelsize=labelsize, textsize=textsize)
            legend_labels.extend([label_umid, label_bmid, label_btm])
            legend_handles.extend([handle_umid, handle_bmid, handle_btm])
        try:
            legend_handles = [handle[0] for handle in legend_handles]
        except:
            pass
        for ax in axes[:-1].ravel():
            ax.set_xticklabels([])
            ax.set_xlabel('')
        for ax in axes[1:-1].ravel():
            ax.tick_params(which='both', top=True, bottom=True, labeltop=False, labelbottom=False, left=True, right=True, labelleft=True, labelright=False)
        axes[0].tick_params(which='both', top=False, bottom=True, labeltop=False, labelbottom=False, left=True, right=True, labelleft=True, labelright=False)
        axes[-1].tick_params(which='both', top=True, bottom=True, labeltop=False, labelbottom=True, left=True, right=True, labelleft=True, labelright=False)
        if denote_cycles is not None:
            if not isinstance(denote_cycles, (tuple, list, np.ndarray)):
                raise ValueError("denote_cycles = type <tuple / list / array>")
            if not isinstance(denote_cycles[0], (tuple, list, np.ndarray)):
                denote_cycles = [denote_cycles]
            for pair_cycles in denote_cycles:
                if pair_cycles in ((23, 24), ('23, 24')):
                    SC = SolarCycle()
                    sc23_dts = SC.sc_23[None]['search_values']
                    sc24_dts = SC.sc_24[None]['search_values']
                    middle = (sc24_dts[0] - sc23_dts[-1]) /2 + sc23_dts[-1]
                    for ith_row, ax in enumerate(axes.ravel()):
                        xdt = date2num(middle)
                        ylim = ax.get_ylim()
                        ax.axvline(xdt, ymin=ylim[0], ymax=ylim[-1], color='r', linestyle='--')
                        if ith_row == 0:
                            arrowprops = dict(facecolor='black', arrowstyle="->")
                            # arrowprops = dict(facecolor='black', arrowstyle="simple")
                            yi = 0.8 * ylim[-1]
                            yj = 0.6 * ylim[-1]
                            x23_left = date2num(sc23_dts[0])
                            x23_right = date2num(sc23_dts[-1])
                            x24_left = date2num(sc24_dts[0])
                            x24_right = date2num(sc24_dts[-1])
                            axes[0].annotate(' ', xy=(sc23_dts[0], yj), xytext=(sc23_dts[-1], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                            axes[0].annotate(' ', xy=(sc24_dts[-1], yj), xytext=(sc24_dts[0], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                            axes[0].text(date2num(sc23_dts[-1]) - 200, yi, 'SC ${}$'.format(denote_cycles[0][0]), horizontalalignment='right', verticalalignment='center')
                            axes[0].text(date2num(sc24_dts[0]) + 200, yi, 'SC ${}$'.format(denote_cycles[0][1]), horizontalalignment='left', verticalalignment='center')
                else:
                    raise ValueError("not yet implemented")
        xdt = date2num(np.concatenate(([sunspot_data['datetime object'].tolist() for sunspot_data in sunspot_datas] + [cme_data['data']['original']['datetime object'].tolist() for cme_data in cme_datas]), axis=0))
        for ax in axes.ravel():
            ax.set_xlim([np.min(xdt), np.max(xdt)])
        if add_legend == True:
            nh, nl = len(legend_handles), len(legend_labels)
            if denote_cycles is None:
                fig.subplots_adjust(bottom=0.2, hspace=0.3)
                fig.legend(handles=legend_handles, labels=legend_labels, ncol=2, loc='lower center', fontsize=textsize, mode='expand', fancybox=True)
            else:
                fig.subplots_adjust(bottom=0.225, hspace=0.325)
                leg_left = fig.legend(title=r'SC ${}$'.format(denote_cycles[0][0]), handles=legend_handles[:nh//2], labels=legend_labels[:nl//2], ncol=2, loc='lower left', fontsize=textsize)
                leg_right = fig.legend(title=r'SC ${}$'.format(denote_cycles[0][1]), handles=legend_handles[nh//2:], labels=legend_labels[nl//2:], ncol=2, loc='lower right', fontsize=textsize)
                for leg in (leg_left, leg_right):
                    leg._legend_box.align = "center"
                    frame = leg.get_frame()
                    # frame.set_facecolor('green')
                    frame.set_edgecolor('k')
        fig.align_ylabels()
        fig.suptitle('Comparison of Observed Frequencies\n of Sunspots & CMEs', fontsize=titlesize)
        savepath = self.get_savepath('sunspot_cme_frequency_comparison_{}'.format(desired_timescale), extension='.png')
        self.display_figure(fig, savepath)

class SpeedDistributionViewer(FigureOptions):

    def __init__(self, datas, saveloc=None):
        super().__init__(saveloc)
        if isinstance(datas, (tuple, list, np.ndarray)):
            self.datas = datas
        else:
            self.datas = [datas]
        self.n = len(self.datas)
        self._SDs = []

    @property
    def SDs(self):
        return self._SDs

    def initialize_optimizations(self, n=None, w=None, edges=None, bin_threshold=5, bin_estimator='reduced gtest', alt_estimator='maximum likelihood', search_parameters='speed', search_conditions='greater than', search_values=20, apply_to='all', bias='left'):
        """ """
        for timeseries in self.datas:
            identifier = timeseries['identifier']
            event_data = timeseries['data']['original']
            Searcher = SearchEvents(event_data, string_keys=('',))
            thresholded_data = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            SD = SpeedDistribution(thresholded_data['speed'], identifier)
            SD.initialize_histograms(n, w, edges, bin_threshold, bias)
            SD.initialize_binned_optimization(initial_parameter_guess=None, estimator=bin_estimator)
            SD.initialize_alternate_optimization(initial_parameter_guess=None, estimator=alt_estimator)
            for estimator in (bin_estimator, alt_estimator):
                optim_result = SD.optimizations[estimator]
                optimizer, prms = optim_result['cls'], optim_result['optimized parameters']
                # optimizer.initialize_dim2_space(prms, estimator, x=None, y=None, xfrac=0.5, xn=27, xspace_by='number', yfrac=0.5, yn=27, yspace_by='number')
                # optimizer.initialize_dim2_space(prms, estimator, x=np.arange(3, 8.01, 0.5), y=np.arange(0.2, 1.01, 0.05))
                optimizer.initialize_dim2_space(prms, estimator, x=np.arange(3, 10.01, 0.1), y=np.arange(0.2, 1.21, 0.01))
                # optimizer.initialize_dim2_space(prms, estimator, x=np.arange(3, 8.01, 0.1), y=np.arange(0.2, 1.01, 0.01))
                optimizer.initialize_dim3_space(prms, estimator, f=optimizer.f_err, args=())
                SD.optimizations[estimator]['cls'] = optimizer
            self._SDs.append(SD)

    def view_tail_histogram(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for idx, (ax, timeseries) in enumerate(zip(axes.ravel(), self.datas)):
            SD = self.SDs[idx]
            ax = SD.subview_tail_histogram(ax, ticksize, labelsize, textsize, titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        if self.n == 1:
            fig.suptitle('Tail of CME Speed Distribution', fontsize=titlesize)
        else:
            fig.suptitle('CME Speed Distributions:\n Tail Comparison', fontsize=titlesize)
        fig, _ = self.get_overlay_legend(fig, axes.ravel()[0], n=1, bottom=0.175, textsize=textsize)
        savepath = self.get_savepath('SDV_speed_tail_histogram', extension='.png')
        self.display_figure(fig, savepath)

    def view_probability_density_fits(self, estimator, ticksize=7, labelsize=8, textsize=8, titlesize=10, show_confidence_interval=False, highlight_tail=False, point_to_tail=False, show_statistics=False, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
            wspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
                wspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                wspace = 0.3
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for idx, (ax, timeseries) in enumerate(zip(axes.ravel(), self.datas)):
            SD = self.SDs[idx]
            ax, handles, labels = SD.subview_probability_density_fitted_histogram(ax, estimator, ticksize, labelsize, textsize, titlesize, show_confidence_interval=show_confidence_interval, highlight_tail=highlight_tail, point_to_tail=point_to_tail, show_statistics=show_statistics)
        fig.subplots_adjust(bottom=0.2, wspace=wspace, hspace=hspace)
        leg = fig.legend(handles=handles, labels=labels, ncol=3, loc='lower center', fontsize=textsize, mode='expand')
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        fig.suptitle('CME Speed Distributions:\n Probability Density Fits', fontsize=titlesize)
        savepath = self.get_savepath('SDV_speed_fits_pdensity__{}'.format(estimator.replace(' ', '_')), extension='.png')
        self.display_figure(fig, savepath)

    def view_empirical_fits(self, estimator, ticksize=7, labelsize=8, textsize=8, titlesize=10, show_confidence_interval=False, highlight_tail=False, point_to_tail=False, show_statistics=False, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
            wspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
                wspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                wspace = 0.3
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for idx, (ax, timeseries) in enumerate(zip(axes.ravel(), self.datas)):
            SD = self.SDs[idx]
            ax, handles, labels = SD.subview_empirically_fitted_histogram(ax, estimator, ticksize, labelsize, textsize, titlesize, show_confidence_interval=show_confidence_interval, highlight_tail=highlight_tail, point_to_tail=point_to_tail, show_statistics=show_statistics)
        fig.subplots_adjust(bottom=0.2, wspace=wspace, hspace=hspace)
        leg = fig.legend(handles=handles, labels=labels, ncol=len(labels), loc='lower center', fontsize=textsize, mode='expand')
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        fig.suptitle('CME Speed Distributions:\n Empirical Fits', fontsize=titlesize)
        savepath = self.get_savepath('SDV_speed_fits_empirical__{}'.format(estimator.replace(' ', '_')), extension='.png')
        self.display_figure(fig, savepath)

    def view_fit_comparisons(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, highlight_tail=False, point_to_tail=False, counts_type='observed frequency', **kwargs):
        """ """
        if counts_type == 'observed frequency':
            title_mod = 'Empirical'
        elif counts_type == 'probability density':
            title_mod = 'Probabiliy Density'
        else:
            raise ValueError("invalid counts_type: {}; available counts_type: 'observed frequency' or 'probability density'")
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
            wspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
                wspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                wspace = 0.3
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for idx, (ax, timeseries) in enumerate(zip(axes.ravel(), self.datas)):
            SD = self.SDs[idx]
            ax = SD.subview_fit_comparisons(ax, ticksize, labelsize, textsize, titlesize, highlight_tail=highlight_tail, point_to_tail=point_to_tail, counts_type=counts_type)
        fig.subplots_adjust(bottom=0.2, wspace=wspace, hspace=hspace)
        handles, labels = ax.get_legend_handles_labels()
        leg = fig.legend(handles=handles, labels=labels, ncol=len(labels), loc='lower center', fontsize=textsize, mode='expand')
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        fig.suptitle('CME Speed Distributions:\n {} Fit Comparisons'.format(title_mod), fontsize=titlesize)
        savepath = self.get_savepath('SDV_speed_fits_comparison', extension='.png')
        self.display_figure(fig, savepath)

    def get_xyz_ticks_with_norm(self, estimator):
        """ """
        x, y, z = [], [], []
        for _SD in self.SDs:
            optimizer = _SD.optimizations[estimator]['cls']
            x.append(optimizer.dim2_space[estimator]['x'])
            y.append(optimizer.dim2_space[estimator]['y'])
            z.append(optimizer.dim3_space[estimator]['z'])
        get_pmin = lambda p : np.round(np.max([np.min(_p) for _p in p]), decimals=1)
        get_pmax = lambda p : np.round(np.min([np.max(_p) for _p in p]), decimals=1)
        xmin, xmax = int(np.round(get_pmin(x))), int(np.round(get_pmax(x)))
        ymin, ymax = get_pmin(y), get_pmax(y)
        zmin, zmax = get_pmin(z), get_pmax(z)
        xw, yw = 0.5, 0.05
        xticks = np.arange(xmin, xmax +xw, xw)
        yticks = np.arange(ymin, ymax +yw, yw)
        if estimator == 'maximum likelihood':
            zmin, zmax = -1.5 * 10**6, 0
            zw = 0.5 * 10**6 /2
        elif estimator in ('reduced gtest', 'reduced chi square'):
            zmin, zmax = 0, 80 * 10**3
            zw = 20 * 10**3 /2
        zticks = np.arange(zmin, zmax +zw, zw)
        norm = self.get_normalized_colormap(zmin, zmax, estimator)
        # levels = SD.estimator_levels[estimator]
        # norm = self.get_normalized_colormap(min(levels), max(levels))
        return xticks, yticks, zticks, norm

    def view_errorspace(self, estimator, show_surface=False, show_contours=False, ticksize=7, labelsize=8, textsize=8, titlesize=10, cmap='plasma', elev=None, azim=None, antialiased=False, shade=False, **kwargs):
        """ """
        nrows = 1
        ncols = self.n
        axes = []
        _view_as = np.array([show_surface, show_contours])
        xticks, yticks, zticks, norm = self.get_xyz_ticks_with_norm(estimator)
        fig = plt.figure(**kwargs)
        for idx, timeseries in enumerate(self.datas):
            SD = self.SDs[idx]
            if np.all(_view_as == True):
                ax = fig.add_subplot(nrows, ncols, idx+1, projection='3d')
                # ax = SD.subview_fit_surface_contours(ax, estimator, cmap)
                fig, ax, cbar = SD.subview_fit_surface(fig, ax, estimator, ticksize, labelsize, textsize, titlesize, cmap, elev, azim, antialiased=antialiased, shade=shade, show_contours=True)
                title_mod = 'Surface Contours'
                path_mod = 'surface_contours'
                adj_btm = 0.5
                hspace = 0.4
                wspace = 0.2
            elif np.any(_view_as == True):
                if show_surface == True:
                    ax = fig.add_subplot(nrows, ncols, idx+1, projection='3d')
                    # ax = SD.subview_fit_surface(ax, estimator, xticks, yticks, zticks, ticksize, labelsize, textsize, titlesize, cmap, elev=elev, azim=azim, antialiased=True, shade=True)
                    fig, ax, cbar = SD.subview_fit_surface(fig, ax, estimator, ticksize, labelsize, textsize, titlesize, cmap, elev, azim, antialiased=antialiased, shade=shade, show_contours=False)
                    title_mod = 'Surface Map'
                    path_mod = 'surface'
                    adj_btm = 0.5
                    hspace = 0.4
                    wspace = 0.2
                else:
                    ax = fig.add_subplot(nrows, ncols, idx+1)
                    ax = SD.subview_fit_contours(ax, estimator, norm, xticks, yticks, ticksize, labelsize, textsize, titlesize, cmap)
                    title_mod = 'Contour Map'
                    path_mod = 'contours'
                    adj_btm = 0.5
                    hspace = 0.4
                    wspace = 0.2
            else:
                raise ValueError("show_surface and show_contours cannot both be False")
            axes.append(ax)
        if ((show_surface == False) and (show_contours == True)):
            self.configure_symmetrical_axes(nrows, ncols, np.array(axes), apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        fig.subplots_adjust(bottom=adj_btm, wspace=wspace, hspace=hspace)
        if estimator == 'maximum likelihood':
            ftit = 'CME Speed Distributions:\n Negative Log-Likelihood Error-Space {}'.format(title_mod)
        else:
            ftit = 'CME Speed Distributions:\n {} Error-Space {}'.format(estimator.title(), title_mod)
        fig.suptitle(ftit, fontsize=titlesize, y=1.05)
        fig.tight_layout()
        savepath = self.get_savepath('SDV_speed_fits_{}__{}'.format(path_mod, estimator.replace(' ', '_')), extension='.png')
        self.display_figure(fig, savepath)

class AnalysisViewer(FigureOptions):

    def __init__(self, datas, saveloc=None):
        super().__init__(saveloc)
        if isinstance(datas, (tuple, list, np.ndarray)):
            self.datas = datas
        else:
            self.datas = [datas]
        self.n = len(self.datas)

    def view_power_law_errors(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, ith_resample=0, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            ax = timeseries.unbiased_estimators.subview_power_law_errors(ax, ticksize, labelsize, textsize, ith_resample)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.15, hspace=hspace)
        fig.legend(handles=handles, labels=labels, ncol=2, loc='lower center', mode='expand', fontsize=textsize)
        fig.suptitle('Max-Spectrum Errors', fontsize=titlesize)
        savepath = self.get_savepath('UB__powerlaw_errors', extension='.png')
        self.display_figure(fig, savepath)

    def view_power_law(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, ith_resample=0, **kwargs):
        """ """
        fig = plt.figure(**kwargs)
        axes_lin, axes_log = [], []
        if self.n == 1:
            nrows, ncols = 1, 1
        elif self.n < 4:
            nrows, ncols = 1, self.n
        else:
            nrows, ncols = 2, int(np.ceil(self.n/2))
        for ith_subplot in range(self.n):
            ax_lin = fig.add_subplot(nrows, ncols, ith_subplot +1)
            ax_log = fig.add_subplot(nrows, ncols, ith_subplot +1, frame_on=False)
            axes_lin.append(ax_lin)
            axes_log.append(ax_log)
        axes_lin, axes_log = np.array(axes_lin), np.array(axes_log)
        for ax_lin, ax_log, timeseries in zip(axes_lin, axes_log, self.datas):
            ax_lin = timeseries.unbiased_estimators.subview_linear_power_law(ax_lin, ticksize, labelsize, textsize, ith_resample)
            ax_log = timeseries.unbiased_estimators.subview_logarithmic_power_law(ax_log, ticksize, labelsize, textsize, ith_resample, timeseries.identifier)
            # ax_lin.set_title(timeseries.identifier, fontsize=titlesize)
        if nrows > 1:
            for ax_lin in axes_lin[:-ncols]:
                ax_lin.set_xticklabels([])
                ax_lin.set_xlabel('')
            for ax_log in axes_log[ncols:]:
                ax_log.set_xticklabels([])
                ax_log.set_xlabel('')
            for ax_log in axes_log.reshape((nrows, ncols))[1:, :].ravel():
                ax_log.tick_params(top=False, labeltop=False)
        if ncols > 1:
            for ax_lin in axes_lin[1::ncols]:
                ax_lin.set_yticklabels([])
                ax_lin.set_ylabel('')
            for ax_log in axes_log[::ncols]:
                ax_log.set_yticklabels([])
                ax_log.set_ylabel('')
        if nrows == 2:
            ax_btm_left = axes_lin.reshape((nrows, ncols))[-1, 0]
            ax_btm_left.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=True, right=False, labelleft=True, labelright=False)
            ax_btm_right = axes_lin.reshape((nrows, ncols))[-1, -1]
            ax_btm_right.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=False, right=False, labelleft=False, labelright=False)
            if ncols > 2:
                for ax in axes_lin.reshape((nrows, ncols))[-1, 1:-1].ravel():
                    ax.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=False, right=False, labelleft=False, labelright=False)
        handles, labels = axes_lin[0].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.15) #, hspace=0.5, wspace=0.5)
        fig.legend(handles=handles, labels=labels, ncol=2, loc='lower center', mode='expand', fontsize=textsize)
        fig.suptitle('Power-Law:\nOptimized Fit of Max Spectrum', fontsize=titlesize, y=1.05)
        savepath = self.get_savepath('UB__powerlaw', extension='.png')
        self.display_figure(fig, savepath)

    def view_alpha_hat_histogram(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, counts_type='observed frequency', **kwargs):
        """ """
        fig, ax = plt.subplots(**kwargs)
        for timeseries in self.datas:
            ax = timeseries.unbiased_estimators.subview_alpha_hat_histogram(ax, ticksize, labelsize, textsize, counts_type, transparency=1/self.n, identifier=timeseries.identifier)
        ax.set_title(r'Histogram of $\hat\alpha$' + '\n via ${}$ Resamples'.format(self.datas[0].unbiased_estimators.nresamples), fontsize=titlesize)
        if self.n >= 4:
            modify_singular = True
        else:
            modify_singular = False
        self.configure_symmetrical_axes(nrows=1, ncols=1, axes=[ax], autoconfigure=True, modify_singular=modify_singular)
        fig, ax = self.get_overlay_legend(fig, ax, self.n, textsize=textsize)
        savepath = self.get_savepath('UB__alpha_histogram', extension='.png')
        self.display_figure(fig, savepath)

    def view_theta_hat_histogram(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, counts_type='observed frequency', **kwargs):
        """ """
        fig, ax = plt.subplots(**kwargs)
        for timeseries in self.datas:
            ax = timeseries.unbiased_estimators.subview_theta_hat_histogram(ax, ticksize, labelsize, textsize, counts_type, transparency=1/self.n, identifier=timeseries.identifier)
        ax.set_title(r'Histogram of $\hat\theta$' + '\n via ${}$ Resamples'.format(self.datas[0].unbiased_estimators.nresamples), fontsize=titlesize)
        if self.n >= 4:
            modify_singular = True
        else:
            modify_singular = False
        self.configure_symmetrical_axes(nrows=1, ncols=1, axes=[ax], autoconfigure=True, modify_singular=modify_singular)
        fig, ax = self.get_overlay_legend(fig, ax, self.n, textsize=textsize)
        savepath = self.get_savepath('UB__theta_histogram', extension='.png')
        self.display_figure(fig, savepath)

    def view_point_estimators(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            ax = timeseries.unbiased_estimators.subview_point_estimators(ax, ticksize, labelsize, textsize)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.15, hspace=hspace)
        fig.legend(handles=handles, labels=labels, ncol=len(labels), loc='lower center', mode='expand', fontsize=textsize)
        fig.suptitle('Extremal Index\nPoint Estimators', fontsize=titlesize)
        savepath = self.get_savepath('UB__point_estimators', extension='.png')
        self.display_figure(fig, savepath)

    def view_inter_exceedance_distribution_histogram(self, speed_threshold, ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            ax = timeseries.inter_exceedances[speed_threshold].distribution.subview_histogram_comparison(ax, ticksize, labelsize, textsize)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, modify_columns=True)
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.15, hspace=hspace)
        fig.legend(handles=handles, labels=labels, ncol=2, loc='lower center', mode='expand', fontsize=textsize)
        fig.suptitle('Distribution of Inter-Exceedance Times', fontsize=titlesize)
        savepath = self.get_savepath('IE__v{}_inter_exceedance_distribution'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_cluster_parameterizations(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        biases = np.array(['threshold', 'first-order', 'baseline'])
        nrows = 2
        ncols = self.n
        if self.n == 1:
            fig, (ax_top, ax_btm) = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
            top_axes = np.array([ax_top])
            btm_axes = np.array([ax_btm])
            axes = np.array([ax_top, ax_btm])
            wspace = 0.2
        else:
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
            top_axes = axes[0, :]
            btm_axes = axes[1, :]
        for ax_top, ax_btm, timeseries in zip(top_axes.ravel(), btm_axes.ravel(), self.datas):
            m1, m2, m3 = [], [], []
            t1, t2, t3 = [], [], []
            speed_thresholds = sorted(list(timeseries.cluster_parameters.keys()))
            for speed_threshold in speed_thresholds:
                cluster_parameters = timeseries.cluster_parameters[speed_threshold]
                for bias, m, t in zip(biases[:-1], (m1, m2), (t1, t2)):
                    m.append(cluster_parameters.moment_estimators[bias])
                    t.append(cluster_parameters.time_thresholds[bias])
                if cluster_parameters.moment_estimators[biases[-1]] is not None:
                    m3.append(cluster_parameters.moment_estimators[biases[-1]])
                if cluster_parameters.time_thresholds[biases[-1]] is not None:
                    t3.append(cluster_parameters.time_thresholds[biases[-1]])
            if len(m3) == 0:
                moment_estimators = {bias : m for bias, m in zip(biases[:-1], (m1, m2))}
            else:
                moment_estimators = {bias : m for bias, m in zip(biases, (m1, m2, m3))}
            if len(t3) == 0:
                time_thresholds = {bias : t for bias, t in zip(biases[:-1], (t1, t2))}
            else:
                time_thresholds = {bias : t for bias, t in zip(biases, (t1, t2, t3))}
            ax_top = cluster_parameters.subview_moment_estimators(ax_top, speed_thresholds, moment_estimators, ticksize, labelsize)
            ax_btm = cluster_parameters.subview_time_thresholds(ax_btm, speed_thresholds, time_thresholds, ticksize, labelsize)
            ax_top.set_title(timeseries.identifier, fontsize=titlesize)
            ax_top.set_xticklabels([])
            ax_top.set_xlabel('')
        if self.n > 1:
            wspace = 0.3
            if self.n > 2:
                wspace = 0.325
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, apply_ticklabels=True, autoconfigure=True, tops_and_bottoms=True)
        handles, labels = top_axes.ravel()[0].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.15, wspace=wspace)
        try:
            if labels[2] in (labels[0], labels[1]):
                fig.legend(handles=[handles[0], handles[1]], labels=[labels[0], labels[1]], ncol=2, loc='lower center', mode='expand', fontsize=textsize)
            else:
                fig.legend(handles=[handles[0], handles[1], handles[2]], labels=[labels[0], labels[1], labels[2]], ncol=3, loc='lower center', mode='expand', fontsize=textsize)
        except:
            fig.legend(handles=handles, labels=labels, ncol=len(labels), loc='lower center', mode='expand', fontsize=textsize)
        title = r'Moment Estimator $\hat\theta$' + '\n' + r'& Corresponding Time Threshold'
        if self.n == 1:
            fig.suptitle(title, fontsize=titlesize, y=1.05)
            axes[-1].tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=True, right=False, labelleft=True, labelright=False)
        else:
            fig.suptitle(title, fontsize=titlesize, y=1.05)
        savepath = self.get_savepath('TC__cluster_parameterizations', extension='.png')
        self.display_figure(fig, savepath)

    def view_intra_times_histogram(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', ticksize=7, labelsize=8, textsize=8, titlesize=10, counts_type='observed frequency', **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions, search_values])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            cluster_parameters = timeseries.cluster_parameters[speed_threshold]
            time_threshold = cluster_parameters.time_thresholds[cluster_parameters.bias]
            CS = ClusteringStatistics(clusters)
            CS.initialize_times_and_durations()
            CS.initialize_intra_times_histogram()
            ax = CS.subview_intra_time_histogram(ax, ticksize, labelsize, textsize, counts_type, time_threshold)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True)
        fig.subplots_adjust(hspace=hspace)
        fig.suptitle('Distribution of Intra-Event Times', fontsize=titlesize)
        savepath = self.get_savepath('TC__v{}_intra_times_histogram'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_intra_durations_histogram(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', ticksize=7, labelsize=8, textsize=8, titlesize=10, counts_type='observed frequency', **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions, search_values])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            cluster_parameters = timeseries.cluster_parameters[speed_threshold]
            time_threshold = cluster_parameters.time_thresholds[cluster_parameters.bias]
            CS = ClusteringStatistics(clusters)
            CS.initialize_times_and_durations()
            CS.initialize_intra_durations_histogram()
            ax = CS.subview_intra_duration_histogram(ax, ticksize, labelsize, textsize, counts_type, time_threshold)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True)
        fig.subplots_adjust(hspace=hspace)
        fig.suptitle('Distribution of Intra-Cluster Durations', fontsize=titlesize)
        savepath = self.get_savepath('TC__v{}_intra_durations_histogram'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_inter_durations_histogram(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', ticksize=7, labelsize=8, textsize=8, titlesize=10, counts_type='observed frequency', **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.3
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions, search_values])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            cluster_parameters = timeseries.cluster_parameters[speed_threshold]
            time_threshold = cluster_parameters.time_thresholds[cluster_parameters.bias]
            CS = ClusteringStatistics(clusters)
            CS.initialize_times_and_durations()
            CS.initialize_inter_durations_histogram()
            ax = CS.subview_inter_duration_histogram(ax, ticksize, labelsize, textsize, counts_type, time_threshold)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True)
        fig.subplots_adjust(hspace=hspace)
        fig.suptitle('Distribution of Inter-Cluster Durations', fontsize=titlesize)
        savepath = self.get_savepath('TC__v{}_inter_durations_histogram'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_relative_cluster_statistics(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        fig = plt.figure(**kwargs)
        nrows = 3 +1
        ncols = self.n
        if self.n == 1:
            compare = False
        else:
            compare = True
        axes = []
        for idx in range(self.n * nrows):
            ax = fig.add_subplot(nrows, self.n, idx+1)
            axes.append(ax)
        axes = np.array(axes).reshape((nrows, self.n))
        for ith_col, timeseries in zip(range(self.n), self.datas):
            (ax_top, ax_mid, ax_btm, ax_txt) = axes[:, ith_col]
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions, search_values])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            CS = ClusteringStatistics(clusters)
            CS.initialize_relative_statistics()
            ax_top = CS.subview_relative_cluster_frequency(ax_top, ticksize, labelsize, textsize, compare)
            ax_mid = CS.subview_relative_ejecta_frequency(ax_mid, ticksize, labelsize, textsize, compare)
            ax_btm = CS.subview_relative_probability(ax_btm, ticksize, labelsize, textsize, compare)
            ax_top.set_title(timeseries.identifier, fontsize=titlesize)
            events_label = '${:,}$ Extreme Events'.format(np.concatenate(clusters['elapsed hour'], axis=0).size)
            clusters_label = '${:,}$ Clusters'.format(clusters['elapsed hour'].size)
            ax_txt.text(0.5, 0.25, '{}\nvia {}'.format(events_label, clusters_label), horizontalalignment='center', verticalalignment='center', fontsize=textsize, transform=ax_txt.transAxes)
            ax_txt.spines['left'].set_visible(False)
            ax_txt.spines['right'].set_visible(False)
            ax_txt.spines['top'].set_visible(False)
            ax_txt.spines['bottom'].set_visible(False)
        if self.n == 1:
            for ax in axes[:-2].ravel():
                ax.set_xticklabels([])
                ax.set_xlabel('')
            for ax in axes[0].ravel():
                ax.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=False, labelleft=True, labelright=False)
            for ax in axes[1:-2].ravel():
                ax.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=False, left=True, right=False, labelleft=True, labelright=False)
            for ax in axes[-2].ravel():
                ax.tick_params(which='both', bottom=True, top=True, labeltop=False, labelbottom=True, left=True, right=False, labelleft=True, labelright=False)
            for ax in axes[-1].ravel():
                ax.tick_params(which='both', bottom=False, top=False, labeltop=False, labelbottom=False, left=False, right=False, labelleft=False, labelright=False)
                ax.patch.set_alpha(0)
        else:
            for ax in axes[0, 1:-1].ravel():
                ax.tick_params(which='both', bottom=True, top=False, labeltop=False, labelbottom=False, left=True, right=True, labelleft=False, labelright=False)
                ax.set_ylabel('')
            self.configure_symmetrical_axes(nrows, ncols, axes[:-1, :], apply_labels=True, apply_ticklabels=True, autoconfigure=True, tops_and_bottoms=True)
            for ax in axes[-1, :].ravel():
                ax.tick_params(which='both', bottom=False, top=False, labeltop=False, labelbottom=False, left=False, right=False, labelleft=False, labelright=False)
                ax.patch.set_alpha(0)
            fig.subplots_adjust(hspace=0.3, wspace=0.4)
        if search_parameters is None:
            fig.suptitle('Clustering Statistics', fontsize=titlesize)
        else:
            fig.suptitle('Clustering Statistics\n{} {} {}\n'.format(search_parameters.title(), search_conditions.title(), search_values), y=1.05, fontsize=titlesize)
        fig.align_ylabels()
        savepath = self.get_savepath('TC__v{}_relative_statistics'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_chronological_cluster_sizes(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
            ydelta = 0
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.45
            ydelta = 5
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions, search_values])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            CS = ClusteringStatistics(clusters)
            CS.subview_chronological_cluster_sizes(ax, ticksize, labelsize, textsize, ydelta=ydelta)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        self.configure_symmetrical_axes(nrows, ncols, axes, apply_labels=True, autoconfigure=True)
        red_handle = self.get_rectangular_handle(fill=True, visible=True, facecolor='r')
        blue_handle = self.get_rectangular_handle(fill=True, visible=True, facecolor='b')
        empty_handle = Rectangle((0, 0), 1, 1, alpha=0)
        handles = [empty_handle, tuple((red_handle, blue_handle)), empty_handle]
        kwargs = {'handler_map' : {tuple: HandlerTuple(None)}}
        fig.subplots_adjust(bottom=0.175, hspace=hspace)
        fig.legend(handles=handles, labels=[' ', 'Consecutive Clusters', ' '], ncol=3, loc='lower center', mode='expand', fontsize=textsize, **kwargs)
        fig.suptitle('Chronological Distribution of Clusters', fontsize=titlesize)
        savepath = self.get_savepath('TC__v{}_chronological_cluster_sizes'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_density_optimized_clusters(self, speed_threshold, show_clusters=False, show_noise=False, label_clusters=False, label_noise=False, zoom_clusters=False, cmap='plasma', ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.45
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            CP = timeseries.cluster_parameters[speed_threshold]
            ax = CP.subview_density_optimized_clusters(ax, cmap, show_clusters, show_noise, label_clusters, label_noise, zoom_clusters, ticksize, labelsize, textsize)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        if nrows > 1:
            if ncols == 1:
                for ax in axes[:-1].ravel():
                    ax.set_xlabel('')
            else:
                for ax in axes[:-1, :].ravel():
                    ax.set_xlabel('')
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        # handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.2, hspace=0.325)
        if len(labels) == 1:
            labels = [' '] + labels + [' ']
        dleg = dict(mode='expand', ncol=len(labels), loc='lower center', scatterpoints=5, fancybox=True, shadow=True, fontsize=textsize)
        fig.legend(handles=handles, labels=labels, handler_map={type(handles[0]) : ScatterHandler()}, **dleg)
        # fig, ax = self.get_overlay_legend(fig, axes.ravel()[0], n=1, bottom=0.215, textsize=textsize)
        fig.suptitle('DBSCAN Clusters', fontsize=titlesize, y=1.05)
        # fig.tight_layout()
        if zoom_clusters == True:
            savepath = self.get_savepath('Optimized_Clusters__v{}_Density_zoom'.format(speed_threshold), extension='.png')
        else:
            savepath = self.get_savepath('Optimized_Clusters__v{}_Density'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_dendrogram(self, speed_threshold, cmap='plasma', label_dendrogram=False, tick_connectors=False, ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        if self.n == 1:
            nrows, ncols = 1, 1
            fig, _ax = plt.subplots(**kwargs)
            axes = np.array([_ax])
            hspace = 0.2
        else:
            if self.n < 4:
                nrows = 1
                ncols = self.n
                hspace = 0.2
            else:
                nrows = 2
                ncols = int(np.ceil(self.n/2))
                hspace = 0.45
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        for ax, timeseries in zip(axes.ravel(), self.datas):
            CP = timeseries.cluster_parameters[speed_threshold]
            facecolors = self.get_discrete_colors_from_colormap(cmap, len(CP.agglomerative_optimizer.levels))
            ax = CP.subview_agglomerative_optimized_dendrogram(ax, facecolors, label_dendrogram, tick_connectors, ticksize, labelsize)
            ax.set_title(timeseries.identifier, fontsize=titlesize)
        if self.n > 1:
            if nrows == 1:
                for ax in axes[1:].ravel():
                    ax.set_ylabel('')
            elif ncols == 1:
                for ax in axes[1:].ravel():
                    ax.set_xlabel('')
            else:
                for ax in axes[:-1, :].ravel():
                    ax.set_xlabel('')
                for ax in axes[:, 1:].ravel():
                    ax.set_ylabel('')
        fig.subplots_adjust(hspace=0.325, wspace=0.3)
        fig.suptitle('Agglomerative Hierarchical Clusters\n{} Linkage'.format(CP.agglomerative_optimizer.linkage_criterion.title()), fontsize=titlesize, y=1.05)
        savepath = self.get_savepath('Optimized_Clusters__v{}_Agglomerative'.format(speed_threshold), extension='.png')
        self.display_figure(fig, savepath)

    def view_tabular_parameter_results(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        rowLabels = ['Median Speed\nThreshold Estimate ' + r'($\frac{km}{s}$)', 'Mean Speed\nThreshold Estimate ' + r'($\frac{km}{s}$)']
        rowLabels += ['Extremal Index\nMoment Estimator ' + r'$\hat\theta$', 'Time Threshold\n$T_C$ (hours)', 'Number of Clusters']
        rowLabels += ['Number of Extreme Events', 'Average\nCluster Size', 'Expected Average\nCluster Size ' + r'$\frac{1}{\hat\theta}$']
        rowLabels += ['Mean & \nStandard Deviation\nof Duration (hours)']
        colLabels = []
        cellText = []
        for timeseries in self.datas:
            colLabels.append(timeseries.identifier)
            UB = timeseries.unbiased_estimators
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            CS = ClusteringStatistics(clusters)
            CS.initialize_relative_statistics()
            CS.initialize_times_and_durations()
            CS.initialize_duration_statistics()
            cluster_parameters = timeseries.cluster_parameters[speed_threshold]
            time_threshold = cluster_parameters.time_thresholds[cluster_parameters.bias]
            theta_hat = cluster_parameters.moment_estimators[cluster_parameters.bias]
            vmedian = UB.speed_thresholds['median']
            vmean = UB.speed_thresholds['mean']
            nclusters = clusters['speed'].size
            nevents = np.concatenate(clusters['speed'], axis=0).size
            avg_cluster_size = np.mean(np.concatenate(clusters['cluster size'], axis=0))
            dur_label = '{0:.2f} ± {1:.2f}'.format(np.round(CS.mean_duration[None], decimals=2), np.round(CS.stdev_duration[None], decimals=2))
            data = [vmedian, vmean, '{0:.2f}'.format(np.round(theta_hat, decimals=2)), time_threshold, nclusters, nevents, '{0:.2f}'.format(np.round(avg_cluster_size, decimals=2)), '{0:.2f}'.format(np.round(1/theta_hat, decimals=2)), dur_label]
            cellText.append(data)
        cellText = np.array(cellText).T
        colColours = ['orange' for label in colLabels]
        rowColours = ['orange' for label in rowLabels]
        nrows, ncols = len(rowLabels), len(colLabels)
        colColours = ['orange' for col in range(ncols)]
        rowColours = ['orange' for row in range(nrows)]
        f_even = lambda : ['bisque' if idx % 2 == 0 else 'peachpuff' for idx in range(ncols)]
        f_odd = lambda : ['peachpuff' if idx % 2 == 0 else 'bisque' for idx in range(ncols)]
        cellColours = [f_even() if row % 2 == 0 else f_odd() for row in range(nrows)]
        fig, ax = plt.subplots(**kwargs)
        table = ax.table(cellText=cellText, colLabels=colLabels, rowLabels=rowLabels, colColours=colColours, rowColours=rowColours, cellColours=cellColours, loc='center', cellLoc='center', bbox=(0, 0, 1, 1))
        ax.axis('off')
        ax.axis('tight')
        table.auto_set_font_size(False)
        table.set_fontsize(textsize)
        for key, cell in table.get_celld().items():
            row, col = key
            if row == 0:
                cell.set_text_props(fontproperties=FontProperties(variant='small-caps', weight='semibold', size=labelsize))
        xscale, yscale = (3, 2)
        table.scale(xscale, yscale)
        savepath = self.get_savepath('TABLE__v{}_parameters'.format(speed_threshold), extension='.png')
        fig.tight_layout()
        self.display_figure(fig, savepath)

    def view_tabular_cluster_size_statistics(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        for timeseries in self.datas:
            identifier = timeseries.identifier.replace(':', '')
            identifier = identifier.replace(r'$V_{linear}$', 'vlinear').replace(r'$V_{20 R_{\odot}, i}$', 'v20ri')
            identifier = identifier.replace(r'$V_{20 R_{\odot}, f}$', 'v20rf').replace(r'$V_{20 R_{\odot}}$', 'v20r')
            identifier = identifier.replace('_', '').replace(' ', '_').replace('$_', '_').replace('_$', '')
            cluster_parameters = timeseries.cluster_parameters[speed_threshold]
            time_threshold = cluster_parameters.time_thresholds[cluster_parameters.bias]
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
            CS = ClusteringStatistics(clusters)
            CS.initialize_relative_statistics()
            CS.initialize_times_and_durations()
            CS.initialize_duration_statistics()
            rel_probability = np.array(['{0:.2f}'.format(rp) for rp in np.round(CS.relative_statistics['probability'], decimals=2)])
            cellText = np.array([CS.relative_statistics['cluster size'], CS.relative_statistics['number of clusters'], CS.relative_statistics['number of ejecta'], rel_probability, CS.duration_labels]).T
            colLabels = ['Cluster\nSize', 'Number of\nClusters', 'Number of\nExtreme Events', 'Relative\nProbability', 'Mean & \nStandard Deviation\nof Duration (hours)'] # 'mean duration'
            fig, ax = plt.subplots(**kwargs)
            table = ax.table(cellText=cellText, colLabels=colLabels, loc='center', cellLoc='center', bbox=(0, 0, 1, 1))
            ax.axis('off')
            ax.axis('tight')
            table.auto_set_font_size(False)
            table.set_fontsize(textsize)
            for key, cell in table.get_celld().items():
                row, col = key
                if row == 0:
                    cell.set_text_props(fontproperties=FontProperties(variant='small-caps', weight='semibold', size=labelsize))
                    cell.set_facecolor('orange')
                else:
                    if col % 2 == 0:
                        cell.set_facecolor('#ffe8d8') # 'bisque', 'peachpuff'
                    else:
                        cell.set_facecolor('#fecfaf') # skyblue
            title = 'Cluster Size Statistics ($T_C$ $=$ ${}$ hours)\n{}'.format(time_threshold, timeseries.identifier)
            ax.set_title(title, fontsize=titlesize)
            xscale, yscale = (3, 2)
            table.scale(xscale, yscale)
            savepath = self.get_savepath('TABLE__{}_v{}_clusters'.format(identifier, speed_threshold), extension='.png')
            fig.tight_layout()
            self.display_figure(fig, savepath)

    def save_clusters(self, speed_threshold, search_parameters=None, search_conditions=None, search_values=None, apply_to='all'):
        """ """
        for timeseries in self.datas:
            cluster_state = True
            identifier = timeseries.identifier.replace(':', '')
            identifier = identifier.replace(r'$V_{linear}$', 'vlinear').replace(r'$V_{20 R_{\odot}, i}$', 'v20ri')
            identifier = identifier.replace(r'$V_{20 R_{\odot}, f}$', 'v20rf').replace(r'$V_{20 R_{\odot}}$', 'v20r')
            identifier = identifier.replace('_', '').replace(' ', '_').replace('$_', '_').replace('_$', '')
            cluster_parameters = timeseries.cluster_parameters[speed_threshold]
            time_threshold = cluster_parameters.time_thresholds[cluster_parameters.bias]
            Searcher = timeseries.temporal_clustering[speed_threshold]
            search_inputs = np.array([search_parameters, search_conditions, search_values])
            if np.all(search_inputs == None):
                clusters = Searcher.clusters
            else:
                try:
                    clusters = Searcher.search(search_parameters, search_conditions, search_values, apply_to)
                except:
                    cluster_state = False
            result = '\nCME CLUSTER DATA:\n\t{} Time Threshold = {} hours\n\n'.format(identifier, time_threshold)
            if cluster_state == True:
                for ith_cluster, (dts, elapsed_hours, speeds, cluster_size) in enumerate(zip(clusters['datetime object'], clusters['elapsed hour'], clusters['speed'], clusters['cluster size'])):
                    result = '{}\n ** CLUSTER #{} (size={}):\n'.format(result, ith_cluster +1, cluster_size[0])
                    result = '{}\n .. DATETIMEs:\n{}\n'.format(result, dts)
                    result = '{}\n .. ELAPSED HOURs:\n{}\n'.format(result, elapsed_hours)
                    result = '{}\n .. SPEEDs:\n{}\n\n'.format(result, speeds)
                    result += '-'*10 + '\n'
            else:
                result += '\nWARNING: NO CLUSTERS MATCH SEARCH CONDITIONS'
            savepath = self.get_savepath('v{}__cluster_data__{}'.format(speed_threshold, identifier), extension='.txt')
            self.display_text(result, savepath)

    def save_optimized_clusters(self, speed_threshold):
        """ """
        for timeseries in self.datas:
            result = ''
            identifier = timeseries.identifier.replace(':', '')
            identifier = identifier.replace(r'$V_{linear}$', 'vlinear').replace(r'$V_{20 R_{\odot}, i}$', 'v20ri')
            identifier = identifier.replace(r'$V_{20 R_{\odot}, f}$', 'v20rf').replace(r'$V_{20 R_{\odot}}$', 'v20r')
            identifier = identifier.replace('_', '').replace(' ', '_').replace('$_', '_').replace('_$', '')
            CP = timeseries.cluster_parameters[speed_threshold]
            if CP.density_optimizer is not None:
                result += str(CP.density_optimizer)
                result += '\n'
            if CP.agglomerative_optimizer is not None:
                result += str(CP.agglomerative_optimizer)
            savepath = self.get_savepath('v{}__optimized_cluster_data__{}'.format(speed_threshold, identifier), extension='.txt')
            self.display_text(result, savepath)


class SupplementaryViewer(FigureOptions):

    def __init__(self, saveloc=None):
        super().__init__(saveloc)

    def view_frechet_distribution(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        fig, axes = plt.subplots(nrows=2, ncols=1, **kwargs)
        (ax_top, ax_btm) = axes
        FD = FrechetDistribution()
        ax_top = FD.subview_probability_density(ax_top, ticksize, labelsize, textsize)
        ax_btm = FD.subview_cumulative_density(ax_btm, ticksize, labelsize, textsize)
        ax_top.set_xticklabels([])
        ax_top.set_xlabel('')
        ax_top.set_title(r'Fr$\acute{e}$chet Distribution', fontsize=titlesize)
        handles, labels = ax_top.get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.175)
        fig.legend(handles=handles, labels=labels, ncol=len(labels)//2, loc='lower center', mode='expand', fontsize=textsize)
        savepath = self.get_savepath('EX__frechet_distribution', extension='.png')
        self.display_figure(fig, savepath)

    def view_cluster_example(self, ticksize=7, labelsize=8, textsize=8, titlesize=10, **kwargs):
        """ """
        time_threshold = 5
        clusters = np.array([np.array([2, 3, 5]), np.array([11, 13, 17]), np.array([22, 26])]) - 2
        ncluster = len(clusters)
        xticks = np.unique(np.sum([cluster.copy().tolist() for cluster in clusters]))
        intra_time, intra_duration, inter_duration = [], [], []
        intra_arrowprops = {'arrowstyle': '|-|', 'color' : 'k'}
        inter_arrowprops = {'arrowstyle': '<->', 'color' : 'gray'}
        facecolors = ('darkorange', 'steelblue', 'purple')
        fig, ax = plt.subplots(**kwargs)
        for ith_cluster in range(len(clusters)):
            curr_cluster = clusters[ith_cluster]
            intra_time = np.diff(curr_cluster).tolist()
            intra_duration = curr_cluster[-1] - curr_cluster[0]
            ax.annotate('', xy=(curr_cluster[-1], 0.95), xycoords='data', xytext=(curr_cluster[0], 0.95), textcoords='data', fontsize=textsize, arrowprops=intra_arrowprops)
            ax.text(curr_cluster[0]+1, 0.925, s=r'$\Delta T_{intra}$')
            if ith_cluster < ncluster - 1:
                next_cluster = clusters[ith_cluster+1]
                inter_duration = next_cluster[0] - curr_cluster[-1]
                ax.annotate('', xy=(next_cluster[0], 1.05), xycoords='data', xytext=(curr_cluster[-1], 1.05), textcoords='data', fontsize=textsize, arrowprops=inter_arrowprops)
                ax.text(curr_cluster[-1]+2, 1.075, s=r'$\Delta T_{inter}$')
            ax.scatter(curr_cluster, np.ones(curr_cluster.size), color=facecolors[ith_cluster], label='Cluster #{}'.format(ith_cluster+1))
        ax.plot([np.nan], [np.nan], color='none', label='$T_C = {}$ hours'.format(time_threshold))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, fontsize=ticksize)
        ax.set_xlabel('Consecutive Differences of Elapsed Hours', fontsize=labelsize)
        ax.set_yticks([0.85, 1.15])
        ax.set_yticklabels([])
        ax.set_xlim([-2, 30])
        ax.set_ylim([0.85, 1.15])
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.subplots_adjust(top=0.8)
        fig.suptitle(r'Example of $3$ Clusters', fontsize=titlesize)
        fig.legend(ncol=4, bbox_transform=fig.transFigure, bbox_to_anchor=(0, 0.75, 1, 0.2), fancybox=True, loc='upper center', mode='expand', fontsize=textsize)
        savepath = self.get_savepath('EX__cluster_example', extension='.png')
        self.display_figure(fig, savepath)
