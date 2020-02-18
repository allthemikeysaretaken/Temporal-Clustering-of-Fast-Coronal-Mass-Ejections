from data_processing import *
from optimization_methods import *
from visual_configuration import *

class SpeedSeriesViewer(VisualConfiguration):

    """
    This class is inherited by all classes pertaining to
    timeseries methods. This class inherits methods from
    `VisualConfiguration` to be used as convenience
    functions for plotting.
    """

    def __init__(self, timestep, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize):
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
        super().__init__(directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize)
        self.timestep = timestep
        self.timestep_label = self.get_timestep_label(timestep)
        self.event_type = 'cme'
        self.ref_parameter = 'speed'
        self.event_type_label = 'CME'
        self.unit_label = self.available['label'][self.event_type]['unit'][self.ref_parameter]
        _available = dict()
        _available['distribution model'] = ('normal distribution', 'lognormal distribution')
        self.normal_distribution_model = DistributionModel('normal distribution')
        self.lognormal_distribution_model = DistributionModel('lognormal distribution')
        self.empty_handle = Rectangle((0, 0), 0.1, 0.1, alpha=0)
        self._available.update(_available)

    @staticmethod
    def get_timestep_label(timestep):
        """
        timestep:
            type <str>
        """
        if timestep in ('second', 'minute', 'hour', 'month', 'year'):
            result = '{}ly'.format(timestep)
        elif timestep == 'day':
            result = 'daily'
        else:
            raise ValueError("invalid timestep: {}".format(timestep))
        return result

    def get_search_label(self, events):
        """

        """
        event_type = 'cme'
        identifiers = events['identifier']
        search_kwargs = identifiers['search kwargs']
        if search_kwargs is None:
            return None
        else:
            result = ''
            _parameters, _conditions, _values = search_kwargs['parameters'], search_kwargs['conditions'], search_kwargs['values']
            f = lambda args : list(args) if isinstance(args, (tuple, list, np.ndarray)) else [args]
            parameters, conditions, values = f(_parameters), f(_conditions), f(_values)
            nps, ncs, nvs = len(parameters), len(conditions), len(values)
            is_same_size = (nps == ncs == nvs)
            if is_same_size == True:
                for idx, (parameter, condition, value) in enumerate(zip(parameters, conditions, values)):
                    status = True
                    for non_parameter in ('solar cycle', 'cycle type'):
                        if non_parameter in parameter:
                            status = False
                            break
                    if status == True:
                        dsyms = self.available['label'][event_type]
                        if parameter in list(dsyms['speed type'].keys()):
                            sym_parameter = dsyms['speed type'][parameter]
                        else:
                            sym_parameter = parameter.title()
                        sym_condition = self.available['label']['comparison'][condition]
                        if isinstance(value, int):
                            result += '{} ${}$ ${:,}$'.format(sym_parameter, sym_condition, value)
                        else:
                            result += '{} ${}$ ${}$'.format(sym_parameter, sym_condition, value)
                        try:
                            unit_label = self.available['label'][event_type]['unit'][parameter]
                            result += ' {}'.format(unit_label)
                        except:
                            for pkey in list(self.available['label'][event_type]['unit'].keys()):
                                if pkey in parameter:
                                    unit_label = self.available['label'][event_type]['unit'][pkey]
                                    result += ' {}'.format(unit_label)
                                    break
                        if idx != nps - 1:
                            result += '\n'
            else:
                raise ValueError("{} parameters for {} conditions and {} values".format(nps, ncs, nvs))
            return result

    def get_tail_label(self, extreme_value, extreme_condition):
        """

        """
        sym = self.available['label']['comparison'][extreme_condition]
        unit = self.available['label']['cme']['unit']['speed']
        label = r'Tail $(V_{CME}$ ' + r'${}$ ${}$ {}$)$'.format(sym, extreme_value, unit)
        return label

    @staticmethod
    def get_number_extreme_events_label(events, extreme_value, extreme_condition):
        """

        """
        distribution = events['lognormal distribution']
        x = distribution['speed']
        S = EventSearcher({'x' : x})
        indices = S.search_events(parameters='x', conditions=extreme_condition, values=extreme_value)[1]
        n = np.sum(indices)
        label = '${:,}$ Extreme Events'.format(n)
        return label

    @staticmethod
    def get_ytypes(is_normalized):
        """

        """
        if is_normalized == True:
            counts_type = 'normalized_counts'
            ykey = 'y normalized'
        else:
            counts_type = 'observed_counts'
            ykey = 'y'
        return counts_type, ykey

    def get_subheader(self, events):
        """

        """
        identifiers, parameters = events['identifier'], events['parameter']
        speed_type = identifiers['speed type']
        speed_label = self.available['label']['cme']['speed type'][speed_type]
        _scs = []
        for key in list(parameters.keys()):
            if 'solar cycle' in key:
                _scs.extend(parameters[key].tolist())
        _scs = np.array(_scs)
        delta = np.diff(_scs)
        if np.all(delta == 0) == True:
            sc_label = 'SC ${}$'.format(_scs[0])
            title = '{}: {}'.format(sc_label, speed_label)
        else:
            title = speed_label[:]
        return title

    def autocorrect_legend_kwargs(self, handles, labels):
        """

        """
        nentries = len(handles)
        if nentries == 1:
            ncol = 3
            handles.insert(0, self.empty_handle)
            handles.append(self.empty_handle)
            labels.insert(0, ' ')
            labels.append(' ')
        elif nentries < 4:
            ncol = nentries
        else:
            ncol = None
            for num in np.arange(2, 6, 1).astype(int):
                if nentries % num == 0:
                    ncol = nentries // num
                ncol = nentries // 2
                break
            if ncol is None:
                ncol = nentries // 3
        return handles, labels, ncol

    def autocorrect_optimized_parameters(self, parameters, distribution_model):
        """

        """
        if distribution_model == 'lognormal distribution':
            return parameters # self.lognormal_distribution_model.from_normal_to_lognormal(parameters)
        else:
            return parameters

    def get_path_suffix(self, distribution_model, is_normalized, show_histogram, error_metrics, confidence_metric, stats_metric, tail_metric, point_to_tail, extreme_value):
        """
        distribution_model:
            type <str>

        is_normalized:
            type <str>

        show_histogram:
            type <bool>

        error_metrics:
            type <str / tuple / list / array> or None

        confidence_metric:
            type <str> or None

        stats_metric:
            type <str> or None

        tail_metric:
            type <str> or None

        point_to_tail:
            type <bool>

        extreme_value:
            type <int / float> or None

        """
        result = '{}'.format(distribution_model.replace(' ', '_'))
        if is_normalized == True:
            result += '_norm'
        if show_histogram == True:
            result += '_hist'
        if error_metrics is not None:
            for _error_metric in error_metrics:
                acronym = self.get_acronym(_error_metric)
                result += '_{}'.format(acronym)
        if confidence_metric is not None:
            acronym = self.get_acronym(confidence_metric)
            result += '_confidence_{}'.format(acronym)
        if stats_metric is not None:
            acronym = self.get_acronym(stats_metric)
            result += '_stats_{}'.format(acronym)
        if tail_metric is not None:
            acronym = self.get_acronym(tail_metric)
            result += '_fill_{}'.format(acronym)
        if point_to_tail == True:
            result += '_tailpoint'
        if extreme_value is not None:
            result += '_{}'.format(extreme_value)
        return result

class ErrorSpaceViewer(SpeedSeriesViewer):

    """

    """

    def __init__(self, timestep, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize):
        super().__init__(timestep, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize)
        self.log_base = np.exp(1)

    def get_errorspace(self, events, distribution_model, error_metric):
        """

        """
        distribution = events[distribution_model]
        optimizations = distribution['optimization']
        errors = optimizations[error_metric]
        fun = errors['fun']
        ED = errors['dimensionalization']
        xmin, xmax = np.min(ED.x), np.max(ED.x)
        ymin, ymax = np.min(ED.y), np.max(ED.y)
        _zmin, _zmax = np.min(ED.z), np.max(ED.z)
        if error_metric in ('chi square', 'g-test'):
            zmin = 0.01
            _nmax = len(str(int(_zmax))) - 1
            zmax = self.round_up(_zmax, nearest=10**_nmax)
            # zmax = self.round_up(_zmax, nearest=10)
            norm = LogNorm(vmin=zmin, vmax=zmax)
            is_reduced = errors['is reduced']
            if is_reduced == True:
                zlabel = self.available['label']['optimization'][error_metric]['reduced'][:]
            else:
                zlabel = self.available['label']['optimization'][error_metric]['ordinary'][:]
            fmt = '%.e'
            cticks = ticker.LogLocator(base=self.log_base)
        elif error_metric in ('maximum likelihood estimation',):
            zmin = 0
            _nmax = len(str(int(_zmax))) - 1
            zmax = self.round_up(_zmax, nearest=10**_nmax)
            norm = Normalize(vmin=zmin, vmax=zmax)
            zlabel = self.available['label']['optimization'][error_metric][:]
            fmt = ticker.FuncFormatter(lambda x, p: format(int(x), ','))
            cticks = np.linspace(zmin, zmax, 5)
        else:
            raise ValueError("invalid error_metric: {}".format(error_metric))
        axis_limits = [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
        mu_axis_label = r'$\mu$'
        sigma_axis_label = r'$\sigma$'
        axis_labels = [mu_axis_label, sigma_axis_label, zlabel]
        f = lambda axis_label, prm : '{} $=$ ${:.2f}$'.format(axis_label, np.round(prm, decimals=2))
        mu_legend_label = '{} $(${}$)$'.format(f(mu_axis_label, ED.prms[0]), self.unit_label)
        sigma_legend_label = '{} $(${}$)$'.format(f(sigma_axis_label, ED.prms[1]), self.unit_label)
        error_label = f(zlabel, np.round(errors['fun'], decimals=2))
        legend_labels = [mu_legend_label, sigma_legend_label, error_label]
        return ED, axis_limits, axis_labels, legend_labels, norm, fmt, cticks, fun

    def equalize_errorspace_axes(self, axes, _xmin, _xmax, _ymin, _ymax, _zmin=None, _zmax=None):
        """

        """
        xmin, xmax = np.max(_xmin), np.min(_xmax)
        ymin, ymax = np.max(_ymin), np.min(_ymax)
        if (_zmin is not None) and (_zmax is not None):
            zmin, zmax = np.max(_zmin), np.min(_zmax)
        for ax in axes.ravel():
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.set_xlim([xmin, xmax])
            ax.set_ylim([ymin, ymax])
            if (_zmin is not None) and (_zmax is not None):
                ax.set_zlim([zmin, zmax])
            ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelsize=self.ticksize)
        return axes

    def subview_errorspace_legends(self, fig, nevents, handles, labels, titles, legend_kwargs):
        """

        """
        if nevents == 1:
            leg = fig.legend(handles=handles, labels=labels, handler_map={type(handles[0]) : ScatterHandler()}, **legend_kwargs)
            leg.set_title(titles[0], prop={'size': self.labelsize})
            leg._legend_box.align = "center"
            frame = leg.get_frame()
            frame.set_edgecolor('k')
        elif nevents == 2:
            w, h = 0.3, 0.2
            for handle, label, title, x0, y0, loc in zip(handles, labels, titles, (0, 0.75), (0, 0), ('lower left', 'lower right')):
                leg = fig.legend(handles=[handle], labels=[label], bbox_to_anchor=(x0, y0, w, h), loc=loc, handler_map={type(handle) : ScatterHandler()}, **legend_kwargs)
                leg.set_title(title, prop={'size': self.labelsize})
                leg._legend_box.align = "center"
                frame = leg.get_frame()
                frame.set_edgecolor('k')
        elif nevents == 3:
            w, h = 0.275, 0.2
            for handle, label, title, x0, y0, loc in zip(handles, labels, titles, (0, 0.375, 0.75), (0, 0, 0), ('lower left', 'lower center', 'lower right')):
                leg = fig.legend(handles=[handle], labels=[label], bbox_to_anchor=(x0, y0, w, h), loc=loc, handler_map={type(handle) : ScatterHandler()}, **legend_kwargs)
                leg.set_title(title, prop={'size': self.labelsize})
                leg._legend_box.align = "center"
                frame = leg.get_frame()
                frame.set_edgecolor('k')
        elif nevents == 4:
            w, h = 0.15, 0.2
            for handle, label, title, x0, y0, loc in zip(handles, labels, titles, (0, 0.225, 0.525, 0.75), (0, 0, 0, 0), ('lower left', 'lower left', 'lower right', 'lower right')):
                leg = fig.legend(handles=[handle], labels=[label], bbox_to_anchor=(x0, y0, w, h), loc=loc, handler_map={type(handle) : ScatterHandler()}, **legend_kwargs)
                leg.set_title(title, prop={'size': self.labelsize})
                leg._legend_box.align = "center"
                frame = leg.get_frame()
                frame.set_edgecolor('k')
        else:
            raise ValueError("not yet implemented for more than 4 events; {} were provided".format(nevents))

    def subview_dim2_errorpace_contours(self, ax, ED, error_metric, levels, cmap, central_color, norm, fmt):
        """

        """
        if levels is None:
            args = (ED.X, ED.Y, ED.Z)
        else:
            args = (ED.X, ED.Y, ED.Z, levels)
        if error_metric in ('chi square', 'g-test'):
            fill_handle = ax.contourf(*args, cmap=cmap, norm=norm, locator=ticker.LogLocator(base=self.log_base))
            line_handle = ax.contour(*args, colors=central_color, linewidths=0.5, locator=ticker.LogLocator(base=self.log_base))
        else:
            fill_handle = ax.contourf(*args, cmap=cmap, norm=norm)
            line_handle = ax.contour(*args, colors=central_color, linewidths=0.5)
        optm_handle = ax.scatter(ED.prms[0], ED.prms[1], color=central_color, marker='*', linewidth=0, s=50)
        handles = [fill_handle, line_handle, optm_handle]
        ax.clabel(line_handle, inline=True, fontsize=self.ticksize, fmt=fmt)
        return ax, handles

    def subview_dim3_errorspace_contours(self, ax, ED, levels, cmap, central_color, norm, lw=3):
        """

        """
        if levels is None:
            args = (ED.X, ED.Y, ED.Z)
        else:
            args = (ED.X, ED.Y, ED.Z, levels)
        offset = np.min(ED.Z) * -1
        fill_handle = ax.contour(*args, cmap=cmap, norm=norm, lw=3, linestyles='solid', offset=offset)
        line_handle = ax.contour(*args, colors=central_color, lw=3, linestyles='solid')
        handles = [fill_handle, line_handle]
        return ax, handles

    def subview_dim3_errorspace_surface(self, ax, ED, cmap, norm, rstride=1, cstride=1, shade=False, alpha=0.7, lw=0.5):
        """

        """
        surface_handle = ax.plot_surface(ED.X, ED.Y, ED.Z, rstride=rstride, cstride=cstride, cmap=cmap, norm=norm, shade=shade, alpha=alpha, lw=lw)
        surface_handle._facecolors2d = surface_handle._facecolors3d
        surface_handle._edgecolors2d = surface_handle._edgecolors3d
        return ax, surface_handle

    def subview_errorspace_colorbar(self, fig, ax, layout, handle, norm, fmt, clabel=None, cticks=None):
        """

        """
        if 'vertical' in layout:
            orientation = 'vertical'
        elif 'horizontal' in layout:
            orientation = 'horizontal'
        else:
            raise ValueError("invalid layout: {}".format(layout))
        if cticks is None:
            cbar = fig.colorbar(handle, ax=ax, orientation=orientation, shrink=0.75, pad=0.1, norm=norm, format=fmt, extend='both')
        else:
            cbar = fig.colorbar(handle, ax=ax, orientation=orientation, shrink=0.75, pad=0.1, norm=norm, format=fmt, extend='both', ticks=cticks)
        cbar.ax.tick_params(labelsize=self.ticksize)
        if clabel is not None:
            cbar.ax.set_title(clabel, fontsize=self.labelsize)
        return fig, ax, cbar

    @staticmethod
    def set_viewing_position(ax, azim=None, elev=None, dist=None):
        """
        ax:
            type <matplotlib object>
        azim:
            type <int / float> or None
        elev:
            type <int / float> or None
        dist:
            type <int / float> or None
        """
        if azim is not None:
            ax.azim = azim
        if elev is not None:
            ax.elev = elev
        if dist is not None:
            ax.dist = dist
        return ax

class SpeedSeries(ErrorSpaceViewer):

    """

    """

    def __init__(self, DB, timestep, directory, ticksize=7, labelsize=8, textsize=9, titlesize=11, headersize=20, cellsize=15):
        """
        DB:
            type <custom class>

        ...
        """
        super().__init__(timestep, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize)
        self.DB = DB
        self._events = []
        self._available.update(dict(DB.available))

    @property
    def events(self):
        return self._events

    def load_events(self, speed_type, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, nan_kwargs=None):
        """
        speed_type:
            type <str>

        ...

        nan_kwargs:
            type <dict> or None
        """
        data = getattr(self.DB, 'cme')
        if data is None:
            raise ValueError("DB has not initialized event_type 'cme'")
        if (search_parameters is None) and (search_conditions is None) and (search_values is None):
            search_kwargs = None
        else:
            search_kwargs = dict()
            search_kwargs['parameters'] = search_parameters
            search_kwargs['conditions'] = search_conditions
            search_kwargs['values'] = search_values
            search_kwargs['apply_to'] = apply_to
            search_kwargs['modifiers'] = modifiers
        if nan_kwargs is None:
            nan_kwargs = dict(parameter=None, policy=None, repl=0)
        elif not isinstance(nan_kwargs, dict):
            raise ValueError("invalid type(nan_kwargs): {}".format(type(nan_kwargs)))
        parameters = self.DB.filter_nans(data, **nan_kwargs)
        parameters['elapsed'] = self.DB.get_elapsed_time(parameters, self.timestep)
        if search_kwargs is not  None:
            S = EventSearcher(parameters)
            parameters = S.search_events(**search_kwargs)[0]
        identifiers = dict()
        identifiers['speed type'] = speed_type
        identifiers['search kwargs'] = search_kwargs
        identifiers['nan kwargs'] = nan_kwargs
        normal_distribution = dict()
        normal_distribution['model'] = self.normal_distribution_model
        normal_distribution['speed'] = np.sort(np.log(parameters[speed_type]))
        lognormal_distribution = dict()
        lognormal_distribution['model'] = self.lognormal_distribution_model
        lognormal_distribution['speed'] = np.sort(parameters[speed_type])
        result = {'parameter' : parameters, 'identifier' : identifiers, 'normal distribution' : normal_distribution, 'lognormal distribution' : lognormal_distribution}
        self._events.append(result)

    def load_histograms(self, distribution_model, **kwargs):
        """
        distribution_model:
            type <str>
        """
        for events in self.events:
            distribution = events[distribution_model]
            vs = distribution['speed']
            H = Histogram(vs, distribution_model)
            H.initialize_histogram(**kwargs)
            events[distribution_model]['histogram'] = H

    def load_optimizations(self, error_metrics, distribution_model, parameter_guess=None, scale='local', reduce_statistic=False, minimizer_kwargs=None, **kwargs):
        """

        """
        if isinstance(error_metrics, str):
            error_metrics = [error_metrics]
        elif not isinstance(error_metrics, (tuple, list, np.ndarray)):
            raise ValueError("invalid type(error_metrics): {}".format(type(error_metrics)))
        for events in self.events:
            for error_metric in error_metrics:
                result = dict()
                distribution = events[distribution_model]
                x, model = distribution['speed'], distribution['model']
                if isinstance(parameter_guess, str):
                    try:
                        parameter_guess = distribution['optimization'][parameter_guess]['parameter']
                    except:
                        raise ValueError("could not obtain parameter_guess from {}".format(parameter_guess))
                if error_metric == 'maximum likelihood estimation':
                    MLE = MaximumLikelihoodEstimation(x, model)
                    optimization_result = MLE.fit(parameter_guess, scale, minimizer_kwargs, **kwargs)
                    y = MLE.model.pdf(optimization_result['parameters'], x)
                    result['optimizer'] = MLE
                    result['y normalized'] = y
                    if 'histogram' in list(distribution.keys()):
                        H = distribution['histogram']
                        result['y'] = y * H.normalization_constant
                    result['statistics'] = MLE.model.get_statistics(optimization_result['parameters'])
                elif error_metric in ('chi square', 'g-test'):
                    if 'histogram' not in list(distribution.keys()):
                        raise ValueError("cannot evaluate ")
                    H = distribution['histogram']
                    BS = BinStatistics(H)
                    optimization_result = BS.fit(error_metric, parameter_guess, scale, reduce_statistic, minimizer_kwargs, **kwargs)
                    y = BS.histogram.model.pdf(optimization_result['parameters'], x)
                    result['optimizer'] = BS
                    result['y normalized'] = y
                    result['y'] = y * H.normalization_constant
                    result['is reduced'] = reduce_statistic
                    result['statistics'] = BS.histogram.model.get_statistics(optimization_result['parameters'])
                else:
                    raise ValueError("invalid error_metric: {}".format(error_metric))
                result['parameter'] = optimization_result['parameters']
                # result['parameter'] = self.autocorrect_optimized_parameters(optimization_result['parameters'], distribution_model)
                result['fun'] = optimization_result['fun']
                if 'optimization' not in list(distribution.keys()):
                    distribution['optimization'] = dict()
                distribution['optimization'][error_metric] = result

    def load_error_dimensionalizations(self, error_metrics, distribution_model, x=None, y=None, xfrac=0.5, xn=7, xspace_by='number', yfrac=0.5, yn=7, yspace_by='number'):
        """

        """
        if isinstance(error_metrics, str):
            error_metrics = [error_metrics]
        elif not isinstance(error_metrics, (tuple, list, np.ndarray)):
            raise ValueError("invalid type(error_metrics): {}".format(type(error_metrics)))
        for events in self.events:
            for error_metric in error_metrics:
                distribution = events[distribution_model]
                optimizer = distribution['optimization'][error_metric]['optimizer']
                prms = distribution['optimization'][error_metric]['parameter']
                ED = ErrorDimensionalization(optimizer, prms)
                ED.load_parameter_space(x, y, xfrac, xn, xspace_by, yfrac, yn, yspace_by)
                distribution['optimization'][error_metric]['dimensionalization'] = ED

    def subview_histograms(self, ax, events, distribution_model, counts_type, histogram_colors, as_steps, as_bars, alpha=1):
        """

        """
        if (as_steps == False) and (as_bars == False):
            raise ValueError("as_steps = False and as_bars = False; at least one must be set to True")
        ncolors = len(histogram_colors)
        distribution = events[distribution_model]
        histogram = distribution['histogram']
        y = getattr(histogram, counts_type)
        if (as_steps == True) and (as_bars == True):
            if ncolors > 1:
                _handle = Rectangle((0, 0), 0.1, 0.1, fill=True, visible=True, facecolor=histogram_colors[0], edgecolor=histogram_colors[1], alpha=alpha)
                ax.step(histogram.midpoints, y, color=histogram_colors[1], linewidth=0.5, where='mid', alpha=alpha)
            else:
                _handle = Rectangle((0, 0), 0.1, 0.1, fill=True, visible=True, facecolor=histogram_colors[0], alpha=alpha)
                ax.step(histogram.midpoints, y, color=histogram_colors[0], linewidth=0.5, where='mid', alpha=alpha)
            ax.bar(histogram.midpoints, y, width=histogram.bin_widths, facecolor=histogram_colors[0], alpha=alpha)
        elif as_steps == True:
            _handle = ax.step(histogram.midpoints, y, color=histogram_colors[0], linewidth=0.5, where='mid', alpha=alpha)
        else:
            if ncolors > 1:
                _handle = ax.bar(histogram.midpoints, y, width=histogram.bin_widths, facecolor=histogram_colors[0], edgecolor=histogram_colors[1], alpha=alpha)
            else:
                _handle = ax.bar(histogram.midpoints, y, width=histogram.bin_widths, facecolor=histogram_colors[0], alpha=alpha)
        return ax, _handle

    def subview_optimized_fits(self, ax, events, distribution_model, error_metrics, ykey, fit_colors):
        """

        """
        if error_metrics is None:
            raise ValueError("error_metrics is None")
        distribution = events[distribution_model]
        x = distribution['speed']
        optimizations = distribution['optimization']
        ncolors = len(fit_colors)
        noptims = len(error_metrics)
        if ncolors < noptims:
            raise ValueError("{} colors for {} optimizations".format(ncolors, noptims))
        _handles, _labels = [], []
        for error_metric, facecolor in zip(error_metrics, fit_colors):
            _optimization = optimizations[error_metric]
            y = _optimization[ykey]
            if error_metric == 'maximum likelihood estimation':
                _label = self.available['label']['optimization'][error_metric]
            else:
                bin_key = 'reduced' if _optimization['is reduced'] == True else 'ordinary'
                _label = self.available['label']['optimization'][error_metric][bin_key]
            _handle = ax.plot(x, y, color=facecolor, alpha=1/noptims)
            _handles.extend(_handle)
            _labels.append(_label)
        return ax, _handles, _labels

    def subview_confidence_interval(self, ax, events, distribution_model, confidence_metric, ykey, confidence_color, nintervals):
        """

        """
        if confidence_metric is None:
            raise ValueError("confidence_metric is None")
        if nintervals < 1:
            raise ValueError("nintervals must be greater than or equal to one")
        distribution = events[distribution_model]
        x = distribution['speed']
        optimization = distribution['optimization']
        metric = optimization[confidence_metric]
        y = metric[ykey]
        statistics = metric['statistics']
        mu, sigma = statistics['mean'], statistics['standard deviation'] # mu, sigma = self.lognormal_distribution_model.from_lognormal_to_normal([m, s])
        _handles = []
        _labels = [r'$\mu \pm \sigma$']
        if nintervals > 1:
            for i in range(nintervals):
                j = i+1
                _label = r'$\mu \pm {} \sigma$'.format(j)
                _labels.append(_label)
        n = len(_labels)
        alphas = np.cumsum(np.ones(n) * 0.05)[::-1]
        for i, alpha in enumerate(alphas):
            j = (i+1) * sigma
            _handle = ax.fill_between(x, y, where=((mu - j) <= x) & (x <= (mu + j)), color=confidence_color, alpha=alpha)
            _handles.append(_handle)
        return ax, _handles, _labels

    def subview_statistics(self, ax, events, distribution_model, stats_metric, stats_colors):
        """

        """
        if stats_metric is None:
            raise ValueError("stats_metric is None")
        if isinstance(stats_colors, str):
            stats_colors = [stats_colors]
        ncolors = len(stats_colors)
        distribution = events[distribution_model]
        optimization = distribution['optimization']
        metric = optimization[stats_metric]
        statistics = metric['statistics']
        nstats = len(list(statistics.keys()))
        if ncolors < nstats:
            raise ValueError("{} colors provided for {} statistics".format(ncolors, nstatistics))
        _handles, _labels, partial_labels = [], [], []
        for facecolor, (stats_key, stats_value) in zip(stats_colors, statistics.items()):
            ylim = ax.get_ylim()
            _handle = ax.axvline(stats_value, ymin=0, ymax=ylim[1], color=facecolor, linestyle=':', marker=None, alpha=0.5)
            _label = '{}'.format(stats_key.title())
            partial_label = '{} $=$ ${}$ {}'.format(_label, int(stats_value), self.unit_label)
            _handles.append(_handle)
            _labels.append(_label)
            partial_labels.append(partial_label)
        ncol = 1 if len(partial_labels) == 2 else 2
        ax.legend(handles=_handles, labels=partial_labels, fontsize=self.ticksize, loc='upper right', ncol=ncol) # loc='lower left', bbox_to_anchor=(0.7, 0.7, 0.25, 0.25))
        return ax, _handles, _labels


    # def subview_statistics(self, ax, events, distribution_model, stats_metric, stats_colors):
    #     """
    #
    #     """
    #     if stats_metric is None:
    #         raise ValueError("stats_metric is None")
    #     if isinstance(stats_colors, str):
    #         stats_colors = [stats_colors]
    #     ncolors = len(stats_colors)
    #     distribution = events[distribution_model]
    #     optimization = distribution['optimization']
    #     metric = optimization[stats_metric]
    #     statistics = metric['statistics']
    #     nstats = len(list(statistics.keys()))
    #     if ncolors < nstats:
    #         raise ValueError("{} colors provided for {} statistics".format(ncolors, nstatistics))
    #     _handles, _labels = [], []
    #     for facecolor, (stats_key, stats_value) in zip(stats_colors, statistics.items()):
    #         ylim = ax.get_ylim()
    #         _handle = ax.axvline(stats_value, ymin=0, ymax=ylim[1], color=facecolor, linestyle=':', marker=None, alpha=0.5)
    #         # _label = '{}'.format(stats_key.title())
    #         _label = '{} $=$ ${}$ {}'.format(stats_key.title(), int(stats_value), self.unit_label)
    #         _handles.append(_handle)
    #         _labels.append(_label)
    #     return ax, _handles, _labels

    def subview_filled_tail(self, ax, events, distribution_model, ykey, extreme_value, extreme_condition, tail_metric, fill_color):
        """

        """
        if extreme_value is None:
            raise ValueError("extreme_value is None")
        _label = self.get_tail_label(extreme_value, extreme_condition)
        distribution = events[distribution_model]
        x = distribution['speed']
        optimizations = distribution['optimization']
        optimization = optimizations[tail_metric]
        y = optimization[ykey]
        S = EventSearcher({'x' : x})
        if distribution_model == 'normal distribution':
            indices = S.search_events(parameters='x', conditions=extreme_condition, values=np.log(extreme_value))[1]
        else:
            indices = S.search_events(parameters='x', conditions=extreme_condition, values=extreme_value)[1]
        _handle = ax.fill_between(x, 0, y, where=indices, color='red', alpha=0.875)
        return ax, _handle, _label

    def subview_arrow_at_extreme_threshold(self, ax, events, distribution_model, ykey, extreme_value, extreme_condition):
        """

        """
        if extreme_value is None:
            raise ValueError("extreme_value is None")
        arrowprops = {'arrowstyle': '->', 'color' : 'k'}
        _label = self.get_tail_label(extreme_value, extreme_condition)
        _y = ax.get_ylim()[1] / 3
        if distribution_model == 'normal distribution':
            _x = np.log(extreme_value)
            xy = (_x, _y)
            xytext = (_x * 0.35, _y * 1.25)
        else:
            distribution = events[distribution_model]
            x = distribution['speed']
            loc = np.where(x == extreme_value)[0]
            _x = x[loc][0]
            xy = (_x, _y)
            xytext = (_x * 1.25, _y * 1.25)
        ax.annotate(_label, xy=xy, xytext=xytext, xycoords='data', textcoords='data', fontsize=self.textsize, arrowprops=arrowprops)
        return ax

    def view_distribution(self, distribution_model, layout, error_metrics=None, confidence_metric=None, stats_metric=None, tail_metric=None, extreme_value=None, extreme_condition='greater than or equal', show_histogram=False, point_to_tail=False, is_normalized=False, histogram_colors=('steelblue', 'darkorange'), fit_colors=('r', 'b', 'darkgreen'), confidence_color='b', stats_colors=('darkred', 'mediumpurple', 'darkblue', 'darkmagenta'), nintervals=3, tail_color='darkred', as_steps=False, as_bars=False, save=False, **kwargs):
        """

        """
        nevents = len(self.events)
        counts_type, ykey = self.get_ytypes(is_normalized)
        if isinstance(error_metrics, str):
            error_metrics = [error_metrics]
        xlabel = r'Speed $(${}$)$'.format(self.unit_label)
        ylabel = counts_type.replace('_', ' ').replace('counts', 'frequency').title()
        header = 'Distribution of CME Speeds'
        vecdir = self.get_vecdir(layout)
        if isinstance(histogram_colors, str):
            histogram_colors = [histogram_colors]
        if isinstance(fit_colors, str):
            fit_colors = [fit_colors]
        legend_kwargs = dict(loc='lower center', mode='expand', fontsize=self.labelsize)
        handles, labels = [], []
        fig, axes = self.get_dim2_figure_and_axes(layout, nevents, **kwargs)
        try:
            hspace = 0.3
            for idx, (ax, events) in enumerate(zip(axes.ravel(), self.events)):
                title = self.get_subheader(events)
                search_label = self.get_search_label(events)
                if search_label is not None:
                    title = '{}\n{}'.format(title, search_label)
                    hspace = 0.4
                if show_histogram == True:
                    ax, _handle = self.subview_histograms(ax, events, distribution_model, counts_type, histogram_colors, as_steps, as_bars, alpha=0.7)
                    if idx == 0:
                        handles.append(_handle)
                        labels.append('Histogram')
                if error_metrics is not None:
                    ax, _handles, _labels = self.subview_optimized_fits(ax, events, distribution_model, error_metrics, ykey, fit_colors)
                    if idx == 0:
                        handles.extend(_handles)
                        labels.extend(_labels)
                if confidence_metric is not None:
                    ax, _handles, _labels = self.subview_confidence_interval(ax, events, distribution_model, confidence_metric, ykey, confidence_color, nintervals)
                    if idx == 0:
                        handles.extend(_handles)
                        labels.extend(_labels)
                if stats_metric is not None:
                    ax, _handles, _labels = self.subview_statistics(ax, events, distribution_model, stats_metric, stats_colors)
                    if idx == 0:
                        handles.extend(_handles)
                        labels.extend(_labels)
                if point_to_tail == True:
                    ax = self.subview_arrow_at_extreme_threshold(ax, events, distribution_model, ykey, extreme_value, extreme_condition)
                if tail_metric is not None:
                    ax, _handle, _label = self.subview_filled_tail(ax, events, distribution_model, ykey, extreme_value, extreme_condition, tail_metric, tail_color)
                    if idx == 0:
                        handles.append(_handle)
                        labels.append(_label)
                ax.set_title(title, fontsize=self.labelsize)
            for ax in axes.ravel():
                ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
                ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
                ax.tick_params(axis='both', which='both', labelsize=self.ticksize)
                ax.grid(color='k', linestyle=':', alpha=0.3)
        except:
            for idx, events in enumerate(self.events):
                _label = self.get_subheader(events)
                search_label = self.get_search_label(events)
                if search_label is not None:
                    _label = '{}\n{}'.format(_label, search_label)
                if show_histogram == True:
                    axes, _handle = self.subview_histograms(axes, events, distribution_model, counts_type, histogram_colors, as_steps, as_bars, alpha=1/nevents)
                    handles.append(_handle)
                    labels.append(_label)
                if error_metrics is not None:
                    axes, _handles, _labels = self.subview_optimized_fits(axes, events, distribution_model, error_metrics, ykey, fit_colors)
                    handles.extend(_handles)
                    labels.extend(_labels)
                if confidence_metric is not None:
                    axes, _handles, _labels = self.subview_confidence_interval(axes, events, distribution_model, confidence_metric, ykey, confidence_color, nintervals)
                    handles.extend(_handles)
                    labels.extend(_labels)
                if stats_metric is not None:
                    axes, _handles, _labels = self.subview_statistics(axes, events, distribution_model, stats_metric, stats_colors)
                    handles.extend(_handles)
                    labels.extend(_labels)
                if point_to_tail == True:
                    axes = self.subview_arrow_at_extreme_threshold(axes, events, distribution_model, ykey, extreme_value, extreme_condition)
                if tail_metric is not None:
                    axes, _handle, _label = self.subview_filled_tail(axes, events, distribution_model, ykey, extreme_value, extreme_condition, tail_metric, tail_color)
                    handles.append(_handle)
                    labels.append(_label)
            axes.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes.tick_params(axis='both', which='both', labelsize=self.ticksize)
            axes.grid(color='k', linestyle=':', alpha=0.3)
        fig, axes = self.add_shared_labels(fig, axes, xlabel, ylabel, vecdir=vecdir)
        fig.suptitle(header, fontsize=self.titlesize)
        handles, labels, ncol = self.autocorrect_legend_kwargs(handles, labels)
        legend_kwargs['ncol'] = ncol
        legend_kwargs['handles'] = handles
        legend_kwargs['labels'] = labels
        fig.subplots_adjust(bottom=0.2, wspace=0.3, hspace=hspace)
        fig.legend(**legend_kwargs)
        if save == True:
            suffix = self.get_path_suffix(distribution_model, is_normalized, show_histogram, error_metrics, confidence_metric, stats_metric, tail_metric, point_to_tail, extreme_value)
            savename = 'speed_distribution_{}_{}'.format(suffix, layout)
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def view_histogram_tail(self, layout, extreme_value=None, extreme_condition='greater than or equal', is_normalized=False, histogram_colors=('steelblue', 'darkorange'), tail_color='r', as_steps=False, as_bars=False, save=False, **kwargs):
        """

        """
        if extreme_value is None:
            raise ValueError("extreme_value is None")
        distribution_model = 'lognormal distribution'
        nevents = len(self.events)
        counts_type, ykey = self.get_ytypes(is_normalized)
        xlim, ylim = (625, 15625), (0.5, 3125)
        xlabel = r'Speed $(${}$)$'.format(self.unit_label)
        ylabel = counts_type.replace('_', ' ').replace('counts', 'frequency').title()
        header = 'Distribution Tail of CMEs'
        vecdir = self.get_vecdir(layout)
        if isinstance(histogram_colors, str):
            histogram_colors = [histogram_colors]
        legend_kwargs = dict(loc='lower center', mode='expand', fontsize=self.labelsize)
        bbox = {'facecolor': 'gray', 'alpha': 0.2, 'pad': 2}
        handles, labels = [], []
        fig, axes = self.get_dim2_figure_and_axes(layout, nevents, **kwargs)
        try:
            hspace = 0.3
            for idx, (ax, events) in enumerate(zip(axes.ravel(), self.events)):
                title = self.get_subheader(events)
                search_label = self.get_search_label(events)
                if search_label is not None:
                    title = '{}\n{}'.format(title, search_label)
                    hspace = 0.4
                ax, _handle = self.subview_histograms(ax, events, distribution_model, counts_type, histogram_colors, as_steps, as_bars, alpha=0.7)
                if idx == 0:
                    handles.append(_handle)
                    labels.append('Histogram')
                nevents_label = self.get_number_extreme_events_label(events, extreme_value, extreme_condition)
                ax.text(0.95, 0.95, nevents_label, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, fontsize=self.labelsize, bbox=bbox)
                _handle = ax.axvline(extreme_value, ymin=0, ymax=ylim[-1], color=tail_color, linestyle='--')
                if idx == 0:
                    tail_label = self.get_tail_label(extreme_value, extreme_condition)
                    handles.append(_handle)
                    labels.append(tail_label)
                ax.set_xscale('log', basex=5)
                ax.set_yscale('log', basey=5)
                ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
                ax.xaxis.set_minor_formatter(ticker.NullFormatter())
                ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
                ax.yaxis.set_minor_formatter(ticker.NullFormatter())
                ax.set_title(title, fontsize=self.labelsize)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                ax.tick_params(axis='both', which='both', labelsize=self.ticksize)
                ax.grid(color='k', linestyle=':', alpha=0.3)
        except:
            for idx, events in enumerate(self.events):
                _label = self.get_subheader(events)
                search_label = self.get_search_label(events)
                if search_label is not None:
                    _label = '{}\n{}'.format(_label, search_label)
                axes, _handle = self.subview_histograms(axes, events, distribution_model, counts_type, histogram_colors, as_steps, as_bars, alpha=1/nevents)
                handles.append(_handle)
                labels.append(_label)
            _handle = axes.axvline(extreme_value, ymin=ylim[0], ymax=ylim[-1], color=tail_color, linestyle='--')
            _label = self.get_tail_label(extreme_value, extreme_condition)
            handles.append(_handle)
            labels.append(_label)
            axes.set_xscale('log', basex=5)
            axes.set_yscale('log', basey=5)
            axes.xaxis.set_major_formatter(ticker.ScalarFormatter())
            axes.xaxis.set_minor_formatter(ticker.NullFormatter())
            axes.yaxis.set_major_formatter(ticker.ScalarFormatter())
            axes.yaxis.set_minor_formatter(ticker.NullFormatter())
            axes.set_xlim(xlim)
            axes.set_ylim(ylim)
            axes.tick_params(axis='both', which='both', labelsize=self.ticksize)
            axes.grid(color='k', linestyle=':', alpha=0.3)
        fig, axes = self.add_shared_labels(fig, axes, xlabel, ylabel, vecdir=vecdir)
        fig.suptitle(header, fontsize=self.titlesize)
        handles, labels, ncol = self.autocorrect_legend_kwargs(handles, labels)
        legend_kwargs['ncol'] = ncol
        legend_kwargs['handles'] = handles
        legend_kwargs['labels'] = labels
        fig.subplots_adjust(bottom=0.2, wspace=0.3, hspace=hspace)
        fig.legend(**legend_kwargs)
        if save == True:
            savename = 'speed_distribution_tail_{}_{}_{}'.format(distribution_model, extreme_value, layout)
        else:
            savename = None
        self.display_image(fig, savename=savename)

    def view_dim2_errorspace(self, distribution_model, layout, error_metrics, levels=None, cmap='plasma', central_color='k', include_colorbar=False, save=False, **kwargs):
        """

        """
        if isinstance(error_metrics, str):
            error_metrics = [error_metrics]
        elif not isinstance(error_metrics, (tuple, list, np.ndarray)):
            raise ValueError("invalid type(error_metrics): {}".format(type(error_metrics)))
        if layout == 'overlay':
            raise ValueError("invalid layout; try a different one")
        nevents = len(self.events)
        vecdir = self.get_vecdir(layout)
        legend_kwargs = dict(fontsize=self.labelsize) # loc='lower center', mode='expand',
        hspace = 0.3
        for error_metric in error_metrics:
            xmin, xmax, ymin, ymax, zmin, zmax = [], [], [], [], [], []
            _opt_handles, _opt_labels = [], []
            titles = []
            header = '{} Error-Space'.format(error_metric.title())
            fig, axes = self.get_dim2_figure_and_axes(layout, nevents, **kwargs)
            if not isinstance(axes, np.ndarray):
                axes = np.array([axes])
            for idx, (ax, events) in enumerate(zip(axes.ravel(), self.events)):
                subheader = self.get_subheader(events)
                search_label = self.get_search_label(events)
                if search_label is None:
                    title = subheader
                else:
                    title = '{}\n{}'.format(subheader, search_label)
                    hspace = 0.4
                if '\n' in title:
                    titles.append('{}\n'.format(title))
                else:
                    titles.append(title)
                ED, axis_limits, axis_labels, legend_labels, norm, fmt, cticks, fun = self.get_errorspace(events, distribution_model, error_metric)
                [(_xmin, _xmax), (_ymin, _ymax), (_zmin, _zmax)] = axis_limits
                [xlabel, ylabel, zlabel] = axis_labels
                [mu_label, sigma_label, error_label] = legend_labels
                opt_label = '{}\n{}\n{}'.format(mu_label, sigma_label, error_label)
                _opt_labels.append(opt_label)
                xmin.append(_xmin)
                xmax.append(_xmax)
                ymin.append(_ymin)
                ymax.append(_ymax)
                zmin.append(_zmin)
                zmax.append(_zmax)
                ax.set_title(title, fontsize=self.labelsize, pad=10)
                ax, _handles = self.subview_dim2_errorpace_contours(ax, ED, error_metric, levels, cmap, central_color, norm, fmt)
                [fill_handle, line_handle, optm_handle] = _handles
                _opt_handles.append(optm_handle)
                if include_colorbar == True:
                    fig, ax, cbar = self.subview_errorspace_colorbar(fig, ax, layout, fill_handle, norm, fmt, clabel=zlabel, cticks=cticks)
            axes = self.equalize_errorspace_axes(axes, xmin, xmax, ymin, ymax, _zmin=None, _zmax=None)
            fig, axes = self.add_shared_labels(fig, axes, xlabel, ylabel, vecdir)
            fig.align_ylabels()
            try:
                fig.subplots_adjust(hspace=0.15 * nevents, wspace=0.25, bottom=0.25)
            except:
                fig.subplots_adjust(bottom=0.25)
            fig.suptitle(header, fontsize=self.titlesize)
            self.subview_errorspace_legends(fig, nevents, _opt_handles, _opt_labels, titles, legend_kwargs)
            if save == True:
                acronym = self.get_acronym(error_metric)
                savename = 'speed_dim2_error_{}_{}_{}'.format(distribution_model, acronym, layout)
            else:
                savename = None
            self.display_image(fig, savename=savename)

    def view_dim3_errorspace(self, distribution_model, layout, error_metrics, show_contours=False, show_surface=False, levels=None, cmap='plasma', central_color='k', include_colorbar=False, rstride=1, cstride=1, shade=False, azim=None, elev=None, dist=None, save=False, **kwargs):
        """

        """
        if (show_contours != True) and (show_surface != True):
            raise ValueError("show_contours and/or show_surface should be set to True")
        if isinstance(error_metrics, str):
            error_metrics = [error_metrics]
        elif not isinstance(error_metrics, (tuple, list, np.ndarray)):
            raise ValueError("invalid type(error_metrics): {}".format(type(error_metrics)))
        if layout == 'overlay':
            raise ValueError("invalid layout; try a different one")
        nevents = len(self.events)
        legend_kwargs = dict(fontsize=self.labelsize) # loc='lower center', mode='expand',
        nrows, ncols = self.get_row_column_numbers(layout, nevents)
        rows, cols = list(range(nrows)), list(range(ncols))
        for error_metric in error_metrics:
            _opt_handles, _opt_labels = [], []
            titles = []
            header = '{} Error-Space'.format(error_metric.title())
            fig = plt.figure(**kwargs)
            for idx, events in enumerate(self.events):
                ax = fig.add_subplot(nrows, ncols, idx+1, projection='3d')
                subheader = self.get_subheader(events)
                search_label = self.get_search_label(events)
                if search_label is None:
                    title = subheader
                else:
                    title = '{}\n{}'.format(subheader, search_label)
                if '\n' in title:
                    titles.append('{}\n'.format(title))
                else:
                    titles.append(title)
                ED, axis_limits, axis_labels, legend_labels, norm, fmt, cticks, fun = self.get_errorspace(events, distribution_model, error_metric)
                [xlabel, ylabel, zlabel] = axis_labels
                [mu_label, sigma_label, error_label] = legend_labels
                opt_label = '{}\n{}\n{}'.format(mu_label, sigma_label, error_label)
                _opt_labels.append(opt_label)
                ax.set_title(title, fontsize=self.labelsize, pad=10)
                if show_contours == True:
                    ax, contour_handles = self.subview_dim3_errorspace_contours(ax, ED, levels, cmap, central_color, norm)
                    [fill_handle, line_handle] = contour_handles
                if show_surface == True:
                    ax, surface_handle = self.subview_dim3_errorspace_surface(ax, ED, cmap, norm, rstride=rstride, cstride=cstride, shade=shade)
                optm_handle = ax.scatter(ED.prms[0], ED.prms[1], fun, color=central_color, marker='*', linewidth=0, s=50)
                _opt_handles.append(optm_handle)
                if include_colorbar == True:
                    try:
                        fig, ax, cbar = self.subview_errorspace_colorbar(fig, ax, layout, surface_handle, norm, fmt, clabel=zlabel, cticks=cticks)
                    except:
                        fig, ax, cbar = self.subview_errorspace_colorbar(fig, ax, layout, fill_handle, norm, fmt, clabel=zlabel, cticks=cticks)
                ax = self.set_viewing_position(ax, azim, elev, dist)
                ax.set_xlabel(xlabel, fontsize=self.labelsize)
                ax.set_ylabel(ylabel, fontsize=self.labelsize)
                ax.set_zlabel(zlabel, fontsize=self.labelsize)
                for axis in ('x', 'y', 'z'):
                    ax.tick_params(axis=axis, labelsize=self.ticksize)
            try:
                fig.subplots_adjust(hspace=0.2 * nevents, wspace=0.25, bottom=0.25)
            except:
                try:
                    fig.subplots_adjust(hspace=0.15 * nevents, wspace=0.25, bottom=0.25)
                except:
                    fig.subplots_adjust(bottom=0.25)
            fig.suptitle(header, fontsize=self.titlesize)
            self.subview_errorspace_legends(fig, nevents, _opt_handles, _opt_labels, titles, legend_kwargs)
            if save == True:
                suffix = ''
                if show_contours == True:
                    suffix += 'contours_'
                if show_surface == True:
                    suffix += 'surface_'
                acronym = self.get_acronym(error_metric)
                savename = 'speed_dim3_error_{}_{}_{}{}'.format(distribution_model, acronym, suffix, layout)
            else:
                savename = None
            self.display_image(fig, savename=savename)

##
