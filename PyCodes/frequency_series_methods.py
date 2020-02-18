from data_processing import *
from optimization_methods import *
from visual_configuration import *

class FrequencySeriesViewer(VisualConfiguration):

    """
    This class is inherited by all classes pertaining to
    timeseries methods. This class inherits methods from
    `VisualConfiguration` to be used as convenience
    functions for plotting.
    """

    def __init__(self, timestep, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize):
        """
        timestep:
            type <str>

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

    @staticmethod
    def get_event_type_label(events):
        """
        event_type:
            type <dict>
        """
        identifiers = events['identifier']
        event_type = identifiers['event type']
        if event_type == 'cme':
            result = '{}'.format(event_type.upper())
        else:
            result = '{}'.format(event_type.title())
        return result

    @staticmethod
    def get_event_total_label(events):
        """
        events:
            type <dict>
        """
        parameters = events['parameter']
        return '${:,}$'.format(np.sum(parameters['count']))

    def get_search_label(self, search_kwargs, event_type):
        """

        """
        result = ''
        if search_kwargs is None:
            return result
        else:
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

    def get_legend_label(self, events, include_search_kwargs=False):
        """

        """
        event_total_label = self.get_event_total_label(events)
        event_type_label = self.get_event_type_label(events)
        label = '{} {}s'.format(event_total_label, event_type_label)
        if include_search_kwargs == True:
            identifiers = events['identifier']
            search_kwargs, event_type = identifiers['search kwargs'], identifiers['event type']
            search_label = self.get_search_label(search_kwargs, event_type)
            label = '{}\n{}'.format(label, search_label)
        return label

    def get_ylabel(self, events):
        """

        """
        event_type_label = self.get_event_type_label(events)
        identifiers = events['identifier']
        search_kwargs, event_type = identifiers['search kwargs'], identifiers['event type']
        search_label = self.get_search_label(search_kwargs, event_type)
        return '{}s\n{}'.format(event_type_label, search_label)

class FrequencySeries(FrequencySeriesViewer):

    def __init__(self, DB, timestep, directory, ticksize=7, labelsize=8, textsize=9, titlesize=11, headersize=20, cellsize=15):
        """
        DB:
            type <custom class>

        timestep:
            type <str>

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
        super().__init__(timestep, directory, ticksize, labelsize, textsize, titlesize, headersize, cellsize)
        self.DB = DB
        self._events = []
        self._available.update(dict(DB.available))

    @property
    def events(self):
        return self._events

    def load_events(self, event_type, search_parameters=None, search_conditions=None, search_values=None, apply_to='all', modifiers=None, dt_prefix='', group_index=0):
        """
        event_type:
            type <str>

        ...

        dt_prefix:
            type <str>

        group_index:
            type <int>
        """
        self.DB.verify_event_type(event_type)
        data = getattr(self.DB, event_type)
        if data is None:
            raise ValueError("DB has not initialized event_type '{}'".format(event_type))
        if (search_parameters is None) and (search_conditions is None) and (search_values is None):
            search_kwargs = None
        else:
            search_kwargs = dict()
            search_kwargs['parameters'] = search_parameters
            search_kwargs['conditions'] = search_conditions
            search_kwargs['values'] = search_values
            search_kwargs['apply_to'] = apply_to
            search_kwargs['modifiers'] = modifiers
        S = EventSearcher(data)
        parameters = S.search_events(**search_kwargs)[0]
        if event_type == 'sunspot':
            parameters = self.DB.get_silso_sunspot_counts_per_time(parameters, self.timestep)
        else:
            parameters = self.DB.get_event_counts_per_time(parameters, self.timestep, dt_prefix, group_index)
        identifiers = dict()
        identifiers['search kwargs'] = search_kwargs
        identifiers['event type'] = event_type
        result = {'parameter' : parameters, 'identifier' : identifiers}
        self._events.append(result)

    def get_chronological_solar_cycles(self):
        solar_cycles = []
        for events in self.events:
            parameters = events['parameter']
            for key in list(parameters.keys()):
                if 'solar cycle' in key:
                    solar_cycles += parameters[key].tolist()
        return np.unique(solar_cycles)

    def subview_frequency_distributions(self, ax, events, facecolor, major, minor, fmt, rotation, label=None):
        """
        ax:
            type <matplotlib object>

        events:
            type <dict>

        facecolor:
            type <str>

        major:
            type <str>

        minor:
            type <str>

        fmt:
            type <str>

        rotation:
            type <int / float>

        label:
            type <str> or None
        """
        parameters, identifiers = events['parameter'], events['identifier']
        if self.timestep == 'year':
            handle = ax.scatter(parameters['datetime'], parameters['count'], color=facecolor, marker='x', label=label)
        else:
            handle = ax.plot(parameters['datetime'], parameters['count'], color=facecolor, linestyle='-', linewidth=1, label=label)
        ax = self.transform_x_as_datetime(ax, major, minor, fmt, rotation)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        return ax, handle

    def subview_solar_cycle_separations(self, axes, solar_cycles, line_color, text_color):
        """
        axes:
            type <matplotlib object / array>

        solar_cycles:
            type <array>

        line_color:
            type <str>

        text_color:
            type <str>
        """
        arrowprops = dict(facecolor='black', arrowstyle="->")
        f = lambda cycle : list(cycle['full']) if isinstance(cycle, dict) else list(cycle)
        if solar_cycles.size > 1:
            if np.all(np.diff(solar_cycles) <= 1) == True:
                for icycle, jcycle in zip(solar_cycles[:-1], solar_cycles[1:]):
                    _icycle = getattr(self.DB, 'sc{}'.format(icycle))
                    _jcycle = getattr(self.DB, 'sc{}'.format(jcycle))
                    prev_cycle, curr_cycle = f(_icycle), f(_jcycle)
                    xprev = (prev_cycle[-1] - prev_cycle[0]) * 0.75 + prev_cycle[0]
                    xcurr = (curr_cycle[-1] - curr_cycle[0]) * 0.25 + curr_cycle[0]
                    xloc = (curr_cycle[0] - prev_cycle[-1]) / 2 + prev_cycle[-1]
                    try:
                        for idx, ax in enumerate(axes.ravel()):
                            ylim = ax.get_ylim()
                            ax.axvline(xloc, ymin=0, ymax=ylim[-1], color=line_color, linestyle='--')
                            if idx == 0:
                                yi, yj = 0.8 * ylim[-1], 0.6 * ylim[-1]
                                ax.annotate(' ', xy=(prev_cycle[0], yj), xytext=(prev_cycle[-1], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                                ax.annotate(' ', xy=(curr_cycle[-1], yj), xytext=(curr_cycle[0], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                                ax.text(xprev, yi, 'SC ${}$'.format(icycle), horizontalalignment='center', verticalalignment='center', fontsize=self.textsize)
                                ax.text(xcurr, yi, 'SC ${}$'.format(jcycle), horizontalalignment='center', verticalalignment='center', fontsize=self.textsize)
                    except:
                        ylim = axes.get_ylim()
                        axes.axvline(xloc, ymin=ylim[0], ymax=ylim[-1], color=line_color, linestyle='--')
                        yi, yj = 0.8 * ylim[-1], 0.6 * ylim[-1]
                        axes.annotate(' ', xy=(prev_cycle[0], yj), xytext=(prev_cycle[-1], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                        axes.annotate(' ', xy=(curr_cycle[-1], yj), xytext=(curr_cycle[0], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                        axes.text(xprev, yi, 'SC ${}$'.format(icycle), horizontalalignment='center', verticalalignment='center', fontsize=self.textsize)
                        axes.text(xcurr, yi, 'SC ${}$'.format(jcycle), horizontalalignment='center', verticalalignment='center', fontsize=self.textsize)
            else:
                for solar_cycle in solar_cycles:
                    _cycle = getattr(self.DB, 'sc{}'.format(solar_cycle))
                    cycle_bounds = f(_cycle)
                    xloc = (cycle_bounds[0] - cycle_bounds[-1]) / 2 + cycle_bounds[0]
                    try:
                        ax = axes.ravel()[0]
                        ylim = ax.get_ylim()
                        yi, yj = 0.8 * ylim[-1], 0.6 * ylim[-1]
                        ax.annotate(' ', xy=(cycle_bounds[0], yj), xytext=(cycle_bounds[-1], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                        ax.text(xloc, yi, 'SC ${}$'.format(solar_cycle), horizontalalignment='center', verticalalignment='center', fontsize=self.textsize)
                    except:
                        ylim = axes.get_ylim()
                        yi, yj = 0.8 * ylim[-1], 0.6 * ylim[-1]
                        axes.annotate(' ', xy=(cycle_bounds[0], yj), xytext=(cycle_bounds[-1], yj), horizontalalignment='center', verticalalignment='center', arrowprops=arrowprops)
                        axes.text(xloc, yi, 'SC ${}$'.format(solar_cycle), horizontalalignment='center', verticalalignment='center', fontsize=self.textsize)
        return axes

    def subview_solar_cycle_legends(self, fig, solar_cycles, legend_kwargs):
        """

        """
        if solar_cycles.size == 1:
            leg = fig.legend(**legend_kwargs)
            leg.set_title('SC ${}$'.format(solar_cycles[0]), prop={'size': self.labelsize})
            leg._legend_box.align = "center"
            frame = leg.get_frame()
            frame.set_edgecolor('k')
        elif solar_cycles.size > 3:
            raise ValueError("not yet implemented for more than 3 solar cycles; {} were provided".format(solar_cycles.size))
        else:
            left_labels, mid_labels, right_labels = [], [], []
            for idx, events in enumerate(self.events):
                event_type_label = self.get_event_type_label(events)
                parameters = events['parameter']
                S = EventSearcher(parameters)
                for ith_cycle, solar_cycle in enumerate(solar_cycles):
                    _ps = S.search_events(parameters='solar cycle', conditions='equal', values=solar_cycle)[0]
                    _label = '${:,}$ {}s'.format(np.sum(_ps['count']), event_type_label)
                    if ith_cycle == 0:
                        left_labels.append(_label)
                    elif ith_cycle == 1:
                        if solar_cycles.size == 2:
                            right_labels.append(_label)
                        else:
                            mid_labels.append(_label)
                    else:
                        right_labels.append(_label)
            if solar_cycles.size == 2:
                w, h = 0.3, 0.2
                leg_left = fig.legend(labels=left_labels, bbox_to_anchor=(0, 0, w, h), loc='lower left', **legend_kwargs)
                leg_right = fig.legend(labels=right_labels, bbox_to_anchor=(0.75, 0, w, h), loc='lower right', **legend_kwargs)
                for leg, solar_cycle in zip((leg_left, leg_right), solar_cycles):
                    leg.set_title('SC ${}$'.format(solar_cycle), prop={'size': self.labelsize})
                    leg._legend_box.align = "center"
                    frame = leg.get_frame()
                    frame.set_edgecolor('k')
            elif solar_cycles.size == 3:
                w, h = 0.275, 0.2
                leg_left = fig.legend(labels=left_labels, bbox_to_anchor=(0, 0, w, h), loc='lower left', **legend_kwargs)
                leg_mid = fig.legend(labels=mid_labels, bbox_to_anchor=(0.375, 0, w, h), loc='lower center', **legend_kwargs)
                leg_right = fig.legend(labels=right_labels, bbox_to_anchor=(0.75, 0, w, h), loc='lower right', **legend_kwargs)
                for leg, solar_cycle in zip((leg_left, leg_mid, leg_right), solar_cycles):
                    leg.set_title('SC ${}$'.format(solar_cycle), prop={'size': self.labelsize})
                    leg._legend_box.align = "center"
                    frame = leg.get_frame()
                    frame.set_edgecolor('k')

    def view_frequency_distributions(self, layout, facecolors=('darkorange', 'steelblue', 'green', 'purple'), major='year', minor='month', fmt="%Y-%m-%d", rotation=15, show_solar_cycle_separations=False, show_solar_cycle_legends=False, save=False, **kwargs):
        """
        layout:
            type <str>

        facecolors:
            type <str / tuple / list / array>

        major:
            type <str>

        minor:
            type <str>

        fmt:
            type <str>

        rotation:
            type <int / float>

        show_cycle_labels:
            type <bool>

        show_cycle_legends:
            type <bool>

        save:
            type <bool>
        """
        nevents = len(self.events)
        if isinstance(facecolors, str):
            facecolors = [facecolors]
        ncolors = len(facecolors)
        if ncolors < nevents:
            raise ValueError("{} facecolors for {} events".format(ncolors, nevents))
        if layout not in ('overlay', 'vertical'):
            raise ValueError("invalid layout '{}' for this figure".format(layout))
        solar_cycles = self.get_chronological_solar_cycles()
        legend_kwargs = dict(fontsize=self.labelsize)
        legend_kwargs['ncol'] = nevents //2 if nevents > 3 else nevents
        header = '{} Distribution of Solar Events'.format(self.timestep_label.title())
        _handles, _labels = [], []
        fig, axes = self.get_dim2_figure_and_axes(layout, nevents, **kwargs)
        try:
            for idx, (ax, events, facecolor) in enumerate(zip(axes.ravel(), self.events, facecolors)):
                _label = self.get_legend_label(events, include_search_kwargs=False)
                ylabel = self.get_ylabel(events)
                ax, _handle = self.subview_frequency_distributions(ax, events, facecolor, major, minor, fmt, rotation, _label)
                _handles.append(_handle)
                _labels.append(_label)
                ax.set_ylabel(ylabel, fontsize=self.labelsize)
                ax.tick_params(axis='x', colors=facecolor, which='both')
                ax.tick_params(axis='y', colors=facecolor, which='both')
                ax.grid(color='k', linestyle=':', alpha=0.3)
                ax.xaxis.label.set_color(facecolor)
                ax.yaxis.label.set_color(facecolor)
                if idx == 0:
                    ax.tick_params(which='both', top=False, bottom=True, labeltop=False, left=True, right=True, labelleft=True, labelright=False, labelsize=self.ticksize) # labelbottom=False,
                else:
                    ax.tick_params(which='both', top=True, bottom=True, labeltop=False, left=True, right=True, labelleft=True, labelright=False, labelsize=self.ticksize) # labelbottom=False,
            axes.ravel()[-1].set_xlabel('Date', fontsize=self.labelsize)
        except:
            for events, facecolor in zip(self.events, facecolors):
                _label = self.get_legend_label(events, include_search_kwargs=True)
                axes, _handle = self.subview_frequency_distributions(axes, events, facecolor, major, minor, fmt, rotation, _label)
                _handles.append(_handle)
                _labels.append(_label)
            axes.set_xlabel('Date', fontsize=self.labelsize)
            axes.set_ylabel('Event Frequency', fontsize=self.labelsize)
            axes.tick_params(axis='both', which='both', labelsize=self.ticksize)
            axes.grid(color='k', linestyle=':', alpha=0.3)
        fig.suptitle(header, fontsize=self.titlesize)
        fig.align_ylabels()
        fig.subplots_adjust(bottom=0.2, hspace=0.3)
        suffix = ''
        if show_solar_cycle_legends == True:
            self.subview_solar_cycle_legends(fig, solar_cycles, legend_kwargs)
            suffix += '_leg'
        else:
            leg = fig.legend(loc='lower center', mode='expand', **legend_kwargs)
            leg._legend_box.align = "center"
            frame = leg.get_frame()
            frame.set_edgecolor('k')
            # frame.set_facecolor('green')
        if show_solar_cycle_separations == True:
            axes = self.subview_solar_cycle_separations(axes, solar_cycles, line_color='r', text_color='k')
            suffix += '_sep'
        if save == True:
            savename = 'frequency_distribution_{}_{}{}'.format(self.timestep_label, layout, suffix)
        else:
            savename = None
        self.display_image(fig, savename=savename)
