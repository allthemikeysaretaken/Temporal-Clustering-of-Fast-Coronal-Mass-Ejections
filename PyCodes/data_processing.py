from search_methods import *
import datetime
import numpy.core.defchararray as np_repl
import numpy_indexed as npi

class SolarCycle():

    """
    This class contains property attributes that can be used to select
    data parameters within a full or partial solar cycle.
    """

    def __init__(self):
        super().__init__()
        cycle_types = ('full', 'high-activity', 'early low-activity', 'late low-activity')
        cycle_numbers = np.arange(24).astype(int) + 1
        self._available = {'solar cycle' : cycle_numbers, 'cycle type' : cycle_types}

    @property
    def available(self):
        return self._available

    @property
    def sc1(self):
        return (datetime.datetime(1755, 2, 1, 0, 0, 0), datetime.datetime(1766, 6, 15, 0, 0, 0))

    @property
    def sc2(self):
        return (datetime.datetime(1766, 6, 15, 0, 0, 0), datetime.datetime(1775, 6, 15, 0, 0, 0))

    @property
    def sc3(self):
        return (datetime.datetime(1775, 6, 15, 0, 0, 0), datetime.datetime(1784, 9, 15, 0, 0, 0))

    @property
    def sc4(self):
        return (datetime.datetime(1784, 9, 15, 0, 0, 0), datetime.datetime(1798, 4, 15, 0, 0, 0))

    @property
    def sc5(self):
        return (datetime.datetime(1798, 4, 15, 0, 0, 0), datetime.datetime(1810, 8, 15, 0, 0, 0))

    @property
    def sc6(self):
        return (datetime.datetime(1810, 8, 15, 0, 0, 0), datetime.datetime(1823, 5, 15, 0, 0, 0))

    @property
    def sc7(self):
        return (datetime.datetime(1823, 5, 15, 0, 0, 0), datetime.datetime(1833, 11, 15, 0, 0, 0))

    @property
    def sc8(self):
        return (datetime.datetime(1833, 11, 15, 0, 0, 0), datetime.datetime(1843, 7, 15, 0, 0 ,0))

    @property
    def sc9(self):
        return (datetime.datetime(1843, 7, 15, 0, 0 ,0), datetime.datetime(1855, 12, 15, 0, 0, 0))

    @property
    def sc10(self):
        return (datetime.datetime(1855, 12, 15, 0, 0, 0), datetime.datetime(1867, 3, 15, 0, 0, 0))

    @property
    def sc11(self):
        return (datetime.datetime(1867, 3, 15, 0, 0, 0), datetime.datetime(1878, 12, 15, 0, 0, 0))

    @property
    def sc12(self):
        return (datetime.datetime(1878, 12, 15, 0, 0, 0), datetime.datetime(1890, 3, 15, 0, 0, 0))

    @property
    def sc13(self):
        return (datetime.datetime(1890, 3, 15, 0, 0, 0), datetime.datetime(1902, 1, 15, 0, 0, 0))

    @property
    def sc14(self):
        return (datetime.datetime(1902, 1, 15, 0, 0, 0), datetime.datetime(1913, 7, 15, 0, 0, 0))

    @property
    def sc15(self):
        return (datetime.datetime(1913, 7, 15, 0, 0, 0), datetime.datetime(1923, 8, 15, 0, 0, 0))

    @property
    def sc16(self):
        return (datetime.datetime(1923, 8, 15, 0, 0, 0), datetime.datetime(1933, 9, 15, 0, 0, 0))

    @property
    def sc17(self):
        return (datetime.datetime(1933, 9, 15, 0, 0, 0), datetime.datetime(1944, 2, 15, 0, 0, 0))

    @property
    def sc18(self):
        return (datetime.datetime(1944, 2, 15, 0, 0, 0), datetime.datetime(1954, 4, 15, 0, 0, 0))

    @property
    def sc19(self):
        return (datetime.datetime(1954, 4, 15, 0, 0, 0), datetime.datetime(1964, 10, 15, 0, 0, 0))

    @property
    def sc20(self):
        return (datetime.datetime(1964, 10, 15, 0, 0, 0), datetime.datetime(1976, 3, 15, 0, 0, 0))

    @property
    def sc21(self):
        return (datetime.datetime(1976, 3, 15, 0, 0, 0), datetime.datetime(1986, 9, 15, 0, 0, 0))

    @property
    def sc22(self):
        return (datetime.datetime(1986, 9, 15, 0, 0, 0), datetime.datetime(1996, 8, 15, 0, 0, 0))

    @property
    def sc23(self):
        res = {}
        res['full'] = (datetime.datetime(1996, 8, 15, 0, 0, 0), datetime.datetime(2008, 12, 15, 0, 0, 0))
        res['high-activity'] = (datetime.datetime(1999, 1, 1, 0, 0, 0), datetime.datetime(2006, 12, 31, 23, 59, 59))
        res['early low-activity'] = (datetime.datetime(1996, 8, 15, 0, 0, 0), datetime.datetime(1998, 12, 31, 23, 59, 59))
        res['late low-activity'] = (datetime.datetime(2007, 1, 1, 0, 0, 0), datetime.datetime(2008, 12, 15, 0, 0, 0))
        return res

    @property
    def sc24(self):
        res = {}
        res['full'] = (datetime.datetime(2008, 12, 15, 0, 0, 0), datetime.datetime(2020, 1, 15, 0, 0, 0))
        res['high-activity'] = (datetime.datetime(2010, 1, 1, 0, 0, 0), datetime.datetime(2016, 12, 31, 23, 59, 59))
        res['early low-activity'] = (datetime.datetime(2008, 12, 15, 0, 0, 0), datetime.datetime(2009, 12, 31, 23, 59, 59))
        res['late low-activity'] = (datetime.datetime(2017, 1, 1, 0, 0, 0), datetime.datetime(2020, 1, 15, 0, 0, 0))
        return res

    def match_datetimes(self, data, key='datetime'):
        """
        data:
            type <dict>

        key:
            type <str>
        """
        if key not in list(data.keys()):
            raise ValueError("cannot match datetimes if data does not contain '{}' key".format(key))
        dts = data[key]
        cycle_numbers = np.zeros(dts.size, dtype=int)
        nchars = len(max(self.available['cycle type'], key=len))
        cycle_types = np.array([' '*nchars] * dts.size, dtype=str)
        S = EventSearcher(data)
        parameters = (key, key)
        conditions = ('greater than or equal', 'less than or equal')
        for cycle_number in self.available['solar cycle']:
            cycle_bounds = getattr(self, 'sc{}'.format(cycle_number))
            if cycle_number in (23, 24):
                for cycle_type in self.available['cycle type']:
                    values = cycle_bounds[cycle_type]
                    try:
                        indices = S.search_events(parameters, conditions, values)[1]
                        cycle_numbers[indices] = cycle_number
                        cycle_types[indices] = cycle_type
                    except:
                        pass
            else:
                try:
                    indices = S.search_events(parameters, conditions, cycle_bounds)[1]
                    cycle_numbers[indices] = cycle_number
                    cycle_types[indices] = 'full'
                except:
                    pass
        return cycle_numbers, cycle_types

class DataParametrization(SolarCycle):

    """
    This class contains methods to filter and search data.
    """

    def __init__(self):
        super().__init__()

    @staticmethod
    def filter_nans_from_parameter(data, parameter, policy, repl):
        """
        data:
            type <dict>

        parameter:
            type <str>

        policy:
            type <str>

        repl:
            type <bool / str / int / float>
        """
        arr = data[parameter]
        if policy == 'discard':
            indices = np.where(~np.isnan(arr))[0]
            data = {key : values[indices] for key, values in data.items()}
        elif policy == 'replace':
            indices = np.where(np.isnan(arr))[0]
            arr[indices] = repl
            data[parameter] = arr
        else:
            raise ValueError("invalid policy: {}".format(policy))
        return data

    def filter_nans(self, data, parameter=None, policy=None, repl=None):
        """
        data:
            type <dict>

        parameter:
            type <str / tuple / list / array> or None

        policy:
            type <str / tuple / list / array> or None

        repl:
            type <bool / str / int / float / tuple / list / array> or None
        """
        result = dict(data)
        if parameter is not None:
            if policy is None:
                raise ValueError("nan_policies cannot be None if keys are specified")
            if repl is None:
                raise ValueError("nan_repls cannot be None if keys are specified")
            if isinstance(parameter, str):
                parameter = [parameter]
            if isinstance(policy, str):
                policy = [policy for key in parameter]
            if isinstance(repl, (bool, str, int, float)):
                repl = [repl for key in parameter]
            for nan_key, nan_policy, nan_repl in zip(parameter, policy, repl):
                result = self.filter_nans_from_parameter(result, nan_key, nan_policy, nan_repl)
        return result

    @staticmethod
    def reorganize_parameters(data, ref_parameters, special_parameters=None):
        """
        data:
            type <dict>

        ref_parameters:
            type <str / tuple / list / array>

        special_parameters:
            type <str / tuple / list / array> or None
        """
        if special_parameters is None:
            return data
        else:
            if isinstance(special_parameters, str):
                special_parameters = [special_parameters]
            elif not isinstance(special_parameters, (tuple, list, np.ndarray)):
                raise ValueError("invalid type(special_parameters): {}".format(type(special_parameters)))
            if isinstance(ref_parameters, str):
                ref_parameters = [ref_parameters]
            elif not isinstance(ref_parameters, (tuple, list, np.ndarray)):
                raise ValueError("invalid type(ref_parameters): {}".format(type(ref_parameters)))
            nrefs, nspecs = len(ref_parameters), len(special_parameters)
            if nrefs != nspecs:
                raise ValueError("{} reference parameters for {} special parameters")
            for ref_parameter, special_key in zip(ref_parameters, special_parameters):
                if ref_parameter in special_key:
                    vs = data[special_key]
                    result = {key : arr for key, arr in data.items() if ref_parameter not in key}
                    result[ref_parameter] = vs
            return result

class TemporalRegularization(DataParametrization):

    def __init__(self):
        super().__init__()
        units, denominators = ('second', 'minute', 'hour', 'day'), (1, 60, 3600, 3600*24)
        identifiers = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
        months = (np.arange(len(identifiers)).astype(int) + 1)
        _available = {'elapsed unit' : dict(zip(units, denominators)), 'months of year' : dict(zip(identifiers, months))}
        self._available.update(_available)

    @staticmethod
    def pad_datetimes(data, timestep):
        """
        data:
            type <dict>

        timestep:
            type <str>
        """
        if 'datetime' not in list(data.keys()):
            raise ValueError("data does not contain the key 'datetime'")
        dts = data['datetime']
        lbound, ubound = dts[0], dts[-1]
        kwargs = {'{}s'.format(timestep) : 1}
        step = datetime.timedelta(**kwargs)
        dts, dates, times = [], [], []
        years, months, days = [], [], []
        hours, minutes, seconds = [], [], []
        while lbound < ubound:
            _dt = lbound.strftime('%Y/%m/%d %H:%M:%S')
            dt = datetime.datetime(int(_dt[:4]), int(_dt[5:7]), int(_dt[8:10]), int(_dt[11:13]), int(_dt[14:16]), int(_dt[17:]))
            dts.append(lbound)
            dates.append(dt.date())
            times.append(dt.time())
            years.append(dt.year)
            months.append(dt.month)
            days.append(dt.day)
            hours.append(dt.hour)
            minutes.append(dt.minute)
            seconds.append(dt.second)
            lbound += step
        keys = ('datetime', 'date', 'time', 'year', 'month', 'day', 'hour', 'minute', 'second')
        values = (dts, dates, times, years, months, days, hours, minutes, seconds)
        res = {key : np.array(arr) for key, arr in zip(keys, values)}
        indices = data['elapsed']
        for key, arr in res.items():
            arr[indices] = data[key]
            res[key] = arr
        res['elapsed'] = np.arange(res['datetime'].size).astype(int)
        return res

    def get_elapsed_time(self, data, timestep):
        """
        data:
            type <dict>

        timestep:
            type <str>
        """
        if 'datetime' not in list(data.keys()):
            raise ValueError("data does not contain the key 'datetime'")
        if timestep not in list(self.available['elapsed unit'].keys()):
            raise ValueError("invalid 'elapsed unit': {}".format(timestep))
        dts = data['datetime']
        res = np.array([dt.total_seconds() for dt in (dts - dts[0])])
        return (res / self.available['elapsed unit'][timestep]).astype(int)

    def extract_datetimes(self, _dates, _times, _date_sep='/', _time_sep=':', integral_months=False):
        """
        _dates:
            type <array>

        _times:
            type <array>

        _date_sep:
            type <str>

        _time_sep:
            type <str>

        integral_months:
            type <bool>
        """
        dts, dates, times = [], [], []
        years, months, days = [], [], []
        hours, minutes, seconds = [], [], []
        if integral_months == True:
            for _date, _time in zip(_dates, _times):
                dchars = np.array(_date.split(sep=_date_sep), dtype=int)
                tchars = np.array(_time.split(sep=_time_sep), dtype=int)
                dt = datetime.datetime(*dchars, *tchars)
                dts.append(dt)
                dates.append(dt.date())
                times.append(dt.time())
                years.append(dt.year)
                months.append(dt.month)
                days.append(dt.day)
                hours.append(dt.hour)
                minutes.append(dt.minute)
                seconds.append(dt.second)
        else:
            for _date, _time in zip(_dates, _times):
                dchars = _date.split(sep=_date_sep)
                tchars = np.array(_time.split(sep=_time_sep), dtype=int)
                dchars[1] = self.available['months of year'][dchars[1]]
                dchars = [int(dchars[2]), int(dchars[1]), int(dchars[0])]
                dt = datetime.datetime(*dchars, *tchars)
                dts.append(dt)
                dates.append(dt.date())
                times.append(dt.time())
                years.append(dt.year)
                months.append(dt.month)
                days.append(dt.day)
                hours.append(dt.hour)
                minutes.append(dt.minute)
                seconds.append(dt.second)
        result = dict()
        temporal_keys = ('datetime', 'date', 'time', 'year', 'month', 'day', 'hour', 'minute', 'second')
        temporal_values = (dts, dates, times, years, months, days, hours, minutes, seconds)
        for key, value in zip(temporal_keys, temporal_values):
            result[key] = np.array(value)
        return result

    def extract_peak_parameters_from_datetimes(self, data, start_prefix='start ', peak_prefix='peak ', end_prefix='end '):
        """
        data:
            type <dict>

        start_prefix:
            type <str>

        peak_prefix:
            type <str>

        end_prefix:
            type <str>
        """
        result = {}
        pdts, pdates, ptimes = [], [], []
        pyears, pmonths, pdays = [], [], []
        phours, pminutes, pseconds = [], [], []
        for d_start, t_start, _t_peak, dt_end in zip(data['{}date'.format(start_prefix)], data['{}time'.format(start_prefix)], data['{}time'.format(peak_prefix)], data['{}datetime'.format(end_prefix)]):
            t_peak = datetime.time(*np.array(_t_peak.split(sep=':')).astype(int))
            if t_start > t_peak:
                pdt = datetime.datetime.combine(dt_end.date(), t_peak)
            else:
                pdt = datetime.datetime.combine(d_start, t_peak)
            pdts.append(pdt)
            pdates.append(pdt.date())
            ptimes.append(pdt.time())
            pyears.append(pdt.year)
            pmonths.append(pdt.month)
            pdays.append(pdt.day)
            phours.append(pdt.hour)
            pminutes.append(pdt.minute)
            pseconds.append(pdt.second)
        keys = ('datetime', 'date', 'time', 'year', 'month', 'day', 'hour', 'minute', 'second')
        arrs = (pdts, pdates, ptimes, pyears, pmonths, pdays, phours, pminutes, pseconds)
        for key, arr in zip(keys, arrs):
            result['{}{}'.format(peak_prefix, key)] = np.array(arr)
        key = '{}datetime'.format(peak_prefix)
        sc_cycles, sc_types = self.match_datetimes(result, key)
        result['{}solar cycle'.format(peak_prefix)] = sc_cycles
        result['{}cycle type'.format(peak_prefix)] = sc_types
        return result

    def get_end_parameters_from_duration(self, data, timestep, dur_key='duration', start_prefix='start ', end_prefix='end '):
        """
        data:
            type <dict>

        timestep:
            type <str>

        dur_key:
            type <str>

        start_prefix:
            type <str>

        end_prefix:
            type <str>
        """
        result = {}
        get_delta = lambda dur : datetime.timedelta(**{'{}s'.format(timestep) : float(dur)})
        duration_delta = np.array([get_delta(dur) for dur in data[dur_key]])
        edts = data['{}datetime'.format(start_prefix)] + duration_delta
        result.update({'{}{}'.format(end_prefix, key) : [] for key in ('date', 'time', 'year', 'month', 'day', 'hour', 'minute', 'second')})
        for dt in edts:
            result['{}date'.format(end_prefix)].append(dt.date())
            result['{}time'.format(end_prefix)].append(dt.time())
            result['{}year'.format(end_prefix)].append(dt.year)
            result['{}month'.format(end_prefix)].append(dt.month)
            result['{}day'.format(end_prefix)].append(dt.day)
            result['{}hour'.format(end_prefix)].append(dt.hour)
            result['{}minute'.format(end_prefix)].append(dt.minute)
            result['{}second'.format(end_prefix)].append(dt.second)
        result = {key : np.array(arr) for key, arr in result.items()}
        result['{}datetime'.format(end_prefix)] = edts
        sc_cycles, sc_types = self.match_datetimes(result, key='{}datetime'.format(end_prefix))
        result['{}solar cycle'.format(end_prefix)] = sc_cycles
        result['{}cycle type'.format(end_prefix)] = sc_types
        if '{}solar cycle'.format(start_prefix) not in list(data.keys()):
            sc_cycles, sc_types = self.match_datetimes(data, key='{}datetime'.format(start_prefix))
            result['{}solar cycle'.format(start_prefix)] = sc_cycles
            result['{}cycle type'.format(start_prefix)] = sc_types
        return result

    def group_regular_parameters(self, time_parameters, other_parameters, timestep, indices):
        """
        time_parameters:
            type <dict>

        other_parameters:
            type <dict>

        timestep:
            type <str>

        indices:
            type <array>
        """
        exclusionary_keys = list(time_parameters.keys()) + ['solar cycle', 'cycle type']
        result = {}
        for key, arr in other_parameters.items():
            if key not in exclusionary_keys:
                mask = None
                if key == 'speed':
                    base = np.ones(time_parameters['year'].size, dtype=int)
                else:
                    if arr.dtype in (bool, np.bool_, float, np.float, int, np.int, np.int64):
                        base = np.zeros(time_parameters['year'].size, dtype=arr.dtype)
                    elif type(arr[0]) in (str, np.str_):
                        n = len(max(arr, key=len))
                        base = np.array([' '*n] * time_parameters['year'].size, dtype=str)
                        mask = np.ones(base.size, dtype=int)
                    else:
                        raise ValueError("invalid dtype: {}".format(arr.dtype))
                base[indices] = arr
                if mask is not None:
                    mask[indices] = 0
                    base[np.where(mask.astype(bool))[0]] = 'none'
                result[key] = base
        result.update(time_parameters)
        solar_cycles, cycle_types = self.match_datetimes(time_parameters)
        is_event = np.zeros(time_parameters['elapsed'].size, dtype=bool)
        is_event[indices] = True
        result['solar cycle'] = solar_cycles
        result['cycle type'] = cycle_types
        result['is event'] = is_event
        return result

    @staticmethod
    def get_event_counts_per_time(data, timestep, dt_prefix='', group_index=0):
        """
        data:
            type <dict>

        search_kwargs:
            type <dict> or None

        dt_prefix:
            type <str>

        group_index:
            type <int>
        """
        if timestep == 'hour':
            primary_condition = (np.diff(data['{}hour'.format(dt_prefix)]) != 0)
            secondary_condition = (np.diff(data['{}day'.format(dt_prefix)]) != 0)
            tertiary_condition = (np.diff(data['{}month'.format(dt_prefix)]) != 0)
            indices = np.where(primary_condition | secondary_condition | tertiary_condition)[0] + 1
        elif timestep == 'day':
            primary_condition = (np.diff(data['{}day'.format(dt_prefix)]) != 0)
            secondary_condition = (np.diff(data['{}month'.format(dt_prefix)]) != 0)
            indices = np.where(primary_condition | secondary_condition)[0] + 1
        elif timestep == 'month':
            primary_condition = (np.diff(data['{}month'.format(dt_prefix)]) != 0)
            indices = np.where(primary_condition)[0] + 1
        elif timestep == 'year':
            primary_condition = (np.diff(data['{}year'.format(dt_prefix)]) != 0)
            indices = np.where(primary_condition)[0] + 1
        else:
            raise ValueError("invalid timestep: {}".format(timestep))
        parameters = {key : [] for key in list(data.keys())}
        counts = []
        arbitrary_key = list(data.keys())[0]
        for key, arr in data.items():
            groups = np.split(arr, indices)
            for group in groups:
                parameters[key].append(group[group_index])
                if key == arbitrary_key:
                    counts.append(len(group))
        nchars = len(dt_prefix)
        result = {}
        for key, arr in parameters.items():
            if key[:nchars] == dt_prefix:
                result[key[nchars:]] = np.array(arr)
        result['count'] = np.array(counts)
        return result

    def get_silso_sunspot_counts_per_time(self, data, timestep):
        """
        data:
            type <dict>

        timestep:
            type <str>
        """
        if timestep == 'day':
            result = dict(data)
        elif timestep in ('month', 'year'):
            condition = np.diff(data[timestep])
            indices = np.where(condition)[0] + 1
            groups = {key : np.split(values, indices) for key, values in data.items()}
            dts, dates, times, years, months, days, counts = [], [], [], [], [], [], []
            for dt, yy, mm, dd, cc in zip(groups['datetime'], groups['year'], groups['month'], groups['day'], groups['count']):
                _dt = dt[0]
                dts.append(_dt)
                dates.append(_dt.date())
                times.append(_dt.time())
                years.append(yy[0])
                months.append(mm[0])
                days.append(dd[0])
                counts.append(np.sum(cc))
            keys = ('datetime', 'date', 'time', 'year', 'month', 'day', 'count')
            values = (dts, dates, times, years, months, days, counts)
            result = {key : np.array(values) for key, values in zip(keys, values)}
            result['source'] = np.array(['SILSO'] * result['datetime'].size)
            solar_cycles, cycle_types = self.match_datetimes(result, key='datetime')
            result['solar cycle'] = solar_cycles
            result['cycle type'] = cycle_types
        else:
            raise ValueError("invalid timestep: {}".format(timestep))
        return result

class DataBase(TemporalRegularization):

    """
    This class contains methods to load solar data.
    """

    def __init__(self, directory):
        """
        directory:
            type <str>
        """
        super().__init__()
        self.directory = directory
        self._cme = None
        self._sunspot = None
        self._flare = None
        ev_types = ('cme', 'flare', 'sunspot')
        soho_lasco_speed_types = ('linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed')
        _available = {'event type' : ev_types, 'soho lasco cme speed type' : soho_lasco_speed_types}
        self._available.update(_available)

    @property
    def cme(self):
        return self._cme

    @property
    def sunspot(self):
        return self._sunspot

    @property
    def flare(self):
        return self._flare

    def load_cmes(self, read_soho_lasco=False):
        """
        read_soho_lasco:
            type <bool>
        """
        if self._cme is not None:
            raise ValueError("'cme' data is already initialized")
        self._cme = {}
        if read_soho_lasco == True:
            filepath = '{}univ_all.txt'.format(self.directory)
            keys = ['date', 'time', 'central position angle', 'angular width']
            keys += list(self.available['soho lasco cme speed type'])
            keys += ['acceleration', 'mass', 'kinetic energy', 'mean position angle']
            data = np.loadtxt(filepath, usecols=np.arange(len(keys)).astype(int), skiprows=4, dtype=str)
            result = {}
            for col, key in enumerate(keys):
                if key in ('date', 'time'):
                    result[key] = data[:, col]
                elif key in self.available['soho lasco cme speed type']:
                    result[key] = np_repl.replace(data[:, col], '----', 'nan').astype(float)
                elif key in ('mass', 'kinetic energy'):
                    tmp = np_repl.replace(data[:, col], '*'*8, 'nan').astype(str)
                    tmp = np_repl.replace(tmp, '*', '').astype(str)
                    result[key] = np_repl.replace(tmp, '-------', 'nan').astype(float)
                elif key == 'acceleration':
                    tmp = np_repl.replace(data[:, col], '------', 'nan').astype(str)
                    result[key] = np_repl.replace(tmp, '*', '').astype(float)
                elif key in ('angular width', 'mean position angle'):
                    result[key] = np.array(data[:, col]).astype(int)
                elif key == 'central position angle':
                    result[key] = np_repl.replace(data[:, col], 'Halo', 'nan').astype(float)
                else:
                    raise ValueError("invalid key: {}".format(key))
            time_parameters = self.extract_datetimes(result['date'], result['time'], _date_sep='/', _time_sep=':', integral_months=True)
            result.update(time_parameters)
            solar_cycles, cycle_types = self.match_datetimes(result, key='datetime')
            result['solar cycle'] = solar_cycles
            result['cycle type'] = cycle_types
            result['is halo'] = np.isnan(result['central position angle'])
            result['source'] = np.array(['SOHO LASCO'] * result['datetime'].size)
            result['event type'] = np.array(['cme'] * result['datetime'].size)
            self._cme.update(result)

    def load_sunspots(self, read_silso=False):
        """
        read_silso:
            type <bool>
        """
        if self._sunspot is not None:
            raise ValueError("'sunspot' data is already initialized")
        self._sunspot = {}
        if read_silso == True:
            filepath = '{}SN_d_tot_V2.0.txt'.format(self.directory)
            keys, cols = ('year', 'month', 'day', 'count'), (0, 1, 2, 4)
            data = np.loadtxt(filepath, usecols=cols, dtype=int)
            result = {}
            for col, key in enumerate(keys):
                if key == 'count':
                    counts = np.copy(data[:, col])
                    counts[counts < 0] = 0
                    result[key] = counts
                else:
                    result[key] = data[:, col]
            result['datetime'] = np.array([datetime.datetime(yy, mm, dd, 0, 0, 0) for yy, mm, dd in zip(result['year'], result['month'], result['day'])])
            solar_cycles, cycle_types = self.match_datetimes(result, key='datetime')
            result['solar cycle'] = solar_cycles
            result['cycle type'] = cycle_types
            result['source'] = np.array(['SILSO'] * result['datetime'].size)
            result['event type'] = np.array(['sunspot'] * result['datetime'].size)
            self._sunspot.update(result)

    def load_flares(self, read_hessi=False):
        """
        read_hessi:
            type <bool>
        """
        if self._flare is not None:
            raise ValueError("'flare' data is already initialized")
        self._flare = {}
        if read_hessi == True:
            filepath = '{}hessi_flare_list.txt'.format(self.directory)
            keys = ('id', 'start date', 'start time', 'peak time', 'end time', 'duration', 'counts per second', 'total counts', 'energy', 'x position', 'y position', 'radial position', 'active region')
            data = np.genfromtxt(filepath, usecols=np.arange(len(keys)).astype(int), skip_header=7, skip_footer=38, dtype=str) # encoding='latin-1',
            result, energies = {}, []
            for col, key in enumerate(keys):
                if key in ('counts per second', 'total counts', 'x position', 'y position', 'radial position', 'duration'):
                    result[key] = np.array(data[:, col], dtype=int)
                elif key == 'energy':
                    for arg in data[:, col]:
                        energy = np.array(arg.split(sep='-'), dtype=int)
                        energies.append(energy)
                else:
                    result[key] = np.array(data[:, col])
            energies = np.array(energies)
            result['minimum energy'] = energies[:, 0]
            result['maximum energy'] = energies[:, 1]
            start_parameters = self.extract_datetimes(result['start date'], result['start time'], _date_sep='-', _time_sep=':', integral_months=False)
            start_parameters = {'start {}'.format(key) : arr for key, arr in start_parameters.items()}
            result.update(start_parameters)
            end_parameters = self.get_end_parameters_from_duration(result, timestep='second', dur_key='duration', start_prefix='start ', end_prefix='end ')
            result.update(end_parameters)
            peak_parameters = self.extract_peak_parameters_from_datetimes(result, start_prefix='start ', peak_prefix='peak ', end_prefix='end ')
            result.update(peak_parameters)
            result['source'] = np.array(['RHESSI'] * result['id'].size)
            result['event type'] = np.array(['flare'] * result['id'].size)
            self._flare.update(result)

    def load_data(self, read_soho_lasco=False, read_silso=False, read_hessi=False):
        """
        read_soho_lasco:
            type <bool>

        read_silso:
            type <bool>

        read_hessi:
            type <bool>
        """
        self.load_cmes(read_soho_lasco)
        self.load_sunspots(read_silso)
        self.load_flares(read_hessi)

    def verify_event_type(self, event_type):
        """
        event_type:
            type <str>
        """
        if event_type not in self.available['event type']:
            raise ValueError("invalid event_type: {}".format(event_type))

##
