import numpy as np
import numpy.core.defchararray as np_repl
import re
import datetime
import numpy_indexed as npi

from search_methods import *

class SolarCycle():

    def __init__(self):
        self.subcycles = (None, 'pre-maximum', 'post-maximum', 'high-activity')
        self.search_parameters = tuple(['datetime object']*3 + ['year'])
        self.search_conditions = [['greater than or equal', 'less than or equal'] for subcycle in self.subcycles]

    @property
    def sc_23(self):
        """ """
        res = {}
        search_values = [(datetime.datetime(1996, 8, 1, 0, 0, 0), datetime.datetime(2008, 12, 31, 23, 0, 0))]
        search_values += [(datetime.datetime(1996, 8, 1, 0), datetime.datetime(1998, 12, 31, 23))]
        search_values += [(datetime.datetime(2007, 1, 1, 0), datetime.datetime(2008, 12, 31, 23))]
        search_values += [(1999, 2006)]
        for subcycle, search_parameter, search_condition, search_value in zip(self.subcycles, self.search_parameters, self.search_conditions, search_values):
            res[subcycle] = {'search_parameters' : search_parameter, 'search_conditions' : search_condition, 'search_values' : search_value}
        return res

    @property
    def sc_24(self):
        """ """
        res = {}
        search_values = [(datetime.datetime(2009, 1, 1, 0), datetime.datetime(2017, 10, 31, 23))]
        search_values += [(datetime.datetime(2009, 1, 1, 0), datetime.datetime(2016, 12, 31, 23))]
        search_values += [(datetime.datetime(2017, 1, 1, 0), datetime.datetime(2017, 10, 31, 23))]
        search_values += [2010, 2016]
        for subcycle, search_parameter, search_condition, search_value in zip(self.subcycles, self.search_parameters, self.search_conditions, search_values):
            res[subcycle] = {'search_parameters' : search_parameter, 'search_conditions' : search_condition, 'search_values' : search_value}
        return res

    @staticmethod
    def get_temporal_search_inputs(years, months=None, days=None, dt_objects=None):
        """ """
        search_parameters = []
        search_conditions = []
        search_values = []
        if dt_objects is not None:
            pass
        else:
            if years is not None:
                if len(years) != 2:
                    raise ValueError("years = (initial_year, final_year)")
                search_parameters += ['year', 'year']
                search_conditions += ['greater than or equal', 'less than or equal']
                if isinstance(years, np.ndarray):
                    search_values += years.tolist()
                else:
                    search_values += list(years)
            if months is not None:
                if len(months) != 2:
                    raise ValueError("months = (initial_month, final_month)")
                search_parameters += ['month', 'month']
                search_conditions += ['greater than or equal', 'less than or equal']
                if isinstance(months, np.ndarray):
                    search_values += months.tolist()
                else:
                    search_values += list(months)
            if days is not None:
                if len(days) != 2:
                    raise ValueError("days = (initial_day, final_day)")
                search_parameters += ['day', 'day']
                search_conditions += ['greater than or equal', 'less than or equal']
                if isinstance(days, np.ndarray):
                    search_values += days.tolist()
                else:
                    search_values += list(days)
        return search_parameters, search_conditions, search_values

    def get_search_kwargs(self, solar_cycle=None, subcycle=None, years=None, months=None, days=None, dt_objects=None):
        """ """
        if solar_cycle in (23, '23'):
            return self.sc_23[subcycle]
        elif solar_cycle in (24, '24'):
            return self.sc_24[subcycle]
        else:
            if solar_cycle is not None:
                raise ValueError("solar_cycle = 23, 24, or None")
            search_parameters, search_conditions, search_values = self.get_temporal_search_inputs(years, months, days, dt_objects)
            return {'search_parameters' : search_parameters, 'search_conditions' : search_conditions, 'search_values' : search_values}

class SohoLascoCharacterMap(SolarCycle):

    def __init__(self):
        super().__init__()
        self.speed_types = ('linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed')

    @property
    def speed_type_mapping(self):
        res = {}
        res['linear speed'] = r'$V_{linear}$'
        res['second order initial speed'] = r'$V_{20 R_{\odot}, i}$'
        res['second order final speed'] = r'$V_{20 R_{\odot}, f}$'
        res['second order 20R speed'] = r'$V_{20 R_{\odot}}$'
        return res

    @property
    def temporal_mapping(self):
        res = {}
        res['23'] = r'SC $23$'
        res['24'] = r'SC $24$'
        res[23] = res['23']
        res[24] = res['24']
        res[None] = ''
        return res

    @staticmethod
    def modify_characters(characters, prev_char, repl_char, dtype):
        """
        characters  :   type <array>
        prev_char   :   type <str> or None
        repl_char   :   type <str>
        dtype       :   float, int, or str
        """
        if prev_char is None:
            return characters.astype(dtype)
        else:
            return np_repl.replace(characters, prev_char, repl_char).astype(dtype)

    def autocorrect_speed(self, args):
        """ """
        return self.modify_characters(args, '----', 'nan', dtype=float)

    def autocorrect_acceleration(self, args):
        """ """
        mod_values = self.modify_characters(args, '------', 'nan', dtype=str)
        return self.modify_characters(mod_values, '*', '', dtype=float)

    def autocorrect_physicals(self, args):
        """ """
        mod_values = self.modify_characters(args, '*'*8, 'nan', dtype=str)
        mod_values = self.modify_characters(mod_values, '*', '', dtype=str)
        return self.modify_characters(mod_values, '-------', 'nan', dtype=float)

    @property
    def parse_parameters(self):
        """ """
        res = {}
        res['central position angle'] = lambda args : self.modify_characters(args, 'Halo', 'nan', dtype=float)
        res['angular width'] = lambda args : args.astype(int)
        res['linear speed'] = self.autocorrect_speed
        res['second order initial speed'] = self.autocorrect_speed
        res['second order final speed'] = self.autocorrect_speed
        res['second order 20R speed'] = self.autocorrect_speed
        res['acceleration'] = self.autocorrect_acceleration
        res['mass'] = self.autocorrect_physicals
        res['kinetic energy'] = self.autocorrect_physicals
        res['mean position angle'] = lambda args : args.astype(int)
        return res

    @staticmethod
    def parse_dates_times(dates, times):
        """ """
        keys = ('year', 'month', 'day', 'hour', 'minute', 'second', 'datetime object')
        dtimes = {key : [] for key in keys}
        for date, time in zip(dates, times):
            year = int(date[:4])
            month = int(date[5:7])
            day = int(date[8:10])
            hour = int(time[:2])
            minute = int(time[3:5])
            second = int(time[6:])
            dtimes['year'].append(year)
            dtimes['month'].append(month)
            dtimes['day'].append(day)
            dtimes['hour'].append(hour)
            dtimes['minute'].append(minute)
            dtimes['second'].append(second)
            dt_object = datetime.datetime(year, month, day, hour, minute, second)
            dtimes['datetime object'].append(dt_object)
        return {key : np.array(value) for key, value in dtimes.items()}

    @staticmethod
    def parse_halos(central_position_angles):
        """ """
        base = np.zeros(central_position_angles.size, dtype=bool)
        indices = np.where(np.isnan(central_position_angles))[0]
        base[indices] = True
        return base

    @staticmethod
    def get_elapsed_hours(dt_objects):
        """ """
        return np.array([(dt_objects[idx] - dt_objects[0]).total_seconds() for idx in range(dt_objects.size)]) / 3600

    def parse_data_by_speeds(self, data, speed_type, nan_policy):
        """ """
        if speed_type not in self.speed_types:
            raise ValueError("unknown speed_type: {}; available speed_types = {}".format(speed_type, self.speed_types))
        if nan_policy == 'discard':
            indices = np.where(~np.isnan(data[speed_type]))
            res = {key : data[key][indices] for key in list(data.keys()) if key not in self.speed_types}
            res['speed'] = data[speed_type][indices]
        elif isinstance(nan_policy, (int, float)):
            indices = np.where(np.isnan(data[speed_type]))
            res = {key : data[key] for key in list(data.keys()) if key not in self.speed_types}
            vv = data[speed_type]
            vv[indices] = nan_policy
            res['speed'] = vv
        else:
            raise ValueError("nan_policy = 'discard' to discard all parameter values at indices for which speeds are nan, or type <int / float> to replace all nan-valued speeds")
        return res

class TimeSeriesConversion(SohoLascoCharacterMap):

    def __init__(self):
        """ """
        super().__init__()
        self.daily_cmes = {}
        self.daily_sunspots = {}

    @staticmethod
    def get_padded_temporal_data(data, delta_hours=1):
        """ """
        padded_dates, padded_times, padded_dts = [], [], []
        padded_years, padded_months, padded_days = [], [], []
        padded_hours, padded_minutes, padded_seconds = [], [], []
        initial_dtime = data['datetime object'][0]
        final_dtime = data['datetime object'][-1]
        step = datetime.timedelta(hours=delta_hours)
        while initial_dtime < final_dtime + step:
            pad_date = initial_dtime.strftime('%Y/%m/%d')
            pad_time = initial_dtime.strftime('%H:%M:%S')
            padded_dates.append(pad_date)
            padded_times.append(pad_time)
            parameters = initial_dtime.strftime('%Y/%m/%d %H:%M:%S')
            pad_year = int(parameters[:4])
            pad_month = int(parameters[5:7])
            pad_day = int(parameters[8:10])
            pad_hour = int(parameters[11:13])
            pad_minute = int(parameters[14:16])
            pad_second = int(parameters[17:])
            pad_dt = datetime.datetime(pad_year, pad_month, pad_day, pad_hour, pad_minute, pad_second)
            padded_dts.append(pad_dt)
            padded_years.append(pad_year)
            padded_months.append(pad_month)
            padded_days.append(pad_day)
            padded_hours.append(pad_hour)
            padded_minutes.append(pad_minute)
            padded_seconds.append(pad_second)
            initial_dtime += step
        result = {'datetime object' : np.array(padded_dts), 'date' : np.array(padded_dates), 'time' : np.array(padded_times)}
        result['year'] = np.array(padded_years)
        result['month'] = np.array(padded_months)
        result['day'] = np.array(padded_days)
        result['hour'] = np.array(padded_hours)
        result['minute'] = np.array(padded_minutes)
        result['second'] = np.array(padded_seconds)
        indices = data['elapsed hour']
        for key in ('datetime object', 'time', 'minute', 'second'):
            result[key][indices] = data[key]
        result['elapsed hour'] = np.arange(result['datetime object'].size)
        return indices, result

    @staticmethod
    def get_padded_speeds(data, indices, size, pad_speed=1):
        """ """
        try:
            base = np.ones(size).astype(int) * pad_speed
        except:
            base = np.zeros(size).astype(int)
        base[indices] = data['speed']
        return base

    @staticmethod
    def get_other_padded_parameters(data, indices, size, exclusive_keys, pad_value=0, pad_booleans=False):
        """ """
        res = {}
        for key in list(data.keys()):
            if key not in exclusive_keys:
                args = data[key]
                if args.dtype == 'str':
                    largest_string = max(args, key=len)
                    n = len(largest_string)
                    base = np.tile(['-' * n], n)
                elif args.dtype == 'bool':
                    base = np.tile([pad_booleans], size)
                elif args.dtype in (float, np.float, int, np.int, np.int64):
                    if isinstance(pad_value, (float, np.float)) or isinstance(args.dtype, (float, np.float)):
                        base = np.ones(size).astype(float) * pad_value
                    else:
                        try:
                            base = np.ones(size).astype(int) * pad_value
                        except:
                            base = np.zeros(size).astype(int)
                else:
                    raise ValueError("invalid dtype={} for data[{}]".format(key, args.dtype))
                base[indices] = data[key]
                res[key] = base
        return res

    def get_padded_data(self, data, delta_hours=1, pad_speed=1, pad_value=0, pad_booleans=False):
        """ """
        indices, padded_data = self.get_padded_temporal_data(data, delta_hours)
        padded_data['speed'] = self.get_padded_speeds(data, indices, padded_data['datetime object'].size, pad_speed)
        other_padded_data = self.get_other_padded_parameters(data, indices, padded_data['datetime object'].size, list(padded_data.keys()), pad_value, pad_booleans)
        padded_data.update(other_padded_data)
        return padded_data

    def count_daily_cmes(self, data):
        """ """
        primary_condition = (np.diff(data['day']) != 0)
        secondary_condition = (np.diff(data['month']) != 0)
        # secondary_condition = (np.diff(data['month']) > 0)
        indices = np.where(primary_condition | secondary_condition)[0] + 1
        gkeys = ('datetime object', 'year', 'month', 'day', 'count')
        res = {key : [] for key in gkeys}
        groups = {key : np.split(value, indices) for key, value in data.items() if key in gkeys}
        for dt_objects, years, months, days in zip(groups['datetime object'], groups['year'], groups['month'], groups['day']):
            res['datetime object'].append(dt_objects[0])
            res['year'].append(years[0])
            res['month'].append(months[0])
            res['day'].append(days[0])
            res['count'].append(dt_objects.size)
        for key, value in res.items():
            self.daily_cmes[key] = np.array(value)
        return self.daily_cmes

    def update_daily_sunspots(self, data):
        """ """
        self.daily_sunspots.update(data)

    @staticmethod
    def count_events_subroutine(daily_data, desired_timescale):
        """ """
        gkeys = ('datetime object', 'year', 'month', 'day', 'count')
        if len(daily_data.keys()) == 0:
            raise ValueError("daily_data is not initialized; daily_data should contain the following keys: {}".format(gkeys))
        timescale_mapping = {'daily' : 'day', 'monthly' : 'month', 'yearly' : 'year'}
        key = timescale_mapping[desired_timescale]
        condition = (np.diff(daily_data[key]) != 0)
        indices = np.where(condition)[0] + 1
        res = {key : [] for key in gkeys}
        groups = {key : np.split(value, indices) for key, value in daily_data.items() if key in gkeys}
        for dt_objects, years, months, days, counts in zip(groups['datetime object'], groups['year'], groups['month'], groups['day'], groups['count']):
            res['datetime object'].append(dt_objects[0])
            res['year'].append(years[0])
            res['month'].append(months[0])
            res['day'].append(days[0])
            res['count'].append(np.sum(counts))
        return {key : np.array(value) for key, value in res.items()}

    def count_monthly_from_daily_cmes(self):
        return self.count_events_subroutine(self.daily_cmes, desired_timescale='monthly')

    def count_yearly_from_daily_cmes(self):
        return self.count_events_subroutine(self.daily_cmes, desired_timescale='yearly')

    def count_monthly_from_daily_sunspots(self):
        return self.count_events_subroutine(self.daily_sunspots, desired_timescale='monthly')

    def count_yearly_from_daily_sunspots(self):
        return self.count_events_subroutine(self.daily_sunspots, desired_timescale='yearly')

class DataProcessing(TimeSeriesConversion):

    def __init__(self):
        """ """
        super().__init__()
        self.cmes = {}
        self.sunspots = {}
        self.available_sources = ('SOHO LASCO', 'SILSO')

    def read_from_file(self, path, source):
        """ """
        if source == 'SOHO LASCO':
            kwargs = {'skiprows' : 4, 'dtype' : str}
            kwargs['usecols'] = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
        elif source == 'SILSO':
            kwargs = {'skiprows' : 0, 'dtype' : int, 'unpack' : True}
            kwargs['usecols'] = (0, 1, 2, 4)
        else:
            raise ValueError("invalid source: {}; available source: {}".format(source, self.available_sources))
        data = np.loadtxt(path, **kwargs)
        return data

    def add_to_storage(self, data, source):
        """ """
        if source == 'SOHO LASCO':
            keys = ['date', 'time', 'central position angle', 'angular width']
            keys += ['linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed']
            keys += ['acceleration', 'mass', 'kinetic energy', 'mean position angle']
            cols = np.arange(len(keys)).astype(int)
            for key, col in zip(keys, cols):
                if key in ('date', 'time'):
                    self.cmes[key] = data[:, col]
                else:
                    f = self.parse_parameters[key]
                    self.cmes[key] = f(data[:, col])
            dtimes = self.parse_dates_times(self.cmes['date'], self.cmes['time'])
            self.cmes.update(dtimes)
            self.cmes['halo'] = self.parse_halos(self.cmes['central position angle'])
        elif source == 'SILSO':
            keys = ['year', 'month', 'day', 'count']
            for key, arr in zip(keys, data):
                if key == 'count':
                    arr[arr < 0] = 0
                self.sunspots[key] = arr
            # self.sunspots['datetime object'] = np.array([datetime.date(year, month, day) for year, month, day in zip(data[0], data[1], data[2])])
            self.sunspots['datetime object'] = np.array([datetime.datetime(year, month, day, 0, 0, 0) for year, month, day in zip(data[0], data[1], data[2])])
        else:
            raise ValueError("invalid source: {}; available source: {}".format(source, self.available_sources))

    def store_data(self, path_soho_lasco=None, path_silso=None):
        """ """
        if path_soho_lasco is not None:
            data = self.read_from_file(path_soho_lasco, 'SOHO LASCO')
            self.add_to_storage(data, 'SOHO LASCO')
        if path_silso is not None:
            data = self.read_from_file(path_silso, 'SILSO')
            self.add_to_storage(data, 'SILSO')

    def dispatch_cme_timeseries(self, speed_type, nan_policy=1, solar_cycle=None, subcycle=None, years=None, months=None, days=None, dt_objects=None, delta_hours=1, pad_speed=1, pad_value=0, pad_booleans=False):
        """ """
        data = self.parse_data_by_speeds(self.cmes, speed_type, nan_policy)
        search_kwargs = self.get_search_kwargs(solar_cycle, subcycle, years, months, days, dt_objects)
        Searcher = SearchEvents(data)
        data = Searcher.search(**search_kwargs)
        # data['elapsed hour'] = self.get_elapsed_hours(data['datetime object']).astype(int)
        # data['elapsed hour'] = np.floor(self.get_elapsed_hours(data['datetime object'])).astype(int)
        data['elapsed hour'] = np.round(self.get_elapsed_hours(data['datetime object'])).astype(int)
        indices = npi.group_by(data['elapsed hour']).argmax(data['speed'])[1]
        unique_data = {key : value[indices] for key, value in data.items()}
        padded_data = self.get_padded_data(unique_data, delta_hours, pad_speed, pad_value, pad_booleans)
        speed_title = self.speed_type_mapping[speed_type]
        cycle_title = self.temporal_mapping[solar_cycle]
        identifier = '{}: {}'.format(cycle_title, speed_title)
        dres = {'original' : data, 'unique' : unique_data, 'padded' : padded_data}
        return {'identifier' : identifier, 'data' : dres}

    def dispatch_sunspot_timeseries(self, solar_cycle=None, subcycle=None, years=None, months=None, days=None, dt_objects=None):
        """ """
        search_kwargs = self.get_search_kwargs(solar_cycle, subcycle, years, months, days, dt_objects)
        Searcher = SearchEvents(self.sunspots)
        daily_data = Searcher.search(**search_kwargs)
        return daily_data
