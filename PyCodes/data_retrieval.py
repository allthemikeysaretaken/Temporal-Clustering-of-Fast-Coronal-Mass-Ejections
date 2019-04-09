import numpy as np
import numpy.core.defchararray as np_repl
import re
import datetime
import numpy_indexed as npi

from search_methods import *

class DataInitializer():

    def __init__(self):
        """ """
        self.data = {}
        self.keys = []
        self.speed_types = ('linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed')

    @staticmethod
    def replace_string_characters(values, prev_char, repl_char, dtype):
        """
        values      :   type <array>
        prev_char   :   type <str> or None
        repl_char   :   type <str>
        dtype       :   float, int, or str
        """
        if prev_char is None:
            if dtype == str:
                return values
            else:
                return values.astype(dtype)
        else:
            return np_repl.replace(values, prev_char, repl_char).astype(dtype)

    @staticmethod
    def create_dt_object(yr, mt, dy, hr, mn, sc):
        """ """
        return datetime.datetime(yr, mt, dy, hr, mn, sc)

    @staticmethod
    def get_relative_elapsed_hours(dt_objs):
        """ """
        return np.array([(dt_objs[idx] - dt_objs[0]).total_seconds() for idx in range(dt_objs.size)]) / 3600

    @staticmethod
    def autocorrect_nans(data, nan_policy, nan_key, nan_repl):
        """ """
        if nan_policy is not None:
            if nan_policy == 'replace':
                indices = np.where(np.isnan(data[nan_key]))
                data[nan_key][indices] = nan_repl
            elif nan_policy == 'discard':
                indices = np.where(~np.isnan(data[nan_key]))
                data = {key : data[key][indices] for key in list(data.keys())}
            else:
                raise ValueError("nan_policy = 'replace' or 'discard'")
        return data

    def read_soho_lasco_catalog(self, path):
        """ """
        kwargs = {'skiprows' : 4, 'dtype' : str}
        kwargs['usecols'] = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
        data = np.loadtxt(path, **kwargs)
        return data

    @staticmethod
    def check_halos(arr):
        """ """
        condition = (arr == 'Halo')
        return condition

    def separate_columns(self, data, repl_char):
        """ """
        keys = ('date', 'time', 'central position angle', 'angular width', 'linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed', 'acceleration', 'mass', 'kinetic energy', 'mean position angle')
        self.data['halo'] = np.array(self.check_halos(data[:, 2]))
        for idx, key in enumerate(keys):
            if key in ('mass', 'kinetic energy'):
                mod_values = self.replace_string_characters(data[:, idx], '*'*8, 'nan', dtype=str)
                mod_values = self.replace_string_characters(mod_values, '*', '', dtype=str)
                self.data[key] = self.replace_string_characters(mod_values, '-------', repl_char, dtype=float)
            elif key == 'acceleration':
                mod_values = self.replace_string_characters(data[:, idx], '------', 'nan', dtype=str)
                self.data[key] = self.replace_string_characters(mod_values, '*', '', dtype=float)
            else:
                if key in ('date', 'time'):
                    prev_char = None
                    dtype = str
                elif key == 'central position angle':
                    prev_char = 'Halo'
                    dtype = float
                elif key == 'angular width':
                    prev_char = None
                    dtype = int
                elif key in ('linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed'):
                    prev_char = '----'
                    dtype = float
                elif key == 'mean position angle':
                    prev_char = None
                    dtype = int
                self.data[key] = self.replace_string_characters(data[:, idx], prev_char, repl_char, dtype)
            self.keys.append(key)

    def store_datetime_objects(self, objs):
        """ """
        self.data['dt object'] = np.array(objs)
        self.keys.append('dt object')

    def store_elapsed_hours(self):
        """ """
        elapsed = self.get_relative_elapsed_hours(self.data['dt object'])
        self.data['elapsed'] = np.round(elapsed).astype(int)
        self.keys.append('elapsed')

    def organize_temporal_columns(self):
        """ """
        keys = ('datetime', 'year', 'month', 'day', 'hour', 'minute', 'second')
        dtypes = tuple([str] + [int for kdx in range(len(keys)-1)])
        dts, ys, ms, ds, hs, mins, secs, objs = [], [], [], [], [], [], [], []
        for iter_date, iter_time in zip(self.data['date'], self.data['time']):
            dt = '{} {}'.format(iter_date, iter_time)
            yr = int(iter_date[:4])
            mt = int(iter_date[5:7])
            dy = int(iter_date[8:10])
            hr = int(iter_time[:2])
            mn = int(iter_time[3:5])
            sc = int(iter_time[6:])
            obj = self.create_dt_object(yr, mt, dy, hr, mn, sc)
            dts.append(dt)
            ys.append(yr)
            ms.append(mt)
            ds.append(dy)
            hs.append(hr)
            mins.append(mn)
            secs.append(sc)
            objs.append(obj)
        values = (dts, ys, ms, ds, hs, mins, secs)
        for key, value, dtype in zip(keys, values, dtypes):
            self.data[key] = np.array(value, dtype=dtype)
            self.keys.append(key)
        return objs

    def execute(self, path, repl_char='nan'):
        """ """
        data = self.read_soho_lasco_catalog(path)
        self.separate_columns(data, repl_char)
        objs = self.organize_temporal_columns()
        self.store_datetime_objects(objs)
        self.store_elapsed_hours()

    def dispatch(self, speed_type, nan_policy='replace', nan_repl=1):
        """ """
        if speed_type not in self.speed_types:
            raise ValueError("unknown speed_type: {}; available speed_types = {}".format(speed_type, self.speed_types))
        if nan_policy is None:
            return self.data
        else:
            res = self.autocorrect_nans(self.data.copy(), nan_policy, speed_type, nan_repl)
            return res

class SolarCycleSelector():

    def __init__(self, data):
        """ """
        self.data = data
        self.Searcher = SearchEjecta(data)
        self.available = {}
        self.available['solar cycle'] = (23, 24)
        self.available['cycle type'] = (None, 'pre-maximum', 'post-maximum', 'high-activity')

    @staticmethod
    def get_search_kwargs(sc, ctype):
        """ """
        res = {}
        if sc == 23:
            if ctype is None:
                initial = datetime.datetime(1996, 8, 1, 0, 0, 0)
                final = datetime.datetime(2008, 12, 31, 23, 0, 0)
                res['search_parameters'] = 'dt object'
            elif ctype == 'pre-maximum':
                initial = datetime.datetime(1996, 8, 1, 0)
                final = datetime.datetime(1998, 12, 31, 23)
                res['search_parameters'] = 'dt object'
            elif ctype == 'post-maximum':
                initial = datetime.datetime(2007, 1, 1, 0)
                final = datetime.datetime(2008, 12, 31, 23)
                res['search_parameters'] = 'dt object'
            elif ctype == 'high-activity':
                initial = 1999
                final = 2006
                res['search_parameters'] = 'year'
        elif sc == 24:
            if ctype is None:
                initial = datetime.datetime(2009, 1, 1, 0)
                final = datetime.datetime(2017, 10, 31, 23)
                res['search_parameters'] = 'dt object'
            elif ctype == 'pre-maximum':
                initial = datetime.datetime(2009, 1, 1, 0)
                final = datetime.datetime(2016, 12, 31, 23)
                res['search_parameters'] = 'dt object'
            elif ctype == 'post-maximum':
                initial = datetime.datetime(2017, 1, 1, 0)
                final = datetime.datetime(2017, 10, 31, 23)
                res['search_parameters'] = 'dt object'
            elif ctype == 'high-activity':
                initial = 2010
                final = 2016
                res['search_parameters'] = 'year'
        res['search_values'] = (initial, final)
        if res['search_parameters'] == 'dt object':
            res['search_conditions'] = ('greater than or equal', 'less than or equal')
        elif res['search_parameters'] == 'year':
            res['search_conditions'] = ('greater than or equal', 'less than or equal')
        return res

    def dispatch(self, sc, ctype=None):
        """ """
        if sc is None:
            return self.data
        else:
            if sc not in self.available['solar cycle']:
                raise ValueError("unknown sc: {}; available solar cycles: {}".format(sc, self.available['solar cycle']))
            if ctype not in self.available['cycle type']:
                raise ValueError("unknown ctype: {}; available cycle types: {}".format(ctype, self.available['cycle type']))
            search_kwargs = self.get_search_kwargs(sc, ctype)
            return self.Searcher.search(**search_kwargs)

class TimeSeriesInitializer():

    def __init__(self, data, solar_cycle=23, cycle_type='high-activity'):
        """ """
        self.data = data
        self.solar_cycle = solar_cycle
        self.cycle_type = cycle_type
        self.speed_types = ('linear speed', 'second order initial speed', 'second order final speed', 'second order 20R speed')
        self.result = {}

    @staticmethod
    def override_values_at_indices(unique_data, key, numerical_mask):
        """ """
        indices = unique_data['elapsed']
        numerical_mask[indices] = unique_data[key]
        return numerical_mask

    def pad_datetimes(self, unique_data, nelapsed):
        """ """
        dt_init, dt_fin = unique_data['dt object'][0], unique_data['dt object'][-1]
        step = datetime.timedelta(hours=nelapsed)
        res = []
        while dt_init < dt_fin + step:
            filler = dt_init.strftime('%Y/%m/%d %H:%M:%S')
            res.append(filler)
            dt_init += step
        return np.array(res)

    def pad_speeds(self, unique_data, speed_type, size):
        """ """
        base_mask = np.ones(size, dtype=int)
        return self.override_values_at_indices(unique_data, speed_type, base_mask)

    def pad_datetime_objects_and_parameters(self, dts, sym=('/|:| ')):
        """ """
        keys = ('year', 'month', 'day', 'hour', 'minute', 'second')
        partials, dt_objs = [], []
        for dt in dts:
            partial = np.array(re.split(sym, dt), dtype=int)
            dt_obj = datetime.datetime(*partial)
            partials.append(partial)
            dt_objs.append(dt_obj)
        partials = np.array(partials).T
        res = {key : partial for key, partial in zip(keys, partials)}
        res['dt object'] = np.array(dt_objs)
        return res

    def get_padded_nondt_parameters(self, unique_data, size):
        """ """
        res = {}
        numerical_mask = np.zeros(size)
        speed_mask = np.ones(size, dtype=int)
        base_keys = ('central position angle', 'angular width', 'mean position angle', 'acceleration', 'mass', 'kinetic energy', 'halo')
        for key in base_keys:
            if key == 'halo':
                base_mask = numerical_mask.copy().astype(bool)
            else:
                base_mask = numerical_mask.copy().astype(unique_data[key].dtype)
            res[key] = self.override_values_at_indices(unique_data, key, base_mask)
        return res

    def get_original_data(self):
        """ """
        res = SolarCycleSelector(self.data).dispatch(self.solar_cycle, self.cycle_type)
        res['elapsed'] = res['elapsed'] - res['elapsed'][0]
        return res

    def get_unique_data(self, original_data, speed_type):
        """ """
        indices = npi.group_by(original_data['elapsed']).argmax(original_data[speed_type])[1]
        return {key : original_data[key][indices] for key in list(original_data.keys())}

    def get_padded_data(self, unique_data, speed_type, nelapsed=1):
        """ """
        res = {}
        dts = self.pad_datetimes(unique_data, nelapsed)
        res[speed_type] = self.pad_speeds(unique_data, speed_type, dts.size)
        dt_kwargs = self.pad_datetime_objects_and_parameters(dts)
        other_kwargs = self.get_padded_nondt_parameters(unique_data, dts.size)
        res.update(dt_kwargs)
        res.update(other_kwargs)
        res['datetime'] = dts
        indices = unique_data['elapsed']
        for key in ('datetime', 'dt object', 'minute', 'second'): # ('year', 'month', 'day', 'hour', 'minute', 'second')
            res[key][indices] = unique_data[key]
        res['elapsed'] = np.array(list(range(dts.size)), dtype=int)
        return res

    def execute(self, speed_type):
        """ """
        original_data = self.get_original_data()
        unique_data = self.get_unique_data(original_data, speed_type)
        padded_data = self.get_padded_data(unique_data, speed_type)
        self.result['original data'] = original_data
        self.result['unique data'] = unique_data
        self.result['padded data'] = padded_data
        for key, item in self.result.items():
            self.result[key]['speed'] = self.result[key][speed_type]

class Data():

    def __init__(self, path):
        """ """
        self.path = path
        self.storage = {}

    def initialize(self):
        """ """
        DI = DataInitializer()
        DI.execute(self.path)
        self.storage['initial data'] = DI

    def dispatch(self, speed_type, solar_cycle=23, cycle_type='high-activity', nan_policy='replace', nan_repl=1):
        """ """
        DI = self.storage['initial data']
        init_data = DI.dispatch(speed_type, nan_policy, nan_repl)
        TimeSeriesData = TimeSeriesInitializer(init_data, solar_cycle, cycle_type)
        TimeSeriesData.execute(speed_type)
        TimeSeriesData.result['identity'] = {'solar cycle' : solar_cycle, 'cycle type' : cycle_type, 'nan policy' : nan_policy, 'nan replacement' : nan_repl, 'speed type' : speed_type}
        return TimeSeriesData.result
