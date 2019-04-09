import numpy as np
from scipy.stats import sem
import operator

class IndexMapper():

    def __init__(self):
        self.comparisons = {}
        self.comparisons['equal'] = operator.eq
        self.comparisons['equality'] = operator.eq
        self.comparisons['exact match'] = operator.eq
        self.comparisons['greater than'] = operator.gt
        self.comparisons['greater than or equal'] = operator.ge
        self.comparisons['less than'] = operator.lt
        self.comparisons['less than or equal'] = operator.le
        self.comparisons['lesser than'] = operator.lt
        self.comparisons['lesser than or equal'] = operator.le
        self.comparisons['not equal'] = operator.ne
        self.collective_types = (tuple, list, np.ndarray)
        self.numerical_types = (float, int, np.float, np.int)
        self.element_types = (str, float, int, np.float, np.int, np.int64)

    @staticmethod
    def from_nearest(data, value):
        """
        data                :   type <array>
        value               :   type <int / float>
        """
        delta = np.abs(data - value)
        loc = np.where(delta == np.min(delta))[0]
        res = np.array([False for i in range(len(data))])
        res[loc] = True
        return res

    @staticmethod
    def from_nearest_forward(data, value):
        """
        data                :   type <array>
        value               :   type <int / float>
        """
        delta = data - value
        try:
            loc = np.where(delta == np.min(delta[delta >= 0]))
        except:
            raise ValueError("no forward-nearest match exists")
        res = np.array([False for i in range(len(data))])
        res[loc] = True
        return res

    @staticmethod
    def from_nearest_backward(data, value):
        """
        data                :   type <array>
        value               :   type <int / float>
        """
        delta = value - data
        try:
            loc = np.where(delta == np.min(delta[delta >= 0]))[0]
        except:
            raise ValueError("no backward-nearest match exists")
        res = np.array([False for i in range(len(data))])
        res[loc] = True
        return res

    @property
    def additional_comparisons(self):
        fmap = {}
        fmap['nearest'] = lambda data, value: self.from_nearest(data, value)
        fmap['nearest forward'] = lambda data, value: self.from_nearest_forward(data, value)
        fmap['nearest backward'] = lambda data, value: self.from_nearest_backward(data, value)
        return fmap

    @property
    def available_conditions(self):
        return list(self.comparisons.keys()) + list(self.additional_comparisons.keys())

    def get_f(self, condition):
        """
        condition           :   type <str>
        """
        if condition in list(self.comparisons.keys()):
            f = self.comparisons[condition]
        elif condition in list(self.additional_comparisons.keys()):
            f = self.additional_comparisons[condition]
        else:
            raise ValueError("unknown condition: {}; available conditions: {}".format(condition, self.available_conditions))
        return f

    @property
    def value_mapping(self):
        vmap = {}
        vmap['mean'] = lambda args : np.mean(args)
        vmap['median'] = lambda args : np.median(args)
        vmap['standard deviation'] = lambda args : np.std(args)
        vmap['standard error'] = lambda args : sem(args)
        return vmap

    @property
    def function_mapping(self):
        fmap = {}
        fmap['delta'] = lambda args : np.diff(args)
        fmap['cumulative sum'] = lambda args : np.cumsum(args)
        return fmap

    @property
    def absolute_function_mapping(self):
        fmap = {}
        fmap['delta'] = lambda args : np.abs(np.diff(args))
        fmap['cumulative sum'] = lambda args : np.cumsum(np.abs(args))
        return fmap

    def single_key_subroutine(self, keys, conditions, values):
        """
        key                 :   type <str>
        conditions          :   type <str> or type <tuple / list / array>
        values              :   type <str> or type <tuple / list / array>
        """
        if isinstance(conditions, str):
            if isinstance(values, self.element_types):
                keys = [keys]
                conditions = [conditions]
                values = [values]
            elif isinstance(values, self.collective_types):
                size = len(values)
                keys = [keys for i in range(size)]
                conditions = [conditions for i in range(size)]
            else:
                raise ValueError("values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        elif isinstance(conditions, self.collective_types):
            size = len(conditions)
            if isinstance(values, self.element_types):
                keys = [keys for i in range(size)]
                values = [values for i in range(size)]
            elif isinstance(values, self.collective_types):
                nvalues = len(values)
                if size != nvalues:
                    raise ValueError("{} search conditions with {} search values".format(size, nvalues))
                keys = [keys for i in range(size)]
            else:
                raise ValueError("values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        else:
            raise ValueError("condition = type <str> or type <{}>".format(self.collective_types))
        return np.array(keys), np.array(conditions), np.array(values)

    def multiple_key_subroutine(self, keys, conditions, values):
        """
        key                 :   type <tuple / list / array>
        conditions          :   type <str> or type <tuple / list / array>
        values              :   type <str> or type <tuple / list / array>
        """
        size = len(keys)
        if isinstance(conditions, str):
            if isinstance(values, self.element_types):
                conditions = [conditions for i in range(size)]
                values = [values for i in range(size)]
            elif isinstance(values, self.collective_types):
                nvalues = len(values)
                if size != nvalues:
                    raise ValueError("{} search keys with {} search values".format(size, nvalues))
                conditions = [conditions for i in range(size)]
            else:
                raise ValueError("values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        elif isinstance(conditions, self.collective_types):
            nconditions = len(conditions)
            if size != nconditions:
                raise ValueError("{} search keys with {} search conditions".format(size, nconditions))
            if isinstance(values, self.element_types):
                values = [values for values in values]
            elif isinstance(values, self.collective_types):
                nvalues = len(values)
                if size != nvalues:
                    raise ValueError("{} search keys with {} search values".format(size, nvalues))
            else:
                raise ValueError("values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        else:
            raise ValueError("condition = type <str> or type <{}>".format(self.collective_types))
        return np.array(keys), np.array(conditions), tuple(values)

    def autocorrect_search_inputs(self, keys, conditions, values, fkeys=None, use_abs=False):
        """
        keys                :   type <tuple / list / array> or type <str>
        conditions          :   type <tuple / list / array> or type <str>
        values              :   type <tuple / list / array> or type <int / float> or type <str>
        fkeys               :   None or type <str> or type <tuple / list / array>
        use_abs             :   type <bool> or type <tuple / list / array>
        """
        if isinstance(keys, str):
            keys, conditions, values = self.single_key_subroutine(keys, conditions, values)
        elif isinstance(keys, self.collective_types):
            keys, conditions, values = self.multiple_key_subroutine(keys, conditions, values)
        else:
            raise ValueError("keys = type <str> or type <{}>".format(self.collective_types))
        if ((fkeys is None) or (isinstance(fkeys, self.element_types))):
            fkeys = [fkeys for key in keys]
        elif isinstance(fkeys, self.collective_types):
            n, m = len(fkeys), len(keys)
            if n != m:
                raise ValueError("{} fkeys with {} respective search args".format(n, m))
        else:
            raise ValueError("fkey = None or type <{}> or type <{}>".format(self.element_types, self.collective_types))
        if isinstance(use_abs, bool):
            use_abs = [use_abs for key in keys]
        elif isinstance(use_abs, self.collective_types):
            n, m = len(use_abs), len(keys)
            if n != m:
                raise ValueError("{} use_abs with {} respective search args".format(n, m))
        else:
            raise ValueError("use_abs = type <bool> or type <{}>".format(self.collective_types))
        return keys, conditions, values, fkeys, use_abs

    def get_value_from_string(self, value, data):
        """
        value               :   type <str>
        data                :   type <array>
        """
        if value not in list(self.value_mapping.keys()):
            raise ValueError("unknown value: {}; available type <str> values: {}".format(value, tuple(list(self.value_mapping.keys()))))
        f = self.value_mapping[value]
        value = f(data)
        return value

class SearchEjecta():

    def __init__(self, data, string_keys=('datetime',)):
        """
        data                :   type <dict>
        string_keys         :   type <tuple / list / array>
        """
        self.data = data
        self.string_keys = string_keys
        self.IndexMap = IndexMapper()
        self.result = {}

    def load_data(self, key, value):
        """ """
        keys = list(self.data.keys())
        size = len(self.data[keys[0]])
        nvalues = len(value)
        if size == nvalues:
            if isinstance(value, np.ndarray):
                self.data[key] = value
            else:
                self.data[key] = np.array(value)
        else:
            raise ValueError("value should contain {} elements, not {} elements".format(size, nvalues))

    def subsect(self, indices):
        """
        indices             :   type <array>
        """
        return {key : value[indices] for key, value in self.data.items()}

    def check_halo_condition(self, use_halo):
        """
        use_halo            :   type <bool> or None

            True
                select halo events

            False
                select non-halo events

            None
                select all events
        """
        if use_halo is None:
            return self.data
        elif isinstance(use_halo, bool):
            indices = np.where(self.data['halo'] == use_halo)[0]
            return self.subsect(indices)
        else:
            raise ValueError("use_halo = type <bool> or None")

    def retrieve_data(self, key, fkey, use_abs=False, use_halo=None):
        """
        key                 :   type <str>
        fkey                :   None or type <str>
        use_abs             :   type <bool>
        """
        initial_data = self.check_halo_condition(use_halo)
        data = initial_data[key]
        # data = self.data[key]
        if fkey is not None:
            if isinstance(fkey, str):
                if use_abs is True:
                    fdata = self.IndexMap.absolute_function_mapping[fkey]
                elif use_abs is False:
                    fdata = self.IndexMap.function_mapping[fkey]
                else:
                    raise ValueError("use_abs = True or False")
                data = fdata(data)
                if fkey == 'delta':
                    if key in self.string_keys:
                        raise ValueError("not yet implemented; delta applies only for explicitly numerical data")
                    else:
                        check = np.issubdtype(data.dtype, np.integer)
                        if check is True:
                            data = data.astype(float)
                        data = np.insert(data, 0, np.nan)
            else:
                raise ValueError("unknown fkey: {}; available fkeys: {}".format(fkey, list(self.IndexMap.function_mapping.keys())))
        return data

    @staticmethod
    def modify_delta_indices(indices):
        """
        indices             :   type <array>
        """
        base = np.array([False for i in range(len(indices))])
        loc = np.array(np.where(indices == True)[0])
        if loc.size > 0:
            partials = loc - 1
            res = np.unique(loc.tolist() + partials.tolist())
            base[indices] = True
            return base
        else:
            return indices

    def get_indices_per_condition(self, key, condition, value, fkey=None, use_abs=False, use_halo=None):
        """
        key                 :   type <tuple / list / array> or type <str>
        condition           :   type <tuple / list / array> or type <str>
        value               :   type <tuple / list / array> or type <int / float>
        fkey                :   None or type <tuple / list / array> or type <str>
        use_abs             :   type <tuple / list / array> or type <bool>
        """
        data = self.retrieve_data(key, fkey, use_abs, use_halo)
        if condition in list(self.IndexMap.comparisons.keys()):
            f = self.IndexMap.comparisons[condition]
        elif condition in list(self.IndexMap.additional_comparisons.keys()):
            f = self.IndexMap.additional_comparisons[condition]
        else:
            raise ValueError("unknown condition: {}; available conditions: {}".format(condition, self.IndexMap.available_conditions))
        res = np.array(f(data, value))
        if fkey == 'delta':
            res = self.modify_delta_indices(res)
        return np.array(res)

    def get_indices(self, search_keys, search_conditions, search_values, fkeys=(), use_abs=(), use_halo=None):
        """
        search_keys         :   type <tuple / list / array> or type <str>
        search_conditions   :   type <tuple / list / array> or type <str>
        search_values       :   type <tuple / list / array> or type <int / float>
        fkeys               :   None or type <tuple / list / array> or type <str>
        use_abs             :   type <tuple / list / array> or type <bool>
        """
        res = []
        for key, condition, value, fkey, use_absolute in zip(search_keys, search_conditions, search_values, fkeys, use_abs):
            if isinstance(value, str):
                if key not in self.string_keys:
                    value = self.IndexMap.get_value_from_string(value, self.data[key])
            ith_indices = self.get_indices_per_condition(key, condition, value, fkey, use_absolute, use_halo)
            res.append(ith_indices)
        return np.array(res)

    @staticmethod
    def select_conjunction(indices, apply_to):
        """
        indices             :   type <array>
        apply_to            :   type <str>
        """
        if apply_to == 'all':
            indices = np.all(indices, axis=0)
        elif apply_to == 'any':
            indices = np.any(indices, axis=0)
        return np.array(indices)

    def search(self, search_parameters, search_conditions, search_values, fkeys=None, ret='data', apply_to='all', use_abs=False, use_halo=None):
        """
        search_parameters   :   type <tuple / list / array> or type <str>
        search_conditions   :   type <tuple / list / array> or type <str>
        search_values       :   type <tuple / list / array> or type <int / float>
        fkeys               :   None or type <tuple / list / array> or type <str>
        ret                 :   type <str>
        apply_to            :   type <str>
        use_abs             :   type <tuple / list / array> or type <bool>
        """
        search_parameters, search_conditions, search_values, fkeys, use_abs = self.IndexMap.autocorrect_search_inputs(search_parameters, search_conditions, search_values, fkeys, use_abs)
        all_indices = self.get_indices(search_parameters, search_conditions, search_values, fkeys, use_abs, use_halo) #, apply_to)
        indices = self.select_conjunction(all_indices, apply_to)
        if np.all(indices == False):
            raise ValueError("no matches found")
        if ret == 'indices':
            return indices
        elif ret == 'data':
            return self.subsect(indices)
        else:
            raise ValueError("ret = 'indices' or 'data'")

class SearchClusters():

    def __init__(self, data, indices, string_keys=('datetime',)):
        """
        data                :   type <dict>
        indices             :   type <array>
        string_keys         :   type <tuple / list / array>
        """
        self.data = data
        self.indices = indices
        self.string_keys = string_keys
        self.IndexMap = IndexMapper()
        self.keys = {'data' : tuple(list(self.data.keys())), 'load' : ('ith cluster', 'cluster size')}
        self.result = {}

    @staticmethod
    def flatten_data(data):
        """
        data                :   type <array>
        """
        return np.concatenate(data, axis=0)

    @property
    def ejecta(self):
        return {key : self.flatten_data(value) for key, value in self.data.items()}

    def initialize_additional_parameters(self):
        key = self.keys['data'][0]
        arr = self.data[key]
        iclusters, cluster_sizes = [], []
        for idx, subarr in enumerate(arr):
            ith_cluster = idx + 1
            for element in subarr:
                iclusters.append(ith_cluster)
                cluster_sizes.append(subarr.size)
        self.result['ith cluster'] = np.array(iclusters)
        self.result['cluster size'] = np.array(cluster_sizes)

    @property
    def ith_cluster(self):
        return self.result['ith cluster']

    @property
    def cluster_size(self):
        return self.result['cluster size']

    @property
    def SearchEjections(self):
        Searcher = SearchEjecta(self.ejecta)
        return Searcher

    def load_parameters(self, Searcher, load_keys=None):
        """
        Searcher            :   type <cls>
        load_keys           :   None or type <tuple / list / array>
        """
        if ((load_keys is not None) or (len(load_keys) > 0)):
            self.initialize_additional_parameters()
            for key in load_keys:
                if key not in self.keys['load']:
                    raise ValueError("unknown key: {}; available load keys: {}".format(key, self.keys['load']))
                if key == 'ith cluster':
                    flat_values = self.ith_cluster
                elif key == 'cluster size':
                    flat_values = self.cluster_size
                Searcher.load_data(key, flat_values)
                tmp = np.split(flat_values, self.indices)
                self.data[key] = np.array(tmp)
        return Searcher

    @staticmethod
    def subsect(data, indices, view_keys=None):
        """
        data                :   type <dict>
        indices             :   type <array>
        view_keys           :   None or type <str> or type <tuple / list / array>
        """
        if view_keys is None:
            return {key : value[indices] for key, value in data.items()}
        else:
            if isinstance(view_keys, str):
                view_keys = [view_keys]
            return {key : value[indices] for key, value in data.items() if key in view_keys}

    def search(self, search_parameters, search_conditions, search_values, fkeys=None, load_keys=('ith cluster', 'cluster size'), apply_to='all', view_keys=None, use_abs=False, use_halo=None):
        """
        search_parameters   :   type <tuple / list / array> or type <str>
        search_conditions   :   type <tuple / list / array> or type <str>
        search_values       :   type <tuple / list / array> or type <int / float>
        fkeys               :   None or type <tuple / list / array> or type <str>
        load_keys           :   None or type <tuple / list / array>
        apply_to            :   type <str>
        view_keys           :   None or type <str> or type <tuple / list / array>
        use_abs             :   type <tuple / list / array> or type <bool>
        """
        Searcher = SearchEjecta(self.ejecta, self.string_keys)
        Searcher = self.load_parameters(Searcher, load_keys)
        flat_indices = Searcher.search(search_parameters, search_conditions, search_values, fkeys, ret='indices', apply_to=apply_to, use_abs=use_abs, use_halo=use_halo)
        loc = np.where(flat_indices == True)[0]
        cluster_indices = np.unique(np.digitize(loc, self.indices))
        return self.subsect(self.data, cluster_indices, view_keys)

## IMPLEMENT INTRA/INTER-DELTA/SUM SEARCH METHODS
