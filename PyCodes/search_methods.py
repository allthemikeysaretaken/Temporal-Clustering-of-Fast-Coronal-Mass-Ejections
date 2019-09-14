import numpy as np
from scipy.stats import sem
import operator

class Symbolizer():

    @property
    def condition_to_symbol_map(self):
        res = {}
        res['greater than'] = '$>$'
        res['greater than or equal'] = '$≥$'
        res['less than'] = '$<$'
        res['less than or equal'] = '$≤$'
        return res

    @property
    def parameter_to_symbol_map(self):
        res = {}
        res['speed'] = r'$V_{CME}$'
        return res

    @property
    def unit_to_symbol_map(self):
        res = {}
        res['speed'] = r'$\frac{km}{s}$'
        return res

class SearchAutoCorrection(Symbolizer):

    def __init__(self):
        super().__init__()
        self.collective_types = (tuple, list, np.ndarray)
        self.numerical_types = (float, int, np.float, np.int)
        self.element_types = (str, float, int, np.float, np.int, np.int64, bool)

    def single_parameter_subroutine(self, search_parameters, search_conditions, search_values):
        """ """
        if isinstance(search_conditions, str):
            if isinstance(search_values, self.element_types):
                search_parameters = [search_parameters]
                search_conditions = [search_conditions]
                search_values = [search_values]
            elif isinstance(search_values, self.collective_types):
                size = len(search_values)
                search_parameters = [search_parameters for i in range(size)]
                search_conditions = [search_conditions for i in range(size)]
            else:
                raise ValueError("search_values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        elif isinstance(search_conditions, self.collective_types):
            size = len(search_conditions)
            if isinstance(search_values, self.element_types):
                search_parameters = [search_parameters for i in range(size)]
                search_values = [search_values for i in range(size)]
            elif isinstance(search_values, self.collective_types):
                nvalues = len(search_values)
                if size != nvalues:
                    raise ValueError("{} search_conditions with {} search_values".format(size, nvalues))
                search_parameters = [search_parameters for i in range(size)]
            else:
                # print("\n TYPE(SEARCH_VALUES) = {}\n".format(type(search_values)))
                # print("\n SEARCH_VALUES = {}\n".format(search_values))
                raise ValueError("search_values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        else:
            raise ValueError("search_conditions = type <str> or type <{}>".format(self.collective_types))
        return search_parameters, search_conditions, search_values

    def multiple_parameter_subroutine(self, search_parameters, search_conditions, search_values):
        """ """
        size = len(search_parameters)
        if isinstance(search_conditions, str):
            if isinstance(search_values, self.element_types):
                search_conditions = [search_conditions for i in range(size)]
                search_values = [search_values for i in range(size)]
            elif isinstance(search_values, self.collective_types):
                nvalues = len(search_values)
                if size != nvalues:
                    raise ValueError("{} search_parameters with {} search_values".format(size, nvalues))
                search_conditions = [search_conditions for i in range(size)]
            else:
                raise ValueError("search_values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        elif isinstance(search_conditions, self.collective_types):
            nconditions = len(search_conditions)
            if size != nconditions:
                raise ValueError("{} search_parameters with {} search_conditions".format(size, nconditions))
            if isinstance(search_values, self.element_types):
                search_values = [search_values for search_values in search_values]
            elif isinstance(search_values, self.collective_types):
                nvalues = len(search_values)
                if size != nvalues:
                    raise ValueError("{} search_parameters with {} search_values".format(size, nvalues))
            else:
                raise ValueError("search_values = type <{}> or type <{}>".format(self.element_types, self.collective_types))
        else:
            raise ValueError("search_conditions = type <str> or type <{}>".format(self.collective_types))
        return np.array(search_parameters), np.array(search_conditions), tuple(search_values)

    def autocorrect_search_inputs(self, search_parameters, search_conditions, search_values, modifiers=None):
        """ """
        search_inputs = np.array([search_parameters, search_conditions])
        if np.all(search_inputs == None):
            raise ValueError("invalid search criteria: all search inputs are None")
        elif np.any(search_inputs == None):
            raise ValueError("invalid search criteria: at least one search input is None")
        else:
            if isinstance(search_parameters, str):
                search_parameters, search_conditions, search_values = self.single_parameter_subroutine(search_parameters, search_conditions, search_values)
            elif isinstance(search_parameters, (tuple, list, np.ndarray)):
                search_parameters, search_conditions, search_values = self.multiple_parameter_subroutine(search_parameters, search_conditions, search_values)
            else:
                raise ValueError("search_parameters = type <str / tuple / list / array>")
        if ((modifiers is None) or (isinstance(modifiers, str))):
            modifiers = [modifiers for search_parameter in search_parameters]
        if len(modifiers) != len(search_parameters):
            raise ValueError("invalid modifiers")
        return search_parameters, search_conditions, search_values, modifiers

class SupplementaryConditions(SearchAutoCorrection):

    def __init__(self):
        super().__init__()

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
        fmap['nearest'] = lambda data, value : self.from_nearest(data, value)
        fmap['nearest forward'] = lambda data, value : self.from_nearest_forward(data, value)
        fmap['nearest backward'] = lambda data, value : self.from_nearest_backward(data, value)
        return fmap

    @property
    def statistical_values(self):
        fmap = {}
        fmap['mean'] = lambda args : np.mean(args)
        fmap['median'] = lambda args : np.median(args)
        fmap['standard deviation'] = lambda args : np.std(args)
        fmap['standard error'] = lambda args : sem(args)
        return fmap

    @property
    def vector_modifiers(self):
        fmap = {}
        fmap['delta'] = lambda args : np.diff(args)
        fmap['absolute delta'] = lambda args : np.abs(np.diff(args))
        fmap['cumulative sum'] = lambda args : np.cumsum(args)
        fmap['absolute cumulative sum'] = lambda args : np.cumsum(np.abs(args))
        return fmap

    @staticmethod
    def verify_halo(data, value):
        """ """
        if isinstance(value, bool):
            inverse = np.invert(value)
            base = np.tile([inverse], data.size)
            indices = np.where(data == value)[0]
            base[indices] = value
            if value == False:
                base = np.invert(base)
        else:
            raise ValueError("value = type <bool>")
        return base

class ConditionMap(SupplementaryConditions):

    def __init__(self):
        super().__init__()
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
        self.comparisons.update(self.additional_comparisons)

    def get_indices_per_condition(self, events, search_parameter, search_conditions, search_values, modifier=None):
        """ """
        data = events[search_parameter]
        if modifier is not None:
            if modifier in list(self.vector_modifiers.keys()):
                f = self.vector_modifiers[modifier]
                data = f(data)
        if isinstance(search_values, str):
            if search_values in list(self.statistical_values.keys()):
                f = self.statistical_values[search_values]
                search_values = f(data)
            # if data.dtype == str:
            #     pass
        if search_conditions in list(self.comparisons.keys()):
            f = self.comparisons[search_conditions]
            res = f(data, search_values)
        if modifier is not None:
            if 'delta' in modifier:
                base = np.full(data.shape[0] +1, False, dtype=bool)
                loc = np.where(res == True)[0]
                if len(loc) > 0:
                    base[loc] = True
                    base[loc +1] = True
                res = base
        return res

    def get_indices(self, events, search_parameters, search_conditions, search_values, modifiers):
        """ """
        res = []
        for parameter, condition, value, modifier in zip(search_parameters, search_conditions, search_values, modifiers):
            indices = self.get_indices_per_condition(events, parameter, condition, value, modifier)
            res.append(indices)
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

class SearchEvents(ConditionMap):

    def __init__(self, events):
        super().__init__()
        self.events = events

    def search(self, search_parameters, search_conditions, search_values, apply_to='all', modifiers=None):
        """ """
        search_parameters, search_conditions, search_values, modifiers = self.autocorrect_search_inputs(search_parameters, search_conditions, search_values)
        indices = self.get_indices(self.events, search_parameters, search_conditions, search_values, modifiers)
        indices = self.select_conjunction(indices, apply_to)
        if np.all(indices == False):
            raise ValueError("no matches found")
        if isinstance(indices[0], np.ndarray):
            indices = indices[0]
        return {key : value[indices] for key, value in self.events.items()}

class SearchClusters(ConditionMap):

    def __init__(self, clusters):
        super().__init__()
        self._clusters = clusters
        self._available_keys = list(self._clusters.keys())

    @property
    def available_keys(self):
        return self._available_keys

    @property
    def break_point_indices(self):
        arbitrary_key = self.available_keys[0]
        return np.cumsum([cluster.size for cluster in self.clusters[arbitrary_key]])

    def add_scale_parameters(self):
        cluster_sizes = []
        cluster_sequence = []
        arbitrary_key = self.available_keys[0]
        for ith_cluster, cluster in enumerate(self.clusters[arbitrary_key]):
            for event in cluster:
                cluster_sequence.append(ith_cluster +1)
                cluster_sizes.append(cluster.size)
        if 'cluster size' not in self.available_keys:
            self._clusters['cluster size'] = np.array(np.split(cluster_sizes, self.break_point_indices))
            self._available_keys.append('cluster size')
        if 'ith cluster' not in self.available_keys:
            self._clusters['ith cluster'] = np.array(np.split(cluster_sequence, self.break_point_indices))
            self._available_keys.append('ith cluster')

    @property
    def clusters(self):
        return self._clusters

    @property
    def events(self):
        return {key : np.concatenate(values, axis=0) for key, values in self.clusters.items()}

    def search(self, search_parameters, search_conditions, search_values, apply_to='all', modifiers=None):
        """ """
        search_parameters, search_conditions, search_values, modifiers = self.autocorrect_search_inputs(search_parameters, search_conditions, search_values, modifiers)
        event_indices = self.get_indices(self.events, search_parameters, search_conditions, search_values, modifiers)
        event_indices = self.select_conjunction(event_indices, apply_to)
        if np.all(event_indices == False):
            raise ValueError("no matches found")
        loc = np.where(event_indices == True)[0]
        cluster_indices = np.unique(np.digitize(loc, self.break_point_indices, right=False))
        return {key : values[cluster_indices] for key, values in self.clusters.items()}






# x = np.linspace(1, 10, 10).astype(int)
# y = np.linspace(1, 100, 10)
# z = np.linspace(-10, 10, 10)
# s = np.array(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'])
# t = x**2
# d = {'x' : x, 'y' : y, 'z' : z, 's' : s, 't' : t}
# Searcher = SearchEvents(d)
#
# # res = Searcher.search('x', 'greater than', 'mean')
# # res = Searcher.search('x', 'greater than', 5, modifiers='cumulative sum')
# # res = Searcher.search('t', 'greater than', 9.1, modifiers='delta')
# # res = Searcher.search('t', 'less than', 5.1, modifiers='cumulative sum')
# res = Searcher.search('t', 'greater than', 4.9, modifiers='cumulative sum')
# # res = Searcher.search(('x', 'y'), ('greater than', 'greater than'), ('mean', 'median'))
# # res = Searcher.search(('x', 't'), ('greater than', 'greater than'), (4, 25))
# # res = Searcher.search(('x', 'y', 't'), ('greater than', 'greater than', 'greater than'), ('mean', 'median', 18), (None, None, 'delta'))
#
# print(Searcher.events)
# # for key, value in Searcher.events.items():
# #     print("\n .. {}:\n{}\n".format(key, value))
#
# print(res)
# # for key, value in res.items():
# #     print("\n .. {}:\n{}\n".format(key, value))
