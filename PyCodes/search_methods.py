import numpy as np
from scipy.stats import sem
import operator

class ConditionMapping():

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


    @property
    def types(self):
        res = {}
        res['collective'] = (tuple, list, np.ndarray)
        res['numerical'] = (float, int, np.float, np.int)
        res['element'] = (str, float, int, np.float, np.int, np.int64, bool)
        return res

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

    def autocorrect_single_parameter_inputs(self, parameters, conditions, values):
        """ """
        if isinstance(conditions, str):
            if isinstance(values, self.types['element']):
                parameters = [parameters]
                conditions = [conditions]
                values = [values]
            elif isinstance(values, self.types['collective']):
                nvalues = len(values)
                parameters = [parameters for i in range(nvalues)]
                conditions = [conditions for i in range(nvalues)]
            else:
                raise ValueError("invalid type(values): {}".format(type(values)))
        elif isinstance(conditions, self.types['collective']):
            nconditions = len(conditions)
            if isinstance(values, self.types['element']):
                parameters = [parameters for i in range(nconditions)]
                values = [values for i in range(nconditions)]
            elif isinstance(values, self.types['collective']):
                nvalues = len(values)
                if nconditions != nvalues:
                    raise ValueError("{} search_conditions with {} search_values".format(nconditions, nvalues))
                parameters = [parameters for i in range(nconditions)]
            else:
                raise ValueError("invalid type(values): {}".format(type(values)))
        else:
            raise ValueError("invalid type(conditions): {}".format(type(conditions)))
        return parameters, conditions, values

    def autocorrect_multiple_parameter_inputs(self, parameters, conditions, values):
        """ """
        nparameters = len(parameters)
        if isinstance(conditions, str):
            if isinstance(values, self.types['element']):
                conditions = [conditions for i in range(nparameters)]
                values = [values for i in range(nparameters)]
            elif isinstance(values, self.types['collective']):
                nvalues = len(values)
                if nparameters != nvalues:
                    raise ValueError("{} parameters for {} values".format(nparameters, nvalues))
                conditions = [conditions for i in range(nparameters)]
            else:
                raise ValueError("invalid type(values): {}".format(type(values)))
        elif isinstance(conditions, self.types['collective']):
            nconditions = len(conditions)
            if nparameters != nconditions:
                raise ValueError("{} parameters for {} conditions".format(nparameters, nconditions))
            if isinstance(values, self.types['element']):
                values = [values for value in values]
            elif isinstance(values, self.types['collective']):
                nvalues = len(values)
                if nparameters != nvalues:
                    raise ValueError("{} parameters for {} values".format(nparameters, nvalues))
            else:
                raise ValueError("invalid type(values): {}".format(type(values)))
        else:
            raise ValueError("invalid type(conditions): {}".format(type(conditions)))
        return np.array(parameters), np.array(conditions), tuple(values)

    def autocorrect_search_inputs(self, parameters, conditions, values, modifiers=None):
        """ """
        if isinstance(parameters, str):
            parameters, conditions, values = self.autocorrect_single_parameter_inputs(parameters, conditions, values)
        elif isinstance(parameters, self.types['collective']):
            parameters, conditions, values = self.autocorrect_multiple_parameter_inputs(parameters, conditions, values)
        else:
            raise ValueError("invalid type(parameters) : {}".format(type(parameters)))
        if (modifiers is None) or (isinstance(modifiers, str)):
            modifiers = [modifiers for parameter in parameters]
        nmodifiers, nparameters = len(modifiers), len(parameters)
        if nmodifiers != nparameters:
            raise ValueError("{} modifiers for {} parameters".format(nmodifiers, nparameters))
        return parameters, conditions, values, modifiers

    def get_indices(self, events, parameters, conditions, values, modifiers):
        """ """
        indices = []
        for parameter, condition, value, modifier in zip(parameters, conditions, values, modifiers):
            data = events[parameter]
            if modifier is not None:
                if modifier in list(self.vector_modifiers.keys()):
                    f = self.vector_modifiers[modifier]
                    data = f(data)
                else:
                    raise ValueError("invalid modifier: {}".format(modifier))
            if isinstance(value, str):
                if value in list(self.statistical_values.keys()):
                    f = self.statistical_values[value]
                    try:
                        value = f(data)
                    except:
                        raise ValueError("invalid type(data) = {} for type(value) = {}".format(type(data), type(value)))
            if condition in list(self.comparisons.keys()):
                f = self.comparisons[condition]
                res = f(data, value)
            else:
                raise ValueError("invalid condition: {}".format(condition))
            if modifier is not None:
                if 'delta' in modifier:
                    base = np.array([False] * (data.size+1))
                    loc = np.where(res == True)[0]
                    if len(loc) > 0:
                        base[loc] = True
                        base[loc+1] = True
                    res = base
            indices.append(res)
        return np.array(indices)

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

class EventSearcher(ConditionMapping):

    def __init__(self, events):
        super().__init__()
        self.events = events

    def search_events(self, parameters=None, conditions=None, values=None, apply_to='all', modifiers=None):
        """ """
        search_args = np.array([parameters, conditions, values])
        condition = (search_args == None)
        if np.all(condition == True):
            result = {key : arr for key, arr in self.events.items()}
            key = list(self.clusters.keys())[0]
            indices = np.arange(self.clusters[key].size).astype(int)
            return result, indices
        else:
            parameters, conditions, values, modifiers = self.autocorrect_search_inputs(parameters, conditions, values, modifiers)
            indices = self.get_indices(self.events, parameters, conditions, values, modifiers)
            indices = self.select_conjunction(indices, apply_to)
            if np.all(indices == False):
                raise ValueError("no matches found")
            result = {key : value[indices] for key, value in self.events.items()}
            return result, indices

class ClusterSearcher(ConditionMapping):

    def __init__(self, clusters):
        super().__init__()
        self._clusters = clusters
        arbitrary_key = list(self._clusters.keys())[0]
        break_point_indices = np.cumsum([cluster.size for cluster in self._clusters[arbitrary_key]]).tolist()
        del break_point_indices[-1]
        self.break_point_indices = np.array(break_point_indices)
        self.load_cluster_parameters()
        self.events = {key : np.concatenate(values, axis=0) for key, values in self._clusters.items()}

    @property
    def clusters(self):
        return self._clusters

    def load_cluster_parameters(self):
        cluster_sizes = []
        cluster_sequence = []
        key = list(self._clusters.keys())[0]
        for ith_cluster, cluster in enumerate(self._clusters[key]):
            for event in cluster:
                cluster_sequence.append(ith_cluster +1)
                cluster_sizes.append(cluster.size)
        result = {}
        keys = list(self._clusters.keys())
        if 'cluster size' not in keys:
            result['cluster size'] = np.array(np.split(cluster_sizes, self.break_point_indices))
        if 'ith cluster' not in keys:
            result['ith cluster'] = np.array(np.split(cluster_sequence, self.break_point_indices))
        self._clusters.update(result)

    def search_clusters(self, parameters=None, conditions=None, values=None, apply_to='all', modifiers=None):
        """ """
        search_args = np.array([parameters, conditions, values])
        condition = (search_args == None)
        if np.all(condition == True):
            result = {key : arr for key, arr in self.clusters.items()}
            key = list(self.clusters.keys())[0]
            indices = np.arange(self.clusters[key].size).astype(int)
            return result, indices
        else:
            parameters, conditions, values, modifiers = self.autocorrect_search_inputs(parameters, conditions, values, modifiers)
            indices = self.get_indices(self.events, parameters, conditions, values, modifiers)
            indices = self.select_conjunction(indices, apply_to)
            if isinstance(indices[0], (tuple, list, np.ndarray)):
                indices = indices[0]
            if np.all(indices == False):
                raise ValueError("no matches found")
            loc = np.where(indices == True)[0]
            indices = np.unique(np.digitize(loc, self.break_point_indices, right=False))
            result = {key : values[indices] for key, values in self.clusters.items()}
            return result, indices


##
