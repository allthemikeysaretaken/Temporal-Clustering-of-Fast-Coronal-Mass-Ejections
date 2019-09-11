import numpy as np
from scipy.special import gamma
from visual_configurations import *

class DataPoint():

    def __init__(self, index, position):
        """ """
        self.index = index
        self.position = position

class PositionalSystem():

    def __init__(self, data):
        """ """
        self.data = data

    def extract_positions(self):
        try:
            res = np.array([dpoint.position for dpoint in self.data])
            return res
        except:
            return self.data

    def extract_indices(self):
        try:
            res = np.array([dpoint.index for dpoint in self.data])
            return res
        except:
            err_msg_intro = "elements of data are not instances of the class DataPoint"
            err_msg_outro = "so these elements represent positions instead of indices"
            msg = "{}, {}".format(err_msg_intro, err_msg_outro)
            raise ValueError(msg)

class DistanceMetrics(VisualEditor):

    def __init__(self, data):
        """ """
        super().__init__()
        self.data = data
        try:
            self.ndim = data.shape[0]
        except:
            self.ndim = len(data)
        self._method_type = ''
        self._distance_metric = ''
        self._axis = 0
        self._identification_labels = None

    @property
    def method_type(self):
        return self._method_type

    @property
    def distance_metric(self):
        return self._distance_metric

    def get_centroid_deltas(self, points):
        """ """
        deltas = []
        for dim in range(self.ndim):
            delta = np.array([self.data[dim] - pt[dim] for pt in points])
            deltas.append(delta)
        return np.array(deltas)

    def get_density_deltas(self, points):
        """ """
        return self.data.T - points

    def get_all_possible_differences_per_dim(self, dim):
        """ """
        x = np.reshape(self.data[dim], (len(self.data[dim]), 1))
        return x - x.transpose()

    def get_hierarchical_deltas(self):
        return np.array([self.get_all_possible_differences_per_dim(dim) for dim in range(self.ndim)])

    @property
    def available_method_types(self):
        res = {}
        res['centroid'] = self.get_centroid_deltas
        res['density'] = self.get_density_deltas
        res['hierarchical'] = self.get_hierarchical_deltas
        return res

    def update_axis(self):
        if self._method_type in ('centroid', 'hierarchical'):
            self._axis = 0
        elif self._method_type in ('density',):
            self._axis = 1

    def update_method_type(self, method_type):
        """ """
        available = list(self.available_method_types.keys())
        if method_type in available:
            self._method_type = method_type
        else:
            raise ValueError("unknown method_type: {}; available method_types: {}".format(method_type, available))
        self.update_axis()

    @property
    def f_delta(self):
        return self.available_method_types[self._method_type]

    def get_euclidean_square_distance(self, points):
        """ """
        deltas = self.f_delta(points)
        return np.sum(np.square(deltas), axis=self._axis)

    def get_euclidean_distance(self, points):
        """ """
        return np.sqrt(self.get_euclidean_square_distance(points))

    def get_manhattan_distance(self, points):
        """ """
        deltas = self.f_delta(points)
        return np.sum(np.abs(deltas), axis=self._axis)

    @property
    def available_distance_metrics(self):
        res = {}
        res['euclidean'] = self.get_euclidean_distance
        res['euclidean square'] = self.get_euclidean_square_distance
        res['manhattan'] = self.get_manhattan_distance
        return res

    def update_distance_metric(self, distance_metric):
        """ """
        available = list(self.available_distance_metrics.keys())
        if distance_metric in available:
            self._distance_metric = distance_metric
        else:
            raise ValueError("unknown distance_metric: {}; available distance_metric: {}".format(distance_metric, available))

    @property
    def f_distance(self):
        return self.available_distance_metrics[self._distance_metric]

    def update_identification_labels(self, k, labels=None):
        """ """
        if labels is None:
            res = [None for idx in range(k)]
        elif labels == 'auto':
            res = np.arange(k).astype(str) # res = np.core.defchararray.upper(np.arange(k).astype(str))
        else:
            nidentifiers = len(labels)
            if nidentifiers >= k:
                res = [labels[idx] for idx in range(k)]
            else:
                res = list(labels)
                for idx in range(k - nidentifiers):
                    res.append(None)
        self._identification_labels = res

    @property
    def identification_labels(self):
        return np.array(self._identification_labels)

class DistanceMatrix(DistanceMetrics):

    def __init__(self, data):
        """ """
        super().__init__(data)
        self._original_distance_matrix = None
        self._reducible_distance_matrix = None

    def initialize_distance_matrix(self, distance_metric=None):
        """ """
        available = list(self.available_distance_metrics.keys())
        if distance_metric is not None:
            self.update_distance_metric(distance_metric)
        if self.distance_metric in available:
            deltas = self.get_hierarchical_deltas()
            if self.distance_metric == 'euclidean square':
                res = np.sum(np.square(deltas), axis=self._axis)
            elif self.distance_metric == 'euclidean':
                res = np.sqrt(np.sum(np.square(deltas), axis=self._axis))
            elif self.distance_metric == 'manhattan':
                res = np.sum(np.abs(deltas), axis=self._axis)
            else:
                raise ValueError("not yet implemented")
        else:
            raise ValueError("unknown distance_metric: {}; available distance_metric: {}".format(distance_metric, available))
        self._original_distance_matrix = res
        self._reducible_distance_matrix = res

    def reset_distance_matrix(self):
        self._original_distance_matrix = None
        self._reducible_distance_matrix = None

    @property
    def original_distance_matrix(self):
        return self._original_distance_matrix

    @property
    def reducible_distance_matrix(self):
        return self._reducible_distance_matrix

    def get_masked_distance_matrix(self, mask_values):
        """ """
        return np.ma.masked_equal(self._reducible_distance_matrix, mask_values, copy=False)

    def merge_subdimensions(self, row, values):
        """ """
        self._reducible_distance_matrix[row, :] = values
        self._reducible_distance_matrix[:, row] = values

    def reduce_subdimensions(self, col):
        """ """
        self._reducible_distance_matrix = np.delete(self._reducible_distance_matrix, obj=col, axis=0)
        self._reducible_distance_matrix = np.delete(self._reducible_distance_matrix, obj=col, axis=1)

    @property
    def ordered_distances(self):
        return np.argsort(self._original_distance_matrix[0])

##
