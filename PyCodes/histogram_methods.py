import numpy as np
from scipy.stats import kurtosis, skew
from scipy.integrate import quad

class EdgeOperations():

    def __init__(self):
        self._edges = None

    @property
    def edges(self):
        return self._edges

    @property
    def midpoints(self):
        return (self.edges[1:] + self.edges[:-1])/2

    @property
    def bin_widths(self):
        return np.diff(self.edges)

    @staticmethod
    def get_edges_by_number(data, n):
        """ """
        if isinstance(n, int):
            return np.linspace(min(data), max(data), n)
        else:
            raise ValueError("n = type <int>")

    @staticmethod
    def get_edges_by_width(data, w):
        """ """
        if isinstance(w, (float, int)):
            return np.arange(min(data), max(data) + 0.9 * w, w)
        else:
            raise ValueError("n = type <float / int>")

    def initialize_edges(self, data, n=None, w=None, edges=None):
        """ """
        args = (n, w, edges)
        err_msg = "set n for number of edges, or w for edge-width, or edges for an array of edges"
        # if np.all(args) is None:
        if (n is None) and (w is None) and (edges is None):
            raise ValueError(err_msg)
        if n is not None:
            if np.any((w, edges)) is not None:
                raise ValueError(err_msg)
            self._edges = self.get_edges_by_number(data, n)
        elif w is not None:
            if np.any((n, edges)) is not None:
                raise ValueError(err_msg)
            self._edges = self.get_edges_by_width(data, w)
        elif edges is not None:
            if np.any((n, w)) is not None:
                raise ValueError(err_msg)
            self._edges = edges

class CountOperations(EdgeOperations):

    def __init__(self):
        super().__init__()
        self._observed_counts = None

    @property
    def observed_counts(self):
        return self._observed_counts

    @property
    def normalization_constant(self):
        return np.sum(self.bin_widths * self.observed_counts)

    @property
    def normalized_counts(self):
        return self.observed_counts / self.normalization_constant

    def initialize_observed_counts(self, data, bias='left'):
        """ """
        if bias == 'left':
            counts, edges = np.histogram(data, self.edges)
            self._edges = edges
        elif bias == 'right':
            counts = np.zeros(len(self.edges) - 1, dtype=int)
            for idx, val in zip(*np.unique(np.searchsorted(self.edges, data, side='left'), return_counts=True)):
                counts[idx - 1] = val
        else:
            raise ValueError("bias = 'left' or 'right'")
        self._observed_counts = counts

class BinThresholding(CountOperations):

    def __init__(self):
        """ """
        super().__init__()

    def left_of_peak_subroutine(self, bin_threshold, idx_bound, idx_peak):
        """ """
        mod_counts, mod_edges = [], []
        while idx_bound <= idx_peak:
            tmp_count = self.observed_counts[idx_bound]
            tmp_edge = self.edges[idx_bound]
            while tmp_count < bin_threshold:
                idx_bound += 1
                tmp_count += self.observed_counts[idx_bound]
                tmp_edge = self.edges[idx_bound]
            mod_counts.append(tmp_count)
            mod_edges.append(tmp_edge)
            idx_bound += 1
        return mod_edges, mod_counts

    def right_of_peak_subroutine(self, bin_threshold, idx_bound, idx_peak):
        """ """
        mod_counts, mod_edges = [], []
        while idx_bound > idx_peak:
            tmp_count = self.observed_counts[idx_bound-1]
            tmp_edge = self.edges[idx_bound]
            while tmp_count < bin_threshold:
                idx_bound -= 1
                tmp_count += self.observed_counts[idx_bound-1]
                tmp_edge = self.edges[idx_bound]
            mod_counts.append(tmp_count)
            mod_edges.append(tmp_edge)
            idx_bound -= 1
        del mod_counts[-1]
        return mod_edges, mod_counts

    def update_edges_and_counts_by_threshold(self, bin_threshold):
        """ """
        if bin_threshold <= 0:
            raise ValueError("bin_threshold > 0")
        updated_edges, updated_counts = [], []
        indices = {'peak' : [np.argmax(self.observed_counts)], 'left' : [0], 'right' : [len(self.edges) - 1]}
        for ileft, ipeak, iright in zip(indices['left'], indices['peak'], indices['right']):
            left_edges, left_counts = self.left_of_peak_subroutine(bin_threshold, ileft, ipeak)
            right_edges, right_counts = self.right_of_peak_subroutine(bin_threshold, iright, ipeak)
            updated_edges += left_edges + right_edges[::-1]
            updated_counts += left_counts + right_counts[::-1]
        self._edges = np.array(updated_edges)
        self._observed_counts = np.array(updated_counts)

    def update_edges_and_counts_by_boundaries(self):
        """ """
        condition = (self._observed_counts > 0)
        nonempty_locations = np.where(condition)[0]
        xmin, xmax = np.min(nonempty_locations), np.max(nonempty_locations)
        self._edges = np.array(self._edges[xmin : xmax +1])
        self._observed_counts = np.array(self._observed_counts[xmin : xmax])

class AdaptiveHistogram(BinThresholding):

    def __init__(self):
        """ """
        super().__init__()
        self._data = None
        self._skew = None
        self._kurtosis = None

    def __str__(self):
        edges = '\n .. edges (shape={}):\n{}\n'.format(self.edges.shape, self.edges)
        midpoints = '\n .. midpoints (shape={}):\n{}\n'.format(self.midpoints.shape, self.midpoints)
        bin_widths = '\n .. bin widths (shape={}):\n{}\n'.format(self.bin_widths.shape, self.bin_widths)
        observed_counts = '\n .. observed counts (shape={}):\n{}\n'.format(self.observed_counts.shape, self.observed_counts)
        normalized_counts = '\n .. normalized counts (shape={}):\n{}\n'.format(self.normalized_counts.shape, self.normalized_counts)
        normalization_constant = '\n .. normalization constant:\n{}\n'.format(self.normalization_constant)
        return '{}{}{}{}{}{}'.format(edges, midpoints, bin_widths, observed_counts, normalized_counts, normalization_constant)

    @property
    def data(self):
        return self._data

    @property
    def skew(self):
        return self._skew

    @property
    def kurtosis(self):
        return self._kurtosis

    def initialize_histogram(self, data, n=None, w=None, edges=None, bin_threshold=None, bias='left', exclude_empty_boundaries=False):
        """ """
        self.initialize_edges(data, n, w, edges)
        self.initialize_observed_counts(data, bias)
        if bin_threshold is not None:
            self.update_edges_and_counts_by_threshold(bin_threshold)
        if exclude_empty_boundaries == True:
            self.update_edges_and_counts_by_boundaries()
        self._skew = skew(data)
        self._kurtosis = kurtosis(data)
        self._data = data

    def get_expectation_values(self, f, args=(), **kwargs):
        """ """
        return np.array([quad(f, lbound, ubound, args=args, **kwargs)[0] * self.data.size for lbound, ubound in zip(self.edges[:-1], self.edges[1:])])

##
