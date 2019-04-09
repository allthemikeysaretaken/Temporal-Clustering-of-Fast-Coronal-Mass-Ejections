import numpy as np
from calculus import *

class BinThresholding():

    """
    purpose:
        (1) --  store histogram counts and edges
        (2) --  convenience function to apply threshold
                to histogram counts and edges
        (3) --  obtain normalization constant
    """

    def __init__(self, edges, counts, threshold):
        """
        edges           :   type <tuple / list / array>
        counts          :   type <tuple / list / array>
        threshold       :   type <int / float>
        """
        self.edges = edges
        self.counts = counts
        self.threshold = threshold
        self.indices = {}
        self.result = {}

    def check_threshold(self):
        if self.threshold is None:
            raise ValueError("threshold = type <int / float>")
        if self.threshold <= 0:
            raise ValueError("threshold > 0")

    def left_of_peak(self, idx_bound, idx_peak):
        """
        idx_bound       :   type <int>
        idx_peak        :   type <int>
        """
        mod_counts, mod_edges = [], []
        while idx_bound <= idx_peak:
            tmp_count = self.counts[idx_bound]
            tmp_edge = self.edges[idx_bound]
            while tmp_count < self.threshold:
                idx_bound += 1
                tmp_count += self.counts[idx_bound]
                tmp_edge = self.edges[idx_bound]
            mod_counts.append(tmp_count)
            mod_edges.append(tmp_edge)
            idx_bound += 1
        return mod_edges, mod_counts

    def right_of_peak(self, idx_bound, idx_peak):
        """
        idx_bound       :   type <int>
        idx_peak        :   type <int>
        """
        mod_counts, mod_edges = [], []
        while idx_bound > idx_peak:
            tmp_count = self.counts[idx_bound-1]
            tmp_edge = self.edges[idx_bound]
            while tmp_count < self.threshold:
                idx_bound -= 1
                tmp_count += self.counts[idx_bound-1]
                tmp_edge = self.edges[idx_bound]
            mod_counts.append(tmp_count)
            mod_edges.append(tmp_edge)
            idx_bound -= 1
        del mod_counts[-1]
        return mod_edges, mod_counts

    def initialize_indices(self):
        self.indices['peak'] = [np.argmax(self.counts)]
        self.indices['left'] = [0]
        self.indices['right'] = [len(self.edges) - 1]

    def update_edges_and_counts(self):
        res_edges, res_counts = [], []
        for ileft, ipeak, iright in zip(self.indices['left'], self.indices['peak'], self.indices['right']):
            left_edges, left_counts = self.left_of_peak(ileft, ipeak)
            right_edges, right_counts = self.right_of_peak(iright, ipeak)
            res_edges += left_edges + right_edges[::-1]
            res_counts += left_counts + right_counts[::-1]
        self.result['edges'] = np.array(res_edges)
        self.result['observed counts'] = np.array(res_counts)

class AdaptiveHistogram():

    """
    purpose:
        (1) --  initialize histogram by bin width, number of bins, or
                preset bin edges
        (2) --  apply bin threshold about central peak
        (3) --  retrieve normalization constant
    """

    def __init__(self, data, threshold=None):
        """
        data            :   type <tuple / list / array>
        threshold       :   type <float / int> or None
        """
        self.data = data
        self.threshold = threshold
        self.size = len(data)
        self.result = {}

    def initialize_edges(self, n, edges_by='number', dtype=float):
        """
        n will supersede edges_by if n is type <tuple / list / array>
        n               :   type <float / int> or type <tuple / list / array>
        edges_by        :   type <str>
        dtype           :   type <float / int>
        """
        if 'edges' not in list(self.result.keys()):
            if isinstance(edges_by, str):
                if edges_by == 'number':
                    if isinstance(n, int):
                        edges = np.linspace(min(self.data), max(self.data), n, dtype=dtype)
                    else:
                        raise ValueError("n = type <int>")
                elif edges_by == 'width':
                    if isinstance(n, (float, int)):
                        edges = np.arange(min(self.data), max(self.data) + 0.9 * n, n, dtype=dtype)
                    else:
                        raise ValueError("n = type <float / int>")
                elif edges_by == 'custom':
                    if isinstance(n, (tuple, list, np.ndarray)):
                        edges = n
                    else:
                        raise ValueError("n = type <tuple / list / array>")
                else:
                    raise ValueError("edges_by = 'number', 'width', or 'custom'")
                self.result['edges'] = edges
        if not isinstance(self.result['edges'], (tuple, list, np.ndarray)):
            raise ValueError("edges = type <tuple / list / array>")

    def simplify_edges(self, intify=False):
        """
        intify          :   type <bool>
        """
        if intify is True:
            self.result['edges'] = np.unique(np.array(self.result['edges'], dtype=int))
        elif intify is not False:
            raise ValueError("intify = True or False")

    def initialize_observed_counts(self, bias='left'):
        """
        bias            :   'left' or 'right'
        """
        if bias == 'left':
            counts, edges = np.histogram(self.data, self.result['edges'])
            self.result['edges'] = edges
        elif bias == 'right':
            counts = np.zeros(len(self.result['edges']) - 1, dtype=int)
            for idx, val in zip(*np.unique(np.searchsorted(self.result['edges'], self.data, side='left'), return_counts=True)):
                counts[idx - 1] = val
        self.result['observed counts'] = counts

    def apply_bin_threshold(self):
        if self.threshold is not None:
            BT = BinThresholding(self.result['edges'], self.result['observed counts'], self.threshold)
            BT.check_threshold()
            BT.initialize_indices()
            BT.update_edges_and_counts()
            self.result['edges'] = BT.result['edges']
            self.result['observed counts'] = BT.result['observed counts']

    def initialize_midpoints(self):
        if 'midpoints' not in list(self.result.keys()):
            edges = self.result['edges']
            self.result['midpoints'] = (edges[1:] + edges[:-1])/2

    def initialize_bin_widths(self):
        self.result['bin width'] = np.diff(self.result['edges'])

    def initialize_normalization(self):
        cnorm = np.sum(self.result['observed counts'] * self.result['bin width'])
        # print("\n{} CNORM:\n{}\n{} OBS COUNTS:\n{}\n{} BIN WIDTHS:\n{}".format(len(cnorm), cnorm, len(self.result['observed counts']), self.result['observed counts'], len(self.result['bin width']), self.result['bin width']))
        self.result['normalized counts'] = self.result['observed counts'] / cnorm
        self.result['normalization constant'] = cnorm

    def execute(self, edge, edges_by='number', edge_bias='left', intify=True):
        """
        edge            :   type <float / int> or type <tuple / list / array>
        edges_by        :   type <str>
        edge_bias       :   type <str>
        intify          :   type <bool>
        """
        self.initialize_edges(edge, edges_by)
        self.simplify_edges(intify)
        self.initialize_observed_counts(edge_bias)
        self.apply_bin_threshold()
        self.initialize_midpoints()
        self.initialize_bin_widths()
        self.initialize_normalization()

    @property
    def edges(self):
        return self.result['edges']

    @property
    def observed_counts(self):
        return self.result['observed counts']

    @property
    def midpoints(self):
        return self.result['midpoints']

    @property
    def bin_widths(self):
        return self.result['bin width']

    @property
    def normalization_constant(self):
        return self.result['normalization constant']

    @property
    def normalized_counts(self):
        """ """
        return self.result['normalized counts']

    def get_expectation_values(self, f, prms):
        """ """
        return np.array([quad(f, lbound, ubound, args=prms)[0] * self.data.size for lbound, ubound in zip(self.edges[:-1], self.edges[1:])])

    def check_threshold(self, counts=None):
        """
        counts      :   type <array> or None
        axis        :   type <int> or None
        """
        if self.threshold is None:
            print("\nBIN THRESHOLD CONDITION IS NOT ACTIVE\n")
        else:
            if counts is None:
                counts = self.observed_counts
            res = counts < self.threshold
            n = np.abs(np.sum(-1 * res))
            # n = len(counts) - np.sum(res[res > 0])
            print("\n\tBIN COUNT >= THRESHOLD:\n{}\n\n\tNUMBER OF BINS < THRESHOLD:\n{}\n".format(res, n))
