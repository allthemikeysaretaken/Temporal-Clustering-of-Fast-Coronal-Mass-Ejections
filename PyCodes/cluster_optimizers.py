from cluster_configurations import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

class DensityMethods(DistanceMetrics):

    def __init__(self, data):
        """ """
        super().__init__(data)
        self._radius = None
        self._minimum_npoints = None
        self._visited = []
        self._groups = []
        self._noise = {}
        self._clusters = {}
        self._cluster_method = None
        self.update_method_type('density')

    def __str__(self):
        data = np.array(self.data)
        string = '\n** METHOD:\n\t{}\nRADIUS: {}\n\tMINIMUM POINT CONDITION FOR CLUSTER: {}\n'.format(self.cluster_method, self.radius, self.minimum_npoints)
        string = '{}\n** (shape={}) DATA:\n{}\n'.format(string, data.shape, data)
        string = '{}\n .. K={} CLUSTERS:\n'.format(string, self.k)
        for ith_cluster in range(self.k):
            key = str(ith_cluster)
            cluster = self.clusters[key]
            string = '{}\ncluster #{} (shape={}):\n{}\n'.format(string, ith_cluster, cluster.shape, cluster)
        for ith_noise in range(len(self.noise)):
            key = str(ith_noise)
            noise = self.noise[key]
            string = '{}\nnoise #{} (shape={}):\n{}\n'.format(string, ith_noise, noise.shape, noise)
        return string

    @property
    def cluster_method(self):
        return self._cluster_method

    @property
    def dpoints(self):
        return np.array([DataPoint(index, position) for index, position in enumerate(self.data.T)])

    def reset_search_parameters(self):
        self._radius = None
        self._minimum_npoints = None

    def reset_clusters(self):
        self._visited = []
        self._groups = []
        self._noise = {}
        self._clusters = {}

    def initialize_parameters(self, radius, minimum_npoints=None):
        """ """
        self._radius = radius
        if minimum_npoints is None:
            self._minimum_npoints = 2 * self.ndim
        else:
            self._minimum_npoints = minimum_npoints

    @property
    def radius(self):
        return self._radius

    @property
    def minimum_npoints(self):
        return self._minimum_npoints

    @property
    def volume(self):
        return np.pi**(self.ndim / 2) * self._radius**self.ndim / gamma(self.ndim/2 + 1)

    def get_points_within_volume(self, dpoint):
        """ """
        distances = self.f_distance(dpoint.position)
        condition = (distances <= self.radius)
        return self.dpoints[condition]

    def validate_point(self, dpoint):
        """ """
        if dpoint.index in self._visited:
            return True
        else:
            return False

    def validate_neighborhood(self, dpoints):
        """ """
        return np.array([self.validate_point(dpoint) for dpoint in dpoints])

    def collect_group(self, dpoint, group):
        """ """
        if self.validate_point(dpoint) == False:
            self._visited.append(dpoint.index)
            group.append(dpoint)
            neighborhood = self.get_points_within_volume(dpoint)
            condition = self.validate_neighborhood(neighborhood)
            indices = np.invert(condition)
            dpoints = neighborhood[indices]
            for jpoint in dpoints:
                group = self.collect_group(jpoint, group)
        return group

    def initialize_groups(self, method):
        """ """
        if method == 'dbscan':
            for dpoint in self.dpoints:
                condition = self.validate_point(dpoint)
                if condition == False:
                    group = self.collect_group(dpoint, group=[])
                    self._groups.append(group)
        elif method == 'mean-shift':
            raise ValueError("not yet implemented")

    @property
    def groups(self):
        return self._groups

    def classify_groups(self):
        ith_noise, ith_cluster = 0, 0
        for ith_group, group in enumerate(self._groups):
            if len(group) >= self.minimum_npoints:
                key = str(ith_cluster)
                self.clusters[key] = PositionalSystem(group).extract_positions().T
                ith_cluster += 1
            else:
                key = str(ith_noise)
                self._noise[key] = PositionalSystem(group).extract_positions().T
                ith_noise += 1

    @property
    def clusters(self):
        return self._clusters

    @property
    def noise(self):
        return self._noise

    @property
    def k(self):
        return len(self._clusters)

    @property
    def cluster_keys(self):
        return np.arange(self.k).astype(str)

    @property
    def m(self):
        return len(self._noise)

    @property
    def noise_keys(self):
        return np.arange(self.m).astype(str)

    def run_dbscan(self, distance_metric=None, radius=None, minimum_npoints=None):
        """ """
        if distance_metric is not None:
            self.update_distance_metric(distance_metric)
        if radius is not None:
            self.reset_search_parameters()
            self.initialize_parameters(radius, minimum_npoints)
        if self._radius is None:
            raise ValueError("invalid radius: None; choose radius > 0")
        self.reset_clusters()
        try:
            self.initialize_groups(method='dbscan')
            self.classify_groups()
            status = True
            msg = None
        except:
            status = False
            msg = 'dbscan was unsuccessful'
        return status, msg

    def run_mean_shift(self, distance_metric=None, radius=None, minimum_npoints=None):
        """ """
        raise ValueError("not yet implemented")

    @property
    def available_density_methods(self):
        res = {}
        res['dbscan'] = lambda distance_metric=None, radius=None, minimum_npoints=None : self.run_dbscan(distance_metric, radius, minimum_npoints)
        res['mean-shift'] = lambda distance_metric=None, radius=None, minimum_npoints=None : self.run_mean_shift(distance_metric, radius, minimum_npoints)
        return res

    def dispatch(self, method, identifiers=None, *args, **kwargs):
        """ """
        available_keys = list(self.available_density_methods.keys())
        if method in available_keys:
            f = self.available_density_methods[method]
            self._cluster_method = method
        else:
            raise ValueError("unknown method: {}; available methods: {}".format(method, available_keys))
        status, msg = f(*args, **kwargs)
        self.update_identification_labels(self.k, identifiers)
        return status, msg

class DendrogramConfiguration(DistanceMatrix):

    def __init__(self, data):
        """ """
        super().__init__(data)
        self._indices = None
        self._levels = []
        self._locations = []
        self._branch_points = []
        self._connector_points = []
        self._connector_labels = []
        self._row_column_pairs = []

    @property
    def indices(self):
        return self._indices

    @property
    def levels(self):
        return np.array(self._levels)

    @property
    def locations(self):
        return np.array(self._locations)

    @property
    def branch_points(self):
        return np.array(self._branch_points)

    @property
    def connector_points(self):
        return np.array(self._connector_points)

    @property
    def connector_labels(self):
        return np.array(self._connector_labels)

    @property
    def row_column_pairs(self):
        return np.array(self._row_column_pairs)

    def initialize_indices(self):
        self._indices = np.arange(len(self._original_distance_matrix))

    def update_levels(self, row, col):
        """ """
        level = self._reducible_distance_matrix[row, col]
        self._levels.append(level)

    def update_locations(self, row, col):
        """ """
        loc = (self._indices[row], self._indices[col])
        self._locations.append(loc)
        self._indices = np.delete(self._indices, obj=col)

    @staticmethod
    def get_branch_point(xleft, xright, ybottom, level, ytop):
        """ """
        xbranch = [xleft, xleft, xright, xright]
        ybranch = [ybottom, level, level, ytop]
        branch_point = [xbranch, ybranch]
        return branch_point

    @staticmethod
    def get_connector_point(xleft, xright, level):
        """ """
        xconnector = np.mean([xleft, xright])
        connector_point = [xconnector, level]
        return connector_point

    def initialize_dendrogram(self, levels, locations, identifiers):
        """ """
        self.update_identification_labels(len(self.levels)+1, identifiers)
        mapping = {xi : {'x' : xj, 'y' : 0, 'label' : label} for xi, xj, label in zip(np.arange(self.ordered_distances.size), self.ordered_distances, self.identification_labels[self.ordered_distances])}
        for idx, (level, location) in enumerate(zip(levels, locations)):
            row, col = location[0], location[1]
            xleft, xright = mapping[row]['x'], mapping[col]['x']
            ybottom, ytop = mapping[row]['y'], mapping[col]['y']
            branch_point = self.get_branch_point(xleft, xright, ybottom, level, ytop)
            connector_point = self.get_connector_point(xleft, xright, level)
            self._branch_points.append(branch_point)
            self._connector_points.append(connector_point)
            mapping[row]['x'] = connector_point[0]
            mapping[row]['y'] = level
            row_label, col_label = mapping[row]['label'], mapping[col]['label']
            identifiers = np.array([row_label, col_label])
            if np.all(identifiers == None):
                connector_label = None
            # elif np.any(identifiers == None):
            #     connector_label = next(item for item in identifiers if item is not None)
            else:
                # connector_label = '{}/{}'.format(*identifiers)
                connector_label = "/".join(identifiers)
            mapping[row]['label'] = connector_label
            self._connector_labels.append(connector_label)

    def subview_dendrogram(self, ax, facecolors, alpha=0.8, label_dendrogram=False, tick_connectors=False, ticksize=7, labelsize=8):
        """ """
        npoints = len(self.connector_points)
        for idx, (branch_point, connector_point, level, identifier, c) in enumerate(zip(self.branch_points, self.connector_points, self.levels, self.connector_labels, facecolors)):
            ax.plot(*branch_point, color=c, alpha=alpha)
            if idx < npoints-1:
                if level == self.levels[idx+1]:
                    label = None
                else:
                    label = identifier
            else:
                label = identifier
            ax.scatter(*connector_point, facecolor=c, label=label, alpha=alpha)
            if label_dendrogram == True:
                ax.text(*connector_point, identifier, color=c, fontsize=labelsize)
        if tick_connectors == True:
            xticks = np.arange(len(self.connector_labels))
            ax.set_xticks(xticks)
            ax.set_xticklabels(self.connector_labels)
            ax.set_xlabel('Ordered Labels')
        else:
            ax.set_xticks([])
            ax.set_xlabel('')
        ax.set_ylabel('Levels')
        return ax

class AgglomerativeHierarchicalMethods(DendrogramConfiguration):

    def __init__(self, data):
        """ """
        super().__init__(data)
        self._mapping = {}
        self._linkage_criterion = ''
        self._dendrogram = {}
        self.update_method_type('hierarchical')

    def __str__(self):
        data = np.array(self.data)
        string = '\n** AGGLOMERATIVE {} METHOD\n\tLINKAGE CRITERION: {}\n'.format(self.method_type, self.linkage_criterion)
        string = '{}\n** (shape={}) DATA:\n{}\n'.format(string, data.shape, data)
        string = '{}\n .. ORIGINAL DISTANCE MATRIX (shape={}):\n{}\n'.format(string, self.original_distance_matrix.shape, self.original_distance_matrix)
        string = '{}\n .. {} DENDROGRAM- ORDERED DISTANCES:\n{}\n'.format(string, len(self.ordered_distances), self.ordered_distances)
        string = '{}\n .. {} DENDROGRAM- LEVELS:\n{}\n'.format(string, len(self.levels), self.levels)
        string = '{}\n .. DENDROGRAM- LOCATIONS (shape={}):\n{}\n'.format(string, self.locations.shape, self.locations)
        string = '{}\n .. DENDROGRAM- BRANCH POINTS (shape={}):\n{}\n'.format(string, self.branch_points.shape, self.branch_points)
        string = '{}\n .. DENDROGRAM- CONNECTOR POINTS (shape={}):\n{}\n'.format(string, self.connector_points.shape, self.connector_points)
        return string

    @property
    def dendrogram(self):
        return self._dendrogram

    @property
    def available_linkage_criteria(self):
        res = {}
        res['single'] = np.min
        res['complete'] = np.max
        return res

    @property
    def f_linkage(self):
        return self.available_linkage_criteria[self.linkage_criterion]

    def update_linkage_criterion(self, linkage_criterion=None):
        """ """
        if linkage_criterion is not None:
            available_criteria = list(self.available_linkage_criteria.keys())
            if linkage_criterion in available_criteria:
                self._linkage_criterion = linkage_criterion
            else:
                raise ValueError("invalid linkage_criterion: {}; available criteria: {}".format(linkage_criterion, available_criteria))

    @property
    def linkage_criterion(self):
        return self._linkage_criterion

    def run_single_linkage(self, distance_metric):
        """ """
        self.update_linkage_criterion('single')
        self.reset_distance_matrix()
        self.initialize_distance_matrix(distance_metric)
        self.initialize_indices()
        status = True
        msg = None
        # try:
        while self._reducible_distance_matrix.size >= 4: ## 2x2 array
            mask = self.get_masked_distance_matrix(mask_values=0)
            locations = np.where(mask == self.f_linkage(mask))
            reduced_locations = locations[0]
            if len(reduced_locations) > 2:
                msg = "WARNING: ties are present"
            row, col = reduced_locations[0], reduced_locations[1]
            self.update_levels(row, col)
            self.update_locations(row, col)
            data = np.array([mask[row], mask[col]])
            values = self.f_linkage(data, axis=0)
            self.merge_subdimensions(row, values)
            self.reduce_subdimensions(col)
        # except:
        #     status = False
        #     msg = 'cross-check not yet implemented'
        #     raise ValueError(msg)
        return status, msg

    def run_complete_linkage(self, distance_metric=None):
        """ """
        raise ValueError("not yet implemented")
        # self.update_linkage_criterion('complete')
        # self.reset_distance_matrix()
        # self.initialize_distance_matrix(distance_metric)
        # self.initialize_indices()
        # status = True
        # msg = None
        # return status, msg

    @property
    def available_linkage_methods(self):
        res = {}
        res['single'] = lambda distance_metric=None : self.run_single_linkage(distance_metric=distance_metric)
        res['complete'] = lambda distance_metric=None : self.run_complete_linkage(distance_metric=distance_metric)
        return res

    def dispatch(self, method, identifiers=None, *args, **kwargs):
        """ """
        available_keys = list(self.available_linkage_methods.keys())
        if method in available_keys:
            f = self.available_linkage_methods[method]
        else:
            raise ValueError("unknown method: {}; available methods: {}".format(method, available_keys))
        status, msg = f(*args, **kwargs)
        self.initialize_dendrogram(self.levels, self.locations, identifiers)
        return status, msg




















##
