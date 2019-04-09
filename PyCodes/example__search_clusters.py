from search_methods import *

def view_result(data, res):
    """
    This function prints the original data and the data that satisfies
    the search condition(s).

    data    :   type <dict>
    res     :   type <dict>
    """
    rkeys = list(res.keys())
    dkeys = list(data.keys())
    for key in dkeys:
        if key not in rkeys:
            raise ValueError("unknown key: {}; available keys: {}".format(key, rkeys))
        rvalues = res[key]
        dvalues = data[key]
        print("-"*20)
        print(".. KEY = {}:".format(key))
        print("="*20)
        print("\n .. ORIGINAL DATA:\n{}".format(dvalues))
        print("\n .. SEARCHED DATA:\n{}\n".format(rvalues))

## SAMPLE DATA
x = np.linspace(-5, 5, 11)
y = x**2
z = np.zeros(x.size)
s = np.array(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'ijk', 'lmnop', 'q'])
data = {'x' : x, 'y' : y, 'z' : z, 's' : s}

## CREATE SAMPLE CLUSTERS BY SPLITTING DATA AT INDICES
indices = np.array([1, 4, 5, 8, 9])
data = {key : np.array(np.split(value, indices)) for key, value in data.items()}

## INSTANTIATE CLASS METHOD
Searcher = SearchClusters(data, indices, string_keys=('s',))

## RETRIEVE SEARCH RESULT
res = Searcher.search('x', 'greater than or equal', 0)
# res = Searcher.search('x', 'nearest', 4.5)
# res = Searcher.search(('s', 'y'), ('exact match', 'greater than'), ('b', 10), apply_to='all')
# res = Searcher.search(('s', 'y'), ('exact match', 'greater than'), ('a', 18), apply_to='any')
# res = Searcher.search(('cluster size', 'y'), ('exact match', 'greater than'), (2, 18), apply_to='all')
# res = Searcher.search(('cluster size', 'y'), ('exact match', 'less than'), (2, 15), apply_to='any')
# res = Searcher.search(('cluster size', 'y'), ('exact match', 'greater than'), (2, 15), apply_to='any')
# res = Searcher.search(('cluster size', 'x'), ('exact match', 'less than'), (1, -3), apply_to='all')
# res = Searcher.search('y', 'nearest', 20.5)
# res = Searcher.search('z', 'not equal', 0)

## VIEW RESULT
view_result(data, res)
