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

## INSTANTIATE CLASS METHOD
Searcher = SearchEjecta(data, string_keys=('s',))

## RETRIEVE SEARCH RESULT
res = Searcher.search('y', 'greater than', 0)
# res = Searcher.search('y', 'greater than', 50, fkeys='cumulative sum')
# res = Searcher.search('x', 'greater than', 5, fkeys='delta')
# res = Searcher.search('y', 'greater than', 5, fkeys='delta')
# res = Searcher.search('y', 'greater than', 5, fkeys='delta', use_abs=True)
# res = Searcher.search('s', 'exact match', 'c')
# res = Searcher.search(('x', 'y', 's'), ('less than', 'greater than', 'not equal'), (0, 1, 'c'))
# res = Searcher.search(('x', 'y'), ('less than', 'greater than'), (0, 1))
# res = Searcher.search('s', 'exact match', 'c')
# res = Searcher.search(('s', 'x'), ('exact match', 'greater than or equal'), ('c', -3))
# res = Searcher.search('x', 'nearest', 0.5)
# res = Searcher.search('x', 'nearest forward', 0.5)
# res = Searcher.search('x', 'nearest forward', 26)
# res = Searcher.search('x', 'nearest backward', 0.5)
# res = Searcher.search('x', 'nearest backward', -26)
# res = Searcher.search('z', 'not equal', 0)

## VIEW RESULT
view_result(data, res)
