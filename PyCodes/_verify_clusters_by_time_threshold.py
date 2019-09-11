import numpy as np

"""
'If the time interval between two fast CMEs is less than tc,
then these CMEs can be grouped into a cluster.' - Ruzmaikin
"""


tc = 5
# elapsed_hours = np.array([1, 4, 5, 6, 8, 11, 15, 20, 21, 26, 30, 35, 41, 45, 50])
elapsed_hours = np.array([1, 6, 8, 11, 15, 20, 21, 26, 30, 35, 41, 45])
inter_exceedances = np.diff(elapsed_hours)
# condition = (tc < inter_exceedances) ## wrong
condition = (tc <= inter_exceedances) ## correct
indices = np.where(condition)[0] + 1
clusters = np.array(np.split(elapsed_hours, indices))

print("\n ** TIME THRESHOLD:\n{} (hours)\n".format(tc))
for cluster in clusters:
    print("\n .. CLUSTER (elapsed hours):\n{}\n".format(cluster))
