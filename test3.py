import numpy as np
size = 10
d = np.array([np.array(range(size)) + x for x in range(0, size ** 2, size)])
print(d)
d = d[np.triu_indices(size, k=1)]
i1 = 6
i2 = 9
n = i1 * (size - 1) - np.sum([x for x in range(i1+1)]) + i2 - 1
print(d[n])
