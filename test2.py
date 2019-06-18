import numpy as np
import matplotlib.pylab as plt

size = 2048
ind = lambda i1, i2: i1 * (size - 1) - np.sum([x for x in range(i1 + 1)]) + i2 - 1


def showhist(a):
    a = np.sort(a)
    a, bins = np.histogram(a, bins=50)
    plt.hist(bins[:-1], bins, weights=a)
    plt.show()


distances = np.load('distances/structure_zero.npy')
pairs = np.load('pairs.npy')
c_o = []
dists = []
for i in pairs:
    d = distances[ind(i[0], i[1])]
    if d < 2:
        c_o.append([i, d])
        dists.append(d)

print(len(c_o))
showhist(dists)
