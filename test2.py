import numpy as np
import matplotlib.pylab as plt

size = 2048
ind = lambda i1, i2: i1 * (size - 1) - np.sum([x for x in range(i1 + 1)]) + i2 - 1


def showhist(a):
    a = np.sort(a)[:3000]
    a, bins = np.histogram(a, bins=50)
    plt.hist(bins[:-1], bins, weights=a)
    plt.show()


a = np.load('distances/structure_0.npy')
# showhist(a)
