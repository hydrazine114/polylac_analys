import numpy as np
import matplotlib.pylab as plt

a = np.load('distances/structure_0.npy')
a = np.sort(a)[:3000]
a, bins = np.histogram(a, bins=50)
plt.hist(bins[:-1], bins, weights=a)
plt.show()

