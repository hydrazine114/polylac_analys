import numpy as np


a = np.array(range(1, 10)).reshape((3, 3))
b = a * a
b = (b % 2 == 0) * b
print(b)
