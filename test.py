import numpy as np
import matplotlib.pyplot as plt

oxygen = np.load('oxygen.npy')
carbo = np.load('carbo.npy')
print(len(carbo))
# pairs = []
# for o in oxygen:
#     for c in carbo:
#         pairs.append([o, c])
# pairs = np.sort(pairs[:])
# np.save('pairs.npy', np.array(pairs))
# a = np.load('pairs.npy')

