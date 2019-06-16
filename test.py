import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

xyz = []
atoms = []
count = 0
struc_num = 1
with open('E:\\lammpcalc\\polylac\\1_100\\1_100_5.dump') as file:
    for line in file:
        if 'Time' not in line and '2087' not in line:
            a = line.split(' ')
            xyz.append([float(a[1]), float(a[2]), float(a[3])])
            atoms.append(a[0])
        else:
            count += 1
        if count / 2 - 1 == struc_num:
            break

box = 28.5
xyz = np.array(xyz)
x = xyz[:, np.newaxis, :] - xyz[np.newaxis, :, :]
x = np.abs(x)
x1 = (x <= (box/2)) * x + (x > (box/2)) * (box - x)
distances = np.sum(x1 ** 2, axis=-1) ** 0.5
distances = distances[np.triu_indices(len(distances), k=1)]
distances = np.sort(distances)[:3000]
a, bins = np.histogram(distances, bins=50)
plt.hist(bins[:-1], bins, weights=a)
plt.show()
"""
[0.76632543 0.77249653 0.84609613 0.87694655 0.9043592  0.90549355
 0.91056293 0.91292166 0.91644093 0.91819996 0.92612455 0.92743041
 0.93257349 0.93990736 0.94169567 0.94703149 0.94741242 0.94802703
 0.95159751 0.95201035]
"""