import numpy as np
xyz = []
atoms = []
with open('test.dump') as file:
    for line in file:
        if 'Time' not in line and '2087' not in line:
            a = line.split(' ')
            xyz.append([float(a[1]), float(a[2]), float(a[3])])
            atoms.append(a[0])
xyz = np.array(xyz)
distances = np.sum((xyz[:, np.newaxis, :] - xyz[np.newaxis, :, :]) ** 2, axis=-1)
distances = distances[np.triu_indices(len(distances), k=1)]
print(np.triu_indices(10))
