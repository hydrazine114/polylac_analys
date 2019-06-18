import numpy as np
import matplotlib.pyplot as plt
import pickle

box = 28.5
period = [1.35, 1.45]


def calc_dist(system, n):
    system = np.array(system)
    x = system[:, np.newaxis, :] - system[np.newaxis, :, :]
    x = np.abs(x)
    x1 = (x <= (box / 2)) * x + (x > (box / 2)) * (box - x)
    distances = np.sum(x1 ** 2, axis=-1) ** 0.5
    distances = distances[np.triu_indices(len(distances), k=1)]
    # np.save('distances/structure_' + str(n), distances)
    return np.sum((distances > period[0]) * (distances < period[1]))


def read():
    with open('E:\\lammpcalc\\polylac\\1_100\\1_100_5.dump') as file:
        xyz = []
        atoms = []
        count = 0
        struc_num = 0
        count2 = 0
        step = 10
        for line in file:
            if 'Time' in line:
                struc_num += 1
            if 'Time' not in line and '2087' not in line and struc_num % step == 0:
                a = line.split(' ')
                xyz.append([float(a[1]), float(a[2]), float(a[3])])
            else:
                count += 1
            if len(xyz) == 2087:
                calc_dist(xyz, count2)
                xyz = []
                atoms = []
                count2 += 1
                print('{:7.3f}%'.format(struc_num / 200))


with open('1_100_5.gro') as file:
    system = []


