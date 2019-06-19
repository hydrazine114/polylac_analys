import numpy as np
import matplotlib.pyplot as plt
from coolfuncs import read_gro

box = 28.5
period = [1.35, 1.45]


def showhist(a):
    a = np.sort(a)
    a = a[:500]
    a, bins = np.histogram(a, bins=50)
    plt.hist(bins[:-1], bins, weights=a)
    plt.show()


def calc_dist(system, n=0):
    system = np.array(system)
    x = system[:, np.newaxis, :] - system[np.newaxis, :, :]
    x = np.abs(x)
    x1 = (x <= (box / 2)) * x + (x > (box / 2)) * (box - x)
    distances = np.sum(x1 ** 2, axis=-1) ** 0.5
    distances = distances[np.triu_indices(len(distances), k=1)]
    # np.save('distances/structure_' + str(n), distances)
    distances = distances.reshape(len(distances))
    return distances


def choose_co(system):
    new_system = []
    for line in system:
        if line[2][0] == 'C' or line[2][0] == 'O':
            new_system.append(line[4:])
    return np.array(new_system)


def choose_co2(system):
    new_system = []
    for line in system:
        if line[0] == 1 or line[0] == 2:
            new_system.append(line[1:])
    return np.array(new_system)


def read_dump(input):
    system = []
    count = 0
    with open(input) as file:
        for line in file:
            if count > 1:
                coords = line.split(' ')
                coords[0] = int(coords[0])
                for i in range(1, 4):
                    coords[i] = float(coords[i])
                system.append(coords)
            count += 1
            if count > 921:
                return system


per = [0, 0.105, 0.115, 0.125, 0.138, 0.147, 0.16]
if True:
    system = read_gro('AFTEROPTwat.gro')[:920]
    coords = choose_co(system)
else:
    system = read_dump('E:\\lammpcalc\\polylac\\1_100\\1_100_5.dump')
    coords = choose_co2(system)

    for i, v in enumerate(per):
        per[i] = v * 10

d = calc_dist(coords)
nums = lambda x1, x2: np.sum((d > x1) * (d < x2))

for i in range(len(per) - 1):
    print(per[i], '-', per[i + 1], ': ', nums(per[i], per[i + 1]))
showhist(d)

"""
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
"""
