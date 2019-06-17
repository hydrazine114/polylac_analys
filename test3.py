import numpy as np
# C O H
with open('E:\\lammpcalc\\polylac\\1_100\\1_100_5.dump') as file:
    count = 0
    atom = []
    for line in file:
        if 'Time' not in line and '2087' not in line:
            atom.append(int(line[0]))
            count += 1
            if count == 2087:
                break
oxygen = []
carbo = []
for index, value in enumerate(atom):
    if value == 1:
        carbo.append(index)
    elif value == 2:
        oxygen.append(index)
oxygen = np.array(oxygen)
carbo = np.array(carbo)
np.save('oxygen.npy', oxygen)
np.save('carbo.npy', carbo)
