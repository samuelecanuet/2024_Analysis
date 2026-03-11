import matplotlib.pyplot as plt


import numpy as np

list_guess = [0]*100
list_out = [0]*100
# reding file with n x y
with open("../out_centers_distorted.txt", "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        if len(parts) != 3:
            continue
        n, x, y = map(float, parts)
        if x == 0 or y == 0:
            list_out[int(n)] = (0, 0)
        plt.scatter(x, y, color='red')
        list_guess[int(n)] = (x, y)


# reding file with n x y
with open("../out_centers_2024.txt", "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        if len(parts) != 3:
            continue
        n, x, y = map(float, parts)
        if x == 0 or y == 0:
            list_out[int(n)] = (0, 0)
        plt.scatter(x, y, color='blue')
        list_out[int(n)] = (x, y)

####### PARAMETERS
with open("../fit_params_2024.txt", "r") as f:
    fit_params = []
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        fit_params.append([float(part) for part in parts])

fX_par = fit_params[0]
fY_par = fit_params[1]


def fit_points(x, par):
    NparX = 3  # Assuming NparX is 3 based on the original code
    NparY = 3  # Assuming NparY is 3 based on the original code
    res = 0.0
    for i in range(NparX + 1):
        for j in range(NparY + 1):
            ij = i * (NparY + 1) + j
            res += par[ij] * (x[0] ** i) * (x[1] ** j)
    return res

### cut all the element in list_guess and list_out that not tuple
list_guess_np = np.array([x for x in list_guess if isinstance(x, tuple)])
list_guess_np_x = [list_guess_np[i][0] for i in range(len(list_guess_np))]
list_out_np = np.array([x for x in list_out if isinstance(x, tuple)])
list_out_np_x = [list_out_np[i][0] for i in range(len(list_out_np))]

list_out_calib_x = [fit_points(x, fX_par) for x in list_out_np]
list_out_calib_y = [fit_points(y, fY_par) for y in list_out_np]

list_out_calib = [(x, y) for x, y in zip(list_out_calib_x, list_out_calib_y)]

plt.scatter(list_out_calib_x, list_out_calib_y, color='green', label='Out Centers')

#################################################################
cmap = plt.get_cmap('plasma')
# Assign colors based on the list_out_calib_y VALUE (not index)
colors = [cmap(list_out_calib_y[i] / max(list_out_calib_y)) for i in range(len(list_out_calib_y))]
list_guess_np_x = [list_guess_np_x[i] * list_out_calib_x[i]/max(list_out_calib_x) for i in range(len(list_out_calib_x))]
fig, ax = plt.subplots(figsize=(9, 8))
for i in range(len(list_out_calib)):
    if list_guess_np_x[i] == 0 or list_out_calib_x[i] == 0:
        continue
    plt.scatter(list_out_calib_x[i]/max(list_out_calib_x), list_guess_np_x[i]/max(list_guess_np_x)-list_out_calib_x[i]/max(list_out_calib_x), color=colors[i], label='Guess X')

plt.show()