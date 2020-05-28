import matplotlib.pyplot as plt


def read_atom(file_name):
    sphere_points_x = []
    sphere_points_y = []
    sphere_points_z = []
    sphere_points_r = []
    sphere_points_c = []

    atoms_colors = {
        'C': "red",
        'N': "green",
        'O': "blue",
        'S': "yellow",
    }

    with open(file_name, "r") as pdb_file:
        for line in pdb_file:
            location = line.split(',')
            sphere_points_x.append(float(location[1]))
            sphere_points_y.append(float(location[2]))
            sphere_points_z.append(float(location[3]))
            sphere_points_r.append(float(location[4]) * 20)
            sphere_points_c.append(atoms_colors.get(location[0]))

    return sphere_points_x, sphere_points_y, sphere_points_z, sphere_points_r, sphere_points_c


""" Main
"""

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
atoms = read_atom("./outputs/atoms.csv")
ax.scatter(atoms[0], atoms[1], atoms[2], 'z', atoms[3], atoms[4])

plt.title('Atoms')
plt.legend(loc=2)
plt.show()
