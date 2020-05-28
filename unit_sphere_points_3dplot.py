import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_sphere(file_name):
    sphere_points_x = []
    sphere_points_y = []
    sphere_points_z = []

    with open(file_name, "r") as pdb_file:
        for line in pdb_file:
            location = line.split(',')
            sphere_points_x.append(float(location[0]))
            sphere_points_y.append(float(location[1]))
            sphere_points_z.append(float(location[2]))

    return sphere_points_x, sphere_points_y, sphere_points_z


""" Main
"""

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sphere = read_sphere("./outputs/unit_sphere_points.csv")
ax.scatter(sphere[0], sphere[1], sphere[2], 'z', 5, "blue")

plt.title('Unit Sphere Points')
plt.legend(loc=2)
plt.show()
