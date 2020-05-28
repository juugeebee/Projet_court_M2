import math

"""The number of points to compute on an atom's sphere"""
SPHERE_POINTS_COUNT = 100

"""The number of closest atoms to take into account for solvant accessibility (0 to take all atoms into account)"""
CLOSEST_ATOMS_COUNT = 50


"""
Math functions
"""


def calculate_surface_area(radius):
    return 4 * math.pi * (radius ** 2)


def calculate_unit_sphere():
    """Saff Kuijlaars algorithm: calculates evenly distributed points on the unit_sphere
    https://web.archive.org/web/20120107030109/http://cgafaq.info/wiki/Evenly_distributed_points_on_sphere#Spirals
    """
    sphere_points = []

    s = 3.6 / math.sqrt(SPHERE_POINTS_COUNT)
    dz = 2.0 / SPHERE_POINTS_COUNT
    long = 0
    z = 1 - dz / 2
    for k in range(SPHERE_POINTS_COUNT):
        r = math.sqrt(1 - z * z)
        location = (math.cos(long) * r, math.sin(long) * r, z)
        sphere_points.append(location)
        z = z - dz
        long = long + s / r

    # save unit sphere to external file to make 3d plot
    save_sphere_to_file(sphere_points, './outputs/unit_sphere_points.csv')

    return sphere_points


def scale_and_move_sphere(unit_sphere_points, radius_w, location):
    """Scales the unit sphere to the atom's radius and moves it to the atom's location
    """
    sphere_points = []
    for unit_sphere_point in unit_sphere_points:
        sphere_point = (
            unit_sphere_point[0] * radius_w + location[0],
            unit_sphere_point[1] * radius_w + location[1],
            unit_sphere_point[2] * radius_w + location[2])
        sphere_points.append(sphere_point)
    return sphere_points


def sort_by_closest_atoms(atom, all_atoms, number_of_closest_atoms):
    """Returns a new list containing all_atoms sorted from closest to farthest from atom
    """
    sorted_all_atoms = sorted(all_atoms, key=lambda another_atom: calculate_distance(atom.location, another_atom.location))
    return sorted_all_atoms[1:number_of_closest_atoms + 1]


def calculate_distance(point, other_point):
    return math.sqrt((other_point[0] - point[0]) ** 2 + (other_point[1] - point[1]) ** 2 + (other_point[2] - point[2]) ** 2)


"""
Utility functions
"""


def save_sphere_to_file(sphere_points, file_name):
    with open(file_name, 'w') as sphere_points_file:
        for sphere_point in sphere_points:
            sphere_points_file.write(str(sphere_point[0]) + ',' + str(sphere_point[1]) + ',' + str(sphere_point[2]) + '\n')


def save_atoms_to_file(atoms, file_name):
    with open(file_name, 'w') as sphere_points_file:
        for atom in atoms:
            sphere_points_file.write(atom.name + ',' + str(atom.location[0]) + ',' + str(atom.location[1]) + ',' + str(
                atom.location[2]) + ',' + str(atom.radius_w) + '\n')


"""
Classes
"""


class Atom:
    """ Atom class containing one atom's following properties:
    - index
    - name
    - residue_name
    - location
    - radius_w (van der walls atom radius)
    - surface_area
    - atom_sphere
    - atom_plus_solvant_sphere
    - accessible_points
    - accessible_points_count
    - accessible_surface_area
    - accessible_surface_area_pct
    """

    index = 1

    atoms_radius_w = {
        'H': 1.2,
        'C': 1.7,
        'N': 1.55,
        #TODO: wikipedia says 1.52A for Oxygen
        #'O': 1.52,
        'O': 1.4,
        'F': 1.47,
        'P': 1.8,
        'S': 1.8,
        'CL': 1.75,
        'CU': 1.4
    }
    solvant_radius_w = atoms_radius_w['O']
    unit_sphere_points = calculate_unit_sphere()

    def __init__(self, pdb_atom_line):
        self.index = Atom.index
        Atom.index = Atom.index + 1

        self.name = pdb_atom_line[12:14].strip()
        self.residue_name = pdb_atom_line[17:20]
        self.location = (float(pdb_atom_line[30:38].strip()), float(pdb_atom_line[38:46].strip()), float(pdb_atom_line[46:54].strip()))
        self.radius_w = Atom.atoms_radius_w[self.name]
        self.surface_area = calculate_surface_area(self.radius_w)
        self.atom_sphere = scale_and_move_sphere(Atom.unit_sphere_points, self.radius_w, self.location)
        self.atom_plus_solvant_sphere = scale_and_move_sphere(Atom.unit_sphere_points, self.radius_w + Atom.solvant_radius_w, self.location)
        self.accessible_points = [True] * len(self.atom_plus_solvant_sphere)
        self.accessible_points_count = 0
        self.accessible_surface_area_pct = 0
        self.accessible_surface_area = 0

    @staticmethod
    def read_pdb(file_name):
        """Parse pdb file and return the list of Atom objects
        """
        atoms = []
        with open(file_name, "r") as pdb_file:
            for line in pdb_file:
                if line[0:6].strip() == 'ATOM':
                    # Remove end of line
                    atoms.append(Atom(line[:-1]))
        return atoms

    def calculate_accessible_points_and_surface(self, all_atoms):
        self.calculate_accessible_points(all_atoms)
        self.calculate_accessible_surface()

        points_count = len(self.accessible_points)
        accessibility_pct = (self.accessible_points_count / points_count) * 100

        print('Résultat atome #' + str(self.index) + '/' + str(len(all_atoms)) + ':',
              'points accessibles=' + str(self.accessible_points_count) + '/' + str(points_count) + ',',
              'surface accessible=' + f'{self.accessible_surface_area:.1f}' + 'A²/' + f'{self.surface_area:.1f}' + 'A²,',
              'pourcentage accessible=' + f'{accessibility_pct:.1f}' + '%')

    def calculate_accessible_points(self, all_atoms):
        """For every point from the atom's atom_plus_solvant_sphere, try to
        see if there is another atom's sphere point at a distance smaller than
        solvant radius.
        If there is one, then mark the point as not accessible
        """
        other_atoms = all_atoms
        if CLOSEST_ATOMS_COUNT > 0:
            other_atoms = sort_by_closest_atoms(self, all_atoms, CLOSEST_ATOMS_COUNT)
        for i in range(len(self.atom_plus_solvant_sphere)):
            point = self.atom_plus_solvant_sphere[i]
            for other_atom in other_atoms:
                if not self.accessible_points[i]:
                    # this point is not accessible: no need to calculate anything
                    break
                if self.index == other_atom.index:
                    # don't compare with itself
                    continue
                for other_atom_point in other_atom.atom_sphere:
                    distance = calculate_distance(point, other_atom_point)
                    if distance < Atom.solvant_radius_w:
                        self.accessible_points[i] = False
                        break

    def calculate_accessible_surface(self):
        """Calculates the accessible surface area, based on the number of
        accessible points
        """
        accessible_points_count = 0
        for is_accessible_point in self.accessible_points:
            if is_accessible_point:
                accessible_points_count = accessible_points_count + 1

        self.accessible_points_count = accessible_points_count
        self.accessible_surface_area_pct = (accessible_points_count / len(atom.accessible_points))
        self.accessible_surface_area = atom.surface_area * atom.accessible_surface_area_pct


"""
Main Program
"""


print('\n-START Génération des atomes...')
atoms = Atom.read_pdb('./inputs/CD59_2J8B.pdb')
# save atoms to external file to make 3d plot
save_atoms_to_file(atoms, './outputs/atoms.csv')
atoms_count = len(atoms)
print('-END Génération des atomes:', atoms_count, 'atomes générés')


print('\n-START Calcul des points et surfaces accessibles...')
for atom in atoms:
    atom.calculate_accessible_points_and_surface(atoms)
print('-END Calcul des points et surfaces accessibles')


print('\n-START Sauvegarde des resultats...')
results_file_name = './outputs/results_' + str(SPHERE_POINTS_COUNT) + 'points_'
if CLOSEST_ATOMS_COUNT > 0:
    results_file_name = results_file_name + str(CLOSEST_ATOMS_COUNT) + 'closestAtoms'
else:
    results_file_name = results_file_name + 'allAtoms'
results_file_name = results_file_name + '.csv'

total_surface_area = 0
total_accessible_surface_area = 0
with open(results_file_name, 'w') as results_file:
    results_file.write('atome,residu,X,Y,Z,surface,surface_accessible,accessibilite\n')
    for atom in atoms:
        points_count = len(atom.accessible_points)
        accessible_points_count = atom.accessible_points_count
        surface_area = atom.surface_area
        accessible_surface_area = atom.accessible_surface_area
        accessibility_pct = (accessible_points_count/points_count) * 100

        total_surface_area = total_surface_area + surface_area
        total_accessible_surface_area = total_accessible_surface_area + accessible_surface_area

        results_file.write(atom.name + ','
                           + atom.residue_name + ','
                           + str(atom.location[0]) + ','
                           + str(atom.location[1]) + ','
                           + str(atom.location[2]) + ','
                           + str(surface_area) + ','
                           + str(accessible_surface_area) + ','
                           + str(accessibility_pct) + '%'
                           + '\n')

    total_accessibility_pct = (total_accessible_surface_area / total_surface_area) * 100
    results_file.write('\nsurfaceTotale,surfaceTotaleAccessible,accessibiliteTotale\n')
    results_file.write(str(total_surface_area) + ',' + str(total_accessible_surface_area) + ',' + str(total_accessibility_pct) + '%\n')

    print('\nResultats:')
    print('Surface totale:', f'{total_surface_area:.1f}' + 'A²')
    print('Surface totale accessible par le solvant:', f'{total_accessible_surface_area:.1f}' + 'A²')
    print("Pourcentace d'accessibilité au solvant:", f'{total_accessibility_pct:.1f}' + '%')
print('\n-END Sauvegarde des resultats')

print('\nProgramme terminé!\n')
