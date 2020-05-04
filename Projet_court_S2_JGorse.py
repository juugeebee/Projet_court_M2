"""Auteur: Julie Bogoin
"""

"""Projet court M2-BIB:  Calcul de la surface exposée au solvant d'une protéine.
"""

import math
import pandas

#############
# FONCTIONS #
#############

def extraction_coordonnes(filename):
    """Extraction des lignes atoms dans une liste
    """

    coord_list = []
    with open(filename,"r") as pdb_file:
        for line in pdb_file:
            if line[0:6].strip() == 'ATOM':
                    coord_list.append(line[:-1])
    return coord_list

def dico_atomes(coord_list):
    """Création d'une listee de dictionnaires regroupant les noms de résidus, des atomes et leurs coodonnées
    """ 
    
    coordonnees = []
    longueur = len(coord_list)
    for i in range(longueur):
        dico = {}
        dico['residue_name'] = coord_list[i][17:20]
        dico['atom_name'] = coord_list[i][12:16]
        dico['X'] = coord_list[i][30:38]
        dico['Y'] = coord_list[i][38:46]
        dico['Z'] = coord_list[i][46:54]
        coordonnees.append(dico)  
    return coordonnees

#P
def creation_sphere(nbre_points):
    """Placer des points uniformement sur une sphere avec l'algorithme Golden Section Spiral: rayon = 1
    """
    
    points_list = []
    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / nbre_points
    for k in range (0, nbre_points):
        y = k * off - 1 + (off / 2)
        phi = k * inc
        points_list.append([math.cos(phi), y, math.sin(phi)])
    return points_list

def transposition_sphere (points_list, coordonnees_atome):
    """Transposer la sphere sur un atome
    """
    
    for i in range(len(points_list)):
        points_list[i][0] = points_list[i][0] + float(coordonnees_atome['X'])
        points_list[i][1] = points_list[i][1] + float(coordonnees_atome['Y'])
        points_list[i][2] = points_list[i][2] + float(coordonnees_atome['Z'])
    return points_list

def calcul_distance (points_list, coordonnees_atome):
    """Calcule la distance entre tous les points d'une sphère et un atome 
    """
    
    distances_list = []
    distances_atomes = []
    #calcul de la distance dans une fenetre de 5 atomes proches
    for i in range(len(points_list)):           
        distance = math.pow( (float(coordonnees_atome['X']) - points_list[i][0]), 2) \
                    + math.pow( (float(coordonnees_atome['Y']) - points_list[i][1]), 2) \
                    + math.pow( (float(coordonnees_atome['Z']) - points_list[i][2]), 2)
        distances_list.append(distance)
        distances_list.sort(reverse=True)
        distances_atomes.append(distances_list)
    return distances_atomes


def calcul_surface(coordonnees_atome) :
    """Calcul de la surface d'une sphère
    """
    
    rayon = 0
    dico_rayon = {'H':1.2,'C':1.8,'N':1.5,'O':1.4,'F':1.47,'S':1.8}
    for i, (atom, r) in enumerate(dico_rayon.items()):
        if(coordonnees_atome['atom_name'] == atom) :
            rayon = r
    surface = 4*math.pi*((rayon)**2)
    return surface

def calcul_accessibilite(distances_atomes, coordonnees_atome, distance_max):
    """ Calcul de la surface exposée au solvant d'un atome
    """

    # Calcul de la surface accessible
    surface = calcul_surface(coordonnees_atome)

    # Calcul du ratio exposé
    compteur = 0
    seuil = 2.8
    for i in range(len(distances_atomes)):
        for j in range(len(distances_atomes[i])):
            if((distances_atomes[i][j]) >= seuil) and ((distances_atomes[i][j]) < distance_max):
                compteur = compteur+1
    ratio = compteur/(len(distances_atomes))

    #Calcul de la surface exposée
    return surface*ratio

def accessibilite_totale(coordonnees, nbre_points, distance_max):
    """ Calcul de la surface exposée au solvant pour tous les atomes
    """

    dico = {}
    
    liste_surfaces = []
    for atome in range(len(coordonnees)) :
        sphere = creation_sphere(nbre_points)
        transposition = transposition_sphere (sphere, coordonnees[atome])
        distance = calcul_distance(transposition, coordonnees[atome])
        accessibilite = calcul_accessibilite(distance, coordonnees[atome], distance_max)
        surface = calcul_surface(coordonnees[atome])
        liste_surfaces.append(surface)

    liste_accessibilites = []
    dico[atome] = accessibilite
    for i , j in enumerate(dico.items()) :
        access =  j[1]
        liste_accessibilites.append(access)

    return liste_surfaces, liste_accessibilites, 

def tableau_final(coordonnees, nbre_points, distance_max):

    df = pandas.DataFrame(columns=['atome','X','Y','Z','residu','accessibilite','surface'])

    acces = accessibilite_totale(coordonnees, nbre_points, distance_max)
    df['accessibilite']=pandas.Series(acces[1])
    df['surface']=pandas.Series(acces[0])

    for atome in range(len(coordonnees)) :    
        df['atome']=pandas.Series(coordonnees[atome]['atom_name'])
        df['residu']=pandas.Series(coordonnees[atome]['residue_name'])
        df['X']=pandas.Series(coordonnees[atome]['X'])
        df['Y']=pandas.Series(coordonnees[atome]['Y'])
        df['Z']=pandas.Series(coordonnees[atome]['Z'])

    return df


###########################
### Programme principal ###
###########################

print('\nCalcul de la surface accessible au solvant en cours...\n')

nbre_points = 1000
distance_max = 3.0

# Extraction des coordonnees de chaque atome du pdb
coord_pdb = extraction_coordonnes('CD59_2J8B.pdb')
coordonnees = dico_atomes(coord_pdb)

df = tableau_final(coordonnees, nbre_points, distance_max)

df.to_csv('results.csv', index='False')

print('\nFichier de résultats généré.')
print('Travail terminé!\n')







