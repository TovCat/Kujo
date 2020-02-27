import reader as r
import numpy as np
import os
import math

# TODO symmetry operations should be a method of the molecule, not the atom
# TODO Atom coordinates should be stored as numpy array
# TODO transformation method in molecule class
# TODO write in file method in molecule class

cluster = []  # list of molecules
size_a = 1  # size of the cluster
size_b = 1
size_c = 1

molecule1 = r.Molecule()  # in our CIF file there're 2 unique molecules
molecule2 = r.Molecule()
file = r.CifFile(os.path.abspath(os.getcwd()) + "\KES48.cif")
file.read_cell_parameters()  # reading cell parameters
file.generate_transform()  # generating a transformation matrix from fractional to cartesian coordinates
file.read_asymm_unit(molecule1)
molecule1.complete()  # making a complete molecule out of an asymmetric unit of primitive cell: it's easier to replicate molecules in such a way
for i in range(len(molecule1.atoms)):
    molecule2.append_atom(molecule1.atoms[i])
for i in range(len(molecule2.atoms)):  # do some symmetry stuff to create molecule 2
    molecule2.atoms[i].rotate_2("b", 1, 1, 0.25)
    molecule2.atoms[i].translate(0, -0.5, 0)
cluster.append(molecule1)
cluster.append(molecule2)
for i1 in range(size_a + 1):  # down below are molecule translations
    for i2 in range(size_b + 1):
        for i3 in range(size_c + 1):
            if i1 == 1 and i2 == 1 and i3 == 0:
                pass
            else:
                new_molecule = r.Molecule()
                for i in range(len(molecule1.atoms)):
                    new_molecule.append_atom(molecule1.atoms[i])
                for i in range(len(new_molecule.atoms)):
                    new_molecule.atoms[i].translate(i1 - 1, i2 - 1, i3)
                cluster.append(new_molecule)
                del new_molecule
for i1 in range(size_a + 1):
    for i2 in range(size_b):
        for i3 in range(size_c):
            if i1 == 1 and i2 == 0 and i3 == 0:
                pass
            else:
                new_molecule = r.Molecule()
                for i in range(len(molecule2.atoms)):
                    new_molecule.append_atom(molecule2.atoms[i])
                for i in range(len(new_molecule.atoms)):
                    new_molecule.atoms[i].translate(i1 - 1, i2, i3)
                new_molecule.type2 = True
                cluster.append(new_molecule)
                del new_molecule
for i in range(len(cluster)):  # transform fractional coordinates to cartesian
    for i1 in range(len(cluster[i].atoms)):
        vector = np.zeros((3, 1))
        vector[0, 0] = cluster[i].atoms[i1].x
        vector[1, 0] = cluster[i].atoms[i1].y
        vector[2, 0] = cluster[i].atoms[i1].z
        vector = np.matmul(file.transform, vector)
        cluster[i].atoms[i1].x = vector[0, 0]
        cluster[i].atoms[i1].y = vector[1, 0]
        cluster[i].atoms[i1].z = vector[2, 0]

mu1 = np.zeros((3, 1))
mu2 = np.zeros((3, 1))
mu1[0, 0] = 5.09777
mu1[1, 0] = 0.07312
mu1[2, 0] = 0.04931
mu2 = mu1
rev_transform = np.zeros((3, 3))
rev_transform = np.linalg.inv(file.transform)
mu2 = np.matmul(rev_transform, mu2)
mu2[0, 0] = -1 * mu2[0, 0]
mu2[2, 0] = -1 * mu2[2, 0]
mu2 = np.matmul(file.transform, mu2)
mu1 = mu1.transpose()
mu2 = mu2.transpose()

d2eA = 0.20819434  # Debye to eÅ
bohr3 = 0.529177210  # bohr in Å
Eh2icm = 219500  # Hartree to cm-1
A = Eh2icm*bohr3*d2eA**2

H = np.zeros((len(cluster), len(cluster)))
for n in range(len(cluster)):
    for m in range(n + 1, len(cluster)):
        v1 = cluster[n].mass_center()
        v2 = cluster[m].mass_center()
        dv = v1 - v2
        dv = dv.transpose()
        r = np.linalg.norm(dv)
        r3 = math.pow(r, 3)
        r5 = math.pow(r, 5)
        if not cluster[n].type2:
            mu_n = mu1
        else:
            mu_n = mu2
        if not cluster[m].type2:
            mu_m = mu1
        else:
            mu_m = mu2
        J = np.inner(mu_n, mu_m)/r3 + (np.inner(mu_n, dv)*np.inner(dv, mu_m))/r5
        H[n, m] = J*A
        H[m, n] = H[n, m]


f = open("cluster.xyz", "w")  #  print out the resulting molecular cluster
num_atoms = len(cluster[0].atoms) * len(cluster)
f.write(repr(num_atoms) + "\n")
f.write("Molecular cluster of FP5\n")
for i in range(len(cluster)):
    for i1 in range(len(cluster[i].atoms)):
        f.write(cluster[i].atoms[i1].label + " " + repr(cluster[i].atoms[i1].x) + " " + repr(cluster[i].atoms[i1].y) + " " + repr(cluster[i].atoms[i1].z) + "\n")
f.close()

f = open("H.txt", "w")
for n in range(len(cluster)):
    for m in range(len(cluster)):
        f.write(repr(H[n, m]) + " ")
    f.write("\n")
f.close()

