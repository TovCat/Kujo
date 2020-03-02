import numpy as np
import math

precision = 1
cluster = []
nucl = ["H", "C", "O"]
nucl_mass = [1, 12, 16]


def periodic_to_float(number, count=0):
    if number.find("(") == -1:
        return float(number)
    else:
        temp = float(number[0:number.find("(")])
        i = 0
        while i < count:
            l1 = len(number[number.find(".") + 1: number.find("(")])
            l2 = len(number[number.find("(") + 1: number.find(")")])
            temp = temp + float(number[number.find("(") + 1:number.find(")")]) * math.pow(10, -1 * (l1 + (i + 1) * l2))
            i = i + 1
        return temp


class Molecule:

    def __init__(self, N: int):
        self.num_atoms = N
        self.atom_label = []
        self.atom_coord = np.zeros((N, 3))

    def inside(self):
        for n in range(self.num_atoms):
            if (self.atom_coord[n, 0] <= 1.1) and (self.atom_coord[n, 1] <= 1.1) and \
                    (self.atom_coord[n, 2] <= 1.1) and (self.atom_coord[n, 0] >= -0.1) and \
                    (self.atom_coord[n, 1] >= -0.1) and (self.atom_coord[n, 2] >= -0.1):
                r = True
                break
            else:
                r = False
        return r

    def inverse(self, x_i, y_i, z_i):
        for i in range(self.num_atoms):
            self.atom_coord[i, 0] = self.atom_coord[i, 0] + 2 * (x_i - self.atom_coord[i, 0])
            self.atom_coord[i, 1] = self.atom_coord[i, 1] + 2 * (y_i - self.atom_coord[i, 1])
            self.atom_coord[i, 2] = self.atom_coord[i, 2] + 2 * (z_i - self.atom_coord[i, 2])


class CifFile:

    def __init__(self, init_path=""):
        # ----------------------------------- #
        # -----------OPEN CIF FILE----------- #
        # ----------------------------------- #
        self.path = init_path
        file = open(self.path, "r")
        self.contents = file.readlines()
        file.close()
        # ----------------------------------- #
        # --------READ CELL PARAMETERS------- #
        # ----------------------------------- #
        for x in self.contents:
            x = x.replace(" ", "")
            x = x.replace("(", "")
            x = x.replace(")", "")
            if "_cell_length_" in x:
                if x[13] == "a":
                    self.cell_a = periodic_to_float(x[14:], precision)
                elif x[13] == "b":
                    self.cell_b = periodic_to_float(x[14:], precision)
                elif x[13] == "c":
                    self.cell_c = periodic_to_float(x[14:], precision)
            if "_cell_angle_" in x:
                if x[12] == "a":
                    self.cell_alpha = periodic_to_float(x[17:], precision)
                elif x[12] == "b":
                    self.cell_beta = periodic_to_float(x[16:], precision)
                elif x[12] == "g":
                    self.cell_gamma = periodic_to_float(x[17:], precision)
        # ----------------------------------- #
        # --GENERATE TRANSFORMATION MATRICES- #
        # ----------------------------------- #
        cosa = np.cos(np.deg2rad(self.cell_alpha))
        cosb = np.cos(np.deg2rad(self.cell_beta))
        cosg = np.cos(np.deg2rad(self.cell_gamma))
        sing = np.sin(np.deg2rad(self.cell_gamma))
        volume = 1.0 - cosa ** 2.0 - cosb ** 2.0 - cosg ** 2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)
        self.transform = np.zeros((3, 3))
        self.transform[0, 0] = self.cell_a
        self.transform[0, 1] = self.cell_b * cosg
        self.transform[0, 2] = self.cell_c * cosb
        self.transform[1, 1] = self.cell_b * sing
        self.transform[1, 2] = self.cell_c * (cosa - cosb * cosg) / sing
        self.transform[2, 2] = self.cell_c * volume / sing
        self.rev_transform = np.linalg.inv(self.transform)

    def clear_contents(self):
        del self.contents

    def read_asymm_unit(self, m: Molecule):
        for x in self.contents:
            index = self.contents.index(x)
            if 149 < index < 174:
                x.strip()
                words = x.split(" ")
                m.atom_coord[index - 150, 0] = periodic_to_float(words[2], precision)
                m.atom_coord[index - 150, 1] = periodic_to_float(words[3], precision)
                m.atom_coord[index - 150, 2] = periodic_to_float(words[4], precision)
                m.atom_label.append(words[1])


def copy_molecule(m1, m2: Molecule):
    for n in range(m1.num_atoms):
        m2.atom_label.append(m1.atom_label[n])
        m2.atom_coord[n, 0] = m1.atom_coord[n, 0]
        m2.atom_coord[n, 1] = m1.atom_coord[n, 1]
        m2.atom_coord[n, 2] = m1.atom_coord[n, 2]


def molecule_coincide(m1, m2: Molecule):
    r = False
    for i in range(m1.num_atoms):
        if (m1.atom_coord[i, 0] == m2.atom_coord[i, 0]) and (m1.atom_coord[i, 1] == m2.atom_coord[i, 1]) and \
                (m1.atom_coord[i, 2] == m2.atom_coord[i, 2]):
            r = True
            break
    return r


inverses = np.array([[0, 0, 0], [0, 0, 0.5], [0, 0, 1.0], [0, 0.5, 0], [0, 0.5, 0.5], [0, 0.5, 1.0],
                     [0, 1.0, 0], [0, 1.0, 0.5], [0, 1.0, 1.0], [0.5, 0, 0], [0.5, 0, 0.5], [0.5, 0, 1.0],
                     [0.5, 0.5, 0],
                     [0.5, 0.5, 0.5], [0.5, 0.5, 1.0], [0.5, 1.0, 0], [0.5, 1.0, 0.5], [0.5, 1.0, 1.0],
                     [1.0, 0, 0], [1.0, 0, 0.5], [1.0, 0, 1.0], [1.0, 0.5, 0], [1.0, 0.5, 0.5], [1.0, 0.5, 1.0],
                     [1.0, 1.0, 0], [1.0, 1.0, 0.5], [1.0, 1.0, 1.0]])
# axis (a - 0, b - 1, c - 2); order; coord1; coord2
screws = np.array([[2, 2, 0, 0.25], [2, 2, 0, 0.75], [2, 2, 0.5, 0.25], [2, 2, 0.5, 0.75], [2, 2, 1.0, 0.25],
                   [2, 2, 1.0, 0.75]])

file = CifFile("D:\[work]\Kujo\KES48.cif")
cluster.append(Molecule(23))
file.read_asymm_unit(cluster[0])
i = 0
while i < len(cluster):
    new_molecule = Molecule(23)
    copy_molecule(cluster[i], new_molecule)
    for i1 in range(inverses.shape[0]):
        new_molecule1 = Molecule(23)
        copy_molecule(new_molecule, new_molecule1)
        new_molecule1.inverse(inverses[i1, 0], inverses[i1, 1], inverses[i1, 2])
        coincide = False
        for i3 in range(len(cluster)):
            if molecule_coincide(cluster[i3], new_molecule1):
                coincide = True
        if new_molecule1.inside() and not coincide:
            cluster.append(new_molecule1)
        del new_molecule1
    del new_molecule
    i = i + 1

for i1 in range(len(cluster)):
   for i in range(23):
       vector = np.zeros((3, 1))
       vector[0, 0] = cluster[i1].atom_coord[i, 0]
       vector[1, 0] = cluster[i1].atom_coord[i, 1]
       vector[2, 0] = cluster[i1].atom_coord[i, 2]
       vector = np.matmul(file.transform, vector)
       cluster[i1].atom_coord[i, 0] = vector[0, 0]
       cluster[i1].atom_coord[i, 1] = vector[1, 0]
       cluster[i1].atom_coord[i, 2] = vector[2, 0]
       print(cluster[i1].atom_label[i], cluster[i1].atom_coord[i, 0], cluster[i1].atom_coord[i, 1],
             cluster[i1].atom_coord[i, 2])
