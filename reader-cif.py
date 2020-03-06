#   reader-cif.py
#       simple .cif file reader that extracts all necessary data from the file
#   Python 3.8 + math, numpy
#
#   Written by Igor Koskin
#   e-mail: osingran@yandex.ru
#
#   last update: 06.03.2020
#   version: 0.01


import math
import numpy as np

nucl = ["H", "C", "O"]
nucl_mass = [1, 12, 16]
nucl_vdw = [0.120, 0.170, 0.152]
precision = 1


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


def hall_parse(s=""):
    if s == "-P 2ybc":
        return 14


def hm_parse(s=""):
    if s == "P12(1)/c1":
        return 14


def hm_alt_parse(s=""):
    if s == "P 21/c":
        return 14


class Molecule:

    def __init__(self, n: int):
        self.num_atoms = n
        self.atom_label = []
        self.atom_coord = np.zeros((n, 3))
        self.inertia_tensor = np.zeros((3, 3))
        self.inertia_eig_val = np.zeros((3, 1))
        self.inertia_eig_vec = np.zeros((3, 3))
        self.threshold = 0.0

    def inside(self):
        a1 = 1 + self.threshold
        a2 = 0 - self.threshold
        r = False
        for n in range(self.num_atoms):
            if (self.atom_coord[n, 0] <= a1) and (self.atom_coord[n, 1] <= a1) and \
                    (self.atom_coord[n, 2] <= a1) and (self.atom_coord[n, 0] >= a2) and \
                    (self.atom_coord[n, 1] >= a2) and (self.atom_coord[n, 2] >= a2):
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

    def translate(self, x_t, y_t, z_t):
        for i in range(self.num_atoms):
            self.atom_coord[i, 0] = self.atom_coord[i, 0] + x_t
            self.atom_coord[i, 1] = self.atom_coord[i, 1] + y_t
            self.atom_coord[i, 2] = self.atom_coord[i, 2] + z_t

    def mirror(self, axis, c1):
        for i in range(self.num_atoms):
            self.atom_coord[i, axis] = self.atom_coord[i, axis] + 2 * (c1 - self.atom_coord[i, axis])

    def rotate(self, axis, order, c1, c2):
        angle = np.deg2rad(360 / order)
        cosa = np.cos(angle)
        sina = np.sin(angle)
        rel1 = 0
        rel2 = 0
        rel3 = 0
        if axis == 0:
            matrix = np.array([[1, 0, 0], [0, cosa, -sina], [0, sina, cosa]])
        elif axis == 1:
            matrix = np.array([[cosa, 0, sina], [0, 1, 0], [-sina, 0, cosa]])
        elif axis == 2:
            matrix = np.array([[cosa, -sina, 0], [sina, cosa, 0], [0, 0, 1]])
        for i in range(self.num_atoms):
            if axis == 0:
                rel1 = self.atom_coord[i, 0]
                rel2 = self.atom_coord[i, 1] - c1
                rel3 = self.atom_coord[i, 2] - c2
            elif axis == 1:
                rel1 = self.atom_coord[i, 0] - c1
                rel2 = self.atom_coord[i, 1]
                rel3 = self.atom_coord[i, 2] - c2
            elif axis == 2:
                rel1 = self.atom_coord[i, 0] - c1
                rel2 = self.atom_coord[i, 1] - c2
                rel3 = self.atom_coord[i, 2]
            vector = np.zeros((1, 3))
            vector[0, 0] = rel1
            vector[0, 1] = rel2
            vector[0, 2] = rel3
            vector = np.matmul(vector, matrix)
            if axis == 0:
                rel1 = vector[0, 0]
                rel2 = vector[0, 1] + c1
                rel3 = vector[0, 2] + c2
            elif axis == 1:
                rel1 = vector[0, 0] + c1
                rel2 = vector[0, 1]
                rel3 = vector[0, 2] + c2
            elif axis == 2:
                rel1 = vector[0, 0] + c1
                rel2 = vector[0, 1] + c2
                rel3 = vector[0, 2]
            self.atom_coord[i, 0] = rel1
            self.atom_coord[i, 1] = rel2
            self.atom_coord[i, 2] = rel3

    def screw(self, axis, order, step, c1, c2):
        self.rotate(axis, order, c1, c2)
        if axis == 0:
            self.translate(step, 0, 0)
        elif axis == 1:
            self.translate(0, step, 0)
        elif axis == 2:
            self.translate(0, 0, step)

    def glide(self, axis, step, c1):
        self.mirror(axis, c1)
        if axis == 0:
            self.translate(step, 0, 0)
        elif axis == 1:
            self.translate(0, step, 0)
        elif axis == 2:
            self.translate(0, 0, step)

    def rotoinversion(self, axis, order, x, y, z):
        if axis == 0:
            self.rotate(axis, order, y, z)
        elif axis == 1:
            self.rotate(axis, order, x, z)
        elif axis == 2:
            self.rotate(axis, order, x, y)
        self.inverse(x, y, z)

    def mass_center(self):
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        mass = 0.0
        for i in range(self.num_atoms):
            sum_x = sum_x + self.atom_coord[i, 0] * nucl_mass[nucl.index(self.atom_label[i])]
            sum_y = sum_y + self.atom_coord[i, 1] * nucl_mass[nucl.index(self.atom_label[i])]
            sum_z = sum_z + self.atom_coord[i, 2] * nucl_mass[nucl.index(self.atom_label[i])]
            mass = mass + nucl_mass[nucl.index(self.atom_label[i])]
        sum_x = sum_x / mass
        sum_y = sum_y / mass
        sum_z = sum_z / mass
        r = np.zeros((3, 1))
        r[0, 0] = sum_x
        r[1, 0] = sum_y
        r[2, 0] = sum_z
        return r

    def inertia(self):
        xx = 0.0
        yy = 0.0
        zz = 0.0
        xy = 0.0
        yz = 0.0
        xz = 0.0
        for i in range(self.num_atoms):
            xx = xx + nucl_mass[nucl.index(self.atom_label[i])] * (self.atom_coord[i, 1] * self.atom_coord[i, 1] +
                                                                   self.atom_coord[i, 2] * self.atom_coord[i, 2])
            yy = yy + nucl_mass[nucl.index(self.atom_label[i])] * (self.atom_coord[i, 0] * self.atom_coord[i, 0] +
                                                                   self.atom_coord[i, 2] * self.atom_coord[i, 2])
            zz = zz + nucl_mass[nucl.index(self.atom_label[i])] * (self.atom_coord[i, 0] * self.atom_coord[i, 0] +
                                                                   self.atom_coord[i, 1] * self.atom_coord[i, 1])
            xy = xy + nucl_mass[nucl.index(self.atom_label[i])] * self.atom_coord[i, 0] * self.atom_coord[i, 1]
            yz = yz + nucl_mass[nucl.index(self.atom_label[i])] * self.atom_coord[i, 1] * self.atom_coord[i, 2]
            xz = xz + nucl_mass[nucl.index(self.atom_label[i])] * self.atom_coord[i, 0] * self.atom_coord[i, 2]
        self.inertia_tensor[0, 0] = xx
        self.inertia_tensor[0, 1] = xy
        self.inertia_tensor[0, 2] = xz
        self.inertia_tensor[1, 0] = xy
        self.inertia_tensor[1, 1] = yy
        self.inertia_tensor[1, 2] = yz
        self.inertia_tensor[2, 0] = xz
        self.inertia_tensor[2, 1] = yz
        self.inertia_tensor[2, 2] = zz
        self.inertia_eig_val, self.inertia_eig_vec = np.linalg.eig(self.inertia_tensor)


def copy_molecule(m1, m2: Molecule):
    for n in range(m1.num_atoms):
        m2.atom_label.append(m1.atom_label[n])
        m2.atom_coord[n, 0] = m1.atom_coord[n, 0]
        m2.atom_coord[n, 1] = m1.atom_coord[n, 1]
        m2.atom_coord[n, 2] = m1.atom_coord[n, 2]


def molecule_coincide(m1, m2: Molecule):
    diff = 0.0
    for i in range(m1.num_atoms):
        diff = (m1.atom_coord[i, 0] - m2.atom_coord[i, 0]) ** 2 + (m1.atom_coord[i, 1] - m2.atom_coord[i, 1]) ** 2 + \
               (m1.atom_coord[i, 2] - m2.atom_coord[i, 2]) ** 2
    if diff < 0.01:
        return True
    else:
        return False


# axis (a - 0, b - 1, c - 2); order; coord1; coord2 - rotations, screws
# axis (a - 0, b - 1, c - 2); coord1 - mirrors, glides
class CellSymmetry:
    def __init__(self, table_number: int):
        if table_number == 1:
            self.no_inverses = True
            self.no_screws = True
            self.no_glides = True
            self.no_rotations = True
            self.no_rotoinversions = True
            self.no_mirrors = True
        if table_number == 2:
            self.inverses = np.array([[0, 0, 0], [0, 0, 0.5], [0, 0, 1.0], [0, 0.5, 0], [0, 0.5, 0.5], [0, 0.5, 1.0],
                                      [0, 1.0, 0], [0, 1.0, 0.5], [0, 1.0, 1.0], [0.5, 0, 0], [0.5, 0, 0.5],
                                      [0.5, 0, 1.0], [0.5, 0.5, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 1.0], [0.5, 1.0, 0],
                                      [0.5, 1.0, 0.5], [0.5, 1.0, 1.0], [1.0, 0, 0], [1.0, 0, 0.5], [1.0, 0, 1.0],
                                      [1.0, 0.5, 0], [1.0, 0.5, 0.5], [1.0, 0.5, 1.0], [1.0, 1.0, 0], [1.0, 1.0, 0.5],
                                      [1.0, 1.0, 1.0]])
            self.no_inverses = False
            self.no_screws = True
            self.no_glides = True
            self.no_rotations = True
            self.no_rotoinversions = True
            self.no_mirrors = True
        if table_number == 3:
            self.rotations = np.array([[1, 2, 0, 0], [1, 2, 0, 0.5], [1, 2, 0, 1], [1, 2, 0.5, 0], [1, 2, 0.5, 0.5],
                                       [1, 2, 0.5, 1], [1, 2, 1, 0], [1, 2, 1, 0.5], [1, 2, 1, 1]])
            self.no_inverses = True
            self.no_screws = True
            self.no_glides = True
            self.no_rotations = False
            self.no_rotoinversions = True
            self.no_mirrors = True
        if table_number == 4:
            self.screws = np.array([[1, 2, 0, 0], [1, 2, 0, 0.5], [1, 2, 0, 1], [1, 2, 0.5, 0], [1, 2, 0.5, 0.5],
                                    [1, 2, 0.5, 1], [1, 2, 1, 0], [1, 2, 1, 0.5], [1, 2, 1, 1]])
            self.no_inverses = True
            self.no_screws = False
            self.no_glides = True
            self.no_rotations = True
            self.no_rotoinversions = True
            self.no_mirrors = True
        if table_number == 5:
            self.rotations = np.array([[1, 2, 0, 0], [1, 2, 0, 0.5], [1, 2, 0, 1], [1, 2, 0.5, 0], [1, 2, 0.5, 0.5],
                                       [1, 2, 0.5, 1], [1, 2, 1, 0], [1, 2, 1, 0.5], [1, 2, 1, 1]])
            self.screws = np.array([[1, 2, 0.25, 0], [1, 2, 0.25, 0.5], [1, 2, 0.25, 1], [1, 2, 0.75, 0],
                                    [1, 2, 0.75, 0.5], [1, 2, 0.75, 1]])
            self.no_inverses = True
            self.no_screws = False
            self.no_glides = True
            self.no_rotations = False
            self.no_rotoinversions = True
            self.no_mirrors = True
        if table_number == 6:
            self.mirrors = np.array([[1, 0.5]])
            self.no_inverses = True
            self.no_screws = True
            self.no_glides = True
            self.no_rotations = True
            self.no_rotoinversions = True
            self.no_mirrors = False

        if table_number == 14:
            self.inverses = np.array([[0, 0, 0], [0, 0, 0.5], [0, 0, 1.0], [0, 0.5, 0], [0, 0.5, 0.5], [0, 0.5, 1.0],
                                      [0, 1.0, 0], [0, 1.0, 0.5], [0, 1.0, 1.0], [0.5, 0, 0], [0.5, 0, 0.5],
                                      [0.5, 0, 1.0], [0.5, 0.5, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 1.0], [0.5, 1.0, 0],
                                      [0.5, 1.0, 0.5], [0.5, 1.0, 1.0], [1.0, 0, 0], [1.0, 0, 0.5], [1.0, 0, 1.0],
                                      [1.0, 0.5, 0], [1.0, 0.5, 0.5], [1.0, 0.5, 1.0], [1.0, 1.0, 0], [1.0, 1.0, 0.5],
                                      [1.0, 1.0, 1.0]])
            self.screws = np.array([[1, 2, 0, 0.25], [1, 2, 0, 0.75], [1, 2, 0.5, 0.25], [1, 2, 0.5, 0.75],
                                    [1, 2, 1.0, 0.25], [1, 2, 1.0, 0.75]])
            self.glides = np.array([[1, 0.25], [1, 0.75]])

            self.no_inverses = False
            self.no_screws = False
            self.no_glides = False
            self.no_rotations = True
            self.no_rotoinversions = True
            self.no_mirrors = True


class CifFile:

    def __init__(self, init_path=""):
        file = open(init_path, "r")
        contents = file.readlines()
        file.close()
        self.it_number = 0
        index = -1
        symmetry_detected = False
        for x in contents:
            index = index + 1
            words = x.split()
            if len(words) != 0:
                for i in range(len(words)):
                    words[i].strip()
                if words[0][0] == "_":
                    if words[0] == "_cell_length_a":
                        self.cell_a = periodic_to_float(words[1], precision)
                    if words[0] == "_cell_length_b":
                        self.cell_b = periodic_to_float(words[1], precision)
                    if words[0] == "_cell_length_c":
                        self.cell_c = periodic_to_float(words[1], precision)
                    if words[0] == "_cell_angle_alpha":
                        self.cell_alpha = periodic_to_float(words[1], precision)
                    if words[0] == "_cell_angle_beta":
                        self.cell_beta = periodic_to_float(words[1], precision)
                    if words[0] == "_cell_angle_gamma":
                        self.cell_gamma = periodic_to_float(words[1], precision)
                    if words[0] == "_cell_angle_volume":
                        self.cell_volume = periodic_to_float(words[1], precision)
                    if words[0] == "_space_group_IT_number" and self.it_number != 0:
                        self.it_number = int(words[1])
                        symmetry_detected = True
                    if words[0] == "_symmetry_space_group_name_H-M_alt" or words[0] == "_symmetry_space_group_name_Hall" or \
                            words[0] == "_symmetry_space_group_name_H-M":
                        words[1].replace("'", "")
                        if words[0] == "_symmetry_space_group_name_H-M_alt" and symmetry_detected is False:
                            self.it_number = hm_alt_parse(words[1])
                        elif words[0] == "_symmetry_space_group_name_H-M" and symmetry_detected is False:
                            self.it_number = hm_parse(words[1])
                        elif words[0] == "_symmetry_space_group_name_Hall" and symmetry_detected is False:
                            self.it_number = hall_parse(words[1])
                if words[0] == "loop_":
                    index_loop = index + 1
                    s = contents[index_loop]
                    loop_list = []
                    while s != " \n":
                        loop_list.append(s)
                        index_loop = index_loop + 1
                        s = contents[index_loop]
                    type_pos = 0
                    fract_x_pos = 0
                    fract_y_pos = 0
                    fract_z_pos = 0
                    loop_end_found = False
                    loop_end = 0
                    for y in loop_list:
                        words_loop = y.split()
                        if words_loop[0] == "_atom_site_type_symbol":
                            type_pos = loop_list.index(y)
                        if words_loop[0] == "_atom_site_fract_x":
                            fract_x_pos = loop_list.index(y)
                        if words_loop[0] == "_atom_site_fract_y":
                            fract_y_pos = loop_list.index(y)
                        if words_loop[0] == "_atom_site_fract_z":
                            fract_z_pos = loop_list.index(y)
                        if loop_end_found is False and words_loop[0][0] != "_":
                            loop_end = loop_list.index(y)
                            loop_end_found = True
                    if type_pos != 0 and fract_x_pos != 0 and fract_y_pos != 0 and fract_z_pos != 0:
                        self.asym_unit = Molecule(len(loop_list) - loop_end)
                        for i in range(loop_end, len(loop_list)):
                            words_coord = loop_list[i].split()
                            self.asym_unit.atom_label.append(words_coord[type_pos])
                            self.asym_unit.atom_coord[i - loop_end, 0] = periodic_to_float(words_coord[fract_x_pos])
                            self.asym_unit.atom_coord[i - loop_end, 1] = periodic_to_float(words_coord[fract_y_pos])
                            self.asym_unit.atom_coord[i - loop_end, 2] = periodic_to_float(words_coord[fract_z_pos])
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
