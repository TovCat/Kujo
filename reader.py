import numpy as np
import math

nucl = ["H", "C", "O"]
nucl_mass = [1, 12, 16]
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


class Atom:
    def __init__(self, label_init=" ", x_init=0.0, y_init=0.0, z_init=0.0):
        self.label = label_init
        self.x = x_init
        self.y = y_init
        self.z = z_init

    def inverse(self, x_i, y_i, z_i):
        self.x = self.x + 2 * (x_i - self.x)
        self.y = self.y + 2 * (y_i - self.y)
        self.z = self.z + 2 * (z_i - self.z)

    def translate(self, x_t, y_t, z_t):
        self.x = self.x + x_t
        self.y = self.y + y_t
        self.z = self.z + z_t

    def rotate_2(self, axis, x_a, y_a, z_a):
        if axis == "a":
            self.y = self.y + 2 * (y_a - self.y)
            self.z = self.z + 2 * (z_a - self.z)
        elif axis == "b":
            self.x = self.x + 2 * (x_a - self.x)
            self.z = self.z + 2 * (z_a - self.z)
        elif axis == "c":
            self.y = self.y + 2 * (y_a - self.y)
            self.z = self.z + 2 * (z_a - self.z)

    def mirror(self, axis, coord):
        vector = np.zeros((3, 1))
        trs = np.zeros((3, 3))
        if axis == "a":
            trs[0, 0] = -1
            trs[1, 1] = 1
            trs[2, 2] = 1
            vector[1, 0] = self.x - coord
            vector[2, 0] = self.y
            vector[3, 0] = self.z
            vector = np.matmul(trs, vector)
            vector[1, 0] = vector[1, 0] + coord
            return vector
        elif axis == "b":
            trs[0, 0] = 1
            trs[1, 1] = -1
            trs[2, 2] = 1
            vector[1, 0] = self.x
            vector[2, 0] = self.y - coord
            vector[3, 0] = self.z
            vector = np.matmul(trs, vector)
            vector[1, 0] = vector[2, 0] + coord
            return vector
        elif axis == "c":
            trs[0, 0] = 1
            trs[1, 1] = 1
            trs[2, 2] = -1
            vector[1, 0] = self.x
            vector[2, 0] = self.y
            vector[3, 0] = self.z - coord
            vector = np.matmul(trs, vector)
            vector[1, 0] = vector[2, 0] + coord
            return vector


class Molecule:

    def __init__(self):
        self.atoms = []
        self.type2 = False

    def append_atom(self, atom: Atom):
        self.atoms.append(Atom(atom.label, atom.x, atom.y, atom.z))

    def mass_center(self):
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        mass = 0.0
        for i in range(len(self.atoms)):
            sum_x = sum_x + self.atoms[i].x * nucl_mass[nucl.index(self.atoms[i].label)]
            sum_y = sum_y + self.atoms[i].y * nucl_mass[nucl.index(self.atoms[i].label)]
            sum_z = sum_z + self.atoms[i].z * nucl_mass[nucl.index(self.atoms[i].label)]
            mass = mass + nucl_mass[nucl.index(self.atoms[i].label)]
        sum_x = sum_x / mass
        sum_y = sum_y / mass
        sum_z = sum_z / mass
        r = np.zeros((3, 1))
        r[0, 0] = sum_x
        r[1, 0] = sum_y
        r[2, 0] = sum_z
        return r

    def complete(self):
        initial_length = len(self.atoms)
        i = 0
        while i < initial_length:
            new_atom = Atom(self.atoms[i].label, self.atoms[i].x, self.atoms[i].y, self.atoms[i].z)
            new_atom.inverse(1, 1, 0)
            self.atoms.append(new_atom)
            i = i + 1


class CifFile:
    #################
    # CELL PARAMETERS
    cell_a = 0.0
    cell_b = 0.0
    cell_c = 0.0
    cell_alpha = 0.0
    cell_beta = 0.0
    cell_gamma = 0.0
    #################
    contents = []  # CIF-file contents
    transform = np.zeros((3, 3))  # fractional to cartesian transformation matrix

    def __init__(self, init_path=""):
        self.path = init_path
        file = open(self.path, "r")
        self.contents = file.readlines()
        file.close()

    def read_cell_parameters(self):
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

    def read_asymm_unit(self, m: Molecule):
        for x in self.contents:
            index = self.contents.index(x)
            if 149 < index < 174:
                x.strip()
                words = x.split(" ")
                x1 = periodic_to_float(words[2], precision)
                y1 = periodic_to_float(words[3], precision)
                z1 = periodic_to_float(words[4], precision)
                m.atoms.append(Atom(words[1], x1, y1, z1))

    def generate_transform(self):
        cosa = np.cos(np.deg2rad(self.cell_alpha))
        cosb = np.cos(np.deg2rad(self.cell_beta))
        cosg = np.cos(np.deg2rad(self.cell_gamma))
        sing = np.sin(np.deg2rad(self.cell_gamma))
        volume = 1.0 - cosa ** 2.0 - cosb ** 2.0 - cosg ** 2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)
        self.transform[0, 0] = self.cell_a
        self.transform[0, 1] = self.cell_b * cosg
        self.transform[0, 2] = self.cell_c * cosb
        self.transform[1, 1] = self.cell_b * sing
        self.transform[1, 2] = self.cell_c * (cosa - cosb * cosg) / sing
        self.transform[2, 2] = self.cell_c * volume / sing
