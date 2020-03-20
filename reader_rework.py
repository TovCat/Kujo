import numpy as np
import math
import dictionaries as dict


def periodic_to_float(number, count=1):
    if number.find("(") == -1:
        return float(number)
    else:
        temp = float(number[0:number.find("(")])
        i = 0
        while i < count:
            l1 = len(number[number.find(".") + 1: number.find("(")])
            l2 = len(number[number.find("(") + 1: number.find(")")])
            temp = temp + float(number[number.find("(") + 1:number.find(")")]) * 10 ** (-1 * (l1 + (i + 1) * l2))
            i = i + 1
        return temp


def transform(m: np.array, trans: np.array):
    for i in range(m.shape[0]):
        m[i:i + 1] = np.matmul(m[i:i + 1], trans)


def parse_eq_xyz(eq: list):
    l = eq
    mult = np.zeros((3, 1))
    add = np.zeros((3, 1))
    for i in range(len(l)):
        if l[i][0] == "-":
            mult[i, 0] = -1
            l[i] = l[i][2:]
        else:
            mult[i, 0] = 1
            l[i] = l[i][1:]
    for i in range(len(l)):
        if l[i][0] == "+":
            l[i] = l[i][1:]
            w = l[i].split("/")
            number = int(w[0]) / int(w[1])
            add[i, 0] = number
        elif l[i][0] == "-":
            l[i] = l[i][1:]
            w = l[i].split("/")
            number = int(w[0]) / int(w[1])
            add[i, 0] = -1 * number


class Molecule:

    def __init__(self, n: int):
        self.num_atoms = n
        self.atom_label = []
        self.atom_coord = np.zeros((n, 3))
        self.threshold = 0.1

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

    def mass_center(self):
        r = np.zeros((1, 3))
        mass = 0.0
        for i in range(self.num_atoms):
            r = r + self.atom_coord[i:i+1] * dict.element_weight[self.atom_label[i]]
            mass = mass + dict.element_weight[self.atom_label[i]]
        r = r / mass
        return r


def copy_molecule(m1, m2: Molecule):
    for n in range(m1.num_atoms):
        m2.atom_label.append(m1.atom_label[n])
        m2.atom_coord[n:n + 1] = m1.atom_coord[n:n + 1]


def molecule_coincide(m1, m2: Molecule):
    diff = 0.0
    for i in range(m1.num_atoms):
        diff = np.linalg.norm(m1.atom_coord[i:i+1] - m2.atom_coord[i:i+1])
    if diff < 0.01:
        return True
    else:
        return False


class CifFile:

    def __init__(self, path=""):
        # trying to open CIF file and read its contents
        try:
            file = open(path, "r")
        except OSError:
            print("Could not open the CIF file at: ", path)
            exit(-1)
        contents = file.readlines()
        file.close()
        # init cif dictionary
        data = {
            "_cell_length_a": "",
            "_cell_length_b": "",
            "_cell_length_c": "",
            "_cell_angle_alpha": "",
            "_cell_angle_beta": "",
            "_cell_angle_gamma": "",
            "_symmetry_space_group_name_H-M": "",
            "_symmetry_space_group_name_Hall": "",
            "loop_": []
        }
        # first survey through contents
        index = 0
        for x in contents:
            words = x.split()
            if len(words) > 0:
                if words[0] in data:
                    if words[0] != "loop_":
                        data[words[0]] = words[1]
                    else:
                        data["loop_"].append(index)
            index = index + 1
        # read the asymmetric unit
        for x in data["loop_"]:
            index_loop = x + 1
            s = contents[index_loop]
            loop_list = []
            while s != "" and s != "loop_":
                loop_list.append(s)
                index_loop = index_loop + 1
                s = contents[index_loop]
                s = s.strip()
            positions = {
                "_atom_site_type_symbol": 0,
                "_atom_site_fract_x": 0,
                "_atom_site_fract_y": 0,
                "_atom_site_fract_z": 0
            }
            loop_end = 0
            for y in loop_list:
                words_loop = y.split()
                if words_loop[0] in positions:
                    positions[words_loop[0]] = loop_list.index(y)
                if words_loop[0][0] != "_":
                    loop_end = loop_list.index(y)
                    break
            if not (0 in positions.values()):
                self.asym_unit = Molecule(len(loop_list) - loop_end)
                for i in range(loop_end, len(loop_list)):
                    words_coord = loop_list[i].split()
                    self.asym_unit.atom_label.append(words_coord[positions["_atom_site_type_symbol"]])
                    self.asym_unit.atom_coord[i - loop_end, 0] = \
                        periodic_to_float(words_coord[positions["_atom_site_fract_x"]])
                    self.asym_unit.atom_coord[i - loop_end, 1] = \
                        periodic_to_float(words_coord[positions["_atom_site_fract_y"]])
                    self.asym_unit.atom_coord[i - loop_end, 2] = \
                        periodic_to_float(words_coord[positions["_atom_site_fract_z"]])
            # convert cell parameters
            self.cell_a = periodic_to_float(data["_cell_length_a"])
            self.cell_b = periodic_to_float(data["_cell_length_b"])
            self.cell_c = periodic_to_float(data["_cell_length_c"])
            self.cell_alpha = periodic_to_float(data["_cell_angle_alpha"])
            self.cell_beta = periodic_to_float(data["_cell_angle_beta"])
            self.cell_gamma = periodic_to_float(data["_cell_angle_gamma"])
            # generate transformation matrices
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


c = CifFile("D:\[work]\Kujo\KES48.cif")
