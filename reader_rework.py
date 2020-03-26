import numpy as np
import math
import dictionaries as dict
import multiprocessing as mp


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
    l = eq.copy()
    mult = np.zeros((1, 3))
    add = np.zeros((1, 3))
    for i in range(len(l)):
        if l[i][0] == "-":
            mult[0, i] = -1
            l[i] = l[i][2:]
        else:
            mult[0, i] = 1
            l[i] = l[i][1:]
    for i in range(len(l)):
        if l[i] != "" and l[i][0] == "+":
            l[i] = l[i][1:]
            w = l[i].split("/")
            number = int(w[0]) / int(w[1])
            add[0, i] = number
        elif l[i] != "" and l[i][0] == "-":
            l[i] = l[i][1:]
            w = l[i].split("/")
            number = int(w[0]) / int(w[1])
            add[0, i] = -1 * number
    return mult, add


class Molecule:

    def __init__(self, n: int):
        self.num_atoms = n
        self.atom_label = []
        self.atom_coord = np.zeros((n, 3))
        self.threshold = 0.1
        self.inertia_tensor = np.zeros((3, 3))
        self.inertia_eig_val = np.zeros((3, 1))
        self.inertia_eig_vec = np.zeros((3, 3))

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
            r = r + self.atom_coord[i:i + 1] * dict.element_weight[self.atom_label[i]]
            mass = mass + dict.element_weight[self.atom_label[i]]
        r = r / mass
        return r

    def inertia(self):
        mass_vector = self.mass_center()
        self.internal_coord = self.atom_coord - mass_vector
        for i in range(self.num_atoms):
            mass = dict.element_weight[self.atom_label[i]]
            self.inertia_tensor[0, 0] = self.inertia_tensor[0, 0] + \
                                        mass * (self.internal_coord[i, 1] ** 2 + self.internal_coord[i, 2] ** 2)
            self.inertia_tensor[1, 1] = self.inertia_tensor[1, 1] + \
                                        mass * (self.internal_coord[i, 0] ** 2 + self.internal_coord[i, 2] ** 2)
            self.inertia_tensor[2, 2] = self.inertia_tensor[2, 2] + \
                                        mass * (self.internal_coord[i, 0] ** 2 + self.internal_coord[i, 1] ** 2)
            self.inertia_tensor[0, 1] = self.inertia_tensor[0, 1] + \
                                        mass * self.internal_coord[i, 0] * self.internal_coord[i, 1]
            self.inertia_tensor[1, 0] = self.inertia_tensor[0, 1]
            self.inertia_tensor[1, 2] = self.inertia_tensor[1, 2] + \
                                        mass * self.internal_coord[i, 1] * self.internal_coord[i, 2]
            self.inertia_tensor[2, 1] = self.inertia_tensor[1, 2]
            self.inertia_tensor[0, 2] = self.inertia_tensor[0, 2] + \
                                        mass * self.internal_coord[i, 0] * self.internal_coord[i, 2]
            self.inertia_tensor[2, 0] = self.inertia_tensor[0, 2]
        self.inertia_eig_val, self.inertia_eig_vec = np.linalg.eig(self.inertia_tensor)
        for n in range(3):
            self.inertia_eig_vec[n:n+1] = self.inertia_eig_vec[n:n+1] / np.linalg.norm(self.inertia_eig_vec[n:n+1])


def copy_molecule(m1, m2: Molecule):
    for n in range(m1.num_atoms):
        m2.atom_label.append(m1.atom_label[n])
        m2.atom_coord[n:n + 1] = m1.atom_coord[n:n + 1]


def molecule_coincide(m1, m2: Molecule):
    for i in range(m1.num_atoms):
        diff = np.linalg.norm(m1.atom_coord[i:i + 1] - m2.atom_coord[i:i + 1])
        if diff < 0.05:
            return True
    return False


class CifFile:

    def __init__(self, path=""):
        # trying to open CIF file and read its contents
        try:
            file = open(path, "r")
            contents = file.readlines()
            file.close()
        except OSError:
            print("Could not open the CIF file at: ", path)
            exit(-1)
        # init cif dictionary
        self.xyz = []
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
            s = contents[index_loop].strip()
            loop_list = []
            while s != "" and s != "loop_":
                loop_list.append(s)
                index_loop = index_loop + 1
                s = contents[index_loop].strip()
            if loop_list[0] == "_symmetry_equiv_pos_as_xyz":
                for i in range(1, len(loop_list)):
                    self.xyz.append(loop_list[i].replace("'", "").split(", "))
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
        # extract xyz eq positions if they're not yet extracted
        if len(self.xyz) == 0:
            if data["_symmetry_space_group_name_Hall"] != "":
                self.xyz = dict.SymOpsHall[data["_symmetry_space_group_name_Hall"]]
            elif data["_symmetry_space_group_name_H-M"] != "":
                data["_symmetry_space_group_name_H-M"] = data["_symmetry_space_group_name_H-M"].replace("(", "")
                data["_symmetry_space_group_name_H-M"] = data["_symmetry_space_group_name_H-M"].replace(")", "")
                self.xyz = dict.SymOpsHall[dict.HM2Hall[data["_symmetry_space_group_name_H-M"]]]
            else:
                print("No symmetry detected in CIF file!")
                exit(-1)


class Cluster:

    def build(self):
        for i1 in range(len(self.cif.xyz)):
            mult, add = parse_eq_xyz(self.cif.xyz[i1])
            for x in range(-3, 3):
                for y in range(-3, 3):
                    for z in range(-3, 3):
                        new_molecule = Molecule(self.cif.asym_unit.num_atoms)
                        copy_molecule(self.cif.asym_unit, new_molecule)
                        new_molecule.atom_coord = new_molecule.atom_coord * mult
                        new_molecule.atom_coord = new_molecule.atom_coord + add
                        new_molecule.atom_coord = new_molecule.atom_coord + x * np.array([1, 0, 0])
                        new_molecule.atom_coord = new_molecule.atom_coord + y * np.array([0, 1, 0])
                        new_molecule.atom_coord = new_molecule.atom_coord + z * np.array([0, 0, 1])
                        coincide = False
                        for i2 in range(len(self.pre_molecules)):
                            if molecule_coincide(self.pre_molecules[i2], new_molecule):
                                coincide = True
                        if new_molecule.inside() and not coincide:
                            self.pre_molecules.append(new_molecule)

    def connectivity(self, line: int):
        for k in range(line+1, len(self.pre_molecules) * self.pre_molecules[0].num_atoms):
            m1 = line // self.pre_molecules[0].num_atoms
            m1n = line % self.pre_molecules[0].num_atoms
            m2 = k // self.pre_molecules[0].num_atoms
            m2n = k % self.pre_molecules[0].num_atoms
            v1 = self.pre_molecules[m1].atom_coord[m1n:m1n + 1]
            v2 = self.pre_molecules[m2].atom_coord[m2n:m2n + 1]
            dist = np.linalg.norm(np.matrix.transpose(v1) - np.matrix.transpose(v2))
            limit_dist = dict.covalent_radius[self.pre_molecules[m1].atom_label[m1n]] + \
                         dict.covalent_radius[self.pre_molecules[m2].atom_label[m2n]]
            if dist <= limit_dist:
                self.bonds[line, k] = 1
                self.bonds[k, line] = 1

    def simplify(self):
        self.mass_centers = np.zeros((len(self.molecules), 3))
        for i in range(len(self.molecules)):
            v = self.molecules[i].mass_center()
            self.mass_centers[i:i + 1] = self.molecules[i].mass_center()
            self.molecules[i].inertia()

    def multiply(self, a: int, b: int, c: int):
        self.mass_centers = np.zeros((len(self.molecules), 3))
        for i in range(len(self.molecules)):
            v = self.molecules[i].mass_center()
            self.mass_centers[i:i + 1] = self.molecules[i].mass_center()
        mass_centers_fract = np.zeros((len(self.molecules), 3))
        for i in range(len(self.molecules)):
            mass_centers_fract[i:i + 1] = np.matmul(self.cif.rev_transform, self.mass_centers[i:i + 1])
        for i in range(len(self.molecules)):
            for x in range(a + 1):
                for y in range(b + 1):
                    for z in range(c + 1):
                        point = mass_centers_fract[i:i + 1] + np.array([x, y, z])
                        if (0.0 <= point[0, 0] <= a) and (0.0 <= point[1, 0] <= b) and (0.0 <= point[2, 0] <= c):
                            new_molecule = Molecule(self.molecules[i].num_atoms)
                            copy_molecule(self.molecules[i], new_molecule)
                            new_molecule.atom_coord = new_molecule.atom_coord + np.array([x, y, z])
                            coincide = False
                            for i2 in range(len(self.molecules)):
                                if molecule_coincide(self.molecules[i2], new_molecule):
                                    coincide = True
                            if not coincide:
                                self.molecules.append(new_molecule)

    def rebuild(self):
        checked = np.zeros((len(self.pre_molecules) * self.pre_molecules[0].num_atoms, 1))
        mol = []
        flag = False
        while not flag:
            for i in range(len(self.pre_molecules) * self.pre_molecules[0].num_atoms):
                if checked[i, 0] == 0:
                    mol.append(i)
                    checked[i, 0] = 1
                    break
            inner_flag = False
            while not inner_flag:
                old_len = len(mol)
                for i1 in range(len(mol)):
                    for i2 in range(len(self.pre_molecules) * self.pre_molecules[0].num_atoms):
                        if self.bonds[i1, i2] == 1:
                            mol.append(i2)
                            checked[i2, 0] = 1
                if old_len == len(mol):
                    inner_flag = True
            flag = True
            for i in range(len(self.pre_molecules) * self.pre_molecules[0].num_atoms):
                if checked[i, 0] != 1:
                    flag = False
            new_molecule = Molecule(len(mol))
            for i in range(len(mol)):
                m = mol[i] // self.pre_molecules[0].num_atoms
                n = mol[i] % self.pre_molecules[0].num_atoms
                new_molecule.atom_label[i] = self.pre_molecules[m].atom_label[n]
                new_molecule.atom_coord[i:i + 1] = self.pre_molecules[m].atom_coord[n:n + 1]
            self.molecules.append(new_molecule)
            mol.clear()

    def clean_up(self):

    def to_cartesian(self, mol):
        for i1 in range(len(mol)):
            for i in range(23):
                vector = np.zeros((3, 1))
                vector[0, 0] = mol[i1].atom_coord[i, 0]
                vector[1, 0] = mol[i1].atom_coord[i, 1]
                vector[2, 0] = mol[i1].atom_coord[i, 2]
                vector = np.matmul(self.cif.transform, vector)
                mol[i1].atom_coord[i, 0] = vector[0, 0]
                mol[i1].atom_coord[i, 1] = vector[1, 0]
                mol[i1].atom_coord[i, 2] = vector[2, 0]

    def print_to_file(self, mol, path):
        try:
            file = open(path, "w")
        except OSError:
            print("Could not write to file at: ", path)
            exit(-1)
        file.write(str(len(mol) * mol[0].num_atoms) + "\n")
        file.write("xyz\n")
        for i1 in range(len(mol)):
            for i2 in range(mol[i1].num_atoms):
                file.write(mol[i1].atom_label[i2] + " ")
                file.write(repr(mol[i1].atom_coord[i2, 0]) + " ")
                file.write(repr(mol[i1].atom_coord[i2, 1]) + " ")
                file.write(repr(mol[i1].atom_coord[i2, 2]) + "\n")
        file.close()

    def __init__(self, a: int, b: int, c: int, path: str):
        self.pre_molecules = []
        self.cif = CifFile(path)
        self.pre_molecules.append(self.cif.asym_unit)
        self.molecules = []
        self.build()
        self.bonds = np.zeros((len(self.pre_molecules) * self.pre_molecules[0].num_atoms, len(self.pre_molecules) *
                               self.pre_molecules[0].num_atoms))
        self.to_cartesian(self.pre_molecules)
        for i in range(len(self.pre_molecules) * self.pre_molecules[0].num_atoms):
            self.connectivity(i)
        #print("Parallel")
        #if __name__ == '__main__':
        #    pool = mp.Pool(mp.cpu_count())
        #    pool.map(self.connectivity, [i for i in range(len(self.pre_molecules) * self.pre_molecules[0].num_atoms)])
        #    pool.close()
        self.print_to_file(self.pre_molecules, "D:\[work]\cluster.xyz")


cl = Cluster(0, 0, 0, "D:\[work]\Kujo\KES48.cif")
