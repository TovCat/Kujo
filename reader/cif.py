from os import getcwd
import numpy as np
import utility.dictionaries
from copy import deepcopy


def periodic_to_float(number, count=1):
    """
    Round periodic number to regular float with accuracy count (1 by default).
    """
    if number.find("(") == -1:
        return float(number)
    else:
        temp = float(number[0:number.find("(")])
        for i in range(count):
            l1 = len(number[number.find(".") + 1: number.find("(")])
            l2 = len(number[number.find("(") + 1: number.find(")")])
            temp += float(number[number.find("(") + 1:number.find(")")]) * 10 ** (-1 * (l1 + (i + 1) * l2))
        return temp


def transform(m: np.array, trans: np.array):
    """
    Perform coordinate transformation.
    """
    for i in range(m.shape[0]):
        m[i:i + 1] = np.matmul(m[i:i + 1], trans)


def rotation_matrix(alpha, beta, gamma):
    """
    Generate rotation matrix for alpha, beta and gamma angles (x, y, z).
    """
    x_rotation = np.array([[1, 0, 0], [0, np.cos(alpha), -1 * np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
    y_rotation = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-1 * np.sin(beta), 0, np.cos(beta)]])
    x_rotation = np.array([[np.cos(gamma), -1 * np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])
    return x_rotation, y_rotation, x_rotation


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
        self.internal_coord = np.zeros((n, 3))
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

    def __eq__(self, other):
        for i in range(self.num_atoms):
            diff = np.linalg.norm(self.atom_coord[i:i + 1] - other.atom_coord[i:i + 1])
            if diff < 0.05:
                return True
        return False

    def __ne__(self, other):
        if not self == other:
            return True
        else:
            return False

    def __add__(self, other):
        r = Molecule(self.num_atoms + other.num_atoms)
        for i in range(self.num_atoms):
            r.atom_label.append(self.atom_label[i])
            r.atom_coord[i:i+1] = self.atom_coord[i:i+1]
        for i in range(other.num_atoms):
            r.atom_label.append(other.atom_label[i])
            r.atom_coord[i+self.num_atoms:i+1+self.num_atoms] = other.atom_coord[i:i+1]
        return r

    def mass_center(self):
        r = np.zeros((1, 3))
        mass = 0.0
        for i in range(self.num_atoms):
            r = r + self.atom_coord[i:i + 1] * utility.dictionaries.element_weight[self.atom_label[i]]
            mass = mass + utility.dictionaries.element_weight[self.atom_label[i]]
        r = r / mass
        return r


def match_rotation(original_mol, mol_to_match: Molecule):
    for x in range(round(1000 * np.pi), 1000):
        for y in range(round(1000 * np.pi), 1000):
            for z in range(round(1000 * np.pi), 1000):
                rotated_mol = original_mol
                x_rot, y_rot, z_rot = rotation_matrix(x, y, z)
                transform(rotated_mol.atom_coord, x_rot)
                transform(rotated_mol.atom_coord, y_rot)
                transform(rotated_mol.atom_coord, z_rot)
                if rotated_mol == mol_to_match:
                    return x, y, z


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
        self.vector_a = np.transpose(np.matmul(self.transform, np.transpose([1, 0, 0])))
        self.vector_b = np.transpose(np.matmul(self.transform, np.transpose([0, 1, 0])))
        self.vector_c = np.transpose(np.matmul(self.transform, np.transpose([0, 0, 1])))
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
                        new_molecule = deepcopy(self.cif.asym_unit)
                        new_molecule.atom_coord *= mult
                        new_molecule.atom_coord += add
                        new_molecule.atom_coord += np.array([x, y, z])
                        coincide = False
                        for i2 in range(len(self.pre_molecules)):
                            if self.pre_molecules[i2] == new_molecule:
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
            limit_dist = utility.dictionaries.covalent_radius[self.pre_molecules[m1].atom_label[m1n]] + \
                         utility.dictionaries.covalent_radius[self.pre_molecules[m2].atom_label[m2n]]
            if dist <= limit_dist:
                self.bonds[line, k] = 1
                self.bonds[k, line] = 1

    def range_matrix_mass_centers(self):
        for i in range(len(self.molecules)):
            self.mass_centers.append(self.molecules.mass_center())
        self.range_matrix_mc = np.zeros((len(self.molecules), len(self.molecules)))
        for n in range(len(self.molecules)):
            for m in range(n + 1, len(self.molecules)):
                self.range_matrix_mc[n, m] = np.linalg.norm(self.mass_centers[n] - self.mass_centers[m])
                self.range_matrix_mc[m, n] = self.range_matrix_mc[m, n]

    def range_matrix_periodic(self):
        self.range_matrix_periodic = np.zeros((len(self.molecules), len(self.molecules)))
        for n in range(len(self.molecules)):
            for m in range(n + 1, len(self.molecules)):
                distances = []
                for x in range(-1, 2):
                    for y in range(3):
                        r = self.mass_centers[n] + self.out_translation[y]
                        distances.append(np.linalg.norm(self.mass_centers[n] - r))
                self.range_matrix_periodic[n, m] = min(distances)
                self.range_matrix_periodic[m, n] = self.range_matrix_periodic[n, m]

    def multiply(self, a: int, b: int, c: int):
        self.mass_centers = []
        for i in range(len(self.molecules)):
            self.mass_centers.append(self.molecules[i].mass_center())
        mass_centers_fract = []
        for i in range(len(self.mass_centers)):
            t = np.transpose(self.mass_centers[i])
            t = np.matmul(self.cif.rev_transform, t)
            mass_centers_fract.append(np.transpose(t))
        for i in range(len(self.molecules)):
            for x in range(a + 1):
                for y in range(b + 1):
                    for z in range(c + 1):
                        new_molecule = Molecule(self.molecules[i].num_atoms)
                        new_molecule = deepcopy(self.molecules[i])
                        trans = np.transpose(np.array([x, y, z]))
                        trans = np.matmul(self.cif.transform, trans)
                        trans = np.transpose(trans)
                        new_molecule.atom_coord = new_molecule.atom_coord + trans
                        coincide = False
                        for i2 in range(len(self.molecules)):
                            if self.molecules[i2] == new_molecule:
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
                    break
            for i1 in mol:
                for i2 in range(len(self.pre_molecules) * self.pre_molecules[0].num_atoms):
                    if self.bonds[i1, i2] == 1 and i2 not in mol:
                        mol.append(i2)
                checked[i1, 0] = 1
            flag = True
            for i in range(len(self.pre_molecules) * self.pre_molecules[0].num_atoms):
                if checked[i, 0] != 1:
                    flag = False
            new_molecule = Molecule(len(mol))
            for i in range(len(mol)):
                m = mol[i] // self.pre_molecules[0].num_atoms
                n = mol[i] % self.pre_molecules[0].num_atoms
                new_molecule.atom_label.append(self.pre_molecules[m].atom_label[n])
                new_molecule.atom_coord[i:i + 1] = self.pre_molecules[m].atom_coord[n:n + 1]
            self.molecules.append(new_molecule)
            mol.clear()
        mc = []
        mc_fract = []
        for i in range(len(self.molecules)):
            mc.append(self.molecules[i].mass_center())
            mc_t = np.transpose(mc[i])
            mc_t = np.matmul(self.cif.rev_transform, mc_t)
            mc_fract.append(np.transpose(mc_t))
        up_limit = 1.001
        low_limit = -0.001
        for_deletion = []
        new_mol = []
        for m in self.molecules:
            n = self.molecules.index(m)
            if mc_fract[n][0, 0] >= up_limit or mc_fract[n][0, 1] >= up_limit or mc_fract[n][0, 2] >= up_limit or mc_fract[n][0, 0] <= low_limit or mc_fract[n][0, 1] <= low_limit or mc_fract[n][0, 2] <= low_limit:
                for_deletion.append(n)
        for i in range(len(self.molecules)):
            if i not in for_deletion:
                new_mol.append(self.molecules[i])
        self.molecules = new_mol

    def to_cartesian(self, mol: Molecule):
        for i in range(mol.num_atoms):
            vector = mol.atom_coord[i:i + 1]
            vector = np.transpose(np.matmul(self.cif.transform, np.transpose(vector)))
            mol.atom_coord[i:i + 1] = vector

    def print_cluster_xyz(self):
        full_path = getcwd()
        file = open(full_path + "/cluster.xyz", "w")
        for n in range(len(self.molecules)):
            total_number = total_number + self.molecules[n].num_atoms
        file.write(repr(total_number) + "\n")
        file.write("XYZ file of molecular cluster generated in Kujo\n")
        for n in range(len(self.molecules)):
            for m in range(self.molecules[n].num_atoms):
                l = self.molecules[n].atom_label[m]
                x = repr(round(self.molecules[n].atom_coord[m, 0], 6))
                y = repr(round(self.molecules[n].atom_coord[m, 1], 6))
                z = repr(round(self.molecules[n].atom_coord[m, 2], 6))
                file.write(l + x + y + z + "\n")
        file.close()

    def print_cluster_readable(self):
        full_path = getcwd()
        file = open(full_path + "/cluster.kujo", "w")
        for n in range(len(self.molecules)):
            file.write(repr(self.molecules[n].num_atoms) + "\n")
            for m in range(self.molecules[n].num_atoms):
                l = self.molecules[n].atom_label[m]
                x = repr(round(self.molecules[n].atom_coord[m, 0], 6))
                y = repr(round(self.molecules[n].atom_coord[m, 1], 6))
                z = repr(round(self.molecules[n].atom_coord[m, 2], 6))
                file.write(l + x + y + z + "\n")
            for i1 in range(3):
                for i2 in range(3):
                    file.write(repr(round(self.molecules[n].inertia_eig_vec[i1, i2]), 6) + " ")
                file.write("\n")
        file.close()

    def __init__(self, cif: CifFile, a, b, c):
        self.pre_molecules = []
        self.cif = cif
        self.pre_molecules.append(self.cif.asym_unit)
        self.molecules = []
        self.bonds = np.zeros((len(self.pre_molecules) * self.pre_molecules[0].num_atoms, len(self.pre_molecules) *
                               self.pre_molecules[0].num_atoms))
        self.mass_centers = []
        self.out_translation = [(a + 1) * self.cif.vector_a, (b + 1) * self.cif.vector_b, (c + 1) * self.cif.vector_c]

