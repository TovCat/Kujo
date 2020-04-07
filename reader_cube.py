import numpy as np
import reader_cif as rc
import dictionaries as dic


class Cube:

    def __init__(self, path=""):
        try:
            file = open(path, "r")
        except OSError:
            print("Could not open the CUBE file at: ", path)
            exit(-1)
        contents = file.readlines()
        file.close()
        words = contents[2].split()
        self.num_atoms = int(words[0])
        self.origin = np.array([float(words[1]), float(words[2]), float(words[3])])
        self.steps = np.zeros((3, 1), dtype=int)
        self.volume = np.zeros((3, 3))
        for i in range(3):
            words = contents[i + 3].split()
            self.steps[i, 0] = int(words[0])
            self.volume[i, 0] = float(words[1])
            self.volume[i, 1] = float(words[2])
            self.volume[i, 2] = float(words[3])
        self.molecule = rc.Molecule(self.num_atoms)
        for i in range(6, 6 + self.num_atoms):
            words = contents[i].split()
            self.molecule.atom_label.append(dic.nrelements[int(words[0])])
            self.molecule.atom_coord[i:i + 1] = np.array([float(words[2]), float(words[3]), float(words[4])])
        self.voxels = np.zeros((self.steps[0, 0], self.steps[1, 0], self.steps[2, 0]))
        sep_contents = []
        for i1 in range(6 + self.num_atoms, len(contents)):
            words = contents[i1].split()
            for i2 in range(len(words)):
                sep_contents.append(float(words[i2]))
        for x in range(self.steps[0, 0]):
            for y in range(self.steps[1, 0]):
                for z in range(self.steps[2, 0]):
                    n = z + y * self.steps[1, 0] + x * self.steps[0, 0]
                    self.voxels[x, y, z] = sep_contents[n]
        self.dv = self.volume[0, 0] * self.volume[1, 1] * self.volume[2, 2]

    def integrate(self, low_limit, up_limit, translate: np.array, is_bohr=False):
        J = 0
        if not is_bohr:
            translate = translate / dic.bohr3
        for x1 in range(low_limit[0, 0], up_limit[0, 0]):
            for y1 in range(low_limit[1, 0], up_limit[1, 0]):
                for z1 in range(low_limit[2, 0], up_limit[2, 0]):
                    for x2 in range(self.steps[0, 0]):
                        for y2 in range(self.steps[1, 0]):
                            for z2 in range(self.steps[2, 0]):
                                v1 = np.array([x1 * self.volume[0, 0], y1 * self.volume[1, 0], z1 * self.volume[2, 0]]) + \
                                     self.origin
                                v2 = np.array([x2 * self.volume[0, 0], y2 * self.volume[1, 0], z2 * self.volume[2, 0]]) + \
                                     self.origin + translate
                                r = np.linalg.norm(v1 - v2)
                                J = J + (self.voxels[x1, y1, z1] * self.voxels[x2, y2, z2] / r) * self.dv
        J = J / dic.A
        return J