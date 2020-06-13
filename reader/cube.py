import numpy as np
import reader.cif as rc
import utility.dictionaries as dic
from copy import deepcopy


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
                    n = z + y * self.steps[2, 0] + x * self.steps[1, 0] * self.steps[2, 0]
                    self.voxels[x, y, z] = sep_contents[n]
        self.dv = self.volume[0, 0] * self.volume[1, 1] * self.volume[2, 2]
        self.grid = []
        for x in range(self.steps[0, 0]):
            for y in range(self.steps[1, 0]):
                for z in range(self.steps[2, 0]):
                    temp = self. origin + np.array([x * self.volume[0, 0], y * self.volume[1, 1], z * self.volume[2, 2]])
                    self.grid.append(temp)


def integrate_cubes(l: list):
    mol1 = l[0]
    mol2 = l[1]
    c1 = l[2]
    c2 = l[2]
    n = l[3]
    for x in len(c.grid):
        rc.transform(c1.grid[x], mol1.rotation)
        rc.transform(c2.grid[x]. mol2.rotation)
    for y1 in range(c.steps[1, 0]):
        for z1 in range(c.steps[2, 0]):
            for x2 in range(c.steps[0, 0]):
                for y2 in range(c.steps[1, 0]):
                    for z2 in range(c.steps[2, 0]):
                        i1 = z1 + y1 * c.steps[2, 0] + n * c.steps[2, 0] * c.steps[1, 0]
                        i2 = z2 + y2 * c.steps[2, 0] + x2 * c.steps[2, 0] * c.steps[1, 0]
                        r = np.linalg.norm(c1.grid[i1] - c2.grid[i2])
                        J = J + ((c.voxels[n, y1, z1] * c.dv) * (c.voxels[x2, y2, z2] * c.dv)) / r
    J = J * dic.A
    return J




