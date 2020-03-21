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
        self.steps = np.zeros((3, 1))
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
            self.molecule.atom_label[i - 6] = dic.nrelements[int(words[0])]
            self.molecule.atom_coord[i:i + 1] = np.array([float(words[2]), float(words[3], float(words[4]))])
        n = 0
        self.voxels = np.zeros((self.steps[0, 0], self.steps[1, 0], self.steps[2, 0]))
        for i in range(6 + self.num_atoms, len(contents)):
            words = contents[i].split()
            for i1 in range(len(words)):
                z = n % self.steps[2, 0]
                y = (n - z) % self.steps[1, 0]
                x = n - y * self.steps[1, 0] - z
                self.voxels[x, y, z] = float(words[i1])
                n = n + 1

