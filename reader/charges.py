import reader.cif
import numpy as np


class Charges:

    def __init__(self):
        self_q = np.zeros((3))

    def read(self, path=""):
        file = open(path, "r")
        contents = file.readlines()
        file.close()
        self.mol = reader.cif.Molecule(len(contents))
        self.q = np.zeros((len(contents)))
        for x in contents:
            words = x.split()
            self.mol.atom_label.append(words[0])
            index = contents.index(x)
            self.mol.atom_coord[index, 0] = float(words[1])
            self.mol.atom_coord[index, 1] = float(words[2])
            self.mol.atom_coord[index, 2] = float(words[3])
            self.q[index] = float(words[4])
